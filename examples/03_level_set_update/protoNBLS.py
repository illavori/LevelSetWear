# ======================================================================
# LevelSetWear, verification of the level set update (example 03)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Prototype of a narrow band level set (NBLS) wear update, without Abaqus.
#
# The worn surface of a flat is followed with two representations and the
# same wear model, and the profiles are compared:
#   HEIGHT : reference height function, table y_s(x), exact Hamilton-Jacobi
#            update for a graph, dy = dh * sqrt(1 + y'^2)
#   NBLS   : phi stored on a fixed Cartesian grid, signed distance, material
#            phi < 0; only a narrow band |phi| < BAND*dx is kept. Per step,
#            the four operations of Section 3.4 of the paper:
#              (i)   zero contour by marching squares (segments)
#              (ii)  Archard depth dh at points on the contour (the contact
#                    points), averaged with the hat kernel of radius LE,
#                    Eq. (kernel), at the segment ends
#              (iii) every grid node takes its distance d to the closest
#                    segment and the depth dh_ext at the closest point
#              (iv)  phi = sign(phi) d + dh_ext, Eq. (lsupdate), which also
#                    reinitialises phi at every update
# Wear model (cheap but self-consistent): rigid cylinder of radius R on a
# Winkler foundation, p = K * penetration, vertical position of the cylinder
# found by load balance (P per unit length), uniform slip per cycle (gross
# slip), Archard dh = k * p * s * DN.
#
#   python protoNBLS.py            -> prints the comparison, saves the figure

import time

import numpy as np
from scipy.spatial import cKDTree

# ------------------------------------------------------------ parameters
R, P, K = 6.0, 18.5, 1.67e5          # mm, N/mm, N/mm^3 (a0 ~ 0.1 mm)
KW, SLIP, DN = 2.9e-8, 0.1, 100      # MPa^-1 (flat of McColl et al.), mm per cycle, cycle jump
NSTEP = 180                          # 18 000 cycles
LE = 0.01                            # element size of the fretting model (mm), kernel radius
DX = LE/16                           # grid spacing of the level set (mm), as in the paper
X0, X1, Y0, Y1 = -0.4, 0.4, -0.02, 0.01
BAND = 6                             # narrow band half width (cells)


# ------------------------------------------------------------ wear model
def cylinderGap(c, x):
    """Vertical position of the cylinder bottom at x (c = lowest point)."""
    return c + R - np.sqrt(R*R - x*x)


def pressure(x, y, w):
    """Winkler pressure at surface points (x, y) with weights w (dx)."""
    lo, hi = y.max() - 0.05, y.max() + 0.05
    for _ in range(80):
        c = 0.5*(lo + hi)
        f = np.sum(K*np.maximum(y - cylinderGap(c, x), 0.0)*w)
        if f > P:
            lo = c
        else:
            hi = c
    return K*np.maximum(y - cylinderGap(c, x), 0.0)


# ------------------------------------------------------------ HEIGHT
class Height:
    def __init__(self):
        self.x = np.arange(X0, X1 + 0.5*DX, DX)
        self.y = np.zeros_like(self.x)

    def step(self):
        w = np.full_like(self.x, DX)
        dh = KW*pressure(self.x, self.y, w)*SLIP*DN
        s = np.gradient(self.y, self.x)
        self.y -= dh*np.sqrt(1.0 + s*s)


# ------------------------------------------------------------ NBLS
class NBLS:
    def __init__(self):
        self.x = np.arange(X0, X1 + 0.5*DX, DX)
        self.y = np.arange(Y0, Y1 + 0.5*DX, DX)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')
        self.phi = self.Y.copy()          # flat surface y = 0, material below
        self.clip()

    def clip(self):
        b = BAND*DX
        np.clip(self.phi, -b, b, out=self.phi)

    def contour(self):
        """Marching squares: segments (n, 2, 2) of phi = 0."""
        f = self.phi
        segs = []
        nx, ny = f.shape
        f00, f10 = f[:-1, :-1], f[1:, :-1]
        f01, f11 = f[:-1, 1:], f[1:, 1:]
        cut = ~(((f00 < 0) == (f10 < 0)) & ((f00 < 0) == (f01 < 0))
                & ((f00 < 0) == (f11 < 0)))
        for i, j in zip(*np.nonzero(cut)):
            xa, ya = self.x[i], self.y[j]
            v = (f[i, j], f[i+1, j], f[i+1, j+1], f[i, j+1])
            c = ((xa, ya), (xa+DX, ya), (xa+DX, ya+DX), (xa, ya+DX))
            pts = []
            for e in range(4):
                a, b = v[e], v[(e+1) % 4]
                if (a < 0) != (b < 0):
                    t = a/(a - b)
                    pa, pb = c[e], c[(e+1) % 4]
                    pts.append((pa[0] + t*(pb[0]-pa[0]),
                                pa[1] + t*(pb[1]-pa[1])))
            if len(pts) == 2:
                segs.append(pts)
            elif len(pts) == 4:                 # saddle: pair by centre sign
                if (0.25*sum(v) < 0) == (v[0] < 0):
                    segs += [[pts[0], pts[3]], [pts[1], pts[2]]]
                else:
                    segs += [[pts[0], pts[1]], [pts[2], pts[3]]]
        return np.array(segs)

    @staticmethod
    def resample(segs, n=4):
        """Points on the segments (n per segment), their weights for the load
        balance (|dx|) and the length of surface each one represents."""
        t = (np.arange(n) + 0.5)/n
        p = segs[:, None, 0, :] + t[None, :, None]*(segs[:, None, 1, :]
                                                    - segs[:, None, 0, :])
        w = np.abs(segs[:, 1, 0] - segs[:, 0, 0])[:, None]/n*np.ones((1, n))
        ell = np.linalg.norm(segs[:, 1] - segs[:, 0], axis=1)[:, None]/n*np.ones((1, n))
        return p.reshape(-1, 2), w.ravel(), ell.ravel()

    def step(self, k=None):
        # (i) current surface
        segs = self.contour()
        # (ii) Archard depth at the contact points, then the kernel average
        # of Eq. (kernel) at the segment ends, hat of radius LE
        pts, w, ell = self.resample(segs)
        dh = KW*pressure(pts[:, 0], pts[:, 1], w)*SLIP*DN
        ends = segs.reshape(-1, 2)
        S = cKDTree(ends).sparse_distance_matrix(cKDTree(pts), LE, output_type='coo_matrix')
        om = (1.0 - S.data/LE)*ell[S.col]
        num = np.bincount(S.row, om*dh[S.col], minlength=len(ends))
        den = np.bincount(S.row, om, minlength=len(ends))
        dhe = (num/np.maximum(den, 1e-30)).reshape(-1, 2)    # (nseg, 2)
        # (iii) distance of every node to the closest segment and depth at
        # the closest point (constant along the normals)
        q = np.c_[self.X.ravel(), self.Y.ravel()]
        nc = 6
        _, idx = cKDTree(segs.mean(axis=1)).query(q, k=nc)
        a = segs[idx, 0, :]                  # (nq, nc, 2)
        d = segs[idx, 1, :] - a
        t = np.einsum('ijk,ijk->ij', q[:, None, :] - a, d)
        t = np.clip(t/np.maximum(np.einsum('ijk,ijk->ij', d, d), 1e-30), 0, 1)
        r = a + t[..., None]*d - q[:, None, :]
        dist2 = np.einsum('ijk,ijk->ij', r, r)
        jm = dist2.argmin(axis=1)
        rows = np.arange(len(q))
        dist = np.sqrt(dist2[rows, jm])
        sm, tm = idx[rows, jm], t[rows, jm]
        dhx = (1.0 - tm)*dhe[sm, 0] + tm*dhe[sm, 1]
        # (iv) Eq. (lsupdate): signed distance plus the extended depth
        self.phi = (np.sign(self.phi.ravel())*dist + dhx).reshape(self.phi.shape)
        self.clip()

    def profile(self):
        """Surface height at the grid columns (first crossing from above)."""
        ys = np.full(len(self.x), np.nan)
        for i in range(len(self.x)):
            f = self.phi[i]
            j = np.nonzero((f[:-1] < 0) & (f[1:] >= 0))[0]
            if len(j):
                j = j[-1]
                ys[i] = self.y[j] + DX*f[j]/(f[j] - f[j+1])
        return ys


def main():
    h, ls = Height(), NBLS()
    th = tl = 0.0
    for k in range(NSTEP):
        t = time.perf_counter()
        h.step()
        th += time.perf_counter() - t
        t = time.perf_counter()
        ls.step(k)
        tl += time.perf_counter() - t
    yl = ls.profile()
    d = yl - h.y
    m = np.abs(h.x) < 0.3
    area_h = -np.trapz(h.y, h.x)
    area_l = -np.trapz(yl, ls.x)
    print('cycles %d, grid %d x %d, band +-%d cells, dx = %.2e mm'
          % (NSTEP*DN, len(ls.x), len(ls.y), BAND, DX))
    print('max depth  HEIGHT %.4f um   NBLS %.4f um'
          % (-1e3*h.y.min(), -1e3*np.nanmin(yl)))
    print('worn area  HEIGHT %.5f um*mm NBLS %.5f um*mm (%.2f %%)'
          % (1e3*area_h, 1e3*area_l, 100*(area_l - area_h)/area_h))
    print('profile difference RMS %.4f um, max %.4f um (%.2f %% of depth)'
          % (1e3*np.sqrt(np.nanmean(d[m]**2)), 1e3*np.nanmax(np.abs(d[m])),
             100*np.nanmax(np.abs(d[m]))/(-h.y.min())))
    print('time HEIGHT %.2f s, NBLS %.2f s' % (th, tl))
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax[0].plot(h.x, 1e3*h.y, '-', lw=2, color='#a8472c', label='height function')
        ax[0].plot(ls.x, 1e3*yl, '--', lw=1.4, color='#00a3ad', label='NBLS')
        ax[0].set_ylabel('y (um)')
        ax[0].legend()
        ax[0].set_title('worn flat after %d cycles' % (NSTEP*DN))
        ax[1].plot(h.x, 1e3*d, color='k', lw=1)
        ax[1].set_ylabel('NBLS - height (um)')
        ax[1].set_xlabel('x (mm)')
        ax[1].set_xlim(-0.3, 0.3)
        fig.tight_layout()
        fig.savefig('protoNBLS.png', dpi=150)
        print('saved protoNBLS.png')
    except ImportError:
        pass


if __name__ == '__main__':
    main()
