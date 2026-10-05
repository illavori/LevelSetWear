C ======================================================================
C     LevelSetWear, level set wear user elements for Abaqus/Standard
C     Version 5 (the code used in the article)
C
C     Written by Iñigo Llavori
C
C     I. Llavori, Wear simulation on a fixed mesh with level sets
C     and immersed frictional contact,
C     submitted for publication. Mondragon Unibertsitatea.
C     https://github.com/illavori/LevelSetWear       BSD 3-Clause
C
C     The worn surfaces of two bodies in contact (body 1 = flat = slave,
C     body 2 = cylinder = master) are level sets on a fixed mesh.
C     Units mm, N, MPa. 2D plane strain. Run with cpus=1 (the elements
C     share data through COMMON blocks). See LevelSetWear_Manual.pdf.
C
C     JTYPE = 1 (U1): 4-node plane strain quad cut by the surface of
C                     body IB, fixed NFIX x NFIX cell rule (see UBULK).
C     JTYPE = 2 (U2): immersed contact element, one per interior vertex
C                     of the slave surface (see UCONT).
C     UMAT          : ghost CPE4 mesh for the visualisation of the U1
C                     results (SDV1-4 stresses, SDV5 volume fraction).
C     UEXTERNALDB   : end of every converged increment, Archard wear
C                     update of both surfaces. End of every step,
C                     contact results (contactOut.txt) and surface
C                     profiles (wearProfile.txt).
C
C     Options in v2opts.txt of the job folder ("NAME value" per line,
C     all 0 by default, which gives the point-wise contact):
C       MORTAR 0 : one penalty and Coulomb law per Gauss point
C       MORTAR 1 : one law per vertex, gap weighted with the hat function
C       MORTAR 2 : one law per vertex, gap weighted with the dual
C                  function 3N-1 (dual mortar, used in the article)
C       WSM 1    : wear of the contact points averaged over their vertex
C                  (implied by MORTAR >= 1)
C       DIAG 1   : contact state changes per iteration (v2diag.txt),
C                  point law only (MORTAR 0)
C
C     Surfaces, one table per body, read once:
C       surf1.txt (body 1, flat), surf2.txt (body 2, cylinder)
C       "N NSUB", then N lines "x y_s", x increasing. Every NSUB-th
C       point (1, 1+NSUB, ...) is a surface node of the mesh: that coarse
C       grid carries the wear depth, the other points the true shape.
C     Level set of body IB: phi = SIDE*(y - y_s(x)), material phi < 0.
C
C     Narrow band level set: if lsgrid.txt exists in the job folder,
C     the surface of each body is a signed distance function on a fixed
C     Cartesian grid, kept in a narrow band (general level set, no height
C     function). This is the mode used for the wear analyses. One line
C     per body:
C       X0 X1 Y0 Y1 DX NBAND NREINIT HKER
C     (grid box, spacing, band half width in cells, redistancing every
C     NREINIT updates, radius of the kernel that carries the Archard
C     depth of the contact points to the surface). phi is initialised
C     from the tables. Outside the grid box the tables are still used
C     and they are not updated, so the box must contain all the wear.
C     Without lsgrid.txt the surfaces are the height tables themselves,
C     updated vertically (used for the contact verification).
C       PHIF  : phi at a point (bilinear on the grid)
C       DSURF : vector from a point to the closest surface point
C       LSUPD : wear update (contour, extension, redistancing)
C       YLS   : surface height from phi (output of the profiles)
C ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     &     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     &     SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     &     DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     &     JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     &     PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      LOGICAL LSTIF, LRES
      INTEGER*8 KT0, KT1, KTR
      DOUBLE PRECISION TBU, TCO, TEX, TST
      INTEGER NBU, NCO, NEX
      COMMON /V2TIM/ TBU, TCO, TEX, TST, NBU, NCO, NEX
C
      LSTIF = LFLAGS(3).EQ.1 .OR. LFLAGS(3).EQ.2
      LRES  = LFLAGS(3).EQ.1 .OR. LFLAGS(3).EQ.5
      DO I = 1, NDOFEL
        RHS(I,1) = 0.D0
        DO J = 1, NDOFEL
          AMATRX(I,J) = 0.D0
        END DO
      END DO
      IF (.NOT.(LSTIF .OR. LRES)) RETURN
      CALL SREAD
C
      CALL SYSTEM_CLOCK(KT0,KTR)
      IF (JTYPE.EQ.1) THEN
        CALL UBULK(RHS,AMATRX,SVARS,ENERGY,NDOFEL,MLVARX,PROPS,
     &       COORDS,U,LSTIF,LRES,JELEM)
        CALL SYSTEM_CLOCK(KT1)
        TBU = TBU + DBLE(KT1-KT0)/DBLE(KTR)
        NBU = NBU + 1
      ELSE IF (JTYPE.EQ.2) THEN
        CALL UCONT(RHS,AMATRX,SVARS,NDOFEL,NNODE,MLVARX,PROPS,
     &       NPROPS,COORDS,U,DU,LSTIF,LRES,JELEM,KSTEP,KINC,DTIME)
        CALL SYSTEM_CLOCK(KT1)
        TCO = TCO + DBLE(KT1-KT0)/DBLE(KTR)
        NCO = NCO + 1
      ELSE
        WRITE(7,*) 'UEL: unknown element type ', JTYPE
        CALL XIT
      END IF
      RETURN
      END
C
C ======================================================================
      SUBROUTINE UBULK(RHS,AMATRX,SVARS,ENERGY,NDOFEL,MLVARX,PROPS,
     &     COORDS,U,LSTIF,LRES,JELEM)
C
C     Cut 4-node quad (B-bar as CPE4). Element class by nodal phi:
C       inside : full integration
C       void   : EPSV * full stiffness (2x2)
C       cut    : material part only, + EPSV * full stiffness
C     The material part of a cut element is the polygon of the parent
C     domain on the material side of the straight cut between the two
C     edge crossings. Two integration schemes (PROPS(6) = NFIX):
C       NFIX = 0 : polygon split in triangles, 6-point rule each
C                  (points move with the surface); inside = 2x2
C       NFIX > 0 : FIXED points, the NFIX x NFIX Gauss rule. The rule
C                  splits the parent square in cells whose sides are the
C                  Gauss weights, each containing its point. The weight
C                  of a point is the material area of its cell, so the
C                  stiffness changes continuously as the surface moves,
C                  and a point keeps its history while its cell has
C                  material (no remapping of the history).
C                  Inside = NFIX x NFIX Gauss.
C     PROPS: (1) E, (2) nu, (3) SIDE, (4) body IB, (5) EPSV, (6) NFIX
C     SVARS: (1:4) material-averaged S11 S22 S33 S12, (5) vol. fraction
C            NFIX > 0: (5 + k) material fraction of the cell of point k
C
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MAXP=64, MAXEL=200000, MAXF=8)
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),
     &     ENERGY(8),PROPS(*),COORDS(2,4),U(*)
      DIMENSION XIN(4),ETN(4),PHI(4)
      DIMENSION PXI(MAXP),PET(MAXP),PW(MAXP)
      DIMENSION QXI(6),QET(6),GXI(4),GET(4),GW(4)
      DIMENSION TL1(6),TL2(6),TW(6),SAVG(4),SVOID(4)
      DIMENSION XG(MAXF),WGT(MAXF),CB(0:MAXF)
      LOGICAL LSTIF, LRES
      DOUBLE PRECISION GOUT
      COMMON /GHOST/ GOUT(5,MAXEL)
C
      DATA XIN / -1.D0,  1.D0,  1.D0, -1.D0 /
      DATA ETN / -1.D0, -1.D0,  1.D0,  1.D0 /
      DATA TL1 / 0.445948490915965D0, 0.445948490915965D0,
     &           0.108103018168070D0, 0.091576213509771D0,
     &           0.091576213509771D0, 0.816847572980459D0 /
      DATA TL2 / 0.445948490915965D0, 0.108103018168070D0,
     &           0.445948490915965D0, 0.091576213509771D0,
     &           0.816847572980459D0, 0.091576213509771D0 /
      DATA TW  / 0.223381589678011D0, 0.223381589678011D0,
     &           0.223381589678011D0, 0.109951743655322D0,
     &           0.109951743655322D0, 0.109951743655322D0 /
C
      E    = PROPS(1)
      PNU  = PROPS(2)
      SIDE = PROPS(3)
      IB   = NINT(PROPS(4))
      EPSV = PROPS(5)
      NFIX = NINT(PROPS(6))
      IF (NFIX.GT.MAXF) THEN
        WRITE(7,*) 'UBULK: NFIX > ', MAXF
        CALL XIT
      END IF
C
C     ---- 2x2 Gauss rule (CPE4 order) and the fixed rule
      G = 1.D0/SQRT(3.D0)
      GXI(1) = -G
      GET(1) = -G
      GXI(2) =  G
      GET(2) = -G
      GXI(3) = -G
      GET(3) =  G
      GXI(4) =  G
      GET(4) =  G
      DO K = 1, 4
        GW(K) = 1.D0
      END DO
      IF (NFIX.GT.0) THEN
        CALL GAULEG(NFIX,XG,WGT)
        CB(0) = -1.D0
        DO I = 1, NFIX
          CB(I) = CB(I-1) + WGT(I)
        END DO
      END IF
C
      NNEG = 0
      NPOS = 0
      DO I = 1, 4
        PHI(I) = PHIF(IB,SIDE,COORDS(1,I),COORDS(2,I))
        IF (ABS(PHI(I)).LT.1.D-9) PHI(I) = 0.D0
        IF (PHI(I).LT.0.D0) NNEG = NNEG + 1
        IF (PHI(I).GT.0.D0) NPOS = NPOS + 1
      END DO
C
      ENERGY(2) = 0.D0
      IF (NPOS.EQ.0) THEN
C       ---- inside
        IF (NFIX.EQ.0) THEN
          CALL KBBAR(COORDS,U,4,GXI,GET,GW,E,PNU,1.D0,LSTIF,LRES,
     &         AMATRX,RHS,NDOFEL,MLVARX,ENER,SAVG,VMAT)
        ELSE
          NP = 0
          DO J = 1, NFIX
            DO I = 1, NFIX
              NP = NP + 1
              PXI(NP) = XG(I)
              PET(NP) = XG(J)
              PW(NP)  = WGT(I)*WGT(J)
              SVARS(5+NP) = 1.D0
            END DO
          END DO
          CALL KBBAR(COORDS,U,NP,PXI,PET,PW,E,PNU,1.D0,LSTIF,LRES,
     &         AMATRX,RHS,NDOFEL,MLVARX,ENER,SAVG,VMAT)
        END IF
        FRAC = 1.D0
      ELSE IF (NNEG.EQ.0) THEN
C       ---- void
        CALL KBBAR(COORDS,U,4,GXI,GET,GW,E,PNU,EPSV,LSTIF,LRES,
     &       AMATRX,RHS,NDOFEL,MLVARX,ENER,SAVG,VMAT)
        DO I = 1, 4
          SAVG(I) = 0.D0
        END DO
        DO K = 1, NFIX*NFIX
          SVARS(5+K) = 0.D0
        END DO
        FRAC = 0.D0
      ELSE
C       ---- cut: material polygon in the parent domain
        NV = 0
        DO I = 1, 4
          J = MOD(I,4) + 1
          IF (PHI(I).LE.0.D0) THEN
            NV = NV + 1
            QXI(NV) = XIN(I)
            QET(NV) = ETN(I)
          END IF
          IF ((PHI(I).LT.0.D0 .AND. PHI(J).GT.0.D0) .OR.
     &        (PHI(I).GT.0.D0 .AND. PHI(J).LT.0.D0)) THEN
            T = PHI(I)/(PHI(I)-PHI(J))
            NV = NV + 1
            QXI(NV) = XIN(I) + T*(XIN(J)-XIN(I))
            QET(NV) = ETN(I) + T*(ETN(J)-ETN(I))
          END IF
        END DO
        NP = 0
        IF (NFIX.EQ.0) THEN
C         moving points: triangles of the polygon, 6-point rule
          DO IT = 2, NV-1
            X1 = QXI(1)
            Y1 = QET(1)
            X2 = QXI(IT)
            Y2 = QET(IT)
            X3 = QXI(IT+1)
            Y3 = QET(IT+1)
            AT = 0.5D0*ABS((X2-X1)*(Y3-Y1) - (X3-X1)*(Y2-Y1))
            DO K = 1, 6
              T3 = 1.D0 - TL1(K) - TL2(K)
              NP = NP + 1
              PXI(NP) = TL1(K)*X1 + TL2(K)*X2 + T3*X3
              PET(NP) = TL1(K)*Y1 + TL2(K)*Y2 + T3*Y3
              PW(NP)  = TW(K)*AT
            END DO
          END DO
        ELSE
C         fixed points: weight = material area of the point's cell
          K = 0
          DO J = 1, NFIX
            DO I = 1, NFIX
              K = K + 1
              AC = CLIPA(NV,QXI,QET,CB(I-1),CB(I),CB(J-1),CB(J))
              SVARS(5+K) = AC/(WGT(I)*WGT(J))
              IF (AC.GT.1.D-14) THEN
                NP = NP + 1
                PXI(NP) = XG(I)
                PET(NP) = XG(J)
                PW(NP)  = AC
              END IF
            END DO
          END DO
        END IF
        IF (NP.GT.0) THEN
          CALL KBBAR(COORDS,U,NP,PXI,PET,PW,E,PNU,1.D0,LSTIF,LRES,
     &         AMATRX,RHS,NDOFEL,MLVARX,ENER,SAVG,VMAT)
        ELSE
          ENER = 0.D0
          VMAT = 0.D0
          DO I = 1, 4
            SAVG(I) = 0.D0
          END DO
        END IF
        ENERGY(2) = ENER
        CALL KBBAR(COORDS,U,4,GXI,GET,GW,E,PNU,EPSV,LSTIF,LRES,
     &       AMATRX,RHS,NDOFEL,MLVARX,ENER,SVOID,VTOT)
        FRAC = VMAT/VTOT
      END IF
      ENERGY(2) = ENERGY(2) + ENER
C
      DO I = 1, 4
        SVARS(I) = SAVG(I)
      END DO
      SVARS(5) = FRAC
      IF (JELEM.GT.MAXEL) THEN
        WRITE(7,*) 'UEL: element ', JELEM, ' > MAXEL = ', MAXEL
        CALL XIT
      END IF
      DO I = 1, 5
        GOUT(I,JELEM) = SVARS(I)
      END DO
      RETURN
      END
C
C ======================================================================
      DOUBLE PRECISION FUNCTION CLIPA(NV,QX,QY,X0,X1,Y0,Y1)
C     Area of the convex polygon (QX, QY) clipped by the rectangle
C     [X0, X1] x [Y0, Y1] (Sutherland-Hodgman, then shoelace).
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MV=16)
      DIMENSION QX(*),QY(*),AX(MV),AY(MV),BX(MV),BY(MV)
C
      N = NV
      DO I = 1, N
        AX(I) = QX(I)
        AY(I) = QY(I)
      END DO
      DO IE = 1, 4
        M = 0
        DO I = 1, N
          J = MOD(I,N) + 1
          IF (IE.EQ.1) THEN
            DI = AX(I) - X0
            DJ = AX(J) - X0
          ELSE IF (IE.EQ.2) THEN
            DI = X1 - AX(I)
            DJ = X1 - AX(J)
          ELSE IF (IE.EQ.3) THEN
            DI = AY(I) - Y0
            DJ = AY(J) - Y0
          ELSE
            DI = Y1 - AY(I)
            DJ = Y1 - AY(J)
          END IF
          IF (DI.GE.0.D0) THEN
            M = M + 1
            BX(M) = AX(I)
            BY(M) = AY(I)
          END IF
          IF ((DI.GE.0.D0 .AND. DJ.LT.0.D0) .OR.
     &        (DI.LT.0.D0 .AND. DJ.GE.0.D0)) THEN
            T = DI/(DI-DJ)
            M = M + 1
            BX(M) = AX(I) + T*(AX(J)-AX(I))
            BY(M) = AY(I) + T*(AY(J)-AY(I))
          END IF
        END DO
        N = M
        IF (N.LT.3) THEN
          CLIPA = 0.D0
          RETURN
        END IF
        DO I = 1, N
          AX(I) = BX(I)
          AY(I) = BY(I)
        END DO
      END DO
      S = 0.D0
      DO I = 1, N
        J = MOD(I,N) + 1
        S = S + AX(I)*AY(J) - AX(J)*AY(I)
      END DO
      CLIPA = 0.5D0*ABS(S)
      RETURN
      END
C
C ======================================================================
      SUBROUTINE GAULEG(N,X,W)
C     Gauss-Legendre points and weights on [-1, 1], increasing x.
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION X(N),W(N)
      PI = 4.D0*ATAN(1.D0)
      M = (N+1)/2
      DO I = 1, M
        Z = COS(PI*(DBLE(I)-0.25D0)/(DBLE(N)+0.5D0))
        DO IT = 1, 100
          P1 = 1.D0
          P2 = 0.D0
          DO J = 1, N
            P3 = P2
            P2 = P1
            P1 = ((2.D0*J-1.D0)*Z*P2 - (J-1.D0)*P3)/J
          END DO
          PP = N*(Z*P1 - P2)/(Z*Z - 1.D0)
          Z1 = Z
          Z = Z1 - P1/PP
          IF (ABS(Z-Z1).LT.1.D-15) GOTO 10
        END DO
   10   X(I) = -Z
        X(N+1-I) = Z
        W(I) = 2.D0/((1.D0-Z*Z)*PP*PP)
        W(N+1-I) = W(I)
      END DO
      RETURN
      END
C
C ======================================================================
      SUBROUTINE UCONT(RHS,AMATRX,SVARS,NDOFEL,NNODE,MLVARX,PROPS,
     &     NPROPS,COORDS,U,DU,LSTIF,LRES,JELEM,KSTEP,KINC,DTIME)
C
C     Embedded contact element, one per vertex j of the slave surface
C     polyline (body 1, flat, material below, SIDE = +1) against the
C     master surface (body 2, cylinder, material above, SIDE = -1).
C
C     The element integrates the two slave segments that meet at j
C     (3 Gauss points each) weighted with the hat function N_j, so the
C     sum over all elements integrates every segment once, and it
C     reports vertex (N_j-weighted average) values like a nodal output.
C     At each Gauss point:
C       - pairing with the closest master segment in the configuration
C         at the START of the increment (U - DU), fixed during the
C         increment (small sliding, updated every increment);
C       - gap to the master chord, corrected with the surface tables
C         (true surfaces y_s(x) vs element chords, both bodies): this
C         removes the faceting of the discretised surfaces;
C       - normal: p = EPSN * <-g>;
C       - tangential: Coulomb, return mapping from the converged
C         traction of the point.
C
C     Connectivity (K = PROPS(4) element rows, M = PROPS(5) lines):
C       slave  node (r, s), r = 0..K from the surface inwards,
C              s = 0, 1, 2 node lines (vertex j on line 1)
C                                           -> local 3*r + s + 1
C       master node (r, c), r = 0..K from the surface inwards,
C              c = 0..M-1 left to right     -> 3*(K+1) + c*(K+1) + r + 1
C     PROPS: (1) EPSN, (2) EPST, (3) mu, (4) K, (5) M, (6) element
C            number offset of the contact elements
C            wear (optional, NPROPS >= 11): (7) k slave, (8) k master,
C            (9) cycle jump (fixed), (10) first wear step, (11) HCRIT
C            (surface change per increment above which a warning is
C            written; no cutback: a PNEWDT request in the iterations
C            reacts to transient Newton iterates and stalls the step)
C            (12, optional) imposed sliding per unit step time (mm). If
C            given and positive, the wear uses p * (12) * DTIME instead of
C            the computed slip increment: continuous unidirectional
C            sliding (block-on-ring) where every point of the contact
C            slides the same known distance, and only the normal contact
C            is solved (Cruzado et al. 2010)
C     For the wear update (UEXTERNALDB) every point stores in /WEARC/
C     its reference x (slave and master), p * slip increment, weight.
C     SVARS: (1:6) vertex x_ref, p, tau, accumulated slip, gap, status
C            (0 open, 1 stick, 2 slip), N_j-weighted averages
C            (6+3k-2, 6+3k-1, 6+3k) tau, accumulated slip and p of point
C            k = 1..6 (24 in total)
C            On entry: values at the start of the increment.
C
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MXN=128, MXC=30, MAXCE=5000, NGS=3)
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),
     &     PROPS(*),COORDS(2,*),U(*),DU(MLVARX,*)
      DIMENSION XS0(2,MXN),GV(2*MXN),TV(2*MXN)
      DIMENSION IQS(4),IQM(4),X4(2,4),SN(4),SM(4),TG(NGS),WG(NGS)
      DIMENSION SEGS(4),SEGM(4,MXC),IROWM(MXC),IOKM(MXC)
      DIMENSION PA(2),PB(2),PS(2),PR(2),D(2),EN(2),ET(2)
      DIMENSION DSV(2),DMV(2)
      LOGICAL LSTIF, LRES, LFND, LEDGE
      DOUBLE PRECISION CO
      INTEGER NCE
      COMMON /CONTO/ CO(6,MAXCE), NCE
      PARAMETER (NPW=6)
      DOUBLE PRECISION DTA
      INTEGER IACT, KEYA
      COMMON /ACTS/ DTA(MAXCE), IACT(NPW,MAXCE), KEYA(MAXCE),
     &     NCALL(MAXCE)
      INTEGER NCALL
      LOGICAL LACT, LTEN
      DOUBLE PRECISION WXS, WXM, WPS, WWT, AKW, CJMP, DYLAST, HCRW
      INTEGER KWEAR, KSLAST
      COMMON /WEARC/ WXS(NPW,MAXCE), WXM(NPW,MAXCE), WPS(NPW,MAXCE),
     &     WWT(NPW,MAXCE), AKW(2), CJMP, DYLAST, HCRW, KWEAR, KSLAST
      DOUBLE PRECISION WYS, WYM
      COMMON /WEARY/ WYS(NPW,MAXCE), WYM(NPW,MAXCE)
      INTEGER IWSM, IDIAG, ISTL, NFL, NFS, NCLO
      DOUBLE PRECISION XFMIN, XFMAX
      COMMON /V2OPT/ XFMIN, XFMAX, IWSM, IDIAG, ISTL(NPW,MAXCE),
     &     NFL(8), NFS(8), NCLO, IMOR
      DIMENSION GVW(2*MXN), TVW(2*MXN)
      LOGICAL LHOLD
C     warning (once per element) when a surface reaches the last element
C     row of the contact element or has left its rows: the worn surface
C     is about to leave the region that the element can see (increase
C     the rows K, property 4 of U2)
      INTEGER IWROW(MAXCE)
      SAVE IWROW
      DATA IWROW /MAXCE*0/
C
      VSL = 0.D0
      IF (NPROPS.GE.12) VSL = PROPS(12)
      IF (NPROPS.GE.11) THEN
        AKW(1) = PROPS(7)
        AKW(2) = PROPS(8)
        CJMP   = PROPS(9)
        KWEAR  = NINT(PROPS(10))
        HCRW   = PROPS(11)
      ELSE
        KWEAR  = 0
      END IF
      LWEAR = 0
      IF (KWEAR.GT.0 .AND. KSTEP.GE.KWEAR) LWEAR = 1
C
      EPSN = PROPS(1)
      EPST = PROPS(2)
      FMU  = PROPS(3)
      KR   = NINT(PROPS(4))
      MM   = NINT(PROPS(5))
      IC   = JELEM - NINT(PROPS(6))
      NSL  = 3*(KR+1)
      IF (NNODE.GT.MXN .OR. MM.GT.MXC .OR. IC.LT.1 .OR.
     &    IC.GT.MAXCE) THEN
        WRITE(7,*) 'UCONT: array limits exceeded, element ', JELEM
        CALL XIT
      END IF
      NCE = MAX(NCE, IC)
C     Active set, as the Lagrange contact of Abaqus: in the first two
C     calls (iterations) of every attempt, the points that were in
C     contact at the end of the previous increment are treated as
C     bilateral (g = 0 enforced by the penalty, tension allowed, no
C     friction). After a wear update the surfaces are apart by the worn
C     depth; with a prescribed load a plain penalty would leave the
C     bodies without support (or, with stiffness only, close a small
C     fraction of the gap). From the third call on: plain unilateral
C     contact, so the converged state has no tension. A new attempt
C     (new KSTEP/KINC, or a cutback: new DTIME) restarts from SVARS.
      KEY = 100000*KSTEP + KINC
      IF (KEYA(IC).NE.KEY .OR. DTA(IC).NE.DTIME) THEN
        DO KP = 1, NPW
          IACT(KP,IC) = 0
          IF (SVARS(6+3*KP).GT.0.D0) IACT(KP,IC) = 1
        END DO
        KEYA(IC) = KEY
        DTA(IC) = DTIME
        NCALL(IC) = 0
      END IF
      IF (LRES) NCALL(IC) = NCALL(IC) + 1
      LTEN = NCALL(IC).LE.2
C     3-point Gauss rule on [0, 1]
      TG(1) = 0.5D0 - 0.5D0*SQRT(0.6D0)
      TG(2) = 0.5D0
      TG(3) = 0.5D0 + 0.5D0*SQRT(0.6D0)
      WG(1) = 5.D0/18.D0
      WG(2) = 8.D0/18.D0
      WG(3) = 5.D0/18.D0
C
C     ---- positions at the start of the increment
      DO IA = 1, NNODE
        XS0(1,IA) = COORDS(1,IA) + U(2*IA-1) - DU(2*IA-1,1)
        XS0(2,IA) = COORDS(2,IA) + U(2*IA)   - DU(2*IA,1)
      END DO
C
C     ---- master surface segments of the patch
      DO IC2 = 0, MM-2
        IOKM(IC2+1) = 0
        DO IR = 0, KR-1
          IF (IOKM(IC2+1).EQ.0) THEN
            CALL MQUAD(KR,NSL,IR,IC2,IQM)
            DO I = 1, 4
              X4(1,I) = COORDS(1,IQM(I))
              X4(2,I) = COORDS(2,IQM(I))
            END DO
            CALL SURFSEG(X4,2,-1.D0,LFND,IVOID,SEGM(1,IC2+1))
            IF (LFND) THEN
              IOKM(IC2+1) = 1
              IROWM(IC2+1) = IR
            ELSE IF (IVOID.EQ.0) THEN
              IOKM(IC2+1) = -1
            END IF
          END IF
        END DO
      END DO
      IF (IWROW(IC).EQ.0) THEN
        DO IC2 = 1, MM-1
          IF ((IOKM(IC2).EQ.1 .AND. IROWM(IC2).EQ.KR-1) .OR.
     &        IOKM(IC2).EQ.0) THEN
            IF (IWROW(IC).EQ.0) WRITE(7,'(A,I7,A,I3,A)')
     &        ' UCONT: WARNING master surface of contact element',
     &        JELEM, ' in or beyond its last row (', KR,
     &        ' rows): increase the rows of the contact elements'
            IWROW(IC) = 1
          END IF
        END DO
      END IF
C
C     ---- Gauss points of the two slave segments
      AJ = 0.D0
      PJ = 0.D0
      TJ = 0.D0
      SJ = 0.D0
      GJ = 0.D0
      XV = 0.D0
      NLED = 0
C     v3, MORTAR = 1: accumulators of the vertex law
      AWM = 0.D0
      GWM = 0.D0
      DGWM = 0.D0
      IF (IMOR.GE.1) THEN
        DO I = 1, NDOFEL
          GVW(I) = 0.D0
          TVW(I) = 0.D0
        END DO
      END IF
      DO ICOL = 0, 1
        LFND = .FALSE.
        IRS = -1
        DO IR = 0, KR-1
          IF (.NOT.LFND) THEN
            IQS(1) = 3*(IR+1) + ICOL + 1
            IQS(2) = 3*(IR+1) + ICOL + 2
            IQS(3) = 3*IR + ICOL + 2
            IQS(4) = 3*IR + ICOL + 1
            DO I = 1, 4
              X4(1,I) = COORDS(1,IQS(I))
              X4(2,I) = COORDS(2,IQS(I))
            END DO
            CALL SURFSEG(X4,1,1.D0,LFND,IVOID,SEGS)
            IF (LFND) IRS = IR
            IF (.NOT.LFND .AND. IVOID.EQ.0) GOTO 60
          END IF
        END DO
   60   CONTINUE
        IF (IWROW(IC).EQ.0 .AND. ((LFND .AND. IRS.EQ.KR-1) .OR.
     &      (.NOT.LFND .AND. IVOID.EQ.1))) THEN
          IWROW(IC) = 1
          WRITE(7,'(A,I7,A,I3,A)') ' UCONT: WARNING slave surface of'
     &      //' contact element', JELEM, ' in or beyond its last row (',
     &      KR, ' rows): increase the rows of the contact elements'
        END IF
        IF (.NOT.LFND) THEN
          DO K = 1, NGS
            KP = ICOL*NGS + K
            SVARS(6+3*KP-2) = 0.D0
            SVARS(6+3*KP-1) = 0.D0
            SVARS(6+3*KP)   = 0.D0
            IACT(KP,IC) = 0
            WWT(KP,IC) = 0.D0
            WPS(KP,IC) = 0.D0
          END DO
          GOTO 90
        END IF
        CALL SHP4(SEGS(1),SEGS(2),SN)
        CALL XAT(SN,IQS,COORDS,PA)
        CALL SHP4(SEGS(3),SEGS(4),SN)
        CALL XAT(SN,IQS,COORDS,PB)
        SLEN = SQRT((PB(1)-PA(1))**2 + (PB(2)-PA(2))**2)
C       vertex j: right end of the left segment, left end of the right
        IF (ICOL.EQ.0) THEN
          XV = PB(1)
        ELSE IF (XV.EQ.0.D0) THEN
          XV = PA(1)
        END IF
        DO K = 1, NGS
          KP = ICOL*NGS + K
          TAUN  = SVARS(6+3*KP-2)
          SLIPN = SVARS(6+3*KP-1)
          SLIP0 = SLIPN
          WWT(KP,IC) = 0.D0
          WPS(KP,IC) = 0.D0
          XIS = SEGS(1) + TG(K)*(SEGS(3)-SEGS(1))
          ETS = SEGS(2) + TG(K)*(SEGS(4)-SEGS(2))
          IF (ICOL.EQ.0) THEN
            HAT = TG(K)
          ELSE
            HAT = 1.D0 - TG(K)
          END IF
          W = WG(K)*SLEN*HAT
          CALL SHP4(XIS,ETS,SN)
          CALL XAT(SN,IQS,XS0,PS)
C         -- pairing: closest master segment (start of increment)
          BEST = 1.D30
          ICB = 0
          DO IC2 = 1, MM-1
            IF (IOKM(IC2).EQ.1) THEN
              CALL MQUAD(KR,NSL,IROWM(IC2),IC2-1,IQM)
              CALL SHP4(SEGM(1,IC2),SEGM(2,IC2),SM)
              CALL XAT(SM,IQM,XS0,PA)
              CALL SHP4(SEGM(3,IC2),SEGM(4,IC2),SM)
              CALL XAT(SM,IQM,XS0,PB)
              D(1) = PB(1) - PA(1)
              D(2) = PB(2) - PA(2)
              DL2 = D(1)**2 + D(2)**2
              T = ((PS(1)-PA(1))*D(1) + (PS(2)-PA(2))*D(2))/DL2
              TC = MIN(1.D0, MAX(0.D0, T))
              DIST = (PA(1)+TC*D(1)-PS(1))**2
     &             + (PA(2)+TC*D(2)-PS(2))**2
              IF (DIST.LT.BEST) THEN
                BEST = DIST
                ICB = IC2
                TB = TC
                LEDGE = (T.LT.-1.D-6 .AND. IC2.EQ.1) .OR.
     &                  (T.GT.1.D0+1.D-6 .AND. IC2.EQ.MM-1)
                DLEN = SQRT(DL2)
                ET(1) = D(1)/DLEN
                ET(2) = D(2)/DLEN
              END IF
            END IF
          END DO
          P = 0.D0
          TAU = 0.D0
          GK = 1.D0
          IF (ICB.GT.0) THEN
            IF (LEDGE) NLED = NLED + 1
            EN(1) = -ET(2)
            EN(2) =  ET(1)
            CALL MQUAD(KR,NSL,IROWM(ICB),ICB-1,IQM)
            XIM = SEGM(1,ICB) + TB*(SEGM(3,ICB)-SEGM(1,ICB))
            ETM = SEGM(2,ICB) + TB*(SEGM(4,ICB)-SEGM(2,ICB))
            CALL SHP4(XIM,ETM,SM)
            DO I = 1, NDOFEL
              GV(I) = 0.D0
              TV(I) = 0.D0
            END DO
            DO I = 1, 4
              IM = IQM(I)
              GV(2*IM-1) = GV(2*IM-1) + SM(I)*EN(1)
              GV(2*IM)   = GV(2*IM)   + SM(I)*EN(2)
              TV(2*IM-1) = TV(2*IM-1) + SM(I)*ET(1)
              TV(2*IM)   = TV(2*IM)   + SM(I)*ET(2)
              IS = IQS(I)
              GV(2*IS-1) = GV(2*IS-1) - SN(I)*EN(1)
              GV(2*IS)   = GV(2*IS)   - SN(I)*EN(2)
              TV(2*IS-1) = TV(2*IS-1) - SN(I)*ET(1)
              TV(2*IS)   = TV(2*IS)   - SN(I)*ET(2)
            END DO
C           gap to the chords with total displacements, slip with DU
            GK = 0.D0
            DGK = 0.D0
            DO IA = 1, NNODE
              DO J = 1, 2
                IDF = 2*(IA-1) + J
                GK = GK + GV(IDF)*(COORDS(J,IA) + U(IDF))
                DGK = DGK + TV(IDF)*DU(IDF,1)
              END DO
            END DO
C           faceting correction: true surfaces vs chords (reference);
C           height tables: vertical offset, NBLS: closest surface point
            CALL XAT(SN,IQS,COORDS,PR)
            CALL DSURF(1,1.D0,PR(1),PR(2),DSV)
            WXS(KP,IC) = PR(1) + DSV(1)
            WYS(KP,IC) = PR(2) + DSV(2)
            CALL XAT(SM,IQM,COORDS,PR)
            CALL DSURF(2,-1.D0,PR(1),PR(2),DMV)
            WXM(KP,IC) = PR(1) + DMV(1)
            WYM(KP,IC) = PR(2) + DMV(2)
            GK = GK + EN(1)*(DMV(1) - DSV(1)) + EN(2)*(DMV(2) - DSV(2))
            IF (IMOR.GE.1) THEN
C             v5: MORTAR 2 weights gap and slip with the dual function
C             Phi_j = 3 N_j - 1 (HAT = N_j), D_j = sum of W as before
              WD = W
              IF (IMOR.EQ.2) WD = WG(K)*SLEN*(3.D0*HAT - 1.D0)
              AWM = AWM + W
              GWM = GWM + WD*GK
              DGWM = DGWM + WD*DGK
              DO I = 1, NDOFEL
                GVW(I) = GVW(I) + WD*GV(I)
                TVW(I) = TVW(I) + WD*TV(I)
              END DO
              WWT(KP,IC) = W
              GOTO 94
            END IF
            LACT = IACT(KP,IC).EQ.1 .AND. LTEN
            IF (GK.LT.0.D0) THEN
              P = -EPSN*GK
              TTR = TAUN + EPST*DGK
              IF (ABS(TTR).LE.FMU*P) THEN
                TAU = TTR
                IST = 1
              ELSE
                SGN = SIGN(1.D0, TTR)
                TAU = FMU*P*SGN
                SLIPN = SLIPN + (ABS(TTR) - FMU*P)/EPST
                IST = 2
              END IF
              IF (LRES) THEN
                DO I = 1, NDOFEL
                  RHS(I,1) = RHS(I,1) + W*(P*GV(I) - TAU*TV(I))
                END DO
              END IF
              IF (LSTIF) THEN
                DO I = 1, NDOFEL
                  DO J = 1, NDOFEL
                    AK = EPSN*GV(I)*GV(J)
                    IF (IST.EQ.1) THEN
                      AK = AK + EPST*TV(I)*TV(J)
                    ELSE
                      AK = AK - FMU*SGN*EPSN*TV(I)*GV(J)
                    END IF
                    AMATRX(I,J) = AMATRX(I,J) + W*AK
                  END DO
                END DO
              END IF
            ELSE IF (LACT) THEN
C             bilateral in the first iterations: tension, no friction
              PT = -EPSN*GK
              IF (LRES) THEN
                DO I = 1, NDOFEL
                  RHS(I,1) = RHS(I,1) + W*PT*GV(I)
                END DO
              END IF
              IF (LSTIF) THEN
                DO I = 1, NDOFEL
                  DO J = 1, NDOFEL
                    AMATRX(I,J) = AMATRX(I,J) + W*EPSN*GV(I)*GV(J)
                  END DO
                END DO
              END IF
            END IF
C           v2, DIAG = 1: state of the point after the contact law,
C           0 open, 1 stick, 2 slip, 3 bilateral, compared between calls
            IF (IDIAG.EQ.1 .AND. LRES) THEN
              JS = 0
              IF (GK.LT.0.D0) THEN
                JS = IST
              ELSE IF (LACT) THEN
                JS = 3
              END IF
              IF (NCALL(IC).GT.1 .AND. JS.NE.ISTL(KP,IC)) THEN
                JC = MIN(NCALL(IC),8)
                NFL(JC) = NFL(JC) + 1
                IF (JS*ISTL(KP,IC).EQ.2) NFS(JC) = NFS(JC) + 1
                XFMIN = MIN(XFMIN, ABS(PS(1)))
                XFMAX = MAX(XFMAX, ABS(PS(1)))
              END IF
              ISTL(KP,IC) = JS
            END IF
C           wear data of the point (Archard: p times slip increment)
            WWT(KP,IC) = W
            WPS(KP,IC) = P*(SLIPN - SLIP0)
            IF (VSL.GT.0.D0) WPS(KP,IC) = P*VSL*DTIME
            AJ = AJ + W
            PJ = PJ + W*P
            TJ = TJ + W*TAU
            SJ = SJ + W*SLIPN
            GJ = GJ + W*GK
          END IF
          SVARS(6+3*KP-2) = TAU
          SVARS(6+3*KP-1) = SLIPN
          SVARS(6+3*KP)   = P
   94     CONTINUE
        END DO
   90   CONTINUE
      END DO
C     v3, MORTAR = 1: contact law of the vertex
      IF (IMOR.GE.1) THEN
        AJ = 0.D0
        PJ = 0.D0
        TJ = 0.D0
        SJ = 0.D0
        GJ = 0.D0
        PVM = 0.D0
        TAUM = SVARS(7)
        SLPM = SVARS(8)
        SLP0M = SLPM
        IF (AWM.GT.0.D0) THEN
          GBM = GWM/AWM
          LHOLD = LTEN .AND. SVARS(6).GT.0.D0
          TAUM = 0.D0
          IF (GBM.LT.0.D0 .OR. LHOLD) THEN
            PVM = -EPSN*GBM
            PBM = FMU*MAX(PVM, 0.D0)
            IF (LHOLD .AND. GBM.GE.0.D0) PBM = FMU*MAX(SVARS(2), 0.D0)
            TTR = SVARS(7) + EPST*DGWM/AWM
            IF (ABS(TTR).LE.PBM) THEN
              TAUM = TTR
              ISTM = 1
            ELSE
              SGN = SIGN(1.D0, TTR)
              TAUM = PBM*SGN
              SLPM = SLPM + (ABS(TTR) - PBM)/EPST
              ISTM = 2
            END IF
            IF (LRES) THEN
              DO I = 1, NDOFEL
                RHS(I,1) = RHS(I,1) + PVM*GVW(I) - TAUM*TVW(I)
              END DO
            END IF
            IF (LSTIF) THEN
              DO I = 1, NDOFEL
                DO J = 1, NDOFEL
                  AK = EPSN*GVW(I)*GVW(J)
                  IF (ISTM.EQ.1) THEN
                    AK = AK + EPST*TVW(I)*TVW(J)
                  ELSE IF (PVM.GT.0.D0 .AND. .NOT.(LHOLD .AND.
     &                     GBM.GE.0.D0)) THEN
                    AK = AK - FMU*SGN*EPSN*TVW(I)*GVW(J)
                  END IF
                  AMATRX(I,J) = AMATRX(I,J) + AK/AWM
                END DO
              END DO
            END IF
          END IF
          AJ = AWM
          PJ = AWM*PVM
          TJ = AWM*TAUM
          SJ = AWM*SLPM
          GJ = GWM
        END IF
        DSM = SLPM - SLP0M
        IF (VSL.GT.0.D0) DSM = VSL*DTIME
        DO KP = 1, NPW
          SVARS(6+3*KP-2) = TAUM
          SVARS(6+3*KP-1) = SLPM
          SVARS(6+3*KP)   = PVM
          IF (WWT(KP,IC).GT.0.D0) WPS(KP,IC) = MAX(PVM, 0.D0)*DSM
        END DO
      END IF
C     v2, WSM = 1: every point of the vertex takes the N_j-weighted
C     average of p*ds over the points of the vertex (same total wear)
      IF (IWSM.EQ.1 .AND. AJ.GT.0.D0) THEN
        PSJ = 0.D0
        DO KP = 1, NPW
          PSJ = PSJ + WWT(KP,IC)*WPS(KP,IC)
        END DO
        DO KP = 1, NPW
          IF (WWT(KP,IC).GT.0.D0) WPS(KP,IC) = PSJ/AJ
        END DO
      END IF
      IF (NLED.GT.0) WRITE(7,*) 'UCONT: pairing at patch edge, ',
     &     'element ', JELEM
C
C     ---- vertex output (N_j-weighted averages)
      SVARS(1) = XV
      DO I = 2, 6
        SVARS(I) = 0.D0
      END DO
      SVARS(5) = 1.D0
      IF (AJ.GT.0.D0) THEN
        SVARS(2) = PJ/AJ
        SVARS(3) = TJ/AJ
        SVARS(4) = SJ/AJ
        SVARS(5) = GJ/AJ
        IF (PJ.GT.0.D0) THEN
          SVARS(6) = 1.D0
          IF (ABS(TJ).GE.0.999D0*FMU*PJ) SVARS(6) = 2.D0
        END IF
      END IF
      DO I = 1, 6
        CO(I,IC) = SVARS(I)
      END DO
      RETURN
      END
C
C ======================================================================
      SUBROUTINE MQUAD(KR,NSL,IR,IC,IQ)
C     Local node numbers (counterclockwise) of the master element in
C     row IR (0 = surface) between master lines IC and IC+1. The master
C     body lies above its surface.
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION IQ(4)
      IQ(1) = NSL + IC*(KR+1) + IR + 1
      IQ(2) = NSL + (IC+1)*(KR+1) + IR + 1
      IQ(3) = NSL + (IC+1)*(KR+1) + IR + 2
      IQ(4) = NSL + IC*(KR+1) + IR + 2
      RETURN
      END
C
C ======================================================================
      SUBROUTINE SURFSEG(X4,IB,SIDE,LFND,IVOID,SEG)
C     Surface segment (phi = 0) of a quad in parent coordinates,
C     SEG = (xi1, eta1, xi2, eta2) ordered by increasing reference x.
C     LFND  : the element contains a surface segment
C     IVOID : 1 if the element has no material
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION X4(2,4),SEG(4),PHI(4),XIN(4),ETN(4),QX(6),QE(6),SN(4)
      LOGICAL LFND
      DATA XIN / -1.D0,  1.D0,  1.D0, -1.D0 /
      DATA ETN / -1.D0, -1.D0,  1.D0,  1.D0 /
C
      TOLP = 1.D-9
      NNEG = 0
      DO I = 1, 4
        PHI(I) = PHIF(IB,SIDE,X4(1,I),X4(2,I))
        IF (ABS(PHI(I)).LT.TOLP) PHI(I) = 0.D0
        IF (PHI(I).LT.0.D0) NNEG = NNEG + 1
      END DO
      LFND = .FALSE.
      IVOID = 0
      IF (NNEG.EQ.0) THEN
        IVOID = 1
        RETURN
      END IF
      NQ = 0
      DO I = 1, 4
        J = MOD(I,4) + 1
        IF (PHI(I).EQ.0.D0) THEN
          NQ = NQ + 1
          QX(NQ) = XIN(I)
          QE(NQ) = ETN(I)
        END IF
        IF (PHI(I)*PHI(J).LT.0.D0) THEN
          T = PHI(I)/(PHI(I)-PHI(J))
          NQ = NQ + 1
          QX(NQ) = XIN(I) + T*(XIN(J)-XIN(I))
          QE(NQ) = ETN(I) + T*(ETN(J)-ETN(I))
        END IF
      END DO
      IF (NQ.LT.2) RETURN
      IF (NQ.GT.2) WRITE(7,*) 'SURFSEG: ', NQ, ' surface points'
      LFND = .TRUE.
      CALL SHP4(QX(1),QE(1),SN)
      XA = 0.D0
      DO I = 1, 4
        XA = XA + SN(I)*X4(1,I)
      END DO
      CALL SHP4(QX(2),QE(2),SN)
      XB = 0.D0
      DO I = 1, 4
        XB = XB + SN(I)*X4(1,I)
      END DO
      IF (XA.LE.XB) THEN
        SEG(1) = QX(1)
        SEG(2) = QE(1)
        SEG(3) = QX(2)
        SEG(4) = QE(2)
      ELSE
        SEG(1) = QX(2)
        SEG(2) = QE(2)
        SEG(3) = QX(1)
        SEG(4) = QE(1)
      END IF
      RETURN
      END
C
C ======================================================================
      SUBROUTINE SHP4(XI,ET,SN)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SN(4)
      SN(1) = 0.25D0*(1.D0-XI)*(1.D0-ET)
      SN(2) = 0.25D0*(1.D0+XI)*(1.D0-ET)
      SN(3) = 0.25D0*(1.D0+XI)*(1.D0+ET)
      SN(4) = 0.25D0*(1.D0-XI)*(1.D0+ET)
      RETURN
      END
C
C ======================================================================
      SUBROUTINE XAT(SN,IQ,XN,P)
C     P = sum SN(i) * XN(:, IQ(i))
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SN(4),IQ(4),XN(2,*),P(2)
      P(1) = 0.D0
      P(2) = 0.D0
      DO I = 1, 4
        P(1) = P(1) + SN(I)*XN(1,IQ(I))
        P(2) = P(2) + SN(I)*XN(2,IQ(I))
      END DO
      RETURN
      END
C
C ======================================================================
      SUBROUTINE KBBAR(COORDS,U,NP,PXI,PET,PW,E,PNU,FACT,LSTIF,LRES,
     &     AMATRX,RHS,NDOFEL,MLVARX,ENER,SAVG,VOL)
C
C     Adds FACT * (B-bar stiffness and internal force) integrated over
C     the NP parent-domain points with weights PW (no Jacobian).
C     In-plane volumetric strain averaged over the points (as CPE4).
C
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MAXP=64, NS=4, ND=8)
      DIMENSION COORDS(2,4),U(*),PXI(*),PET(*),PW(*)
      DIMENSION AMATRX(NDOFEL,NDOFEL),RHS(MLVARX,*),SAVG(NS)
      DIMENSION XIN(4),ETN(4),DNX(4,MAXP),DNY(4,MAXP),DV(MAXP)
      DIMENSION BV(ND),B(NS,ND),DMAT(NS,NS),DB(NS,ND),EPS(NS),SIG(NS)
      LOGICAL LSTIF, LRES
      DATA XIN / -1.D0,  1.D0,  1.D0, -1.D0 /
      DATA ETN / -1.D0, -1.D0,  1.D0,  1.D0 /
C
      FAC = E/((1.D0+PNU)*(1.D0-2.D0*PNU))
      DO I = 1, NS
        SAVG(I) = 0.D0
        DO J = 1, NS
          DMAT(I,J) = 0.D0
        END DO
      END DO
      DO I = 1, 3
        DO J = 1, 3
          DMAT(I,J) = FAC*PNU
        END DO
        DMAT(I,I) = FAC*(1.D0-PNU)
      END DO
      DMAT(4,4) = FAC*(1.D0-2.D0*PNU)/2.D0
C
      VOL = 0.D0
      DO K = 1, NP
        XJ11 = 0.D0
        XJ12 = 0.D0
        XJ21 = 0.D0
        XJ22 = 0.D0
        DO I = 1, 4
          DXI = 0.25D0*XIN(I)*(1.D0+ETN(I)*PET(K))
          DET = 0.25D0*ETN(I)*(1.D0+XIN(I)*PXI(K))
          XJ11 = XJ11 + DXI*COORDS(1,I)
          XJ12 = XJ12 + DXI*COORDS(2,I)
          XJ21 = XJ21 + DET*COORDS(1,I)
          XJ22 = XJ22 + DET*COORDS(2,I)
        END DO
        DETJ = XJ11*XJ22 - XJ12*XJ21
        IF (DETJ.LE.0.D0) THEN
          WRITE(7,*) 'UEL: non-positive Jacobian'
          CALL XIT
        END IF
        DO I = 1, 4
          DXI = 0.25D0*XIN(I)*(1.D0+ETN(I)*PET(K))
          DET = 0.25D0*ETN(I)*(1.D0+XIN(I)*PXI(K))
          DNX(I,K) = ( XJ22*DXI - XJ12*DET)/DETJ
          DNY(I,K) = (-XJ21*DXI + XJ11*DET)/DETJ
        END DO
        DV(K) = PW(K)*DETJ
        VOL = VOL + DV(K)
      END DO
C
      DO I = 1, 4
        BV(2*I-1) = 0.D0
        BV(2*I)   = 0.D0
        DO K = 1, NP
          BV(2*I-1) = BV(2*I-1) + DNX(I,K)*DV(K)/VOL
          BV(2*I)   = BV(2*I)   + DNY(I,K)*DV(K)/VOL
        END DO
      END DO
C
      ENER = 0.D0
      DO K = 1, NP
        DO J = 1, ND
          DO I = 1, NS
            B(I,J) = 0.D0
          END DO
        END DO
        DO I = 1, 4
          B(1,2*I-1) = DNX(I,K)
          B(2,2*I)   = DNY(I,K)
          B(4,2*I-1) = DNY(I,K)
          B(4,2*I)   = DNX(I,K)
        END DO
        DO J = 1, ND
          IF (MOD(J,2).EQ.1) THEN
            BVOLK = DNX((J+1)/2,K)
          ELSE
            BVOLK = DNY(J/2,K)
          END IF
          CORR = (BV(J)-BVOLK)/2.D0
          B(1,J) = B(1,J) + CORR
          B(2,J) = B(2,J) + CORR
        END DO
C
        DO I = 1, NS
          EPS(I) = 0.D0
          DO J = 1, ND
            EPS(I) = EPS(I) + B(I,J)*U(J)
          END DO
        END DO
        DO I = 1, NS
          SIG(I) = 0.D0
          DO J = 1, NS
            SIG(I) = SIG(I) + DMAT(I,J)*EPS(J)
          END DO
          SAVG(I) = SAVG(I) + SIG(I)*DV(K)/VOL
          ENER = ENER + 0.5D0*FACT*SIG(I)*EPS(I)*DV(K)
        END DO
C
        W = FACT*DV(K)
        IF (LRES) THEN
          DO J = 1, ND
            DO I = 1, NS
              RHS(J,1) = RHS(J,1) - B(I,J)*SIG(I)*W
            END DO
          END DO
        END IF
        IF (LSTIF) THEN
          DO J = 1, ND
            DO I = 1, NS
              DB(I,J) = 0.D0
              DO L = 1, NS
                DB(I,J) = DB(I,J) + DMAT(I,L)*B(L,J)
              END DO
            END DO
          END DO
          DO I = 1, ND
            DO J = 1, ND
              DO L = 1, NS
                AMATRX(I,J) = AMATRX(I,J) + B(L,I)*DB(L,J)*W
              END DO
            END DO
          END DO
        END IF
      END DO
      RETURN
      END
C
C ======================================================================
      SUBROUTINE SREAD
C     Reads the surface tables of both bodies once per job.
C     (INQUIRE first: inside Abaqus a failed OPEN aborts the run.)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (NTM=5001)
      DOUBLE PRECISION XT, YT
      INTEGER NT, ITRD, NSB
      COMMON /SURFT/ XT(NTM,2), YT(NTM,2), NT(2), ITRD, NSB(2)
      CHARACTER*256 OUTDIR
      CHARACTER*9 FNAME(2)
      INTEGER LENOUT
      LOGICAL LEX
      DATA FNAME / 'surf1.txt', 'surf2.txt' /
C
      IF (ITRD.EQ.1) RETURN
      CALL GETOUTDIR(OUTDIR,LENOUT)
      DO IB = 1, 2
        INQUIRE(FILE=OUTDIR(1:LENOUT)//'\'//FNAME(IB),EXIST=LEX)
        IF (.NOT.LEX) THEN
          WRITE(7,*) 'UEL: missing surface table ', FNAME(IB)
          CALL XIT
        END IF
        OPEN(69,FILE=OUTDIR(1:LENOUT)//'\'//FNAME(IB),STATUS='OLD')
        READ(69,*) NT(IB), NSB(IB)
        IF (NT(IB).GT.NTM .OR. NT(IB).LT.2) THEN
          WRITE(7,*) 'UEL: bad number of points in ', FNAME(IB)
          CALL XIT
        END IF
        DO I = 1, NT(IB)
          READ(69,*) XT(I,IB), YT(I,IB)
        END DO
        CLOSE(69)
        WRITE(7,*) 'UEL: surface table ', FNAME(IB), NT(IB)
      END DO
      CALL LSINIT(OUTDIR,LENOUT)
      ITRD = 1
      RETURN
      END
C
C ======================================================================
      DOUBLE PRECISION FUNCTION YSURF(IB,X)
C     Surface height of body IB at x: linear between table points,
C     constant beyond the ends.
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (NTM=5001)
      DOUBLE PRECISION XT, YT
      INTEGER NT, ITRD, NSB
      COMMON /SURFT/ XT(NTM,2), YT(NTM,2), NT(2), ITRD, NSB(2)
C
      N = NT(IB)
      IF (X.LE.XT(1,IB)) THEN
        YSURF = YT(1,IB)
      ELSE IF (X.GE.XT(N,IB)) THEN
        YSURF = YT(N,IB)
      ELSE
        IL = 1
        IR = N
   10   IF (IR-IL.GT.1) THEN
          IM = (IL+IR)/2
          IF (XT(IM,IB).LE.X) THEN
            IL = IM
          ELSE
            IR = IM
          END IF
          GOTO 10
        END IF
        T = (X-XT(IL,IB))/(XT(IR,IB)-XT(IL,IB))
        YSURF = YT(IL,IB) + T*(YT(IR,IB)-YT(IL,IB))
      END IF
      RETURN
      END
C
C ======================================================================
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
C     LOP = 0 : open contactOut.txt and wearProfile.txt
C     LOP = 2 : end of every converged increment. In the wear steps
C               (KSTEP >= KWEAR) Archard wear depth at every contact
C               point (slave x for body 1, paired master x for body 2)
C                 dh = k_IB * p * |ds| * cycle jump
C               With lsgrid.txt (narrow band level set) the depths are
C               carried to the zero contour with a kernel of radius
C               HKER, extended to the band and added to phi (LSUPD).
C               Without it they are projected on the coarse grid of body
C               IB (its surface nodes) with hat weights, interpolated to
C               the table points and applied along the normal (the flat
C               surface moves down, the cylinder surface up).
C     LOP = 1 : start of increment; when the step changes, output of
C               the previous step (contact results and profiles)
C     LOP = 3 : end of the analysis, output of the last step
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION TIME(2)
      PARAMETER (MAXCE=5000, NTM=5001, NPW=6)
      DOUBLE PRECISION CO, XT, YT, WXS, WXM, WPS, WWT, AKW, CJMP, DYLAST
      DOUBLE PRECISION HCRW
      INTEGER NCE, NT, ITRD, NSB, KWEAR, KSLAST
      COMMON /CONTO/ CO(6,MAXCE), NCE
      COMMON /SURFT/ XT(NTM,2), YT(NTM,2), NT(2), ITRD, NSB(2)
      COMMON /WEARC/ WXS(NPW,MAXCE), WXM(NPW,MAXCE), WPS(NPW,MAXCE),
     &     WWT(NPW,MAXCE), AKW(2), CJMP, DYLAST, HCRW, KWEAR, KSLAST
      DIMENSION AG(NTM), BG(NTM), XG(NTM), DY(NTM)
      CHARACTER*256 OUTDIR
      INTEGER LENOUT
      SAVE TLAST
      PARAMETER (MGX=8001, MGY=241)
      DOUBLE PRECISION PG, GX0, GY0, GDX, HKER
      INTEGER NGX, NGY, NBND, NREI, NUPD, LSON
      COMMON /LSG/ PG(MGX,MGY,2), GX0(2), GY0(2), GDX(2), HKER(2),
     &     NGX(2), NGY(2), NBND(2), NREI(2), NUPD(2), LSON
      DOUBLE PRECISION WYS, WYM
      COMMON /WEARY/ WYS(NPW,MAXCE), WYM(NPW,MAXCE)
      LOGICAL LDUMP
      SAVE LDUMP
      INTEGER IWSM, IDIAG, ISTL, NFL, NFS, NCLO
      DOUBLE PRECISION XFMIN, XFMAX
      COMMON /V2OPT/ XFMIN, XFMAX, IWSM, IDIAG, ISTL(NPW,MAXCE),
     &     NFL(8), NFS(8), NCLO, IMOR
      INTEGER*8 KT0, KT1, KTR
      DOUBLE PRECISION TBU, TCO, TEX, TST
      INTEGER NBU, NCO, NEX
      COMMON /V2TIM/ TBU, TCO, TEX, TST, NBU, NCO, NEX
      CHARACTER*16 ONAME
      LOGICAL LOPT
C
      IF (LOP.EQ.0) THEN
        CALL GETOUTDIR(OUTDIR,LENOUT)
C       v2: options, timers and diagnostic
        IWSM = 0
        IDIAG = 0
        IMOR = 0
        INQUIRE(FILE=OUTDIR(1:LENOUT)//'\v2opts.txt',EXIST=LOPT)
        IF (LOPT) THEN
          OPEN(75,FILE=OUTDIR(1:LENOUT)//'\v2opts.txt',STATUS='OLD')
   5      READ(75,*,END=6) ONAME, IVAL
          IF (ONAME(1:3).EQ.'WSM') IWSM = IVAL
          IF (ONAME(1:4).EQ.'DIAG') IDIAG = IVAL
          IF (ONAME(1:6).EQ.'MORTAR') IMOR = IVAL
          GOTO 5
   6      CLOSE(75)
        END IF
        WRITE(7,'(A,I2,A,I2,A,I2)') ' UEL version 5: WSM =', IWSM,
     &    ', DIAG =', IDIAG, ', MORTAR =', IMOR
        TBU = 0.D0
        TCO = 0.D0
        TEX = 0.D0
        NBU = 0
        NCO = 0
        NEX = 0
        CALL SYSTEM_CLOCK(KT0,KTR)
        TST = DBLE(KT0)/DBLE(KTR)
        DO I = 1, 8
          NFL(I) = 0
          NFS(I) = 0
        END DO
        XFMIN = 1.D30
        XFMAX = 0.D0
        IF (IDIAG.EQ.1) THEN
          OPEN(76,FILE=OUTDIR(1:LENOUT)//'\v2diag.txt',STATUS='UNKNOWN')
          WRITE(76,'(A)') '# step inc closed chg(2..8+) ss(2..8+) '//
     &      '|x|min |x|max'
        END IF
        OPEN(71,FILE=OUTDIR(1:LENOUT)//'\contactOut.txt',
     &       STATUS='UNKNOWN')
        OPEN(72,FILE=OUTDIR(1:LENOUT)//'\wearProfile.txt',
     &       STATUS='UNKNOWN')
C       diagnostic: contact points of every converged wear increment
C       (x, y of slave and master, weight, p * slip), if wdump.on exists
        INQUIRE(FILE=OUTDIR(1:LENOUT)//'\wdump.on',EXIST=LDUMP)
        IF (LDUMP) OPEN(74,FILE=OUTDIR(1:LENOUT)//'\wearDump.txt',
     &       STATUS='UNKNOWN')
        KSLAST = 0
        TLAST = 0.D0
        DYLAST = 0.D0
      ELSE IF (LOP.EQ.1) THEN
        IF (KSTEP.NE.KSLAST .AND. KSLAST.GT.0) THEN
          WRITE(7,'(A,I4,3(A,F9.2,A,I9))') ' TIMERS end of step',
     &      KSLAST, ' | UBULK', TBU, ' s, calls', NBU, ' | UCONT', TCO,
     &      ' s, calls', NCO, ' | wear update', TEX, ' s, calls', NEX
          CALL SYSTEM_CLOCK(KT0,KTR)
          WRITE(7,'(A,I4,A,F10.2,A)') ' WALL end of step', KSLAST,
     &      ':', DBLE(KT0)/DBLE(KTR) - TST, ' s since the start'
          CALL WOUT(KSLAST,TLAST)
        END IF
        KSLAST = KSTEP
      ELSE IF (LOP.EQ.2) THEN
        CALL SYSTEM_CLOCK(KT0,KTR)
        TLAST = TIME(2)
        IF (IDIAG.EQ.1) THEN
          NCLO = 0
          DO IC = 1, NCE
            DO K = 1, NPW
              IF (ISTL(K,IC).EQ.1 .OR. ISTL(K,IC).EQ.2) NCLO = NCLO+1
            END DO
          END DO
          IF (XFMIN.GT.1.D29) XFMIN = -1.D0
          WRITE(76,'(2I6,I6,14I6,2F10.4)') KSTEP, KINC, NCLO,
     &      (NFL(I), I = 2, 8), (NFS(I), I = 2, 8), XFMIN, XFMAX
          DO I = 1, 8
            NFL(I) = 0
            NFS(I) = 0
          END DO
          XFMIN = 1.D30
          XFMAX = 0.D0
        END IF
        IF (LDUMP .AND. KWEAR.GT.0 .AND. KSTEP.GE.KWEAR) THEN
          NQ = 0
          DO IC = 1, NCE
            DO K = 1, NPW
              IF (WWT(K,IC).GT.0.D0) NQ = NQ + 1
            END DO
          END DO
          WRITE(74,'(A,3I8,2ES14.6)') 'INC', KSTEP, KINC, NQ,
     &      AKW(1)*CJMP, AKW(2)*CJMP
          DO IC = 1, NCE
            DO K = 1, NPW
              IF (WWT(K,IC).GT.0.D0) WRITE(74,'(6ES20.12)')
     &          WXS(K,IC), WYS(K,IC), WXM(K,IC), WYM(K,IC),
     &          WWT(K,IC), WPS(K,IC)
            END DO
          END DO
        END IF
        IF (KWEAR.GT.0 .AND. KSTEP.GE.KWEAR .AND. ITRD.EQ.1 .AND.
     &      LSON.EQ.1) THEN
C         NBLS: update of the grid level sets
          DYLAST = 0.D0
          DO IB = 1, 2
            SIDE = 1.D0
            IF (IB.EQ.2) SIDE = -1.D0
            CALL LSUPD(IB,SIDE,DYM)
            DYLAST = DYLAST + DYM
            IF (DYM.GT.HCRW) WRITE(7,'(A,I2,A,ES11.3,A,I4,A,I6)')
     &        'UEXTERNALDB: WARNING body', IB, ' surface change', DYM,
     &        ' > HCRIT, step', KSTEP, ' inc', KINC
          END DO
        ELSE IF (KWEAR.GT.0 .AND. KSTEP.GE.KWEAR .AND. ITRD.EQ.1) THEN
          DYLAST = 0.D0
          DO IB = 1, 2
            NS = MAX(NSB(IB),1)
            NG = (NT(IB)-1)/NS + 1
            DO M = 1, NG
              XG(M) = XT(1+(M-1)*NS,IB)
              AG(M) = 0.D0
              BG(M) = 0.D0
            END DO
C           wear depth of every point, projected on the coarse grid
            DO IC = 1, NCE
              DO K = 1, NPW
                W = WWT(K,IC)
                IF (W.GT.0.D0) THEN
                  IF (IB.EQ.1) THEN
                    X = WXS(K,IC)
                  ELSE
                    X = WXM(K,IC)
                  END IF
                  DH = AKW(IB)*CJMP*WPS(K,IC)
                  IF (X.GT.XG(1) .AND. X.LT.XG(NG)) THEN
                    IL = 1
                    IR = NG
   10               IF (IR-IL.GT.1) THEN
                      IM = (IL+IR)/2
                      IF (XG(IM).LE.X) THEN
                        IL = IM
                      ELSE
                        IR = IM
                      END IF
                      GOTO 10
                    END IF
                    T = (X-XG(IL))/(XG(IR)-XG(IL))
                    AG(IL) = AG(IL) + (1.D0-T)*W*DH
                    BG(IL) = BG(IL) + (1.D0-T)*W
                    AG(IR) = AG(IR) + T*W*DH
                    BG(IR) = BG(IR) + T*W
                  END IF
                END IF
              END DO
            END DO
            DO M = 1, NG
              IF (BG(M).GT.0.D0) THEN
                AG(M) = AG(M)/BG(M)
              ELSE
                AG(M) = 0.D0
              END IF
            END DO
C           table points: interpolated depth, along the normal
            DO I = 1, NT(IB)
              M = (I-1)/NS + 1
              T = DBLE(MOD(I-1,NS))/DBLE(NS)
              DH = AG(M)*(1.D0-T)
              IF (M.LT.NG) DH = DH + AG(M+1)*T
              I1 = MAX(1,I-1)
              I2 = MIN(NT(IB),I+1)
              S = (YT(I2,IB)-YT(I1,IB))/(XT(I2,IB)-XT(I1,IB))
              DY(I) = DH*SQRT(1.D0+S*S)
            END DO
            DYM = 0.D0
            DO I = 1, NT(IB)
              DYM = MAX(DYM, DY(I))
              IF (IB.EQ.1) THEN
                YT(I,IB) = YT(I,IB) - DY(I)
              ELSE
                YT(I,IB) = YT(I,IB) + DY(I)
              END IF
            END DO
            DYLAST = DYLAST + DYM
            IF (DYM.GT.HCRW) WRITE(7,'(A,I2,A,ES11.3,A,I4,A,I6)')
     &        'UEXTERNALDB: WARNING body', IB, ' surface change', DYM,
     &        ' > HCRIT, step', KSTEP, ' inc', KINC
          END DO
        END IF
        CALL SYSTEM_CLOCK(KT1)
        TEX = TEX + DBLE(KT1-KT0)/DBLE(KTR)
        NEX = NEX + 1
      ELSE IF (LOP.EQ.3) THEN
        WRITE(7,'(A,3(A,F9.2,A,I9))') ' TIMERS end of analysis',
     &    ' | UBULK', TBU, ' s, calls', NBU, ' | UCONT', TCO,
     &    ' s, calls', NCO, ' | wear update', TEX, ' s, calls', NEX
        IF (IDIAG.EQ.1) CLOSE(76)
        CALL WOUT(KSTEP,TIME(2))
        CLOSE(71)
        CLOSE(72)
        IF (LDUMP) CLOSE(74)
      END IF
      RETURN
      END
C
C ======================================================================
      SUBROUTINE WOUT(KS,TT)
C     Output at the end of step KS: contact results of all contact
C     elements (unit 71) and the surface profiles of both bodies on
C     their coarse grids (unit 72).
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MAXCE=5000, NTM=5001)
      DOUBLE PRECISION CO, XT, YT
      INTEGER NCE, NT, ITRD, NSB
      COMMON /CONTO/ CO(6,MAXCE), NCE
      COMMON /SURFT/ XT(NTM,2), YT(NTM,2), NT(2), ITRD, NSB(2)
C
      WRITE(71,'(A,I4,A,ES16.8)') 'STEP', KS, ' TIME', TT
      DO IC = 1, NCE
        WRITE(71,'(I6,5ES16.8,F4.0)') IC, (CO(I,IC), I = 1, 6)
      END DO
      FLUSH(71)
      WRITE(72,'(A,I4,A,ES16.8)') 'STEP', KS, ' TIME', TT
      DO IB = 1, 2
        NS = MAX(NSB(IB),1)
        NG = (NT(IB)-1)/NS + 1
        WRITE(72,'(A,I2,I6)') 'BODY', IB, NG
        SIDE = 1.D0
        IF (IB.EQ.2) SIDE = -1.D0
        DO M = 1, NG
          I = 1 + (M-1)*NS
          WRITE(72,'(2ES20.12)') XT(I,IB), YLS(IB,SIDE,XT(I,IB))
        END DO
      END DO
      FLUSH(72)
      RETURN
      END
C
C ======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &     RPL,DDSDDT,DRPLDE,DRPLDT,
     &     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     &     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     &     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C     Ghost mesh: negligible stiffness, STATEV(1:5) = UBULK results of
C     element NOEL - PROPS(1) (material-averaged S11..S12, fraction).
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     &     DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     &     TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     &     DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
      PARAMETER (MAXEL=200000, EGHOST=1.D-8)
      DOUBLE PRECISION GOUT
      COMMON /GHOST/ GOUT(5,MAXEL)
C
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = 0.D0
        END DO
        DDSDDE(I,I) = EGHOST
        STRESS(I) = 0.D0
      END DO
      IEL = NOEL - NINT(PROPS(1))
      IF (IEL.GE.1 .AND. IEL.LE.MAXEL) THEN
        DO I = 1, MIN(5,NSTATV)
          STATEV(I) = GOUT(I,IEL)
        END DO
      END IF
      RETURN
      END
C
C ======================================================================
C     NBLS: narrow band level set on a fixed Cartesian grid per body.
C     Common block (identical in every routine that uses it):
C       PG(i,j,IB)  phi of body IB at grid node (X0 + (i-1) DX,
C                   Y0 + (j-1) DX), signed distance, material < 0,
C                   clipped to +-NBND*DX outside the band
C       LSON        1 if lsgrid.txt was found (NBLS active)
C ======================================================================
      SUBROUTINE LSINIT(OUTDIR,LENOUT)
C     Reads lsgrid.txt (if present) and initialises phi of both bodies
C     as the signed distance to their surface tables.
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MGX=8001, MGY=241)
      DOUBLE PRECISION PG, GX0, GY0, GDX, HKER
      INTEGER NGX, NGY, NBND, NREI, NUPD, LSON
      COMMON /LSG/ PG(MGX,MGY,2), GX0(2), GY0(2), GDX(2), HKER(2),
     &     NGX(2), NGY(2), NBND(2), NREI(2), NUPD(2), LSON
      CHARACTER*256 OUTDIR
      INTEGER LENOUT
      LOGICAL LEX
C
      LSON = 0
      INQUIRE(FILE=OUTDIR(1:LENOUT)//'\lsgrid.txt',EXIST=LEX)
      IF (.NOT.LEX) RETURN
      OPEN(69,FILE=OUTDIR(1:LENOUT)//'\lsgrid.txt',STATUS='OLD')
      DO IB = 1, 2
        READ(69,*) X0, X1, Y0, Y1, DX, NB, NR, HK
        NGX(IB) = NINT((X1-X0)/DX) + 1
        NGY(IB) = NINT((Y1-Y0)/DX) + 1
        IF (NGX(IB).GT.MGX .OR. NGY(IB).GT.MGY .OR. NB.LT.2) THEN
          WRITE(7,*) 'LSINIT: grid too large or band too narrow'
          CALL XIT
        END IF
        GX0(IB) = X0
        GY0(IB) = Y0
        GDX(IB) = DX
        NBND(IB) = NB
        NREI(IB) = MAX(NR,1)
        HKER(IB) = HK
        NUPD(IB) = 0
        SIDE = 1.D0
        IF (IB.EQ.2) SIDE = -1.D0
        B = NB*DX
        DO J = 1, NGY(IB)
          Y = Y0 + (J-1)*DX
          DO I = 1, NGX(IB)
            X = X0 + (I-1)*DX
            S = SIDE*(Y - YSURF(IB,X))
            IF (ABS(S).GT.4.D0*B) THEN
              PG(I,J,IB) = SIGN(B,S)
            ELSE
              D = MIN(DPLINE(IB,X,Y,2.D0*B), B)
              IF (S.LT.0.D0) D = -D
              PG(I,J,IB) = D
            END IF
          END DO
        END DO
        WRITE(7,'(A,I2,A,2I6,A,ES11.3,A,I3)') 'LSINIT: body', IB,
     &    ' grid', NGX(IB), NGY(IB), ' dx', DX, ' band', NB
      END DO
      CLOSE(69)
      LSON = 1
      RETURN
      END
C
C ======================================================================
      DOUBLE PRECISION FUNCTION DPLINE(IB,X,Y,RMAX)
C     Distance from (X, Y) to the surface table polyline of body IB,
C     searching the segments within RMAX in x (RMAX if none closer).
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (NTM=5001)
      DOUBLE PRECISION XT, YT
      INTEGER NT, ITRD, NSB
      COMMON /SURFT/ XT(NTM,2), YT(NTM,2), NT(2), ITRD, NSB(2)
C
      N = NT(IB)
      DPLINE = RMAX
      XA = X - RMAX
      IL = 1
      IR = N
   10 IF (IR-IL.GT.1) THEN
        IM = (IL+IR)/2
        IF (XT(IM,IB).LE.XA) THEN
          IL = IM
        ELSE
          IR = IM
        END IF
        GOTO 10
      END IF
      DO I = IL, N-1
        IF (XT(I,IB).GT.X+RMAX) GOTO 20
        DX = XT(I+1,IB) - XT(I,IB)
        DY = YT(I+1,IB) - YT(I,IB)
        DL2 = DX*DX + DY*DY
        T = 0.D0
        IF (DL2.GT.0.D0) T = ((X-XT(I,IB))*DX + (Y-YT(I,IB))*DY)/DL2
        T = MIN(1.D0, MAX(0.D0, T))
        D = SQRT((XT(I,IB)+T*DX-X)**2 + (YT(I,IB)+T*DY-Y)**2)
        DPLINE = MIN(DPLINE, D)
      END DO
   20 CONTINUE
      RETURN
      END
C
C ======================================================================
      DOUBLE PRECISION FUNCTION PHIF(IB,SIDE,X,Y)
C     Level set of body IB at (X, Y), material < 0: bilinear on the
C     grid (NBLS) or SIDE*(y - y_s(x)) (tables, or outside the grid).
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MGX=8001, MGY=241)
      DOUBLE PRECISION PG, GX0, GY0, GDX, HKER
      INTEGER NGX, NGY, NBND, NREI, NUPD, LSON
      COMMON /LSG/ PG(MGX,MGY,2), GX0(2), GY0(2), GDX(2), HKER(2),
     &     NGX(2), NGY(2), NBND(2), NREI(2), NUPD(2), LSON
C
      IF (LSON.EQ.1) THEN
        FI = (X-GX0(IB))/GDX(IB)
        FJ = (Y-GY0(IB))/GDX(IB)
        IF (FI.GE.0.D0 .AND. FI.LE.DBLE(NGX(IB)-1) .AND.
     &      FJ.GE.0.D0 .AND. FJ.LE.DBLE(NGY(IB)-1)) THEN
          I = MIN(INT(FI)+1, NGX(IB)-1)
          J = MIN(INT(FJ)+1, NGY(IB)-1)
          TX = FI - DBLE(I-1)
          TY = FJ - DBLE(J-1)
          PHIF = (1.D0-TX)*(1.D0-TY)*PG(I,J,IB)
     &         + TX*(1.D0-TY)*PG(I+1,J,IB)
     &         + TX*TY*PG(I+1,J+1,IB)
     &         + (1.D0-TX)*TY*PG(I,J+1,IB)
          RETURN
        END IF
      END IF
      PHIF = SIDE*(Y - YSURF(IB,X))
      RETURN
      END
C
C ======================================================================
      SUBROUTINE DSURF(IB,SIDE,X,Y,DV)
C     Vector from (X, Y) to the true surface of body IB: vertical
C     offset with the tables, -phi grad(phi)/|grad(phi)|^2 with NBLS.
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MGX=8001, MGY=241)
      DOUBLE PRECISION PG, GX0, GY0, GDX, HKER
      INTEGER NGX, NGY, NBND, NREI, NUPD, LSON
      COMMON /LSG/ PG(MGX,MGY,2), GX0(2), GY0(2), GDX(2), HKER(2),
     &     NGX(2), NGY(2), NBND(2), NREI(2), NUPD(2), LSON
      DIMENSION DV(2)
C
      IF (LSON.EQ.0) THEN
        DV(1) = 0.D0
        DV(2) = YSURF(IB,X) - Y
        RETURN
      END IF
      H = 0.5D0*GDX(IB)
      P0 = PHIF(IB,SIDE,X,Y)
      GX = (PHIF(IB,SIDE,X+H,Y) - PHIF(IB,SIDE,X-H,Y))/(2.D0*H)
      GY = (PHIF(IB,SIDE,X,Y+H) - PHIF(IB,SIDE,X,Y-H))/(2.D0*H)
      G2 = GX*GX + GY*GY
      DV(1) = 0.D0
      DV(2) = 0.D0
      IF (G2.GT.0.D0) THEN
        DV(1) = -P0*GX/G2
        DV(2) = -P0*GY/G2
      END IF
      RETURN
      END
C
C ======================================================================
      DOUBLE PRECISION FUNCTION YLS(IB,SIDE,X)
C     Surface height of body IB at X for the output of the profiles:
C     first zero crossing of phi scanning the grid column from the
C     void side (above the flat, below the cylinder).
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MGX=8001, MGY=241)
      DOUBLE PRECISION PG, GX0, GY0, GDX, HKER
      INTEGER NGX, NGY, NBND, NREI, NUPD, LSON
      COMMON /LSG/ PG(MGX,MGY,2), GX0(2), GY0(2), GDX(2), HKER(2),
     &     NGX(2), NGY(2), NBND(2), NREI(2), NUPD(2), LSON
C
      YLS = YSURF(IB,X)
      IF (LSON.EQ.0) RETURN
      FI = (X-GX0(IB))/GDX(IB)
      IF (FI.LT.0.D0 .OR. FI.GT.DBLE(NGX(IB)-1)) RETURN
      N = NGY(IB)
      DO K = 1, N-1
        IF (SIDE.GT.0.D0) THEN
          J1 = N - K + 1
          J2 = J1 - 1
        ELSE
          J1 = K
          J2 = K + 1
        END IF
        Y1 = GY0(IB) + (J1-1)*GDX(IB)
        Y2 = GY0(IB) + (J2-1)*GDX(IB)
        P1 = PHIF(IB,SIDE,X,Y1)
        P2 = PHIF(IB,SIDE,X,Y2)
        IF (P1.GE.0.D0 .AND. P2.LT.0.D0) THEN
          YLS = Y1 + (Y2-Y1)*P1/(P1-P2)
          RETURN
        END IF
      END DO
      RETURN
      END
C
C ======================================================================
      SUBROUTINE LSUPD(IB,SIDE,DYM)
C     Archard wear update of the level set of body IB (end of a
C     converged increment):
C       1. zero contour by marching squares (segments, grid cells)
C       2. wear depth at the segment ends: kernel average (hat of radius
C          HKER, weights of the points) of dh = k * cycle jump * p*slip
C          of the contact points (true surface positions)
C       3. every node near a segment gets the distance to the closest
C          segment and the depth interpolated there (extension of the
C          speed along the normals)
C       4. phi <- phi + dh_ext, redistanced every NREINIT updates:
C          phi <- sign(phi)*distance + dh_ext; clipped to the band
C     The grid never moves; phi increases (the material recedes).
C     DYM: maximum depth of the update.
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (MGX=8001, MGY=241, MSEG=40000, MQ=30000)
      DOUBLE PRECISION PG, GX0, GY0, GDX, HKER
      INTEGER NGX, NGY, NBND, NREI, NUPD, LSON
      COMMON /LSG/ PG(MGX,MGY,2), GX0(2), GY0(2), GDX(2), HKER(2),
     &     NGX(2), NGY(2), NBND(2), NREI(2), NUPD(2), LSON
      PARAMETER (MAXCE=5000, NPW=6)
      DOUBLE PRECISION CO, WXS, WXM, WPS, WWT, AKW, CJMP, DYLAST, HCRW
      INTEGER NCE, KWEAR, KSLAST
      COMMON /CONTO/ CO(6,MAXCE), NCE
      COMMON /WEARC/ WXS(NPW,MAXCE), WXM(NPW,MAXCE), WPS(NPW,MAXCE),
     &     WWT(NPW,MAXCE), AKW(2), CJMP, DYLAST, HCRW, KWEAR, KSLAST
      DOUBLE PRECISION WYS, WYM
      COMMON /WEARY/ WYS(NPW,MAXCE), WYM(NPW,MAXCE)
      DIMENSION SG(4,MSEG), SH(2,MSEG), QX(MQ), QY(MQ), QH(MQ), QW(MQ)
      DIMENSION DMIN(MGX,MGY), DHE(MGX,MGY)
      DIMENSION V(4), CX(4), CY(4), PX(4), PY(4)
      SAVE SG, SH, QX, QY, QH, QW, DMIN, DHE
C
      DYM = 0.D0
      DX = GDX(IB)
      NX = NGX(IB)
      NY = NGY(IB)
C
C     ---- contact points and their depths
      NQ = 0
      HMAX = 0.D0
      DO IC = 1, NCE
        DO K = 1, NPW
          IF (WWT(K,IC).GT.0.D0 .AND. NQ.LT.MQ) THEN
            NQ = NQ + 1
            IF (IB.EQ.1) THEN
              QX(NQ) = WXS(K,IC)
              QY(NQ) = WYS(K,IC)
            ELSE
              QX(NQ) = WXM(K,IC)
              QY(NQ) = WYM(K,IC)
            END IF
            QH(NQ) = AKW(IB)*CJMP*WPS(K,IC)
            QW(NQ) = WWT(K,IC)
            HMAX = MAX(HMAX, QH(NQ))
          END IF
        END DO
      END DO
      IF (HMAX.LE.0.D0) RETURN
C
C     ---- 1. zero contour (marching squares, material phi < 0)
      NS = 0
      DO J = 1, NY-1
        DO I = 1, NX-1
          V(1) = PG(I,J,IB)
          V(2) = PG(I+1,J,IB)
          V(3) = PG(I+1,J+1,IB)
          V(4) = PG(I,J+1,IB)
          NNEG = 0
          DO K = 1, 4
            IF (V(K).LT.0.D0) NNEG = NNEG + 1
          END DO
          IF (NNEG.GT.0 .AND. NNEG.LT.4) THEN
            CX(1) = GX0(IB) + (I-1)*DX
            CY(1) = GY0(IB) + (J-1)*DX
            CX(2) = CX(1) + DX
            CY(2) = CY(1)
            CX(3) = CX(2)
            CY(3) = CY(1) + DX
            CX(4) = CX(1)
            CY(4) = CY(3)
            NP = 0
            DO K = 1, 4
              K2 = MOD(K,4) + 1
              IF ((V(K).LT.0.D0) .NEQV. (V(K2).LT.0.D0)) THEN
                T = V(K)/(V(K)-V(K2))
                NP = NP + 1
                PX(NP) = CX(K) + T*(CX(K2)-CX(K))
                PY(NP) = CY(K) + T*(CY(K2)-CY(K))
              END IF
            END DO
            IF (NS+2.GT.MSEG) THEN
              WRITE(7,*) 'LSUPD: too many contour segments'
              CALL XIT
            END IF
            IF (NP.EQ.2) THEN
              NS = NS + 1
              SG(1,NS) = PX(1)
              SG(2,NS) = PY(1)
              SG(3,NS) = PX(2)
              SG(4,NS) = PY(2)
            ELSE IF (NP.EQ.4) THEN
C             saddle: pair the crossings by the sign of the centre
              VC = 0.25D0*(V(1)+V(2)+V(3)+V(4))
              IF ((VC.LT.0.D0) .EQV. (V(1).LT.0.D0)) THEN
                NS = NS + 1
                SG(1,NS) = PX(1)
                SG(2,NS) = PY(1)
                SG(3,NS) = PX(4)
                SG(4,NS) = PY(4)
                NS = NS + 1
                SG(1,NS) = PX(2)
                SG(2,NS) = PY(2)
                SG(3,NS) = PX(3)
                SG(4,NS) = PY(3)
              ELSE
                NS = NS + 1
                SG(1,NS) = PX(1)
                SG(2,NS) = PY(1)
                SG(3,NS) = PX(2)
                SG(4,NS) = PY(2)
                NS = NS + 1
                SG(1,NS) = PX(3)
                SG(2,NS) = PY(3)
                SG(3,NS) = PX(4)
                SG(4,NS) = PY(4)
              END IF
            END IF
          END IF
        END DO
      END DO
C
C     ---- 2. depth at the segment ends (kernel average)
      R = HKER(IB)
      DO IS = 1, NS
        DO IE = 1, 2
          XE = SG(2*IE-1,IS)
          YE = SG(2*IE,IS)
          A = 0.D0
          B = 0.D0
          DO IQ = 1, NQ
            IF (ABS(QX(IQ)-XE).LT.R) THEN
              D = SQRT((QX(IQ)-XE)**2 + (QY(IQ)-YE)**2)
              IF (D.LT.R) THEN
                W = (1.D0 - D/R)*QW(IQ)
                A = A + W*QH(IQ)
                B = B + W
              END IF
            END IF
          END DO
          SH(IE,IS) = 0.D0
          IF (B.GT.0.D0) SH(IE,IS) = A/B
          DYM = MAX(DYM, SH(IE,IS))
        END DO
      END DO
C
C     ---- 3. closest segment of the nodes near the contour
      BIG = 1.D30
      DO J = 1, NY
        DO I = 1, NX
          DMIN(I,J) = BIG
        END DO
      END DO
      NW = NBND(IB) + 1
      DO IS = 1, NS
        AX = SG(1,IS)
        AY = SG(2,IS)
        BX = SG(3,IS) - AX
        BY = SG(4,IS) - AY
        BL2 = BX*BX + BY*BY
        I1 = MAX(1, INT((MIN(AX,SG(3,IS))-GX0(IB))/DX) + 1 - NW)
        I2 = MIN(NX, INT((MAX(AX,SG(3,IS))-GX0(IB))/DX) + 2 + NW)
        J1 = MAX(1, INT((MIN(AY,SG(4,IS))-GY0(IB))/DX) + 1 - NW)
        J2 = MIN(NY, INT((MAX(AY,SG(4,IS))-GY0(IB))/DX) + 2 + NW)
        DO J = J1, J2
          Y = GY0(IB) + (J-1)*DX
          DO I = I1, I2
            X = GX0(IB) + (I-1)*DX
            T = 0.D0
            IF (BL2.GT.0.D0) T = ((X-AX)*BX + (Y-AY)*BY)/BL2
            T = MIN(1.D0, MAX(0.D0, T))
            D = SQRT((AX+T*BX-X)**2 + (AY+T*BY-Y)**2)
            IF (D.LT.DMIN(I,J)) THEN
              DMIN(I,J) = D
              DHE(I,J) = SH(1,IS) + T*(SH(2,IS)-SH(1,IS))
            END IF
          END DO
        END DO
      END DO
C
C     ---- 4. update (and redistancing), clipped to the band
      NUPD(IB) = NUPD(IB) + 1
      LRE = 0
      IF (MOD(NUPD(IB),NREI(IB)).EQ.0) LRE = 1
      BND = NBND(IB)*DX
      DO J = 1, NY
        DO I = 1, NX
          IF (DMIN(I,J).LT.BIG) THEN
            P = PG(I,J,IB)
            IF (LRE.EQ.1) THEN
              IF (P.LT.0.D0) THEN
                P = -DMIN(I,J)
              ELSE
                P = DMIN(I,J)
              END IF
            END IF
            P = P + DHE(I,J)
            PG(I,J,IB) = MIN(BND, MAX(-BND, P))
          END IF
        END DO
      END DO
      RETURN
      END
