C ======================================================================
C     LevelSetWear, reference mesh-motion model for Abaqus/Standard
C     (UMESHMOTION, the comparison model of the article)
C
C     Written by Iñigo Llavori
C
C     I. Llavori, Wear simulation on a fixed mesh with level sets
C     and immersed frictional contact,
C     submitted for publication. Mondragon Unibertsitatea.
C     https://github.com/illavori/LevelSetWear       BSD 3-Clause
C
C     Reference mesh-motion model (UMESHMOTION) used for the comparison
C     of the article. Derived from the subroutine published by
C     J. J. Madge in Appendix A1 of his PhD thesis (University of
C     Nottingham, 2008). This 2D version was written by the author
C     during the development of the 3D model of I. Llavori et al.,
C     International Journal of Fatigue 143 (2021) 106012, which did
C     not publish it.
C     Local Archard wear on both contact surfaces (2D). Run with cpus=1
C     (the nodes share data through COMMON blocks). Revised in 2026:
C       - Slip bookkeeping keyed to each slave node and robust to
C         restarts and to abandoned increments, parameters read once
C         per job, double precision throughout.
C       - *ADAPTIVE MESH CONSTRAINT, TYPE=VELOCITY: ULOCAL is a mesh
C         velocity, so the wear depth is divided by DTIME.
C       - Diagnostic files: slaves -> ContactOut_slave.txt,
C         masters -> ContactOut_mast.txt (they were swapped).
C       - Bounds checks on MAXNOD / MAXCNT (stop with XIT).
C       - HCRIT read from file 'hCrit' if present (default 1e-3 mm).
C       - hCritic.txt holds only the maximum accepted wear depth.
C       - JRCD checked for slave nodes (no wear if data missing).
C       - KSTEP=1 guard removed (never reached in restart jobs).
C
C     Block-on-ring option: if a file vSlip exists in the job folder, it
C     holds the imposed sliding per unit step time (mm). The slip
C     increment of every contact node is then vSlip*DTIME instead of
C     the CSLIP increment of Abaqus: continuous unidirectional sliding,
C     where only the normal contact is solved (Cruzado et al., Rev.
C     Metal. 46, 2010). Without the file the behaviour is that of the
C     fretting version.
C ======================================================================
       SUBROUTINE UMESHMOTION(UREF,ULOCAL,NODE,NNDOF,
     &        LNODETYPE,ALOCAL,NDIM,TIME,DTIME,PNEWDT,
     &         KSTEP,KINC,KMESHSWEEP,JMATYP,JGVBLOCK)
C
       INCLUDE 'ABA_PARAM.INC'
C
       PARAMETER (MAXNOD=100000, MAXCNT=2000, NELEMMAX=500)
       PARAMETER (MAXWRN=20)
       DIMENSION ULOCAL(NDIM),JELEMLIST(NELEMMAX)
       DIMENSION ALOCAL(NDIM,*),TIME(2)
       DIMENSION JMATYP(*),JGVBLOCK(*)
       DIMENSION ARRAY(15)
C
       INTEGER MASTER, LFLAG, RFLAG, ISL, CNT1, IOS
       INTEGER JTYP, JRCD, JRCD1, JRCD2, NELEMS, LTRN, JELEMTYPE
       CHARACTER*256 OUTDIR
       INTEGER LENOUTDIR
       LOGICAL LEXIST
C
C     ---- Persistent data (one job = one process)
C     Real data first, integer data in a separate block (alignment).
C       SPRESS  : last contact pressure of each slave node
C       SXCRD   : last x coordinate of each slave node
C       SLIPST  : CSLIP at the end of the last CONVERGED increment
C       SLIPCD  : CSLIP of the current (candidate) increment
C       SINC    : slip increment of each slave in its latest increment
C       SINCC   : slip increment of the last CONVERGED increment
C       SPRC    : contact pressure of the last CONVERGED increment
C       CYCJMP  : cycle jump (read once)
C       HCRIT   : max wear depth per node and increment [mm]
C       WMAX    : max wear depth accepted so far (<= HCRIT)
C       VSLIP   : imposed sliding per unit time (0: computed slip)
       DOUBLE PRECISION SPRESS, SXCRD, SLIPST, SLIPCD, SINC, SINCC
       DOUBLE PRECISION SPRC, CYCJMP, HCRIT, WMAX, VSLIP
       COMMON /WEARR/ SPRESS(MAXCNT), SXCRD(MAXCNT), SLIPST(MAXCNT),
     &   SLIPCD(MAXCNT), SINC(MAXCNT), SINCC(MAXCNT), SPRC(MAXCNT),
     &   CYCJMP, HCRIT, WMAX, VSLIP
C       ISCLK/IMCLK : number of slave/master nodes registered
C       ISREG/IMREG : global node number -> local slave/master index
C       KINCL       : increment in which each slave was last visited
C       INITD       : slave already initialised in this job
C       IPARRD      : parameters already read in this job
C       NCYLND      : number of cylinder nodes (master criterion)
C       NINCR       : increments per fretting cycle
C       NWARN       : JRCD warnings written so far
       INTEGER ISCLK, IMCLK, ISREG, IMREG, KINCL, INITD
       INTEGER IPARRD, NCYLND, NINCR, NWARN
       COMMON /WEARI/ ISCLK, IMCLK, ISREG(MAXNOD), IMREG(MAXNOD),
     &   KINCL(MAXCNT), INITD(MAXCNT), IPARRD, NCYLND, NINCR, NWARN
C
C     ---- Read parameters and open output files once per job
       IF (IPARRD.EQ.0) THEN
         CALL GETOUTDIR(OUTDIR,LENOUTDIR)
         OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\cycleJump')
         READ(68,*) CYCJMP
         CLOSE(68)
         OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\nIncr')
         READ(68,*) NINCR
         CLOSE(68)
         OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\cylinderNodes.txt')
         READ(68,*) NCYLND
         CLOSE(68)
C       Optional file with the critical wear depth
C       (INQUIRE first: inside Abaqus a failed OPEN aborts the run
C        even with IOSTAT)
         HCRIT = 1.0D-3
         INQUIRE(FILE=OUTDIR(1:LENOUTDIR)//'\hCrit',EXIST=LEXIST)
         IF (LEXIST) THEN
           OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\hCrit',STATUS='OLD')
           READ(68,*,IOSTAT=IOS) HCRIT
           IF (IOS.NE.0 .OR. HCRIT.LE.0.0D0) HCRIT = 1.0D-3
           CLOSE(68)
         END IF
         WRITE(7,*) 'UMESHMOTION: HCRIT = ', HCRIT
         VSLIP = 0.0D0
         INQUIRE(FILE=OUTDIR(1:LENOUTDIR)//'\vSlip',EXIST=LEXIST)
         IF (LEXIST) THEN
           OPEN(68,FILE=OUTDIR(1:LENOUTDIR)//'\vSlip',STATUS='OLD')
           READ(68,*,IOSTAT=IOS) VSLIP
           IF (IOS.NE.0) VSLIP = 0.0D0
           CLOSE(68)
         END IF
         WRITE(7,*) 'UMESHMOTION: imposed sliding per unit time = ',
     &        VSLIP
         OPEN(16,FILE=OUTDIR(1:LENOUTDIR)//'\ContactOut_slave.txt',
     &        STATUS='UNKNOWN')
         OPEN(17,FILE=OUTDIR(1:LENOUTDIR)//'\DummyNodes.txt',
     &        STATUS='UNKNOWN')
         OPEN(18,FILE=OUTDIR(1:LENOUTDIR)//'\ContactOut_mast.txt',
     &        STATUS='UNKNOWN')
         OPEN(33,FILE=OUTDIR(1:LENOUTDIR)//'\hCritic.txt',
     &        STATUS='UNKNOWN')
         WMAX = 0.0D0
         WRITE(33,*) WMAX
         IPARRD = 1
       END IF
C
       IF (NODE.GT.MAXNOD) THEN
         WRITE(7,*) 'UMESHMOTION: NODE ', NODE, ' > MAXNOD = ', MAXNOD
         CALL XIT
       END IF
C
C     ---- Nodal contact data from Abaqus
C     JTYP was left undefined in the original version; 0 reproduces the
C     behaviour it had in practice.
       JTYP   = 0
       JRCD   = 0
       NELEMS = NELEMMAX
       CALL GETNODETOELEMCONN(NODE,NELEMS,JELEMLIST,JELEMTYPE,JRCD,
     &      JGVBLOCK)
       CALL GETVRMAVGATNODE(NODE,JTYP,'CSTRESS',ARRAY,
     &      JRCD1,JELEMLIST,NELEMS,JMATYP,JGVBLOCK)
       CPRESS = ARRAY(1)
       CSHEAR = ARRAY(2)
       CALL GETVRMAVGATNODE(NODE,JTYP,'CDISP',ARRAY,
     &      JRCD2,JELEMLIST,NELEMS,JMATYP,JGVBLOCK)
       COPEN  = ARRAY(1)
       CSLIP  = ARRAY(2)
       CALL GETVRN(NODE,'COORD',ARRAY,JRCD,JGVBLOCK,LTRN)
       XCOORD = ARRAY(1)
       YCOORD = ARRAY(2)
C
C     ---- Master (cylinder) or slave (flat)
C     NOTE: relies on the cylinder instance being numbered first.
       MASTER = 0
       IF (NODE.LT.(NCYLND+2)) MASTER = 1
C
       IF (MASTER.EQ.0) THEN
C       ---- Slave node: own contact data
C       Missing contact data: no wear for this node and increment,
C       slip bookkeeping untouched.
         IF (JRCD1.NE.0 .OR. JRCD2.NE.0) THEN
           IF (NWARN.LT.MAXWRN) THEN
             NWARN = NWARN+1
             WRITE(7,*) 'UMESHMOTION: no CSTRESS/CDISP at node ',
     &            NODE, ' KINC ', KINC, ' JRCD ', JRCD1, JRCD2
           END IF
           CPRESS = 0.0D0
           DSLIP  = 0.0D0
           GOTO 100
         END IF
         IF (ISREG(NODE).EQ.0) THEN
           IF (ISCLK.GE.MAXCNT) THEN
             WRITE(7,*) 'UMESHMOTION: more than MAXCNT = ', MAXCNT,
     &            ' slave nodes. Increase MAXCNT.'
             CALL XIT
           END IF
           ISCLK       = ISCLK+1
           ISREG(NODE) = ISCLK
         END IF
         ISL = ISREG(NODE)
         SXCRD(ISL)  = XCOORD
C       Slip reference. First visit in this job: the reference is the
C       current slip (no wear in that first increment). When KINC has
C       advanced, the previous candidate belonged to a converged
C       increment and becomes the reference. A retried (cut back)
C       increment keeps the same KINC, so the reference is untouched.
         IF (INITD(ISL).EQ.0) THEN
           SLIPST(ISL) = CSLIP
           KINCL(ISL)  = KINC
           INITD(ISL)  = 1
         ELSE IF (KINCL(ISL).NE.KINC) THEN
C         Commit the previous (converged) increment for the masters
           SINCC(ISL)  = SLIPCD(ISL) - SLIPST(ISL)
           SPRC(ISL)   = SPRESS(ISL)
           SLIPST(ISL) = SLIPCD(ISL)
           KINCL(ISL)  = KINC
         END IF
         SLIPCD(ISL) = CSLIP
         SPRESS(ISL) = CPRESS
         SINC(ISL)   = CSLIP - SLIPST(ISL)
         DSLIP       = SINC(ISL)
       ELSE
C       ---- Master node: interpolate from the two nearest slaves
C       Uses the slaves' last CONVERGED increment (one-increment lag,
C       independent of node order and never from abandoned attempts).
         IF (IMREG(NODE).EQ.0) THEN
           IMCLK       = IMCLK+1
           IMREG(NODE) = IMCLK
         END IF
         XERRR = 1.0D0
         XERRL = 1.0D0
         LFLAG = 0
         RFLAG = 0
         PR = 0.0D0
         PL = 0.0D0
         SR = 0.0D0
         SL = 0.0D0
         XR = XCOORD
         XL = XCOORD
         DO CNT1 = 1, ISCLK
           DX = ABS(XCOORD-SXCRD(CNT1))
           IF (SXCRD(CNT1).GT.XCOORD) THEN
             IF (DX.LT.XERRR) THEN
               PR = SPRC(CNT1)
               SR = SINCC(CNT1)
               XR = SXCRD(CNT1)
               XERRR = DX
               RFLAG = 1
             END IF
           ELSE
             IF (DX.LT.XERRL) THEN
               PL = SPRC(CNT1)
               SL = SINCC(CNT1)
               XL = SXCRD(CNT1)
               XERRL = DX
               LFLAG = 1
             END IF
           END IF
         END DO
         IF (LFLAG.EQ.1 .AND. RFLAG.EQ.1 .AND. XR.NE.XL) THEN
           CPRESS = PL + (PR-PL)*(XCOORD-XL)/(XR-XL)
           DSLIP  = SL + (SR-SL)*(XCOORD-XL)/(XR-XL)
         ELSE IF (LFLAG.EQ.1) THEN
           CPRESS = PL
           DSLIP  = SL
         ELSE IF (RFLAG.EQ.1) THEN
           CPRESS = PR
           DSLIP  = SR
         ELSE
           CPRESS = 0.0D0
           DSLIP  = 0.0D0
         END IF
       END IF
C
C     ---- Archard wear depth for this node and increment
C     UREF carries the wear coefficient from *ADAPTIVE MESH CONSTRAINT.
 100   CONTINUE
       IF (VSLIP.GT.0.0D0) DSLIP = VSLIP*DTIME
       WDIST = 0.0D0
       IF (CPRESS.GT.0.0D0) THEN
         WDIST = CYCJMP*CPRESS*ABS(DSLIP)*UREF
       END IF
C
C     Move the node inwards along the local surface normal.
C     TYPE=VELOCITY: the mesh moves ULOCAL*DTIME in this increment.
       IF (DTIME.GT.0.0D0) THEN
         ULOCAL(2) = ULOCAL(2) - WDIST/DTIME
       END IF
C
C     ---- Cut back the increment if the wear step is too large
       IF (WDIST.GT.HCRIT) THEN
         PNEWDT0 = PNEWDT
         PNEWDT1 = 0.99D0*HCRIT/WDIST
         IF (PNEWDT1.LT.PNEWDT0) THEN
           PNEWDT = PNEWDT1
           WRITE(7,*) 'CHANGING TIME INCREMENT FROM ', PNEWDT0
           WRITE(7,*) 'TO ', PNEWDT
           WRITE(7,*) 'BASED ON NODE ', NODE
         END IF
       ELSE IF (WDIST.GT.WMAX) THEN
C       Largest accepted wear depth so far: overwrite hCritic.txt
         WMAX = WDIST
         REWIND(33)
         WRITE(33,*) WMAX
       END IF
C
C     ---- Diagnostics
       IF (MASTER.EQ.0) THEN
         WRITE(16,1000) NODE,CPRESS,CSHEAR,MASTER,CSLIP,DSLIP,WDIST,
     &        XCOORD,YCOORD
       ELSE
         WRITE(18,1000) NODE,CPRESS,CSHEAR,MASTER,CSLIP,DSLIP,WDIST,
     &        XCOORD,YCOORD
       END IF
       IF (KINC.EQ.1) THEN
         WRITE(17,*) '0 ','0',XCOORD,YCOORD
       END IF
C
 1000  FORMAT(I6,1X,E13.6,1X,E13.6,1X,I3,1X,E13.6,
     &        1X,E13.6,1X,E13.6,1X,E13.6,1X,E13.6)
C
       RETURN
       END
