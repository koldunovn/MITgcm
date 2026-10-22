#include "SEAICE_OPTIONS.h"
#include "CPP_EEOPTIONS.h"
#ifdef ALLOW_EXCH2
#include "W2_OPTIONS.h"
#endif

#ifdef SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ

CBOP
C     !ROUTINE: SEAICE_FP_EXCH_UV_XY_RL_AD

C     !INTERFACE:
      SUBROUTINE SEAICE_FP_EXCH_UV_XY_RL_AD(
     U                          uPhi, vPhi,
     I                          withSigns, exchAdj, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE SEAICE_FP_EXCH_UV_XY_RL_AD
C     | o Adjoint of EXCH_UV_XY_RL (eesupp/src/exch_uv_xy_rx.template)
C     |   for the EVP fixed-point adjoint
C     |   (SEAICE_EVP_SUBCYCLES_FPADJ_AD): uPhi, vPhi hold the
C     |   adjoint of the exchanged field on entry and the adjoint
C     |   of the field before the exchange on exit.
C     *==========================================================*
C     | exchAdj (SEAICE_fpExchAdj) selects the adjoint:
C     |  1 : the one MITgcm's AD tools use. Both TAF and Tapenade
C     |      differentiate exch2_uv_3d_rl.f (exch2_ad_diff.list)
C     |      and call the hand-written EXCH2_RL2_CUBE_AD for the
C     |      edge exchanges (pkg/exch2/exch2.flow, and
C     |      exch2_rl2_cube_b.F for Tapenade). This is that
C     |      adjoint written out: the transposed corner copies of
C     |      EXCH2_UV_3D_RL, then EXCH2_RL2_CUBE_AD with
C     |      EXCH_UPDATE_CORNERS and with EXCH_IGNORE_CORNERS.
C     |  2 : the same, with SEAICE_FP_EXCH2_RL2_CUBE_ADX instead
C     |      of EXCH2_RL2_CUBE_AD (below; the exact transpose of
C     |      the edge exchange also where the target regions of two
C     |      neighbours overlap at the cube corners). Once
C     |      EXCH2_RX2_CUBE_AD runs its GET-adjoint neighbour loop
C     |      backwards as well, options 1 and 2 are the same.
C     | Only the ALLOW_EXCH2 arm without W2_USE_R1_ONLY is
C     | implemented (cs32x15); the other arms stop.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     uPhi      :: adjoint of the 2 components of a vector field
C     vPhi      ::  with overlap regions
C     withSigns :: Flag controlling whether vector is signed.
C     exchAdj   :: 1 or 2 (above)
C     myThid    :: my Thread Id. number
      _RL uPhi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vPhi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      LOGICAL withSigns
      INTEGER exchAdj
      INTEGER myThid
CEOP

      IF ( exchAdj.NE.1 .AND. exchAdj.NE.2 ) THEN
       STOP 'ABNORMAL END: S/R SEAICE_FP_EXCH_UV_XY_RL_AD: exchAdj'
      ENDIF
#ifdef ALLOW_EXCH2
#ifdef W2_USE_R1_ONLY
      STOP 'S/R SEAICE_FP_EXCH_UV_XY_RL_AD: W2_USE_R1_ONLY not covered'
#else
      CALL SEAICE_FP_EXCH2_UV_3D_RL_AD(
     U                     uPhi, vPhi,
     I                     withSigns, 1, exchAdj, myThid )
#endif
#else /* ALLOW_EXCH2 */
      STOP 'S/R SEAICE_FP_EXCH_UV_XY_RL_AD: needs ALLOW_EXCH2'
#endif /* ALLOW_EXCH2 */

      RETURN
      END

#ifdef ALLOW_EXCH2
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_FP_EXCH2_UV_3D_RL_AD

C     !INTERFACE:
      SUBROUTINE SEAICE_FP_EXCH2_UV_3D_RL_AD(
     U                           uPhi, vPhi,
     I                           withSigns, myNz, exchAdj, myThid )

C     !DESCRIPTION:
C     Adjoint of EXCH2_UV_3D_RL (pkg/exch2/exch2_uv_3d_rx.template,
C     with DO_CORNER_COPY_V2U undefined there, :3): the
C     forward is EXCH2_RL2_CUBE with EXCH_IGNORE_CORNERS (:70-73), with
C     EXCH_UPDATE_CORNERS (:74-77), then the corner copies (:122-222);
C     the adjoint runs them backwards.

C     !USES:
      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"

C     !INPUT/OUTPUT PARAMETERS:
C     uPhi, vPhi :: adjoint of the vector field (in: after, out: before
C                   the exchange)
      INTEGER myNz
      _RL uPhi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNz,nSx,nSy)
      _RL vPhi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNz,nSx,nSy)
      LOGICAL withSigns
      INTEGER exchAdj
      INTEGER myThid

C     !LOCAL VARIABLES:
      INTEGER OLw, OLe, OLn, OLs, exchWidthX, exchWidthY
CEOP

      OLw        = OLx
      OLe        = OLx
      OLn        = OLy
      OLs        = OLy
      exchWidthX = OLx
      exchWidthY = OLy

C     adjoint of the corner copies (:79-227)
      CALL SEAICE_FP_EXCH2_UV_CORNERS_AD(
     U                     uPhi, vPhi,
     I                     withSigns, myNz, myThid )

C     adjoint of :74-77, then of :70-73
      IF ( exchAdj.EQ.1 ) THEN
       CALL EXCH2_RL2_CUBE_AD( uPhi, vPhi, withSigns, 'Cg',
     I            OLw, OLe, OLs, OLn, myNz,
     I            exchWidthX, exchWidthY,
     I            EXCH_UPDATE_CORNERS, myThid )
       CALL EXCH2_RL2_CUBE_AD( uPhi, vPhi, withSigns, 'Cg',
     I            OLw, OLe, OLs, OLn, myNz,
     I            exchWidthX, exchWidthY,
     I            EXCH_IGNORE_CORNERS, myThid )
      ELSE
       CALL SEAICE_FP_EXCH2_RL2_CUBE_ADX( uPhi, vPhi, withSigns, 'Cg',
     I            OLw, OLe, OLs, OLn, myNz,
     I            exchWidthX, exchWidthY,
     I            EXCH_UPDATE_CORNERS, myThid )
       CALL SEAICE_FP_EXCH2_RL2_CUBE_ADX( uPhi, vPhi, withSigns, 'Cg',
     I            OLw, OLe, OLs, OLn, myNz,
     I            exchWidthX, exchWidthY,
     I            EXCH_IGNORE_CORNERS, myThid )
      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_FP_EXCH2_UV_CORNERS_AD

C     !INTERFACE:
      SUBROUTINE SEAICE_FP_EXCH2_UV_CORNERS_AD(
     U                           uPhi, vPhi,
     I                           withSigns, myNz, myThid )

C     !DESCRIPTION:
C     Adjoint of the corner copies of EXCH2_UV_3D_RL
C     (exch2_uv_3d_rx.template:79-227; DO_CORNER_COPY_V2U is
C     undefined at :3). Each copy a = s*b becomes b_ad = b_ad + s*a_ad,
C     a_ad = 0, in the reverse order of the forward statements; the
C     four corners are taken in reverse order. W2_FILL_NULL_REGIONS
C     (constant fill values) is not covered.

C     !USES:
      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"

C     !INPUT/OUTPUT PARAMETERS:
      INTEGER myNz
      _RL uPhi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNz,nSx,nSy)
      _RL vPhi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,myNz,nSx,nSy)
      LOGICAL withSigns
      INTEGER myThid

C     !LOCAL VARIABLES:
      INTEGER bi, bj, myTile, k
CEOP

#ifdef W2_FILL_NULL_REGIONS
      STOP 'S/R SEAICE_FP_EXCH2_UV_CORNERS_AD: W2_FILL_NULL_REGIONS'
#endif

      IF (useCubedSphereExchange) THEN

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        myTile = W2_myTileList(bi,bj)

C     adjoint of :201-222 (East and North edges)
        IF ( exch2_isEedge(myTile) .EQ. 1 .AND.
     &       exch2_isNedge(myTile) .EQ. 1 ) THEN
         DO k=myNz,1,-1
          IF ( OLx.GE.2 .AND. OLy.GE.2 ) THEN
C            vPhi(sNx+1,sNy+2,k,bi,bj)=uPhi(sNx+2,sNy,k,bi,bj)
            uPhi(sNx+2,sNy,k,bi,bj) = uPhi(sNx+2,sNy,k,bi,bj)
     &                              + vPhi(sNx+1,sNy+2,k,bi,bj)
            vPhi(sNx+1,sNy+2,k,bi,bj) = 0. _d 0
C            uPhi(sNx+2,sNy+1,k,bi,bj)=vPhi(sNx,sNy+2,k,bi,bj)
            vPhi(sNx,sNy+2,k,bi,bj) = vPhi(sNx,sNy+2,k,bi,bj)
     &                              + uPhi(sNx+2,sNy+1,k,bi,bj)
            uPhi(sNx+2,sNy+1,k,bi,bj) = 0. _d 0
          ENDIF
         ENDDO
        ENDIF

C     adjoint of :173-199 (East and South edges)
        IF ( exch2_isEedge(myTile) .EQ. 1 .AND.
     &       exch2_isSedge(myTile) .EQ. 1 ) THEN
         DO k=myNz,1,-1
          IF ( OLx.GE.2 .AND. OLy.GE.2 ) THEN
           IF ( withSigns ) THEN
C            vPhi(sNx+1,0,k,bi,bj)=-uPhi(sNx+2,1,k,bi,bj)
            uPhi(sNx+2,1,k,bi,bj) = uPhi(sNx+2,1,k,bi,bj)
     &                            - vPhi(sNx+1,0,k,bi,bj)
            vPhi(sNx+1,0,k,bi,bj) = 0. _d 0
C            uPhi(sNx+2,0,k,bi,bj)=-vPhi(sNx,0,k,bi,bj)
            vPhi(sNx,0,k,bi,bj) = vPhi(sNx,0,k,bi,bj)
     &                          - uPhi(sNx+2,0,k,bi,bj)
            uPhi(sNx+2,0,k,bi,bj) = 0. _d 0
           ELSE
C            vPhi(sNx+1,0,k,bi,bj)= uPhi(sNx+2,1,k,bi,bj)
            uPhi(sNx+2,1,k,bi,bj) = uPhi(sNx+2,1,k,bi,bj)
     &                            + vPhi(sNx+1,0,k,bi,bj)
            vPhi(sNx+1,0,k,bi,bj) = 0. _d 0
C            uPhi(sNx+2,0,k,bi,bj)= vPhi(sNx,0,k,bi,bj)
            vPhi(sNx,0,k,bi,bj) = vPhi(sNx,0,k,bi,bj)
     &                          + uPhi(sNx+2,0,k,bi,bj)
            uPhi(sNx+2,0,k,bi,bj) = 0. _d 0
           ENDIF
          ENDIF
         ENDDO
        ENDIF

C     adjoint of :145-171 (West and North edges)
        IF ( exch2_isWedge(myTile) .EQ. 1 .AND.
     &       exch2_isNedge(myTile) .EQ. 1 ) THEN
         DO k=myNz,1,-1
          IF ( OLx.GE.2 .AND. OLy.GE.2 ) THEN
           IF ( withSigns ) THEN
C            vPhi(0,sNy+2,k,bi,bj)=-uPhi(0,sNy,k,bi,bj)
            uPhi(0,sNy,k,bi,bj) = uPhi(0,sNy,k,bi,bj)
     &                          - vPhi(0,sNy+2,k,bi,bj)
            vPhi(0,sNy+2,k,bi,bj) = 0. _d 0
C            uPhi(0,sNy+1,k,bi,bj)=-vPhi(1,sNy+2,k,bi,bj)
            vPhi(1,sNy+2,k,bi,bj) = vPhi(1,sNy+2,k,bi,bj)
     &                            - uPhi(0,sNy+1,k,bi,bj)
            uPhi(0,sNy+1,k,bi,bj) = 0. _d 0
           ELSE
C            vPhi(0,sNy+2,k,bi,bj)= uPhi(0,sNy,k,bi,bj)
            uPhi(0,sNy,k,bi,bj) = uPhi(0,sNy,k,bi,bj)
     &                          + vPhi(0,sNy+2,k,bi,bj)
            vPhi(0,sNy+2,k,bi,bj) = 0. _d 0
C            uPhi(0,sNy+1,k,bi,bj)= vPhi(1,sNy+2,k,bi,bj)
            vPhi(1,sNy+2,k,bi,bj) = vPhi(1,sNy+2,k,bi,bj)
     &                            + uPhi(0,sNy+1,k,bi,bj)
            uPhi(0,sNy+1,k,bi,bj) = 0. _d 0
           ENDIF
          ENDIF
         ENDDO
        ENDIF

C     adjoint of :122-143 (West and South edges)
        IF ( exch2_isWedge(myTile) .EQ. 1 .AND.
     &       exch2_isSedge(myTile) .EQ. 1 ) THEN
         DO k=myNz,1,-1
          IF ( OLx.GE.2 .AND. OLy.GE.2 ) THEN
C            vPhi(0,0,k,bi,bj)=uPhi(0,1,k,bi,bj)
            uPhi(0,1,k,bi,bj) = uPhi(0,1,k,bi,bj)
     &                        + vPhi(0,0,k,bi,bj)
            vPhi(0,0,k,bi,bj) = 0. _d 0
C            uPhi(0,0,k,bi,bj)=vPhi(1,0,k,bi,bj)
            vPhi(1,0,k,bi,bj) = vPhi(1,0,k,bi,bj)
     &                        + uPhi(0,0,k,bi,bj)
            uPhi(0,0,k,bi,bj) = 0. _d 0
          ENDIF
         ENDDO
        ENDIF

       ENDDO
      ENDDO

      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CBOP
C     !ROUTINE: SEAICE_FP_EXCH2_RL2_CUBE_ADX

C     !INTERFACE:
      SUBROUTINE SEAICE_FP_EXCH2_RL2_CUBE_ADX(
     U            array1, array2,
     I            signOption, fieldCode,
     I            myOLw, myOLe, myOLs, myOLn, myNz,
     I            exchWidthX, exchWidthY,
     I            cornerMode, myThid )

C     !DESCRIPTION:
C     Two components vector field AD-Exchange:
C     Tile-edge overlap-region of a 2 component vector field is added to
C     corresponding near-edge interior data point and then zero out.
C     SEAICE_FP_EXCH2_RL2_CUBE_ADX (EVP fixed-point adjoint, option
C     SEAICE_fpExchAdj = 2): EXCH2_RL2_CUBE_AD (pkg/exch2/
C     exch2_rx2_cube_ad.template, RX -> RL) with one change:
C     the first neighbour loop (AD of the GET from buffer) runs
C     N = nN,..,1 instead of N = 1,..,nN. The forward EXCH2_RL2_CUBE
C     gets the neighbours N = 1,..,nN in this order; where the target
C     regions of two neighbours overlap (corner halo points with
C     EXCH_UPDATE_CORNERS), the last neighbour's value stays, so the
C     adjoint of that point belongs to the last one. The loop of the
C     template gives it to the first one and zeroes it.

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#include "W2_EXCH2_BUFFER.h"

C     !INPUT/OUTPUT PARAMETERS:
C     array1      :: 1rst component array with edges to exchange.
C     array2      :: 2nd  component array with edges to exchange.
C     signOption  :: Flag controlling whether vector is signed.
C     fieldCode   :: field code (position on staggered grid)
C     myOLw,myOLe :: West and East overlap region sizes.
C     myOLs,myOLn :: South and North overlap region sizes.
C     exchWidthX  :: Width of data region exchanged in X.
C     exchWidthY  :: Width of data region exchanged in Y.
C     cornerMode  :: halo-corner-region treatment: update/ignore corner region
C     myThid      :: Thread number of this instance of S/R EXCH...

      INTEGER myOLw, myOLe, myOLs, myOLn, myNz
      _RL array1(1-myOLw:sNx+myOLe,
     &           1-myOLs:sNy+myOLn,
     &           myNz, nSx, nSy)
      _RL array2(1-myOLw:sNx+myOLe,
     &           1-myOLs:sNy+myOLn,
     &           myNz, nSx, nSy)
      LOGICAL signOption
      CHARACTER*2 fieldCode
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER cornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     e2_msgHandles :: Synchronization and coordination data structure used to
C                   :: coordinate access to e2Bufr1_RL or to regulate message
C                   :: buffering. In PUT communication sender will increment
C                   :: handle entry once data is ready in buffer. Receiver will
C                   :: decrement handle once data is consumed from buffer.
C                   :: For MPI MSG communication MPI_Wait uses handle to check
C                   :: Isend has cleared. This is done in routine after receives.
C     note: a) current implementation does not use e2_msgHandles for communication
C              between threads: all-threads barriers are used (see CNH note below).
C              For a 2-threads synchro communication (future version),
C              e2_msgHandles should be shared (in common block, moved to BUFFER.h)
      INTEGER bi, bj
C     Variables for working through W2 topology
      INTEGER e2_msgHandles( 2, W2_maxNeighbours, nSx, nSy )
      INTEGER thisTile, farTile, N, nN, oN
      INTEGER tIlo1, tIhi1, tJlo1, tJhi1, oIs1, oJs1
      INTEGER tIlo2, tIhi2, tJlo2, tJhi2, oIs2, oJs2
      INTEGER tIStride, tJStride
      INTEGER tKlo, tKhi, tKStride
      INTEGER i1Lo, i1Hi, j1Lo, j1Hi, k1Lo, k1Hi
      INTEGER i2Lo, i2Hi, j2Lo, j2Hi, k2Lo, k2Hi
      LOGICAL updateCorners

#ifdef ALLOW_USE_MPI
      INTEGER iBufr1, iBufr2, nri, nrj
C     MPI stuff (should be in a routine call)
      INTEGER mpiStatus(MPI_STATUS_SIZE)
      INTEGER mpiRc
      INTEGER wHandle
#endif
CEOP

      updateCorners = cornerMode .EQ. EXCH_UPDATE_CORNERS
C-    Tile size of arrays to exchange:
      i1Lo  = 1-myOLw
      i1Hi  = sNx+myOLe
      j1Lo  = 1-myOLs
      j1Hi  = sNy+myOLn
      k1Lo  = 1
      k1Hi  = myNz
      i2Lo  = 1-myOLw
      i2Hi  = sNx+myOLe
      j2Lo  = 1-myOLs
      j2Hi  = sNy+myOLn
      k2Lo  = 1
      k2Hi  = myNz

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     Prevent anyone to access shared buffer while an other thread modifies it
      CALL BAR2( myThid )

C--   Extract from buffer (either from level 1 if local exch,
C                     or level 2 if coming from an other Proc)
C  AD: = fill buffer from my-tile-edge overlap-region, level 1 or 2 depending
C  AD:   on local (to this Proc) or remote Proc tile destination
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        thisTile=W2_myTileList(bi,bj)
        nN=exch2_nNeighbours(thisTile)
        DO N=nN,1,-1
          CALL EXCH2_GET_UV_BOUNDS(
     I               fieldCode, exchWidthX, updateCorners,
     I               thisTile, N,
     O               tIlo1, tIhi1, tJlo1, tJhi1,
     O               tIlo2, tIhi2, tJlo2, tJhi2,
     O               tiStride, tjStride,
     O               oIs1, oJs1, oIs2, oJs2,
     I               myThid )
          tKLo=1
          tKHi=myNz
          tKStride=1

C     From buffer, get my points
C     (tIlo:tIhi:tiStride,tJlo:tJhi,tJStride,tKlo:tKhi,tKStride) in "array1,2":
C     Note: when transferring data within a process:
C      o e2Bufr entry to read is entry associated with opposing send record
C      o e2_msgHandle entry to read is entry associated with opposing send record.
          CALL EXCH2_AD_GET_RL2(
     I               tIlo1, tIhi1, tIlo2, tIhi2, tiStride,
     I               tJlo1, tJhi1, tJlo2, tJhi2, tjStride,
     I               tKlo, tKhi, tkStride,
     I               thisTile, N, bi, bj,
     I               e2BufrRecSize, W2_maxNeighbours, nSx, nSy,
     O               iBuf1Filled(N,bi,bj), iBuf2Filled(N,bi,bj),
     O               e2Bufr1_RL, e2Bufr2_RL,
     U               array1(1-myOLw,1-myOLs,1,bi,bj),
     U               array2(1-myOLw,1-myOLs,1,bi,bj),
     I               i1Lo, i1Hi, j1Lo, j1Hi, k1Lo, k1Hi,
     I               i2Lo, i2Hi, j2Lo, j2Hi, k2Lo, k2Hi,
     U               e2_msgHandles,
     I               W2_myCommFlag(N,bi,bj), myThid )
        ENDDO
       ENDDO
      ENDDO

C     Wait until all threads finish filling buffer
      CALL BAR2( myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_USE_MPI
      IF ( usingMPI ) THEN
C  AD: all MPI part is acting on buffer and is identical to forward code,
C  AD: except a) the buffer level: send from lev.2, receive into lev.1
C  AD:        b) the length of transferred buffer (<- match the ad_put/ad_get)

      _BEGIN_MASTER( myThid )

C--   Send my data (in buffer, level 1) to target Process
      DO bj=1,nSy
       DO bi=1,nSx
        thisTile=W2_myTileList(bi,bj)
        nN=exch2_nNeighbours(thisTile)
        DO N=1,nN
C-    Skip the call if this is an internal exchange
         IF ( W2_myCommFlag(N,bi,bj) .EQ. 'M' ) THEN
          CALL EXCH2_SEND_RL2(
     I               thisTile, N,
     I               e2BufrRecSize,
     I               iBuf1Filled(N,bi,bj),    iBuf2Filled(N,bi,bj),
     I               e2Bufr1_RL(1,N,bi,bj,2), e2Bufr2_RL(1,N,bi,bj,2),
     O               e2_msgHandles(1,N,bi,bj),
     I               W2_myCommFlag(N,bi,bj), myThid )
         ENDIF
        ENDDO
       ENDDO
      ENDDO

C--   Receive data (in buffer, level 1) from source Process
      DO bj=1,nSy
       DO bi=1,nSx
        thisTile=W2_myTileList(bi,bj)
        nN=exch2_nNeighbours(thisTile)
        DO N=1,nN
C-    Skip the call if this is an internal exchange
         IF ( W2_myCommFlag(N,bi,bj) .EQ. 'M' ) THEN
          farTile=exch2_neighbourId(N,thisTile)
          oN = exch2_opposingSend(N,thisTile)
          CALL EXCH2_GET_UV_BOUNDS(
     I               fieldCode, exchWidthX, updateCorners,
     I               farTile, oN,
     O               tIlo1, tIhi1, tJlo1, tJhi1,
     O               tIlo2, tIhi2, tJlo2, tJhi2,
     O               tiStride, tjStride,
     O               oIs1, oJs1, oIs2, oJs2,
     I               myThid )
          nri = 1 + (tIhi1-tIlo1)/tiStride
          nrj = 1 + (tJhi1-tJlo1)/tjStride
          iBufr1 = nri*nrj*myNz
          nri = 1 + (tIhi2-tIlo2)/tiStride
          nrj = 1 + (tJhi2-tJlo2)/tjStride
          iBufr2 = nri*nrj*myNz
C       Receive from neighbour N to fill buffer and later on the array
          CALL EXCH2_RECV_RL2(
     I               thisTile, N,
     I               e2BufrRecSize,
     I               iBufr1, iBufr2,
     I               e2Bufr1_RL(1,N,bi,bj,1), e2Bufr2_RL(1,N,bi,bj,1),
     I               W2_myCommFlag(N,bi,bj), myThid )
         ENDIF
        ENDDO
       ENDDO
      ENDDO

C--   Clear message handles/locks
      DO bj=1,nSy
       DO bi=1,nSx
        thisTile=W2_myTileList(bi,bj)
        nN=exch2_nNeighbours(thisTile)
        DO N=1,nN
C     Note: In a between process tile-tile data transport using
C           MPI the sender needs to clear an Isend wait handle here.
C           In a within process tile-tile data transport using true
C           shared address space/or direct transfer through commonly
C           addressable memory blocks the receiver needs to assert
C           that he has consumed the buffer the sender filled here.
         farTile=exch2_neighbourId(N,thisTile)
         IF     ( W2_myCommFlag(N,bi,bj) .EQ. 'M' ) THEN
          wHandle = e2_msgHandles(1,N,bi,bj)
          CALL MPI_Wait( wHandle, mpiStatus, mpiRc )
          wHandle = e2_msgHandles(2,N,bi,bj)
          CALL MPI_Wait( wHandle, mpiStatus, mpiRc )
         ELSEIF ( W2_myCommFlag(N,bi,bj) .EQ. 'P' ) THEN
         ELSE
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      _END_MASTER( myThid )
C     Everyone waits until master-thread finishes receiving
      CALL BAR2( myThid )

      ENDIF
#endif /* ALLOW_USE_MPI */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Post sends into buffer (buffer level 1):
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        thisTile=W2_myTileList(bi,bj)
        nN=exch2_nNeighbours(thisTile)
        DO N=1,nN
          farTile=exch2_neighbourId(N,thisTile)
          oN = exch2_opposingSend(N,thisTile)
          CALL EXCH2_GET_UV_BOUNDS(
     I               fieldCode, exchWidthX, updateCorners,
     I               farTile, oN,
     O               tIlo1, tIhi1, tJlo1, tJhi1,
     O               tIlo2, tIhi2, tJlo2, tJhi2,
     O               tiStride, tjStride,
     O               oIs1, oJs1, oIs2, oJs2,
     I               myThid )
          tKLo=1
          tKHi=myNz
          tKStride=1
C-    Put my points in buffer for neighbour N to fill points
C     (tIlo:tIhi:tiStride,tJlo:tJhi,tJStride,tKlo:tKhi,tKStride)
C     in its copy of "array1" & "array2".
          CALL EXCH2_AD_PUT_RL2(
     I               tIlo1, tIhi1, tIlo2, tIhi2, tiStride,
     I               tJlo1, tJhi1, tJlo2, tJhi2, tjStride,
     I               tKlo, tKhi, tkStride,
     I               oIs1, oJs1, oIs2, oJs2,
     I               thisTile, N,
     I               e2BufrRecSize,
     O               e2Bufr1_RL(1,N,bi,bj,1), e2Bufr2_RL(1,N,bi,bj,1),
     I               array1(1-myOLw,1-myOLs,1,bi,bj),
     I               array2(1-myOLw,1-myOLs,1,bi,bj),
     I               i1Lo, i1Hi, j1Lo, j1Hi, k1Lo, k1Hi,
     I               i2Lo, i2Hi, j2Lo, j2Hi, k2Lo, k2Hi,
     O               e2_msgHandles(1,N,bi,bj),
     I               W2_myCommFlag(N,bi,bj), signOption, myThid )
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

#endif /* ALLOW_EXCH2 */

#endif /* SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ */
