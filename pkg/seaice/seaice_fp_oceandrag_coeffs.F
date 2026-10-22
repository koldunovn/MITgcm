#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

C     SEAICE_FP_OCEANDRAG_COEFFS and its hand-written adjoint
C     SEAICE_FP_OCEANDRAG_COEFFS_AD, for the modified subcycle F_mod of
C     the EVP fixed-point adjoint (seaice_evp_fp_subcycle.F,
C     seaice_evp_fp_subcycle_ad.F). SEAICE_FP_OCEANDRAG_COEFFS is
C     seaice_oceandrag_coeffs.F with uVel, vVel as arguments instead of
C     DYNVARS.h (no field with a derivative passes through COMMON in the
C     modified subcycle): the same statements, and the same file head,
C     so
C     that the same CPP arm is compiled (without ALLOW_OBCS this head
C     does not define OBCS_UVICE_OLD: the arm with maskInW/maskInS).

#ifdef SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ

CBOP
C     !ROUTINE: SEAICE_FP_OCEANDRAG_COEFFS
C     !INTERFACE:
      SUBROUTINE SEAICE_FP_OCEANDRAG_COEFFS(
     I     uIceLoc, vIceLoc, HEFFMLoc, uVel, vVel,
     O     CwatC,
     I     iStep, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_OCEANDRAG_COEFFS
C     | o Compute the drag coefficients for ice-ocean drag,
C     |   so that we can use the same code for different solvers
C     *==========================================================*
C     | written by Martin Losch, Oct 2012
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: my Thread Id. number
C     iStep  :: current sub-time step iterate
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER iStep
C     u/vIceLoc :: local copies of the current ice velocity
      _RL uIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     HEFFMLoc  :: local copy of land-sea masks
      _RL HEFFMLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     uVel, vVel :: ocean velocity (DYNVARS.h fields, as arguments)
      _RL uVel      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVel      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C     CwatC     :: drag coefficients
      _RL CwatC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef SEAICE_CGRID
C     === local variables ===
C     i,j,bi,bj,ksrf :: loop indices
      INTEGER i,j,bi,bj
      INTEGER kSrf
      _RL     tempVar, tempMin
      _RL     dragCoeff
CEOP

      IF ( usingPCoords ) THEN
       kSrf = Nr
      ELSE
       kSrf = 1
      ENDIF
      tempMin = SEAICEdWatMin*SEAICEdWatMin

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
C     non-linear water drag coefficients CwatC (DWATN)
#ifdef OBCS_UVICE_OLD
          tempVar = 0.25 _d 0*(
     &            ( ( uIceLoc(i  ,j,bi,bj)-uVel(i  ,j,kSrf,bi,bj) )
     &             +( uIceLoc(i+1,j,bi,bj)-uVel(i+1,j,kSrf,bi,bj) )
     &            )**2
     &          + ( ( vIceLoc(i, j ,bi,bj)-vVel(i, j ,kSrf,bi,bj) )
     &             +( vIceLoc(i,j+1,bi,bj)-vVel(i,j+1,kSrf,bi,bj) )
     &            )**2 )
#else /* OBCS_UVICE_OLD */
          tempVar = 0.25 _d 0*(
     &            ( ( uIceLoc(i  ,j,bi,bj)-uVel(i  ,j,kSrf,bi,bj) )
     &             *maskInW( i ,j,bi,bj)
     &             +( uIceLoc(i+1,j,bi,bj)-uVel(i+1,j,kSrf,bi,bj) )
     &             *maskInW(i+1,j,bi,bj) )**2
     &          + ( ( vIceLoc(i,j  ,bi,bj)-vVel(i,j  ,kSrf,bi,bj) )
     &             *maskInS(i, j ,bi,bj)
     &             +( vIceLoc(i,j+1,bi,bj)-vVel(i,j+1,kSrf,bi,bj) )
     &             *maskInS(i,j+1,bi,bj) )**2 )
#endif /* OBCS_UVICE_OLD */
          IF ( YC(i,j,bi,bj) .LT. ZERO ) THEN
           dragCoeff = SEAICE_waterDrag_south*rhoConst
          ELSE
           dragCoeff = SEAICE_waterDrag      *rhoConst
          ENDIF
          CwatC(i,j,bi,bj) = SEAICEdWatMin
C     This is the AD-proof version of
C         IF ( dragCoeff * SQRT(tempVar) .GT. SEAICEdWatMin )
C     which was previously implemented as
C         IF ( tempVar .GT. (SEAICEdWatMin/dragCoeff)**2 )
C     I like this better:
          IF ( dragCoeff*dragCoeff * tempVar .GT. tempMin )
     &         CwatC(i,j,bi,bj) = dragCoeff*SQRT(tempVar)
          CwatC(i,j,bi,bj) = CwatC(i,j,bi,bj) * HEFFMLoc(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* SEAICE_CGRID */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_FP_OCEANDRAG_COEFFS_AD
C     !INTERFACE:
      SUBROUTINE SEAICE_FP_OCEANDRAG_COEFFS_AD(
     I     uIceLoc, vIceLoc, HEFFMLoc, uVel, vVel,
     U     adCwatC, aduIceLoc, advIceLoc, aduVel, advVel,
     I     myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_FP_OCEANDRAG_COEFFS_AD
C     | o adjoint of SEAICE_FP_OCEANDRAG_COEFFS (above), written
C     |   after Tapenade's reverse mode of it: adCwatC is
C     |   used and set to zero on the points the forward writes,
C     |   the adjoints of the velocities are incremented; the
C     |   forward values are recomputed (no AD runtime)
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL uIceLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL HEFFMLoc (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uVel     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVel     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL adCwatC  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL aduIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL advIceLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL aduVel   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL advVel   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      INTEGER myThid
CEOP

#ifdef SEAICE_CGRID
C     === local variables ===
C     mW1, mW2, mS1, mS2 :: maskInW/S of the two velocity points (1
C                           with OBCS_UVICE_OLD)
      INTEGER i,j,bi,bj
      INTEGER kSrf
      _RL tempVar, tempMin, dragCoeff
      _RL adtempVar, tb0, tb1
      _RL mW1, mW2, mS1, mS2

      IF ( usingPCoords ) THEN
       kSrf = Nr
      ELSE
       kSrf = 1
      ENDIF
      tempMin = SEAICEdWatMin*SEAICEdWatMin

      DO bj=myByHi(myThid),myByLo(myThid),-1
       DO bi=myBxHi(myThid),myBxLo(myThid),-1
        DO j=sNy+OLy-1,1-OLy,-1
         DO i=sNx+OLx-1,1-OLx,-1
#ifdef OBCS_UVICE_OLD
          mW1 = 1. _d 0
          mW2 = 1. _d 0
          mS1 = 1. _d 0
          mS2 = 1. _d 0
#else /* OBCS_UVICE_OLD */
          mW1 = maskInW( i ,j,bi,bj)
          mW2 = maskInW(i+1,j,bi,bj)
          mS1 = maskInS(i,j  ,bi,bj)
          mS2 = maskInS(i,j+1,bi,bj)
#endif /* OBCS_UVICE_OLD */
C     recompute the forward values (seaice_oceandrag_coeffs.F:74-98)
          tempVar = 0.25 _d 0*(
     &            ( ( uIceLoc(i  ,j,bi,bj)-uVel(i  ,j,kSrf,bi,bj) )*mW1
     &             +( uIceLoc(i+1,j,bi,bj)-uVel(i+1,j,kSrf,bi,bj) )*mW2
     &            )**2
     &          + ( ( vIceLoc(i,j  ,bi,bj)-vVel(i,j  ,kSrf,bi,bj) )*mS1
     &             +( vIceLoc(i,j+1,bi,bj)-vVel(i,j+1,kSrf,bi,bj) )*mS2
     &            )**2 )
          IF ( YC(i,j,bi,bj) .LT. ZERO ) THEN
           dragCoeff = SEAICE_waterDrag_south*rhoConst
          ELSE
           dragCoeff = SEAICE_waterDrag      *rhoConst
          ENDIF
C     reverse of CwatC = CwatC*HEFFMLoc (:107) and of
C     CwatC = dragCoeff*SQRT(tempVar) where it exceeds SEAICEdWatMin
C     (:99-106)
          adCwatC(i,j,bi,bj) = HEFFMLoc(i,j,bi,bj)*adCwatC(i,j,bi,bj)
          adtempVar = 0. _d 0
          IF ( dragCoeff*dragCoeff * tempVar .GT. tempMin
     &         .AND. tempVar .NE. 0. _d 0 )
     &         adtempVar = dragCoeff*adCwatC(i,j,bi,bj)
     &         / ( 2. _d 0*SQRT(tempVar) )
          adCwatC(i,j,bi,bj) = 0. _d 0
C     reverse of tempVar (:74-93)
          tb0 = 2. _d 0*(
     &           ( uIceLoc(i,j,bi,bj)-uVel(i,j,kSrf,bi,bj) )*mW1
     &         + ( uIceLoc(i+1,j,bi,bj)-uVel(i+1,j,kSrf,bi,bj) )*mW2 )
     &         * 0.25 _d 0*adtempVar
          tb1 = 2. _d 0*(
     &           ( vIceLoc(i,j,bi,bj)-vVel(i,j,kSrf,bi,bj) )*mS1
     &         + ( vIceLoc(i,j+1,bi,bj)-vVel(i,j+1,kSrf,bi,bj) )*mS2 )
     &         * 0.25 _d 0*adtempVar
          advIceLoc(i,j  ,bi,bj) = advIceLoc(i,j  ,bi,bj) + mS1*tb1
          advVel(i,j  ,kSrf,bi,bj) = advVel(i,j  ,kSrf,bi,bj) - mS1*tb1
          advIceLoc(i,j+1,bi,bj) = advIceLoc(i,j+1,bi,bj) + mS2*tb1
          advVel(i,j+1,kSrf,bi,bj) = advVel(i,j+1,kSrf,bi,bj) - mS2*tb1
          aduIceLoc(i  ,j,bi,bj) = aduIceLoc(i  ,j,bi,bj) + mW1*tb0
          aduVel(i  ,j,kSrf,bi,bj) = aduVel(i  ,j,kSrf,bi,bj) - mW1*tb0
          aduIceLoc(i+1,j,bi,bj) = aduIceLoc(i+1,j,bi,bj) + mW2*tb0
          aduVel(i+1,j,kSrf,bi,bj) = aduVel(i+1,j,kSrf,bi,bj) - mW2*tb0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* SEAICE_CGRID */

      RETURN
      END

#endif /* SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ */
