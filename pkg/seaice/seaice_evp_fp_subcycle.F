#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#else
# define OBCS_UVICE_OLD
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C     The modified subcycle F_mod of the fixed-point adjoint of the mEVP
C     subcycle loop (SEAICEuseEVPfixedPointAdj), and what it needs:
C       SEAICE_EVP_FP_FROZEN      :: the frozen state at the iterate x
C       SEAICE_EVP_FP_SUBCYCLE    :: F_mod, one modified subcycle
C     (SEAICE_FP_OCEANDRAG_COEFFS, SEAICE_OCEANDRAG_COEFFS with the
C     ocean velocity as an argument, is in seaice_fp_oceandrag_coeffs.F:
C     its CPP arms follow seaice_oceandrag_coeffs.F, whose head does
C     not define OBCS_UVICE_OLD)
C     F_mod is one subcycle of SEAICE_EVP_SUBCYCLES (the loop body,
C     seaice_evp.F:610-1330) without its final EXCH_UV_XY_RL,
C     with the changes of the fixed-point adjoint (the "modified
C     subcycle", whose transpose is iterated). The changes act
C     only on the derivative: at the iterate x with the frozen state of
C     x (SEAICE_EVP_FP_FROZEN), F_mod(x) is the subcycle of x, except
C     where the raise of alpha and beta changes the pseudo time step
C     (fpRaise). The frozen fields are separate arguments, so that the
C     derivative of F_mod with respect to x and the inputs is the
C     derivative of the modified subcycle (the frozen state is held
C     constant).
C     The hand-written adjoint is seaice_evp_fp_subcycle_ad.F.

C     CPP options F_mod does not implement: the bodies are compiled out
C     (a STOP), and SEAICE_CHECK refuses SEAICEuseEVPfixedPointAdj
#if (defined SEAICE_ALLOW_SIDEDRAG || defined SEAICE_ALLOW_MOM_ADVECTION \
  || defined SEAICE_ALLOW_CLIPZETA || defined SEAICE_EVP_ELIMINATE_UNDERFLOWS \
  || defined SEAICE_ITD || defined SEAICE_DELTA_SMOOTHREG \
  || defined SEAICE_DYN_STABLE_ADJOINT || !defined OBCS_UVICE_OLD)
# define SEAICE_FPADJ_NOT_COVERED
#endif

#ifdef SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ

CBOP
C     !ROUTINE: SEAICE_EVP_FP_FROZEN
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_FP_FROZEN(
     I     uIce, vIce, press0, tensileStrFac, seaiceMassC,
     I     recip_ecc2, fpRaise, fpKappa0,
     O     uIceFr, vIceFr, deltaFr, raiseC,
     I     iEVPstep, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_FP_FROZEN
C     | o The frozen state of F_mod at the iterate (uIce, vIce):
C     |   the values that the modified subcycle holds constant
C     |   (they carry no derivative):
C     |   uIceFr, vIceFr :: the ice velocity of the frozen drag
C     |   deltaFr        :: Delta at C points (deltaC of
C     |                     seaice_evp.F:698-708 at
C     |                     the iterate)
C     |   raiseC         :: the raise factor r of alpha at C
C     |                     points (1 where not raised)
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
#include "SEAICE_GRID.h"

C     !INPUT/OUTPUT PARAMETERS:
C     uIce, vIce     :: the iterate's ice velocity
C     press0, tensileStrFac, seaiceMassC :: SEAICE.h fields
C     recip_ecc2     :: 1/e^2 of SEAICE_EVP (seaice_evp.F:156-159)
C     fpRaise        :: raise alpha and beta (else raiseC = 1)
C     fpKappa0       :: kappa0 of the raise; <= 0: 4/c per C point
      _RL uIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tensileStrFac (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL recip_ecc2
      LOGICAL fpRaise
      _RL fpKappa0
      _RL uIceFr        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceFr        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaFr       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL raiseC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER iEVPstep
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef SEAICE_FPADJ_NOT_COVERED
      STOP 'ABNORMAL END: F_mod: a SEAICE CPP option is not covered'
#endif
#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP ) \
  && !defined SEAICE_FPADJ_NOT_COVERED
C     === Local variables ===
C     e11L, e22L, e12L :: strain rates of the iterate
C     zetaFr         :: zeta at C points of the iterate
C     sR, qR, cR     :: s = gamma/(alpha beta), s/kappa0, the C-grid
C                       constant c
      INTEGER i, j, bi, bj
      _RL e11L (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22L (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12L (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ep    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL em    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL e12Csq(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL deltaSq, deltaCreg, tmp, zetaFr, sR, qR, cR

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          uIceFr(i,j,bi,bj) = uIce(i,j,bi,bj)
          vIceFr(i,j,bi,bj) = vIce(i,j,bi,bj)
          e11L(i,j,bi,bj)   = 0. _d 0
          e22L(i,j,bi,bj)   = 0. _d 0
          e12L(i,j,bi,bj)   = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C     the strain rates of the iterate (seaice_evp.F:622-627)
      CALL SEAICE_CALC_STRAINRATES(
     I     uIce, vIce,
     O     e11L, e22L, e12L,
     I     iEVPstep, myTime, myIter, myThid )
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          e12Csq(i,j)        = 0. _d 0
          deltaFr(i,j,bi,bj) = 0. _d 0
          raiseC(i,j,bi,bj)  = 1. _d 0
          ep(i,j) = e11L(i,j,bi,bj) + e22L(i,j,bi,bj)
          em(i,j) = e11L(i,j,bi,bj) - e22L(i,j,bi,bj)
         ENDDO
        ENDDO
C     Delta at C points, the statements of
C     seaice_evp.F:675-708
        IF ( SEAICEetaZmethod .EQ. 0 ) THEN
         DO j=1-OLy+1,sNy+OLy-1
          DO i=1-OLx+1,sNx+OLx-1
           tmp = 0.25 *
     &          ( e12L(i,j  ,bi,bj) + e12L(i+1,j  ,bi,bj)
     &          + e12L(i,j+1,bi,bj) + e12L(i+1,j+1,bi,bj) )
           e12Csq(i,j) = tmp*tmp
          ENDDO
         ENDDO
        ELSEIF ( SEAICEetaZmethod .EQ. 3 ) THEN
         DO j=1-OLy+1,sNy+OLy-1
          DO i=1-OLx+1,sNx+OLx-1
           e12Csq(i,j) = 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &          ( rAz(i  ,j  ,bi,bj)*e12L(i  ,j  ,bi,bj)**2
     &          + rAz(i+1,j  ,bi,bj)*e12L(i+1,j  ,bi,bj)**2
     &          + rAz(i  ,j+1,bi,bj)*e12L(i  ,j+1,bi,bj)**2
     &          + rAz(i+1,j+1,bi,bj)*e12L(i+1,j+1,bi,bj)**2 )
          ENDDO
         ENDDO
        ENDIF
        DO j=0,sNy+1
         DO i=0,sNx+1
          deltaSq = ep(i,j)**2 + recip_ecc2 * em(i,j)**2
     &         + recip_ecc2 * 4. _d 0 * e12Csq(i,j)
#ifdef ALLOW_AUTODIFF
          deltaFr(i,j,bi,bj) = 0. _d 0
          IF ( deltaSq .GT. 0. _d 0 )
     &         deltaFr(i,j,bi,bj) = SQRT(deltaSq)
#else
          deltaFr(i,j,bi,bj) = SQRT(deltaSq)
#endif /* ALLOW_AUTODIFF */
         ENDDO
        ENDDO
C     the raise factor of alpha and beta: r = MAX(1, SQRT(s/kappa0))
C     on ice (seaiceMassC > 0 and HEFFM > 0), 1 elsewhere, with
C     s = zeta deltaTdyn recip_rA / seaiceMassC / (alpha beta) and
C     zeta = zetaC of seaice_evp.F:716-720 at the iterate;
C     kappa0 = 4/c,
C     c = 4 (1 + recip_ecc2) rA (recip_dxF^2 + recip_dyF^2), where
C     fpKappa0 <= 0 (the default)
        IF ( fpRaise ) THEN
         DO j=0,sNy+1
          DO i=0,sNx+1
           deltaCreg = MAX(deltaFr(i,j,bi,bj),SEAICE_deltaMin)
           zetaFr = HALF*( press0(i,j,bi,bj)
     &          * ( 1. _d 0 + tensileStrFac(i,j,bi,bj) )
     &          )/deltaCreg
           IF ( seaiceMassC(i,j,bi,bj) .GT. 0. _d 0 .AND.
     &          HEFFM(i,j,bi,bj) .GT. 0. _d 0 ) THEN
            sR = zetaFr * SEAICE_deltaTdyn * recip_rA(i,j,bi,bj)
     &           / seaiceMassC(i,j,bi,bj)
     &           / ( SEAICE_evpAlpha * SEAICE_evpBeta )
            IF ( fpKappa0 .GT. 0. _d 0 ) THEN
             qR = sR / fpKappa0
            ELSE
             cR = 4. _d 0 * ( 1. _d 0 + recip_ecc2 ) * rA(i,j,bi,bj)
     &          * ( recip_dxF(i,j,bi,bj)*recip_dxF(i,j,bi,bj)
     &            + recip_dyF(i,j,bi,bj)*recip_dyF(i,j,bi,bj) )
             qR = 0. _d 0
             IF ( cR .GT. 0. _d 0 ) qR = sR / ( 4. _d 0 / cR )
            ENDIF
            raiseC(i,j,bi,bj) = MAX( 1. _d 0, SQRT(qR) )
           ENDIF
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_EVP_FP_SUBCYCLE
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_FP_SUBCYCLE(
     U     uIce, vIce,
     U     seaice_sigma1, seaice_sigma2, seaice_sigma12,
     U     e11, e22, e12, zetaZ, deltaC,
     U     stressDivergenceX, stressDivergenceY,
     U     DWATN, CbotC,
     U     FORCEX, FORCEY, press, zeta, eta,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA,
     I     areaW, areaS,
     I     uIceFr, vIceFr, deltaFr, raiseC,
     I     fpTheta, fpThetaRigid, fpRigidDelta,
     I     fpThinArea, fpThinMass, fpRaise, fpPrLive,
     I     kSrf, COSWAT, SINWAT, recip_ecc2,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     iEVPstep, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_FP_SUBCYCLE
C     | o F_mod: one subcycle of the mEVP loop (the loop body of
C     |   SEAICE_EVP_SUBCYCLES, seaice_evp.F:610-1330)
C     |   without its final EXCH_UV_XY_RL (:1271), with the
C     |   changes of the modified subcycle:
C     |   - theta on dzeta/dDelta: zeta from
C     |     Delta_theta = deltaFr + theta*(deltaC - deltaFr),
C     |     theta = fpThetaRigid where deltaFr < fpRigidDelta;
C     |   - the replacement pressure from deltaFr (fpPrLive:
C     |     from deltaC); PRESS0 and tensileStrFac live;
C     |   - alpha, beta raised by raiseC (fpRaise);
C     |   - DWATN from (uIceFr, vIceFr) where AREA < fpThinArea
C     |     or seaiceMassC < fpThinMass, else from (uIce, vIce).
C     |   With the frozen state of x (SEAICE_EVP_FP_FROZEN) and
C     |   fpRaise = .FALSE., F_mod(x) = the subcycle of x.
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
#include "SEAICE_GRID.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     The fields keep the names of SEAICE_EVP_SUBCYCLES; every active
C     field is an argument (none passes through COMMON).
C     carry x (read and written; grid points the subcycle does not write keep
C     their input value): uIce, vIce, seaice_sigma1/2/12, e11, e22,
C     e12, zetaZ, deltaC, stressDivergenceX/Y, DWATN, CbotC, FORCEX,
C     FORCEY, press, zeta, eta
C     denom1, denom2, evpAlphaC/Z, evpBetaU/V :: constant mEVP fields
C                       of SEAICE_EVP (no adaptive EVP)
C     press0 .. areaS :: inputs p (as SEAICE_EVP_SUBCYCLES)
C     uIceFr, vIceFr, deltaFr, raiseC :: the frozen state
C                       (SEAICE_EVP_FP_FROZEN; passive)
C     fpTheta .. fpPrLive :: the settings of the modified subcycle
C     kSrf .. recip_deltaT :: scalars of SEAICE_EVP
C     iEVPstep       :: the subcycle number (passed to the callees)
      _RL uIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetaZ         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DWATN         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL CbotC         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zeta          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL eta           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL denom1        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL denom2        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpAlphaC     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpAlphaZ     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpBetaU      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpBetaV      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tensileStrFac (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassU   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassV   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uIceNm1       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceNm1       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uVel          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVel          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HEFF          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL AREA          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL areaW         (1:sNx,1:sNy,nSx,nSy)
      _RL areaS         (1:sNx,1:sNy,nSx,nSy)
      _RL uIceFr        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceFr        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaFr       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL raiseC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fpTheta, fpThetaRigid, fpRigidDelta
      _RL fpThinArea, fpThinMass
      LOGICAL fpRaise, fpPrLive
      INTEGER kSrf
      _RL COSWAT
      _RS SINWAT
      _RL recip_ecc2
      _RL evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT
      INTEGER iEVPstep
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef SEAICE_FPADJ_NOT_COVERED
      STOP 'ABNORMAL END: F_mod: a SEAICE CPP option is not covered'
#endif
#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP ) \
  && !defined SEAICE_FPADJ_NOT_COVERED

C     === Local variables ===
C     the locals of SEAICE_EVP_SUBCYCLES, and
C     dwatnFr        :: DWATN of the frozen ice velocity
C     aC, aZ, d1, d2, bU, bV :: alpha at C, Z points, denom1, denom2,
C                       beta at U, V points (raised with fpRaise)
C     thetaC, deltaTh :: theta of the point, Delta_theta
C     deltaP, zetaP  :: Delta and zeta of the replacement pressure
      INTEGER i, j, bi, bj
      _RL seaice_div    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL seaice_tension(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL seaice_shear  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig11         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig22         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ep            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL em            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL e12Csq        (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL pressC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL zetaC         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL aC            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL aZ            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL d1            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL d2            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL bU            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL bV            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dwatnFr       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaCreg, deltaSq, tmp
      _RL thetaC, deltaTh, deltaP, zetaP
      _RL sumNorm, denomU, denomV
      _RL locMaskU, locMaskV
      _RL betaFacP1U,  betaFacP1V
      _RL betaFacU,    betaFacV

C     first calculate strain rates and bulk moduli/viscosities
C     (seaice_evp.F:622-627)
      CALL SEAICE_CALC_STRAINRATES(
     I     uIce, vIce,
     O     e11, e22, e12,
     I     iEVPstep, myTime, myIter, myThid )

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          seaice_div    (i,j) = 0. _d 0
          seaice_tension(i,j) = 0. _d 0
          seaice_shear  (i,j) = 0. _d 0
          pressC        (i,j) = 0. _d 0
          e12Csq        (i,j) = 0. _d 0
          zetaC         (i,j) = 0. _d 0
          zetaZ   (i,j,bi,bj) = 0. _d 0
          deltaC  (i,j,bi,bj) = 0. _d 0
C     alpha, denom1/2, beta of SEAICE_EVP (constant without aEVP)
          aC(i,j) = evpAlphaC(i,j,bi,bj)
          aZ(i,j) = evpAlphaZ(i,j,bi,bj)
          d1(i,j) = denom1(i,j,bi,bj)
          d2(i,j) = denom2(i,j,bi,bj)
          bU(i,j) = evpBetaU(i,j,bi,bj)
          bV(i,j) = evpBetaV(i,j,bi,bj)
         ENDDO
        ENDDO
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          ep(i,j) = e11(i,j,bi,bj) + e22(i,j,bi,bj)
          em(i,j) = e11(i,j,bi,bj) - e22(i,j,bi,bj)
         ENDDO
        ENDDO
C     average strain rates to C points (seaice_evp.F:675-697)
        IF ( SEAICEetaZmethod .EQ. 0 ) THEN
         DO j=1-OLy+1,sNy+OLy-1
          DO i=1-OLx+1,sNx+OLx-1
           tmp = 0.25 *
     &          ( e12(i,j  ,bi,bj) + e12(i+1,j  ,bi,bj)
     &          + e12(i,j+1,bi,bj) + e12(i+1,j+1,bi,bj) )
           e12Csq(i,j) = tmp*tmp
          ENDDO
         ENDDO
        ELSEIF ( SEAICEetaZmethod .EQ. 3 ) THEN
         DO j=1-OLy+1,sNy+OLy-1
          DO i=1-OLx+1,sNx+OLx-1
           e12Csq(i,j) = 0.25 _d 0 * recip_rA(i,j,bi,bj) *
     &          ( rAz(i  ,j  ,bi,bj)*e12(i  ,j  ,bi,bj)**2
     &          + rAz(i+1,j  ,bi,bj)*e12(i+1,j  ,bi,bj)**2
     &          + rAz(i  ,j+1,bi,bj)*e12(i  ,j+1,bi,bj)**2
     &          + rAz(i+1,j+1,bi,bj)*e12(i+1,j+1,bi,bj)**2 )
          ENDDO
         ENDDO
        ENDIF
        DO j=0,sNy+1
         DO i=0,sNx+1
          deltaSq = ep(i,j)**2 + recip_ecc2 * em(i,j)**2
     &         + recip_ecc2 * 4. _d 0 * e12Csq(i,j)
#ifdef ALLOW_AUTODIFF
C     avoid sqrt of 0
          deltaC(i,j,bi,bj) = 0. _d 0
          IF ( deltaSq .GT. 0. _d 0 )
     &         deltaC(i,j,bi,bj) = SQRT(deltaSq)
#else
          deltaC(i,j,bi,bj) = SQRT(deltaSq)
#endif /* ALLOW_AUTODIFF */
C     fixed-point adjoint: zeta from Delta_theta, whose value is
C     deltaC's (deltaFr = deltaC at the iterate) and whose derivative
C     is theta times deltaC's; theta = fpThetaRigid at the near-rigid
C     points deltaFr < fpRigidDelta
          thetaC = fpTheta
          IF ( deltaFr(i,j,bi,bj) .LT. fpRigidDelta )
     &         thetaC = fpThetaRigid
          deltaTh = deltaFr(i,j,bi,bj)
     &         + thetaC*( deltaC(i,j,bi,bj) - deltaFr(i,j,bi,bj) )
          deltaCreg  = MAX(deltaTh,SEAICE_deltaMin)
          zetaC(i,j) = HALF*( press0(i,j,bi,bj)
     &         * ( 1. _d 0 + tensileStrFac(i,j,bi,bj) )
     &         )/deltaCreg
C     the replacement pressure from deltaFr (P_r frozen; fpPrLive: live)
          deltaP = deltaFr(i,j,bi,bj)
          IF ( fpPrLive ) deltaP = deltaC(i,j,bi,bj)
          zetaP = HALF*( press0(i,j,bi,bj)
     &         * ( 1. _d 0 + tensileStrFac(i,j,bi,bj) )
     &         )/MAX(deltaP,SEAICE_deltaMin)
C     recompute pressure (seaice_evp.F:782-791, with zetaP,
C     deltaP)
          pressC(i,j) =
     &         ( press0(i,j,bi,bj) * ( 1. _d 0 - SEAICEpressReplFac )
     &         + TWO*zetaP*deltaP*SEAICEpressReplFac
     &           /(1. _d 0 + tensileStrFac(i,j,bi,bj))
     &         ) * (1. _d 0 - tensileStrFac(i,j,bi,bj))
         ENDDO
        ENDDO
C     fixed-point adjoint: raised alpha at C
        IF ( fpRaise ) THEN
         DO j=0,sNy+1
          DO i=0,sNx+1
           aC(i,j) = SEAICE_evpAlpha*raiseC(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDIF
C     compute zetaZ by simple averaging (seaice_evp.F:743-757;
C     deltaZ,
C     read only with SEAICEuseTEM, is not formed)
        DO j=1,sNy+1
         DO i=1,sNx+1
          sumNorm = HEFFM(i,j,  bi,bj)+HEFFM(i-1,j,  bi,bj)
     &            + HEFFM(i,j-1,bi,bj)+HEFFM(i-1,j-1,bi,bj)
          IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
          zetaZ(i,j,bi,bj) = sumNorm *
     &         ( zetaC(i,  j) + zetaC(i-1,j-1)
     &         + zetaC(i-1,j) + zetaC(i,  j-1) )
         ENDDO
        ENDDO
#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
C     save eta, zeta, pressure for diagnostics
C     (seaice_evp.F:792-803)
         DO j=1,sNy
          DO i=1,sNx
           press(i,j,bi,bj) = pressC(i,j)
           zeta (i,j,bi,bj) = zetaC(i,j)
           eta  (i,j,bi,bj) = zetaC(i,j)*recip_ecc2
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */
C     the RHS of the stress equations (seaice_evp.F:855-874,
C     the arm
C     without SEAICEuseTEM)
        DO j=0,sNy
         DO i=0,sNx
          seaice_div    (i,j) =
     &         ( 2. _d 0 *zetaC(i,j)*ep(i,j) - pressC(i,j)
     &         ) * HEFFM(i,j,bi,bj)
          seaice_tension(i,j) = 2. _d 0*zetaC(i,j)
     &         * em(i,j) * HEFFM(i,j,bi,bj)
         ENDDO
        ENDDO
        DO j=1,sNy+1
         DO i=1,sNx+1
          seaice_shear  (i,j) =
     &         2. _d 0*zetaZ(i,j,bi,bj)*e12(i,j,bi,bj)
         ENDDO
        ENDDO
C     first step stress equations; fixed-point adjoint:
C     denom1, denom2 from the raised alpha as
C     in the aEVP arm seaice_evp.F:884-891
        IF ( fpRaise ) THEN
         DO j=0,sNy
          DO i=0,sNx
           d1(i,j) = 1. _d 0 /  aC(i,j)
           d2(i,j) = d1(i,j)
          ENDDO
         ENDDO
        ENDIF
        DO j=0,sNy
         DO i=0,sNx
C     sigma1 and sigma2 are computed on C points
C     (seaice_evp.F:896-926)
          seaice_sigma1 (i,j,bi,bj) = ( seaice_sigma1 (i,j,bi,bj)
     &         * ( aC(i,j) - evpRevFac )
     &         + seaice_div(i,j)
     &         ) * d1(i,j)
     &         *HEFFM(i,j,bi,bj)
          seaice_sigma2 (i,j,bi,bj) = ( seaice_sigma2 (i,j,bi,bj)
     &         * ( aC(i,j) - evpRevFac )
     &         + seaice_tension(i,j)*recip_evpRevFac
     &         ) * d2(i,j)
     &         *HEFFM(i,j,bi,bj)
C     recover sigma11 and sigma22
          sig11(i,j) = 0.5 _d 0 *
     &         ( seaice_sigma1(i,j,bi,bj)+seaice_sigma2(i,j,bi,bj) )
          sig22(i,j) = 0.5 _d 0 *
     &         ( seaice_sigma1(i,j,bi,bj)-seaice_sigma2(i,j,bi,bj) )
         ENDDO
        ENDDO
C     sigma12 is computed on Z points; fixed-point adjoint:
C     alpha at Z points = alpha times the
C     maximum of raiseC over the 4 adjacent C points (the aEVP arm
C     seaice_evp.F:929-938 averages)
        IF ( fpRaise ) THEN
         DO j=1,sNy+1
          DO i=1,sNx+1
           aZ(i,j) = SEAICE_evpAlpha * MAX(
     &          MAX(raiseC(i,  j,bi,bj),raiseC(i-1,j-1,bi,bj)),
     &          MAX(raiseC(i-1,j,bi,bj),raiseC(i,  j-1,bi,bj)) )
           d2(i,j) = 1. _d 0 /  aZ(i,j)
          ENDDO
         ENDDO
        ENDIF
        DO j=1,sNy+1
         DO i=1,sNx+1
          seaice_sigma12(i,j,bi,bj) = ( seaice_sigma12(i,j,bi,bj)
     &         * ( aZ(i,j) - evpRevFac )
     &         + seaice_shear(i,j)*recip_evpRevFac
     &         ) * d2(i,j)
         ENDDO
        ENDDO
C     compute divergence of stress tensor
C     (seaice_evp.F:956-986)
        DO j=1,sNy
         DO i=1,sNx
          stressDivergenceX(i,j,bi,bj) =
     &         ( sig11(i  ,j  ) * _dyF(i  ,j,bi,bj)
     &         - sig11(i-1,j  ) * _dyF(i-1,j,bi,bj)
     &         + seaice_sigma12(i,j+1,bi,bj) * _dxV(i,j+1,bi,bj)
     &         - seaice_sigma12(i,j  ,bi,bj) * _dxV(i,j  ,bi,bj)
     &         ) * recip_rAw(i,j,bi,bj)
C     extra metric terms
     &         + k2AtU(i,j,bi,bj)
     &         * HALF*( seaice_sigma12(i,j+1,bi,bj)
     &                + seaice_sigma12(i,j  ,bi,bj) )
     &         - k1AtU(i,j,bi,bj)
     &         * HALF*( sig22(i,j) + sig22(i-1,j) )
          stressDivergenceY(i,j,bi,bj) =
     &         ( sig22(i,j  ) * _dxF(i,j  ,bi,bj)
     &         - sig22(i,j-1) * _dxF(i,j-1,bi,bj)
     &         + seaice_sigma12(i+1,j,bi,bj) * _dyU(i+1,j,bi,bj)
     &         - seaice_sigma12(i  ,j,bi,bj) * _dyU(i  ,j,bi,bj)
     &         ) * recip_rAs(i,j,bi,bj)
C     extra metric terms
     &         + k1AtV(i,j,bi,bj)
     &         * HALF*( seaice_sigma12(i  ,j,bi,bj)
     &                + seaice_sigma12(i+1,j,bi,bj) )
     &         - k2AtV(i,j,bi,bj)
     &         * HALF*( sig11(i,j) + sig11(i,j-1) )
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C     set up rhs for stepping the velocity field
C     (seaice_evp.F:1051-1085)
C     fixed-point adjoint: DWATN of the frozen ice velocity where the
C     ice is
C     thin (AREA < fpThinArea or seaiceMassC < fpThinMass), of the
C     iterate elsewhere; the value is the same at the iterate
      CALL SEAICE_FP_OCEANDRAG_COEFFS(
     I     uIce, vIce, HEFFM, uVel, vVel,
     O     DWATN,
     I     iEVPstep, myTime, myIter, myThid )
      CALL SEAICE_FP_OCEANDRAG_COEFFS(
     I     uIceFr, vIceFr, HEFFM, uVel, vVel,
     O     dwatnFr,
     I     iEVPstep, myTime, myIter, myThid )
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
          IF ( AREA(i,j,bi,bj) .LT. fpThinArea .OR.
     &         seaiceMassC(i,j,bi,bj) .LT. fpThinMass )
     &         DWATN(i,j,bi,bj) = dwatnFr(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#ifdef SEAICE_ALLOW_BOTTOMDRAG
      CALL SEAICE_BOTTOMDRAG_COEFFS(
     I     uIce, vIce, HEFFM,
     I     HEFF, AREA,
     O     CbotC,
     I     iEVPstep, myTime, myIter, myThid )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          bU(i,j) = evpBetaU(i,j,bi,bj)
          bV(i,j) = evpBetaV(i,j,bi,bj)
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
C     masks of the drag term, anti symmetric drag force, ice ocean
C     stress and Coriolis terms (seaice_evp.F:1086-1148)
          locMaskU = seaiceMassU(i,j,bi,bj)
          locMaskV = seaiceMassV(i,j,bi,bj)
          IF ( locMaskU .NE. 0. _d 0 ) locMaskU = 1. _d 0
          IF ( locMaskV .NE. 0. _d 0 ) locMaskV = 1. _d 0
          FORCEX(i,j,bi,bj)=FORCEX0(i,j,bi,bj)+
     &         ( 0.5 _d 0 * ( DWATN(i,j,bi,bj)+DWATN(i-1,j,bi,bj) ) *
     &         COSWAT * uVel(i,j,kSrf,bi,bj)
     &         - SIGN(SINWAT, _fCori(i,j,bi,bj))* 0.5 _d 0 *
     &         ( DWATN(i  ,j,bi,bj) * 0.5 _d 0 *
     &          (vVel(i  ,j  ,kSrf,bi,bj)-vIce(i  ,j  ,bi,bj)
     &          +vVel(i  ,j+1,kSrf,bi,bj)-vIce(i  ,j+1,bi,bj))
     &         + DWATN(i-1,j,bi,bj) * 0.5 _d 0 *
     &          (vVel(i-1,j  ,kSrf,bi,bj)-vIce(i-1,j  ,bi,bj)
     &          +vVel(i-1,j+1,kSrf,bi,bj)-vIce(i-1,j+1,bi,bj))
     &         )*locMaskU ) * areaW(i,j,bi,bj)
          FORCEY(i,j,bi,bj)=FORCEY0(i,j,bi,bj)+
     &         ( 0.5 _d 0 * ( DWATN(i,j,bi,bj)+DWATN(i,j-1,bi,bj) ) *
     &         COSWAT * vVel(i,j,kSrf,bi,bj)
     &         + SIGN(SINWAT, _fCori(i,j,bi,bj)) * 0.5 _d 0 *
     &         ( DWATN(i,j  ,bi,bj) * 0.5 _d 0 *
     &          (uVel(i  ,j  ,kSrf,bi,bj)-uIce(i  ,j  ,bi,bj)
     &          +uVel(i+1,j  ,kSrf,bi,bj)-uIce(i+1,j  ,bi,bj))
     &         + DWATN(i,j-1,bi,bj) * 0.5 _d 0 *
     &          (uVel(i  ,j-1,kSrf,bi,bj)-uIce(i  ,j-1,bi,bj)
     &          +uVel(i+1,j-1,kSrf,bi,bj)-uIce(i+1,j-1,bi,bj))
     &         )*locMaskV ) * areaS(i,j,bi,bj)
C     coriols terms
          FORCEX(i,j,bi,bj)=FORCEX(i,j,bi,bj) + 0.5 _d 0*(
     &           seaiceMassC(i  ,j,bi,bj) * _fCori(i  ,j,bi,bj)
     &         * 0.5 _d 0*( vIce(i  ,j,bi,bj)+vIce(i  ,j+1,bi,bj) )
     &         + seaiceMassC(i-1,j,bi,bj) * _fCori(i-1,j,bi,bj)
     &         * 0.5 _d 0*( vIce(i-1,j,bi,bj)+vIce(i-1,j+1,bi,bj) )
     &         )
          FORCEY(i,j,bi,bj)=FORCEY(i,j,bi,bj) - 0.5 _d 0*(
     &           seaiceMassC(i,j  ,bi,bj) * _fCori(i,j  ,bi,bj)
     &         * 0.5 _d 0*( uIce(i,j  ,bi,bj)+uIce(i+1,  j,bi,bj) )
     &         + seaiceMassC(i,j-1,bi,bj) * _fCori(i,j-1,bi,bj)
     &         * 0.5 _d 0*( uIce(i,j-1,bi,bj)+uIce(i+1,j-1,bi,bj) )
     &         )
         ENDDO
        ENDDO
C     step momentum equations with ice-ocean stress term treated
C     implicitly (seaice_evp.F:1174-1251); fixed-point adjoint:
C     beta at U, V points = beta times the
C     maximum of raiseC over the 2 adjacent C points (the aEVP arm
C     seaice_evp.F:1177-1188 averages alpha)
        IF ( fpRaise ) THEN
         DO j=1,sNy
          DO i=1,sNx
           bU(i,j) = SEAICE_evpBeta
     &          * MAX(raiseC(i-1,j,bi,bj),raiseC(i,j,bi,bj))
           bV(i,j) = SEAICE_evpBeta
     &          * MAX(raiseC(i,j-1,bi,bj),raiseC(i,j,bi,bj))
          ENDDO
         ENDDO
        ENDIF
        DO j=1,sNy
         DO i=1,sNx
          betaFacU   = bU(i,j)*recip_deltaT
          betaFacV   = bV(i,j)*recip_deltaT
          tmp        = evpStarFac*recip_deltaT
          betaFacP1V = betaFacV + tmp
          betaFacP1U = betaFacU + tmp
C     (SEAICE_evpAreaReg > 0, seaice_evp.F:1202-1210, is
C     refused by
C     SEAICE_CHECK with SEAICEuseEVPfixedPointAdj)
          denomU = seaiceMassU(i,j,bi,bj)*betaFacP1U
     &         + 0.5 _d 0*( DWATN(i,j,bi,bj) + DWATN(i-1,j,bi,bj) )
     &         * COSWAT * areaW(i,j,bi,bj)
          denomV = seaiceMassV(i,j,bi,bj)*betaFacP1V
     &         + 0.5 _d 0*( DWATN(i,j,bi,bj) + DWATN(i,j-1,bi,bj) )
     &         * COSWAT * areaS(i,j,bi,bj)
#ifdef SEAICE_ALLOW_BOTTOMDRAG
          denomU = denomU + areaW(i,j,bi,bj)
     &         * 0.5 _d 0*( CbotC(i,j,bi,bj) + CbotC(i-1,j,bi,bj) )
          denomV = denomV + areaS(i,j,bi,bj)
     &         * 0.5 _d 0*( CbotC(i,j,bi,bj) + CbotC(i,j-1,bi,bj) )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
          IF ( denomU .EQ. 0. _d 0 ) denomU = 1. _d 0
          IF ( denomV .EQ. 0. _d 0 ) denomV = 1. _d 0
          uIce(i,j,bi,bj) = seaiceMaskU(i,j,bi,bj) *
     &         ( seaiceMassU(i,j,bi,bj)*betaFacU
     &         * uIce(i,j,bi,bj)
     &         + seaiceMassU(i,j,bi,bj)*recip_deltaT*evpStarFac
     &         * uIceNm1(i,j,bi,bj)
     &         + FORCEX(i,j,bi,bj)
     &         + stressDivergenceX(i,j,bi,bj) ) / denomU
          vIce(i,j,bi,bj) = seaiceMaskV(i,j,bi,bj) *
     &         ( seaiceMassV(i,j,bi,bj)*betaFacV
     &         * vIce(i,j,bi,bj)
     &         + seaiceMassV(i,j,bi,bj)*recip_deltaT*evpStarFac
     &         * vIceNm1(i,j,bi,bj)
     &         + FORCEY(i,j,bi,bj)
     &         + stressDivergenceY(i,j,bi,bj) ) / denomV
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C     (no CALL EXCH_UV_XY_RL here: the driver applies the adjoint of
C     the exchange separately, SEAICE_FP_EXCH_UV_XY_RL_AD)

#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END

#endif /* SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ */
