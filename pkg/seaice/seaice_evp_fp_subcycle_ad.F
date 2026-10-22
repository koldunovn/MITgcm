#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#else
# define OBCS_UVICE_OLD
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C     The hand-written adjoint F_mod_b of the modified subcycle F_mod
C     (SEAICE_EVP_FP_SUBCYCLE, seaice_evp_fp_subcycle.F) of the
C     fixed-point adjoint of the mEVP subcycle loop:
C       SEAICE_EVP_FP_SUBCYCLE_AD     :: F_mod_b
C       (SEAICE_FP_OCEANDRAG_COEFFS_AD, the adjoint of the drag variant,
C        is in seaice_fp_oceandrag_coeffs.F)
C       SEAICE_FP_STRAINRATES_AD      :: adjoint of
C                                        SEAICE_CALC_STRAINRATES
C     Written by hand, statement by statement after the reverse mode of
C     F_mod that Tapenade generates (used as a reference only; it is
C     not part of this code), with no AD runtime: the forward values
C     the reverse sweep needs
C     are recomputed (F_mod is called on work copies, and the tile-local
C     intermediates are recomputed in the reverse loops) instead of
C     PUSH/POP. Calling convention:
C     - the carry x (uIce .. eta) is input only; its adjoint (aduIce ..
C       adeta) is overwritten: lambda on input, F_x^T lambda on output
C       (grid points F_mod does not write pass lambda through);
C     - the adjoints of the inputs p (adpress0 .. adareaS) are
C       incremented by F_p^T lambda;
C     - the frozen state, the settings and the constant mEVP fields are
C       passive.

C     CPP options F_mod does not implement (seaice_evp_fp_subcycle.F)
#if (defined SEAICE_ALLOW_SIDEDRAG || defined SEAICE_ALLOW_MOM_ADVECTION \
  || defined SEAICE_ALLOW_CLIPZETA || defined SEAICE_EVP_ELIMINATE_UNDERFLOWS \
  || defined SEAICE_ITD || defined SEAICE_DELTA_SMOOTHREG \
  || defined SEAICE_DYN_STABLE_ADJOINT || !defined OBCS_UVICE_OLD)
# define SEAICE_FPADJ_NOT_COVERED
#endif

#ifdef SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ

CBOP
C     !ROUTINE: SEAICE_EVP_FP_SUBCYCLE_AD
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_FP_SUBCYCLE_AD(
     I     uIce, vIce,
     I     seaice_sigma1, seaice_sigma2, seaice_sigma12,
     I     e11, e22, e12, zetaZ, deltaC,
     I     stressDivergenceX, stressDivergenceY,
     I     DWATN, CbotC,
     I     FORCEX, FORCEY, press, zeta, eta,
     U     aduIce, advIce,
     U     adsigma1, adsigma2, adsigma12,
     U     ade11, ade22, ade12, adzetaZ, addeltaC,
     U     adstrDivX, adstrDivY,
     U     adDWATN, adCbotC,
     U     adFORCEX, adFORCEY, adpress, adzeta, adeta,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA,
     I     areaW, areaS,
     U     adpress0, adtensileStrFac,
     U     adseaiceMassC, adseaiceMassU, adseaiceMassV,
     U     adFORCEX0, adFORCEY0, aduIceNm1, advIceNm1,
     U     aduVel, advVel, adHEFF, adAREA,
     U     adareaW, adareaS,
     I     uIceFr, vIceFr, deltaFr, raiseC,
     I     fpTheta, fpThetaRigid, fpRigidDelta,
     I     fpThinArea, fpThinMass, fpRaise, fpPrLive,
     I     kSrf, COSWAT, SINWAT, recip_ecc2,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     iEVPstep, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_FP_SUBCYCLE_AD
C     | o F_mod_b: the adjoint of F_mod at the iterate x with the
C     |   frozen state (uIceFr, vIceFr, deltaFr, raiseC) and the
C     |   settings fpTheta .. fpPrLive (seaice_evp_fp_subcycle.F)
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
C     uIce .. eta    :: the carry x (as SEAICE_EVP_FP_SUBCYCLE; input)
C     aduIce .. adeta :: its adjoint (lambda in, F_x^T lambda out)
C     denom1 .. evpBetaV :: constant mEVP fields of SEAICE_EVP
C     press0 .. areaS :: the inputs p
C     adpress0 .. adareaS :: their adjoints (incremented)
C     uIceFr .. raiseC :: the frozen state (SEAICE_EVP_FP_FROZEN)
C     fpTheta .. fpPrLive :: the settings of F_mod
C     kSrf .. recip_deltaT :: scalars of SEAICE_EVP
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
      _RL aduIce        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL advIce        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adsigma1      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adsigma2      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adsigma12     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade11         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade22         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade12         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adzetaZ       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL addeltaC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adstrDivX     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adstrDivY     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adDWATN       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adCbotC       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEX      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEY      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adpress       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adzeta        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adeta         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
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
      _RL adpress0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adtensileStrFac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaiceMassC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaiceMassU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaiceMassV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEX0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEY0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL aduIceNm1     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL advIceNm1     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL aduVel        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL advVel        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL adHEFF        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adAREA        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adareaW       (1:sNx,1:sNy,nSx,nSy)
      _RL adareaS       (1:sNx,1:sNy,nSx,nSy)
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
      STOP 'ABNORMAL END: F_mod_b: a SEAICE CPP option is not covered'
#endif
#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP ) \
  && !defined SEAICE_FPADJ_NOT_COVERED

C     === Local variables ===
C     y*             :: F_mod(x) (the forward, recomputed on copies)
C     ad*            :: adjoints of the tile-local fields of F_mod
C     dwFrAd, uFrAd, vFrAd :: adjoint of the frozen drag coefficient
C                       and (discarded) of the frozen ice velocity
C     d2C, d2Z       :: denom2 as the sigma2 (C points) and the sigma12
C                       (Z points) updates read it
      INTEGER i, j, bi, bj
      _RL yuIce (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yvIce (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ysig1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ysig2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ysig12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ye11  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ye22  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ye12  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yzetaZ(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ydeltaC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ysdX  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ysdY  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yDWATN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yCbotC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yFX   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yFY   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ypress(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yzeta (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yeta  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dwFrAd(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uFrAd (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vFrAd (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ep            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL em            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL zetaC         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL aC            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL aZ            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL d1            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL d2C           (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL d2Z           (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL bU            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL bV            (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adep          (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adem          (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL ade12Csq      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adzetaC       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adpressC      (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL addiv         (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adtension     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adshear       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adsig11       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL adsig22       (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL deltaSq, deltaTh, deltaCreg, thetaC, deltaP, zetaP, maxP
      _RL tmp, sumNorm, denomU, denomV, numU, numV
      _RL locMaskU, locMaskV, sgnFac, aA, aB
      _RL betaFacP1U, betaFacP1V, betaFacU, betaFacV
      _RL tb0, tb1, tb2, tb3, tb4, tb5, tb6, tb7, tb8
      _RL addenomU, addenomV, adlocMaskU, adlocMaskV
      _RL addeltaSq, addeltaTh, addeltaCreg
      _RL addeltaP, adzetaP, admaxP
      _RL adtmp
      LOGICAL zeroU, zeroV

#ifdef SEAICE_ALLOW_BOTTOMDRAG
      IF ( SEAICEbasalDragK2 .GT. 0. _d 0 ) THEN
C     the bottom drag coefficient with SEAICEbasalDragK2 > 0 has no
C     hand-written adjoint (refused by SEAICE_CHECK)
       STOP 'ABNORMAL END: F_mod_b: SEAICEbasalDragK2 > 0'
      ENDIF
#endif
      IF ( SEAICE_no_slip .AND. SEAICE_2ndOrderBC ) THEN
       STOP 'ABNORMAL END: F_mod_b: SEAICE_2ndOrderBC'
      ENDIF

C--   the forward F_mod(x) on copies of the carry
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          yuIce (i,j,bi,bj) = uIce(i,j,bi,bj)
          yvIce (i,j,bi,bj) = vIce(i,j,bi,bj)
          ysig1 (i,j,bi,bj) = seaice_sigma1(i,j,bi,bj)
          ysig2 (i,j,bi,bj) = seaice_sigma2(i,j,bi,bj)
          ysig12(i,j,bi,bj) = seaice_sigma12(i,j,bi,bj)
          ye11  (i,j,bi,bj) = e11(i,j,bi,bj)
          ye22  (i,j,bi,bj) = e22(i,j,bi,bj)
          ye12  (i,j,bi,bj) = e12(i,j,bi,bj)
          yzetaZ(i,j,bi,bj) = zetaZ(i,j,bi,bj)
          ydeltaC(i,j,bi,bj) = deltaC(i,j,bi,bj)
          ysdX  (i,j,bi,bj) = stressDivergenceX(i,j,bi,bj)
          ysdY  (i,j,bi,bj) = stressDivergenceY(i,j,bi,bj)
          yDWATN(i,j,bi,bj) = DWATN(i,j,bi,bj)
          yCbotC(i,j,bi,bj) = CbotC(i,j,bi,bj)
          yFX   (i,j,bi,bj) = FORCEX(i,j,bi,bj)
          yFY   (i,j,bi,bj) = FORCEY(i,j,bi,bj)
          ypress(i,j,bi,bj) = press(i,j,bi,bj)
          yzeta (i,j,bi,bj) = zeta(i,j,bi,bj)
          yeta  (i,j,bi,bj) = eta(i,j,bi,bj)
          dwFrAd(i,j,bi,bj) = 0. _d 0
          uFrAd (i,j,bi,bj) = 0. _d 0
          vFrAd (i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL SEAICE_EVP_FP_SUBCYCLE(
     U     yuIce, yvIce, ysig1, ysig2, ysig12,
     U     ye11, ye22, ye12, yzetaZ, ydeltaC,
     U     ysdX, ysdY, yDWATN, yCbotC,
     U     yFX, yFY, ypress, yzeta, yeta,
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

C--   reverse of the second tile loop of F_mod: the momentum update
C     (seaice_evp.F:1174-1251) and FORCEX/Y (:1086-1148)
      DO bj=myByHi(myThid),myByLo(myThid),-1
       DO bi=myBxHi(myThid),myBxLo(myThid),-1
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          bU(i,j) = evpBetaU(i,j,bi,bj)
          bV(i,j) = evpBetaV(i,j,bi,bj)
         ENDDO
        ENDDO
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
        DO j=sNy,1,-1
         DO i=sNx,1,-1
C     recompute the coefficients of the update
          betaFacU   = bU(i,j)*recip_deltaT
          betaFacV   = bV(i,j)*recip_deltaT
          tmp        = evpStarFac*recip_deltaT
          betaFacP1V = betaFacV + tmp
          betaFacP1U = betaFacU + tmp
          denomU = seaiceMassU(i,j,bi,bj)*betaFacP1U
     &         + 0.5 _d 0*( yDWATN(i,j,bi,bj) + yDWATN(i-1,j,bi,bj) )
     &         * COSWAT * areaW(i,j,bi,bj)
          denomV = seaiceMassV(i,j,bi,bj)*betaFacP1V
     &         + 0.5 _d 0*( yDWATN(i,j,bi,bj) + yDWATN(i,j-1,bi,bj) )
     &         * COSWAT * areaS(i,j,bi,bj)
#ifdef SEAICE_ALLOW_BOTTOMDRAG
          denomU = denomU + areaW(i,j,bi,bj)
     &         * 0.5 _d 0*( yCbotC(i,j,bi,bj) + yCbotC(i-1,j,bi,bj) )
          denomV = denomV + areaS(i,j,bi,bj)
     &         * 0.5 _d 0*( yCbotC(i,j,bi,bj) + yCbotC(i,j-1,bi,bj) )
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
          zeroU = denomU .EQ. 0. _d 0
          zeroV = denomV .EQ. 0. _d 0
          IF ( zeroU ) denomU = 1. _d 0
          IF ( zeroV ) denomV = 1. _d 0
C     reverse of vIce(i,j) = seaiceMaskV*( massV*betaFacV*vIce
C     + massV*recip_deltaT*evpStarFac*vIceNm1 + FORCEY + strDivY )/denomV
          numV = seaiceMassV(i,j,bi,bj)*betaFacV*vIce(i,j,bi,bj)
     &         + seaiceMassV(i,j,bi,bj)*recip_deltaT*evpStarFac
     &         * vIceNm1(i,j,bi,bj)
     &         + yFY(i,j,bi,bj) + ysdY(i,j,bi,bj)
          tb6 = seaiceMaskV(i,j,bi,bj)*advIce(i,j,bi,bj)/denomV
          advIce(i,j,bi,bj) = seaiceMassV(i,j,bi,bj)*betaFacV*tb6
          tb4 = recip_deltaT*evpStarFac*tb6
          adseaiceMassV(i,j,bi,bj) = adseaiceMassV(i,j,bi,bj)
     &         + vIce(i,j,bi,bj)*betaFacV*tb6
     &         + vIceNm1(i,j,bi,bj)*tb4
          adFORCEY(i,j,bi,bj) = adFORCEY(i,j,bi,bj) + tb6
          adstrDivY(i,j,bi,bj) = adstrDivY(i,j,bi,bj) + tb6
          addenomV = -( numV*tb6/denomV )
          advIceNm1(i,j,bi,bj) = advIceNm1(i,j,bi,bj)
     &         + seaiceMassV(i,j,bi,bj)*tb4
C     reverse of the uIce(i,j) update
          numU = seaiceMassU(i,j,bi,bj)*betaFacU*uIce(i,j,bi,bj)
     &         + seaiceMassU(i,j,bi,bj)*recip_deltaT*evpStarFac
     &         * uIceNm1(i,j,bi,bj)
     &         + yFX(i,j,bi,bj) + ysdX(i,j,bi,bj)
          tb6 = seaiceMaskU(i,j,bi,bj)*aduIce(i,j,bi,bj)/denomU
          aduIce(i,j,bi,bj) = seaiceMassU(i,j,bi,bj)*betaFacU*tb6
          tb4 = recip_deltaT*evpStarFac*tb6
          adseaiceMassU(i,j,bi,bj) = adseaiceMassU(i,j,bi,bj)
     &         + uIce(i,j,bi,bj)*betaFacU*tb6
     &         + uIceNm1(i,j,bi,bj)*tb4
          adFORCEX(i,j,bi,bj) = adFORCEX(i,j,bi,bj) + tb6
          adstrDivX(i,j,bi,bj) = adstrDivX(i,j,bi,bj) + tb6
          addenomU = -( numU*tb6/denomU )
          aduIceNm1(i,j,bi,bj) = aduIceNm1(i,j,bi,bj)
     &         + seaiceMassU(i,j,bi,bj)*tb4
C     reverse of IF ( denomU/V .EQ. 0. ) denomU/V = 1.
          IF ( zeroV ) addenomV = 0. _d 0
          IF ( zeroU ) addenomU = 0. _d 0
C     reverse of denomU, denomV (seaice_evp.F:1211-1224)
#ifdef SEAICE_ALLOW_BOTTOMDRAG
          tb6 = areaS(i,j,bi,bj)*0.5 _d 0*addenomV
          adCbotC(i,j  ,bi,bj) = adCbotC(i,j  ,bi,bj) + tb6
          adCbotC(i,j-1,bi,bj) = adCbotC(i,j-1,bi,bj) + tb6
          tb6 = areaW(i,j,bi,bj)*0.5 _d 0*addenomU
          adCbotC(i  ,j,bi,bj) = adCbotC(i  ,j,bi,bj) + tb6
          adCbotC(i-1,j,bi,bj) = adCbotC(i-1,j,bi,bj) + tb6
          adareaS(i,j,bi,bj) = adareaS(i,j,bi,bj)
     &         + ( yCbotC(i,j,bi,bj) + yCbotC(i,j-1,bi,bj) )
     &         * 0.5 _d 0*addenomV
          adareaW(i,j,bi,bj) = adareaW(i,j,bi,bj)
     &         + ( yCbotC(i,j,bi,bj) + yCbotC(i-1,j,bi,bj) )
     &         * 0.5 _d 0*addenomU
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
          tb6 = COSWAT*0.5 _d 0*addenomV
          adareaS(i,j,bi,bj) = adareaS(i,j,bi,bj)
     &         + ( yDWATN(i,j,bi,bj) + yDWATN(i,j-1,bi,bj) )*tb6
          adseaiceMassV(i,j,bi,bj) = adseaiceMassV(i,j,bi,bj)
     &         + betaFacP1V*addenomV
          tb5 = areaS(i,j,bi,bj)*tb6
          adDWATN(i,j  ,bi,bj) = adDWATN(i,j  ,bi,bj) + tb5
          adDWATN(i,j-1,bi,bj) = adDWATN(i,j-1,bi,bj) + tb5
          adseaiceMassU(i,j,bi,bj) = adseaiceMassU(i,j,bi,bj)
     &         + betaFacP1U*addenomU
          tb6 = COSWAT*0.5 _d 0*addenomU
          adareaW(i,j,bi,bj) = adareaW(i,j,bi,bj)
     &         + ( yDWATN(i,j,bi,bj) + yDWATN(i-1,j,bi,bj) )*tb6
          tb5 = areaW(i,j,bi,bj)*tb6
          adDWATN(i  ,j,bi,bj) = adDWATN(i  ,j,bi,bj) + tb5
          adDWATN(i-1,j,bi,bj) = adDWATN(i-1,j,bi,bj) + tb5
         ENDDO
        ENDDO
        DO j=sNy,1,-1
         DO i=sNx,1,-1
C     reverse of the Coriolis terms (seaice_evp.F:1134-1146)
          tb6 = -( _fCori(i,j  ,bi,bj)*0.5 _d 0**2
     &         * adFORCEY(i,j,bi,bj) )
          tb8 = -( _fCori(i,j-1,bi,bj)*0.5 _d 0**2
     &         * adFORCEY(i,j,bi,bj) )
          adseaiceMassC(i,j-1,bi,bj) = adseaiceMassC(i,j-1,bi,bj)
     &         + ( uIce(i,j-1,bi,bj) + uIce(i+1,j-1,bi,bj) )*tb8
          tb4 = seaiceMassC(i,j-1,bi,bj)*tb8
          aduIce(i  ,j-1,bi,bj) = aduIce(i  ,j-1,bi,bj) + tb4
          aduIce(i+1,j-1,bi,bj) = aduIce(i+1,j-1,bi,bj) + tb4
          adseaiceMassC(i,j,bi,bj) = adseaiceMassC(i,j,bi,bj)
     &         + ( uIce(i,j,bi,bj) + uIce(i+1,j,bi,bj) )*tb6
          tb5 = seaiceMassC(i,j,bi,bj)*tb6
          aduIce(i  ,j,bi,bj) = aduIce(i  ,j,bi,bj) + tb5
          aduIce(i+1,j,bi,bj) = aduIce(i+1,j,bi,bj) + tb5
          tb6 = _fCori(i  ,j,bi,bj)*0.5 _d 0**2*adFORCEX(i,j,bi,bj)
          tb8 = _fCori(i-1,j,bi,bj)*0.5 _d 0**2*adFORCEX(i,j,bi,bj)
          adseaiceMassC(i-1,j,bi,bj) = adseaiceMassC(i-1,j,bi,bj)
     &         + ( vIce(i-1,j,bi,bj) + vIce(i-1,j+1,bi,bj) )*tb8
          tb4 = seaiceMassC(i-1,j,bi,bj)*tb8
          advIce(i-1,j  ,bi,bj) = advIce(i-1,j  ,bi,bj) + tb4
          advIce(i-1,j+1,bi,bj) = advIce(i-1,j+1,bi,bj) + tb4
          adseaiceMassC(i,j,bi,bj) = adseaiceMassC(i,j,bi,bj)
     &         + ( vIce(i,j,bi,bj) + vIce(i,j+1,bi,bj) )*tb6
          tb5 = seaiceMassC(i,j,bi,bj)*tb6
          advIce(i,j  ,bi,bj) = advIce(i,j  ,bi,bj) + tb5
          advIce(i,j+1,bi,bj) = advIce(i,j+1,bi,bj) + tb5
C     recompute the masks of the drag term (:1103-1106)
          locMaskU = seaiceMassU(i,j,bi,bj)
          locMaskV = seaiceMassV(i,j,bi,bj)
          IF ( locMaskU .NE. 0. _d 0 ) locMaskU = 1. _d 0
          IF ( locMaskV .NE. 0. _d 0 ) locMaskV = 1. _d 0
          sgnFac = 0.5 _d 0*SIGN(SINWAT, _fCori(i,j,bi,bj))
C     reverse of FORCEY(i,j) = FORCEY0 + ( 0.5*(DWATN(i,j)+DWATN(i,j-1))
C     *COSWAT*vVel(i,j) + sgnFac*( DWATN(i,j)*0.5*aA
C     + DWATN(i,j-1)*0.5*aB )*locMaskV )*areaS
          aA = uVel(i  ,j  ,kSrf,bi,bj) - uIce(i  ,j  ,bi,bj)
     &       + uVel(i+1,j  ,kSrf,bi,bj) - uIce(i+1,j  ,bi,bj)
          aB = uVel(i  ,j-1,kSrf,bi,bj) - uIce(i  ,j-1,bi,bj)
     &       + uVel(i+1,j-1,kSrf,bi,bj) - uIce(i+1,j-1,bi,bj)
          adFORCEY0(i,j,bi,bj) = adFORCEY0(i,j,bi,bj)
     &         + adFORCEY(i,j,bi,bj)
          tb6 = areaS(i,j,bi,bj)*adFORCEY(i,j,bi,bj)
          adareaS(i,j,bi,bj) = adareaS(i,j,bi,bj)
     &         + ( COSWAT*0.5 _d 0*( ( yDWATN(i,j,bi,bj)
     &         + yDWATN(i,j-1,bi,bj) )*vVel(i,j,kSrf,bi,bj) )
     &         + sgnFac*( ( 0.5 _d 0*yDWATN(i,j,bi,bj)*aA
     &         + 0.5 _d 0*yDWATN(i,j-1,bi,bj)*aB )*locMaskV ) )
     &         * adFORCEY(i,j,bi,bj)
          adFORCEY(i,j,bi,bj) = 0. _d 0
          tb5 = COSWAT*0.5 _d 0*tb6
          tb3 = sgnFac*tb6
          tb0 = locMaskV*tb3
          adlocMaskV = ( 0.5 _d 0*yDWATN(i,j,bi,bj)*aA
     &         + 0.5 _d 0*yDWATN(i,j-1,bi,bj)*aB )*tb3
          adDWATN(i,j  ,bi,bj) = adDWATN(i,j  ,bi,bj)
     &         + aA*0.5 _d 0*tb0 + vVel(i,j,kSrf,bi,bj)*tb5
          adDWATN(i,j-1,bi,bj) = adDWATN(i,j-1,bi,bj)
     &         + aB*0.5 _d 0*tb0 + vVel(i,j,kSrf,bi,bj)*tb5
          advVel(i,j,kSrf,bi,bj) = advVel(i,j,kSrf,bi,bj)
     &         + ( yDWATN(i,j,bi,bj) + yDWATN(i,j-1,bi,bj) )*tb5
          tb7 = yDWATN(i,j,bi,bj)*0.5 _d 0*tb0
          tb1 = yDWATN(i,j-1,bi,bj)*0.5 _d 0*tb0
          aduVel(i  ,j-1,kSrf,bi,bj) = aduVel(i  ,j-1,kSrf,bi,bj) + tb1
          aduVel(i+1,j-1,kSrf,bi,bj) = aduVel(i+1,j-1,kSrf,bi,bj) + tb1
          aduIce(i  ,j-1,bi,bj) = aduIce(i  ,j-1,bi,bj) - tb1
          aduIce(i+1,j-1,bi,bj) = aduIce(i+1,j-1,bi,bj) - tb1
          aduVel(i  ,j  ,kSrf,bi,bj) = aduVel(i  ,j  ,kSrf,bi,bj) + tb7
          aduVel(i+1,j  ,kSrf,bi,bj) = aduVel(i+1,j  ,kSrf,bi,bj) + tb7
          aduIce(i  ,j  ,bi,bj) = aduIce(i  ,j  ,bi,bj) - tb7
          aduIce(i+1,j  ,bi,bj) = aduIce(i+1,j  ,bi,bj) - tb7
C     reverse of FORCEX(i,j) = FORCEX0 + ( 0.5*(DWATN(i,j)+DWATN(i-1,j))
C     *COSWAT*uVel(i,j) - sgnFac*( DWATN(i,j)*0.5*aA
C     + DWATN(i-1,j)*0.5*aB )*locMaskU )*areaW
          aA = vVel(i  ,j  ,kSrf,bi,bj) - vIce(i  ,j  ,bi,bj)
     &       + vVel(i  ,j+1,kSrf,bi,bj) - vIce(i  ,j+1,bi,bj)
          aB = vVel(i-1,j  ,kSrf,bi,bj) - vIce(i-1,j  ,bi,bj)
     &       + vVel(i-1,j+1,kSrf,bi,bj) - vIce(i-1,j+1,bi,bj)
          adFORCEX0(i,j,bi,bj) = adFORCEX0(i,j,bi,bj)
     &         + adFORCEX(i,j,bi,bj)
          tb2 = areaW(i,j,bi,bj)*adFORCEX(i,j,bi,bj)
          tb1 = COSWAT*0.5 _d 0*tb2
          aduVel(i,j,kSrf,bi,bj) = aduVel(i,j,kSrf,bi,bj)
     &         + ( yDWATN(i,j,bi,bj) + yDWATN(i-1,j,bi,bj) )*tb1
          tb0 = -( sgnFac*tb2 )
          tb3 = locMaskU*tb0
          adareaW(i,j,bi,bj) = adareaW(i,j,bi,bj)
     &         + ( COSWAT*0.5 _d 0*( ( yDWATN(i,j,bi,bj)
     &         + yDWATN(i-1,j,bi,bj) )*uVel(i,j,kSrf,bi,bj) )
     &         - sgnFac*( ( 0.5 _d 0*yDWATN(i,j,bi,bj)*aA
     &         + 0.5 _d 0*yDWATN(i-1,j,bi,bj)*aB )*locMaskU ) )
     &         * adFORCEX(i,j,bi,bj)
          adFORCEX(i,j,bi,bj) = 0. _d 0
          adlocMaskU = ( 0.5 _d 0*yDWATN(i,j,bi,bj)*aA
     &         + 0.5 _d 0*yDWATN(i-1,j,bi,bj)*aB )*tb0
          adDWATN(i  ,j,bi,bj) = adDWATN(i  ,j,bi,bj)
     &         + aA*0.5 _d 0*tb3 + uVel(i,j,kSrf,bi,bj)*tb1
          adDWATN(i-1,j,bi,bj) = adDWATN(i-1,j,bi,bj)
     &         + aB*0.5 _d 0*tb3 + uVel(i,j,kSrf,bi,bj)*tb1
          tb5 = yDWATN(i-1,j,bi,bj)*0.5 _d 0*tb3
          advVel(i-1,j  ,kSrf,bi,bj) = advVel(i-1,j  ,kSrf,bi,bj) + tb5
          advVel(i-1,j+1,kSrf,bi,bj) = advVel(i-1,j+1,kSrf,bi,bj) + tb5
          advIce(i-1,j  ,bi,bj) = advIce(i-1,j  ,bi,bj) - tb5
          advIce(i-1,j+1,bi,bj) = advIce(i-1,j+1,bi,bj) - tb5
          tb4 = yDWATN(i,j,bi,bj)*0.5 _d 0*tb3
          advVel(i,j  ,kSrf,bi,bj) = advVel(i,j  ,kSrf,bi,bj) + tb4
          advVel(i,j+1,kSrf,bi,bj) = advVel(i,j+1,kSrf,bi,bj) + tb4
          advIce(i,j  ,bi,bj) = advIce(i,j  ,bi,bj) - tb4
          advIce(i,j+1,bi,bj) = advIce(i,j+1,bi,bj) - tb4
C     reverse of locMaskU/V = seaiceMassU/V, 1 where it is not 0
          IF ( seaiceMassV(i,j,bi,bj) .EQ. 0. _d 0 )
     &         adseaiceMassV(i,j,bi,bj) = adseaiceMassV(i,j,bi,bj)
     &         + adlocMaskV
          IF ( seaiceMassU(i,j,bi,bj) .EQ. 0. _d 0 )
     &         adseaiceMassU(i,j,bi,bj) = adseaiceMassU(i,j,bi,bj)
     &         + adlocMaskU
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   reverse of SEAICE_BOTTOMDRAG_COEFFS: with SEAICEbasalDragK2 = 0
C     (the only case allowed, above) it writes nothing, so CbotC passes
C     through and HEFF, AREA get no adjoint

C--   reverse of the frozen drag on thin ice and of the drag of the
C     iterate (seaice_evp_fp_subcycle.F; the adjoint of the frozen ice
C     velocity, uFrAd/vFrAd, is discarded)
      DO bj=myByHi(myThid),myByLo(myThid),-1
       DO bi=myBxHi(myThid),myBxLo(myThid),-1
        DO j=sNy+OLy-1,1-OLy,-1
         DO i=sNx+OLx-1,1-OLx,-1
          IF ( AREA(i,j,bi,bj) .LT. fpThinArea .OR.
     &         seaiceMassC(i,j,bi,bj) .LT. fpThinMass ) THEN
           dwFrAd(i,j,bi,bj)  = dwFrAd(i,j,bi,bj) + adDWATN(i,j,bi,bj)
           adDWATN(i,j,bi,bj) = 0. _d 0
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL SEAICE_FP_OCEANDRAG_COEFFS_AD(
     I     uIceFr, vIceFr, HEFFM, uVel, vVel,
     U     dwFrAd, uFrAd, vFrAd, aduVel, advVel,
     I     myThid )
      CALL SEAICE_FP_OCEANDRAG_COEFFS_AD(
     I     uIce, vIce, HEFFM, uVel, vVel,
     U     adDWATN, aduIce, advIce, aduVel, advVel,
     I     myThid )

C--   reverse of the first tile loop of F_mod (stresses,
C     seaice_evp.F:635-989)
      DO bj=myByHi(myThid),myByLo(myThid),-1
       DO bi=myBxHi(myThid),myBxLo(myThid),-1
C     recompute the tile-local fields of the forward
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          ep(i,j) = ye11(i,j,bi,bj) + ye22(i,j,bi,bj)
          em(i,j) = ye11(i,j,bi,bj) - ye22(i,j,bi,bj)
          zetaC(i,j) = 0. _d 0
          aC (i,j) = evpAlphaC(i,j,bi,bj)
          aZ (i,j) = evpAlphaZ(i,j,bi,bj)
          d1 (i,j) = denom1(i,j,bi,bj)
          d2C(i,j) = denom2(i,j,bi,bj)
          d2Z(i,j) = denom2(i,j,bi,bj)
          adep     (i,j) = 0. _d 0
          adem     (i,j) = 0. _d 0
          ade12Csq (i,j) = 0. _d 0
          adzetaC  (i,j) = 0. _d 0
          adpressC (i,j) = 0. _d 0
          addiv    (i,j) = 0. _d 0
          adtension(i,j) = 0. _d 0
          adshear  (i,j) = 0. _d 0
          adsig11  (i,j) = 0. _d 0
          adsig22  (i,j) = 0. _d 0
         ENDDO
        ENDDO
        DO j=0,sNy+1
         DO i=0,sNx+1
          thetaC = fpTheta
          IF ( deltaFr(i,j,bi,bj) .LT. fpRigidDelta )
     &         thetaC = fpThetaRigid
          deltaTh = deltaFr(i,j,bi,bj)
     &         + thetaC*( ydeltaC(i,j,bi,bj) - deltaFr(i,j,bi,bj) )
          deltaCreg  = MAX(deltaTh,SEAICE_deltaMin)
          zetaC(i,j) = HALF*( press0(i,j,bi,bj)
     &         * ( 1. _d 0 + tensileStrFac(i,j,bi,bj) )
     &         )/deltaCreg
         ENDDO
        ENDDO
        IF ( fpRaise ) THEN
         DO j=0,sNy+1
          DO i=0,sNx+1
           aC(i,j) = SEAICE_evpAlpha*raiseC(i,j,bi,bj)
          ENDDO
         ENDDO
         DO j=0,sNy
          DO i=0,sNx
           d1(i,j)  = 1. _d 0 /  aC(i,j)
           d2C(i,j) = d1(i,j)
          ENDDO
         ENDDO
         DO j=1,sNy+1
          DO i=1,sNx+1
           aZ(i,j) = SEAICE_evpAlpha * MAX(
     &          MAX(raiseC(i,  j,bi,bj),raiseC(i-1,j-1,bi,bj)),
     &          MAX(raiseC(i-1,j,bi,bj),raiseC(i,  j-1,bi,bj)) )
           d2Z(i,j) = 1. _d 0 /  aZ(i,j)
          ENDDO
         ENDDO
        ENDIF
C     reverse of the divergence of the stress tensor (:956-986)
        DO j=sNy,1,-1
         DO i=sNx,1,-1
          tb2 = recip_rAs(i,j,bi,bj)*adstrDivY(i,j,bi,bj)
          tb1 = k1AtV(i,j,bi,bj)*HALF*adstrDivY(i,j,bi,bj)
          tb0 = -( k2AtV(i,j,bi,bj)*HALF*adstrDivY(i,j,bi,bj) )
          adstrDivY(i,j,bi,bj) = 0. _d 0
          adsig11(i,j-1) = adsig11(i,j-1) + tb0
          adsigma12(i,j,bi,bj) = adsigma12(i,j,bi,bj) + tb1
          adsigma12(i+1,j,bi,bj) = adsigma12(i+1,j,bi,bj) + tb1
     &         + _dyU(i+1,j,bi,bj)*tb2
          adsig22(i,j)   = adsig22(i,j)   + _dxF(i,j  ,bi,bj)*tb2
          adsig22(i,j-1) = adsig22(i,j-1) - _dxF(i,j-1,bi,bj)*tb2
          tb1 = recip_rAw(i,j,bi,bj)*adstrDivX(i,j,bi,bj)
          adsig11(i,j) = adsig11(i,j) + tb0 + _dyF(i,j,bi,bj)*tb1
          tb0 = k2AtU(i,j,bi,bj)*HALF*adstrDivX(i,j,bi,bj)
          adsigma12(i,j,bi,bj) = adsigma12(i,j,bi,bj) + tb0
     &         - _dyU(i,j,bi,bj)*tb2 - _dxV(i,j,bi,bj)*tb1
          tb2 = -( k1AtU(i,j,bi,bj)*HALF*adstrDivX(i,j,bi,bj) )
          adstrDivX(i,j,bi,bj) = 0. _d 0
          adsig22(i  ,j) = adsig22(i  ,j) + tb2
          adsig22(i-1,j) = adsig22(i-1,j) + tb2
          adsigma12(i,j+1,bi,bj) = adsigma12(i,j+1,bi,bj) + tb0
     &         + _dxV(i,j+1,bi,bj)*tb1
          adsig11(i-1,j) = adsig11(i-1,j) - _dyF(i-1,j,bi,bj)*tb1
         ENDDO
        ENDDO
C     reverse of the sigma12 update on Z points (:945-948)
        DO j=sNy+1,1,-1
         DO i=sNx+1,1,-1
          tb1 = d2Z(i,j)*adsigma12(i,j,bi,bj)
          adsigma12(i,j,bi,bj) = ( aZ(i,j) - evpRevFac )*tb1
          adshear(i,j) = adshear(i,j) + recip_evpRevFac*tb1
         ENDDO
        ENDDO
C     reverse of sig11, sig22 and the sigma1, sigma2 updates on C
C     points (:896-924)
        DO j=sNy,0,-1
         DO i=sNx,0,-1
          adsigma1(i,j,bi,bj) = adsigma1(i,j,bi,bj)
     &         + 0.5 _d 0*adsig22(i,j) + 0.5 _d 0*adsig11(i,j)
          adsigma2(i,j,bi,bj) = adsigma2(i,j,bi,bj)
     &         + 0.5 _d 0*adsig11(i,j) - 0.5 _d 0*adsig22(i,j)
          adsig22(i,j) = 0. _d 0
          adsig11(i,j) = 0. _d 0
          tb1 = d2C(i,j)*HEFFM(i,j,bi,bj)*adsigma2(i,j,bi,bj)
          adsigma2(i,j,bi,bj) = ( aC(i,j) - evpRevFac )*tb1
          adtension(i,j) = adtension(i,j) + recip_evpRevFac*tb1
          tb1 = d1(i,j)*HEFFM(i,j,bi,bj)*adsigma1(i,j,bi,bj)
          adsigma1(i,j,bi,bj) = ( aC(i,j) - evpRevFac )*tb1
          addiv(i,j) = addiv(i,j) + tb1
         ENDDO
        ENDDO
C     reverse of seaice_shear (:870-871)
        DO j=sNy+1,1,-1
         DO i=sNx+1,1,-1
          adzetaZ(i,j,bi,bj) = adzetaZ(i,j,bi,bj)
     &         + ye12(i,j,bi,bj)*2. _d 0*adshear(i,j)
          ade12(i,j,bi,bj) = ade12(i,j,bi,bj)
     &         + yzetaZ(i,j,bi,bj)*2. _d 0*adshear(i,j)
          adshear(i,j) = 0. _d 0
         ENDDO
        ENDDO
C     reverse of seaice_div, seaice_tension (:861-865)
        DO j=sNy,0,-1
         DO i=sNx,0,-1
          tb0 = HEFFM(i,j,bi,bj)*addiv(i,j)
          tb1 = HEFFM(i,j,bi,bj)*2. _d 0*adtension(i,j)
          adtension(i,j) = 0. _d 0
          adzetaC(i,j) = adzetaC(i,j) + em(i,j)*tb1
     &         + ep(i,j)*2. _d 0*tb0
          adem(i,j) = adem(i,j) + zetaC(i,j)*tb1
          addiv(i,j) = 0. _d 0
          adep(i,j) = adep(i,j) + zetaC(i,j)*2. _d 0*tb0
          adpressC(i,j) = adpressC(i,j) - tb0
         ENDDO
        ENDDO
#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
C     reverse of the diagnostics copies (:792-803)
         DO j=sNy,1,-1
          DO i=sNx,1,-1
           adzetaC(i,j) = adzetaC(i,j)
     &          + recip_ecc2*adeta(i,j,bi,bj) + adzeta(i,j,bi,bj)
           adeta (i,j,bi,bj) = 0. _d 0
           adzeta(i,j,bi,bj) = 0. _d 0
           adpressC(i,j) = adpressC(i,j) + adpress(i,j,bi,bj)
           adpress(i,j,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */
C     reverse of zetaZ (:743-752)
        DO j=sNy+1,1,-1
         DO i=sNx+1,1,-1
          sumNorm = HEFFM(i,j,  bi,bj)+HEFFM(i-1,j,  bi,bj)
     &            + HEFFM(i,j-1,bi,bj)+HEFFM(i-1,j-1,bi,bj)
          IF ( sumNorm.GT.0. _d 0 ) sumNorm = 1. _d 0 / sumNorm
          tb1 = sumNorm*adzetaZ(i,j,bi,bj)
          adzetaZ(i,j,bi,bj) = 0. _d 0
          adzetaC(i  ,j  ) = adzetaC(i  ,j  ) + tb1
          adzetaC(i-1,j-1) = adzetaC(i-1,j-1) + tb1
          adzetaC(i-1,j  ) = adzetaC(i-1,j  ) + tb1
          adzetaC(i  ,j-1) = adzetaC(i  ,j-1) + tb1
         ENDDO
        ENDDO
C     reverse of pressC, zetaP, zetaC, Delta_theta and deltaC
C     (:698-791 with the fixed-point arms of F_mod)
        DO j=sNy+1,0,-1
         DO i=sNx+1,0,-1
C     recompute the point's forward values
          IF ( SEAICEetaZmethod .EQ. 0 ) THEN
           tmp = 0.25 *
     &          ( ye12(i,j  ,bi,bj) + ye12(i+1,j  ,bi,bj)
     &          + ye12(i,j+1,bi,bj) + ye12(i+1,j+1,bi,bj) )
           deltaSq = ep(i,j)**2 + recip_ecc2 * em(i,j)**2
     &          + recip_ecc2 * 4. _d 0 * ( tmp*tmp )
          ELSE
           deltaSq = ep(i,j)**2 + recip_ecc2 * em(i,j)**2
     &          + recip_ecc2 * 4. _d 0 * ( 0.25 _d 0
     &          * recip_rA(i,j,bi,bj) *
     &          ( rAz(i  ,j  ,bi,bj)*ye12(i  ,j  ,bi,bj)**2
     &          + rAz(i+1,j  ,bi,bj)*ye12(i+1,j  ,bi,bj)**2
     &          + rAz(i  ,j+1,bi,bj)*ye12(i  ,j+1,bi,bj)**2
     &          + rAz(i+1,j+1,bi,bj)*ye12(i+1,j+1,bi,bj)**2 ) )
          ENDIF
          thetaC = fpTheta
          IF ( deltaFr(i,j,bi,bj) .LT. fpRigidDelta )
     &         thetaC = fpThetaRigid
          deltaTh = deltaFr(i,j,bi,bj)
     &         + thetaC*( ydeltaC(i,j,bi,bj) - deltaFr(i,j,bi,bj) )
          deltaCreg = MAX(deltaTh,SEAICE_deltaMin)
          deltaP = deltaFr(i,j,bi,bj)
          IF ( fpPrLive ) deltaP = ydeltaC(i,j,bi,bj)
          maxP  = MAX(deltaP,SEAICE_deltaMin)
          zetaP = HALF*( press0(i,j,bi,bj)
     &         * ( 1. _d 0 + tensileStrFac(i,j,bi,bj) ) )/maxP
C     reverse of pressC = ( press0*(1-SEAICEpressReplFac)
C     + TWO*zetaP*deltaP*SEAICEpressReplFac/(1+tensileStrFac) )
C     * (1-tensileStrFac)
          tb2 = 1. _d 0 + tensileStrFac(i,j,bi,bj)
          tb3 = zetaP*deltaP/tb2
          tb0 = ( 1. _d 0 - tensileStrFac(i,j,bi,bj) )*adpressC(i,j)
          tb1 = TWO*SEAICEpressReplFac*tb0/tb2
          adzetaP  = deltaP*tb1
          addeltaP = zetaP*tb1
          tb4 = press0(i,j,bi,bj)/maxP
          adtensileStrFac(i,j,bi,bj) = adtensileStrFac(i,j,bi,bj)
     &         + tb4*HALF*adzetaP
     &         - ( ( 1. _d 0 - SEAICEpressReplFac )*press0(i,j,bi,bj)
     &         + TWO*SEAICEpressReplFac*tb3 )*adpressC(i,j)
     &         - tb3*tb1
          adpressC(i,j) = 0. _d 0
          tb5 = ( tensileStrFac(i,j,bi,bj) + 1. _d 0 )*HALF*adzetaP
     &         / maxP
          adpress0(i,j,bi,bj) = adpress0(i,j,bi,bj)
     &         + ( 1. _d 0 - SEAICEpressReplFac )*tb0 + tb5
          admaxP = -( tb4*tb5 )
C     reverse of zetaP = HALF*(..)/MAX(deltaP,SEAICE_deltaMin); deltaP
C     = deltaFr is passive, deltaC with fpPrLive
          IF ( .NOT.( deltaP .LT. SEAICE_deltaMin ) )
     &         addeltaP = addeltaP + admaxP
          IF ( fpPrLive ) addeltaC(i,j,bi,bj) = addeltaC(i,j,bi,bj)
     &         + addeltaP
C     reverse of zetaC = HALF*(..)/deltaCreg
          tb4 = press0(i,j,bi,bj)/deltaCreg
          tb5 = ( tensileStrFac(i,j,bi,bj) + 1. _d 0 )*HALF*adzetaC(i,j)
     &         / deltaCreg
          adtensileStrFac(i,j,bi,bj) = adtensileStrFac(i,j,bi,bj)
     &         + tb4*HALF*adzetaC(i,j)
          adzetaC(i,j) = 0. _d 0
          adpress0(i,j,bi,bj) = adpress0(i,j,bi,bj) + tb5
          addeltaCreg = -( tb4*tb5 )
C     reverse of deltaCreg = MAX(deltaTh,SEAICE_deltaMin) and
C     Delta_theta = deltaFr + thetaC*(deltaC - deltaFr)
          addeltaTh = addeltaCreg
          IF ( deltaTh .LT. SEAICE_deltaMin ) addeltaTh = 0. _d 0
          addeltaC(i,j,bi,bj) = addeltaC(i,j,bi,bj) + thetaC*addeltaTh
C     reverse of deltaC = SQRT(deltaSq) where deltaSq > 0, 0 elsewhere
          addeltaSq = 0. _d 0
          IF ( deltaSq .GT. 0. _d 0 )
     &         addeltaSq = addeltaC(i,j,bi,bj)/( 2. _d 0*SQRT(deltaSq) )
          addeltaC(i,j,bi,bj) = 0. _d 0
          adep(i,j) = adep(i,j) + 2. _d 0*ep(i,j)*addeltaSq
          adem(i,j) = adem(i,j) + 2. _d 0*em(i,j)*recip_ecc2*addeltaSq
          ade12Csq(i,j) = ade12Csq(i,j)
     &         + recip_ecc2*4. _d 0*addeltaSq
         ENDDO
        ENDDO
C     reverse of e12Csq (:675-697)
        IF ( SEAICEetaZmethod .EQ. 0 ) THEN
         DO j=sNy+OLy-1,1-OLy+1,-1
          DO i=sNx+OLx-1,1-OLx+1,-1
           tmp = 0.25 *
     &          ( ye12(i,j  ,bi,bj) + ye12(i+1,j  ,bi,bj)
     &          + ye12(i,j+1,bi,bj) + ye12(i+1,j+1,bi,bj) )
           adtmp = 2. _d 0*tmp*ade12Csq(i,j)
           ade12Csq(i,j) = 0. _d 0
           tb0 = 0.25*adtmp
           ade12(i  ,j  ,bi,bj) = ade12(i  ,j  ,bi,bj) + tb0
           ade12(i+1,j  ,bi,bj) = ade12(i+1,j  ,bi,bj) + tb0
           ade12(i  ,j+1,bi,bj) = ade12(i  ,j+1,bi,bj) + tb0
           ade12(i+1,j+1,bi,bj) = ade12(i+1,j+1,bi,bj) + tb0
          ENDDO
         ENDDO
        ELSEIF ( SEAICEetaZmethod .EQ. 3 ) THEN
         DO j=sNy+OLy-1,1-OLy+1,-1
          DO i=sNx+OLx-1,1-OLx+1,-1
           tb0 = recip_rA(i,j,bi,bj)*0.25 _d 0*ade12Csq(i,j)
           ade12Csq(i,j) = 0. _d 0
           ade12(i  ,j  ,bi,bj) = ade12(i  ,j  ,bi,bj)
     &          + 2. _d 0*ye12(i  ,j  ,bi,bj)*rAz(i  ,j  ,bi,bj)*tb0
           ade12(i+1,j  ,bi,bj) = ade12(i+1,j  ,bi,bj)
     &          + 2. _d 0*ye12(i+1,j  ,bi,bj)*rAz(i+1,j  ,bi,bj)*tb0
           ade12(i  ,j+1,bi,bj) = ade12(i  ,j+1,bi,bj)
     &          + 2. _d 0*ye12(i  ,j+1,bi,bj)*rAz(i  ,j+1,bi,bj)*tb0
           ade12(i+1,j+1,bi,bj) = ade12(i+1,j+1,bi,bj)
     &          + 2. _d 0*ye12(i+1,j+1,bi,bj)*rAz(i+1,j+1,bi,bj)*tb0
          ENDDO
         ENDDO
        ENDIF
C     reverse of ep, em (:667-672)
        DO j=sNy+OLy,1-OLy,-1
         DO i=sNx+OLx,1-OLx,-1
          ade11(i,j,bi,bj) = ade11(i,j,bi,bj) + adem(i,j) + adep(i,j)
          ade22(i,j,bi,bj) = ade22(i,j,bi,bj) + adep(i,j) - adem(i,j)
          adem(i,j) = 0. _d 0
          adep(i,j) = 0. _d 0
         ENDDO
        ENDDO
C     zetaZ and deltaC are set to zero on every point first (:663-664)
        DO j=sNy+OLy,1-OLy,-1
         DO i=sNx+OLx,1-OLx,-1
          addeltaC(i,j,bi,bj) = 0. _d 0
          adzetaZ (i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   reverse of SEAICE_CALC_STRAINRATES (:622-627)
      CALL SEAICE_FP_STRAINRATES_AD(
     U     aduIce, advIce, ade11, ade22, ade12,
     I     myThid )

#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_FP_STRAINRATES_AD
C     !INTERFACE:
      SUBROUTINE SEAICE_FP_STRAINRATES_AD(
     U     aduFld, advFld, ade11Loc, ade22Loc, ade12Loc,
     I     myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_FP_STRAINRATES_AD
C     | o adjoint of SEAICE_CALC_STRAINRATES
C     |   (seaice_calc_strainrates.F; linear in the
C     |   velocity) without the SEAICE_2ndOrderBC arm (:158-189,
C     |   refused by SEAICE_CHECK with SEAICEuseEVPfixedPointAdj)
C     |   and without OBCS (OBCS_UVICE_OLD): ade11/22/12Loc are
C     |   used and set to zero on the points the forward writes,
C     |   aduFld, advFld are incremented
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
      _RL aduFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL advFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade11Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade22Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade12Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
CEOP

#ifdef SEAICE_CGRID
C     === local variables ===
      INTEGER i, j, bi, bj
      _RS hFacU, hFacV, noSlipFac
      _RL addudx, addvdy, addudy, addvdx, aduave, advave, tb0

      noSlipFac = 0. _d 0
      IF ( SEAICE_no_slip ) noSlipFac = 1. _d 0

      DO bj=myByHi(myThid),myByLo(myThid),-1
       DO bi=myBxHi(myThid),myBxLo(myThid),-1
C     reverse of the strain rates at Z points (:132-150)
        DO j=sNy+OLy,1-OLy+1,-1
         DO i=sNx+OLx,1-OLx+1,-1
          hFacU = SIMaskU(i,j,bi,bj) - SIMaskU(i,j-1,bi,bj)
          hFacV = SIMaskV(i,j,bi,bj) - SIMaskV(i-1,j,bi,bj)
          tb0 = HEFFM(i,j,bi,bj)*0.5 _d 0*HEFFM(i-1,j-1,bi,bj)
     &        * HEFFM(i-1,j,bi,bj)*HEFFM(i,j-1,bi,bj)
     &        * ade12Loc(i,j,bi,bj)
          aduave = hFacU*2.0 _d 0*_recip_dyU(i,j,bi,bj)*noSlipFac
     &         * ade12Loc(i,j,bi,bj) - k2AtZ(i,j,bi,bj)*tb0
          advave = hFacV*2.0 _d 0*_recip_dxV(i,j,bi,bj)*noSlipFac
     &         * ade12Loc(i,j,bi,bj) - k1AtZ(i,j,bi,bj)*tb0
          ade12Loc(i,j,bi,bj) = 0. _d 0
          addudy = tb0
          addvdx = tb0
C     dvdx = (vFld(i,j)-vFld(i-1,j))*recip_dxV, vave = 0.5*(..+..)
          tb0 = _recip_dxV(i,j,bi,bj)*addvdx
          advFld(i  ,j,bi,bj) = advFld(i  ,j,bi,bj)
     &         + 0.5 _d 0*advave + tb0
          advFld(i-1,j,bi,bj) = advFld(i-1,j,bi,bj)
     &         + 0.5 _d 0*advave - tb0
C     dudy = (uFld(i,j)-uFld(i,j-1))*recip_dyU, uave = 0.5*(..+..)
          tb0 = _recip_dyU(i,j,bi,bj)*addudy
          aduFld(i,j  ,bi,bj) = aduFld(i,j  ,bi,bj)
     &         + 0.5 _d 0*aduave + tb0
          aduFld(i,j-1,bi,bj) = aduFld(i,j-1,bi,bj)
     &         + 0.5 _d 0*aduave - tb0
         ENDDO
        ENDDO
C     reverse of the strain rates at C points (:88-108)
        DO j=sNy+OLy-1,1-OLy,-1
         DO i=sNx+OLx-1,1-OLx,-1
          addudx = ade11Loc(i,j,bi,bj)
          advave = k2AtC(i,j,bi,bj)*ade11Loc(i,j,bi,bj)
          addvdy = ade22Loc(i,j,bi,bj)
          aduave = k1AtC(i,j,bi,bj)*ade22Loc(i,j,bi,bj)
          ade11Loc(i,j,bi,bj) = 0. _d 0
          ade22Loc(i,j,bi,bj) = 0. _d 0
C     dvdy = recip_dyF*(vFld(i,j+1)-vFld(i,j)), vave = 0.5*(..+..)
          tb0 = _recip_dyF(i,j,bi,bj)*addvdy
          advFld(i,j+1,bi,bj) = advFld(i,j+1,bi,bj)
     &         + 0.5 _d 0*advave + tb0
          advFld(i,j  ,bi,bj) = advFld(i,j  ,bi,bj)
     &         + 0.5 _d 0*advave - tb0
C     dudx = recip_dxF*(uFld(i+1,j)-uFld(i,j)), uave = 0.5*(..+..)
          tb0 = _recip_dxF(i,j,bi,bj)*addudx
          aduFld(i+1,j,bi,bj) = aduFld(i+1,j,bi,bj)
     &         + 0.5 _d 0*aduave + tb0
          aduFld(i  ,j,bi,bj) = aduFld(i  ,j,bi,bj)
     &         + 0.5 _d 0*aduave - tb0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* SEAICE_CGRID */

      RETURN
      END

#endif /* SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ */
