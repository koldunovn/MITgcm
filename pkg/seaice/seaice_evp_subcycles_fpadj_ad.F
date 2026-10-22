#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#else
# define OBCS_UVICE_OLD
#endif
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

C     The fixed-point adjoint of the mEVP subcycle loop
C     (SEAICEuseEVPfixedPointAdj; see doc/phys_pkgs/seaice.rst):
C       SEAICE_EVP_SUBCYCLES_FPADJ_AD :: the hand-written adjoint of
C                                        SEAICE_EVP_SUBCYCLES_FPADJ
C       SEAICE_EVP_FPADJ_CORE         :: its work
C       SEAICE_EVP_FPADJ_RERUN        :: reruns subcycles from x_0
C       SEAICE_FPADJ_COPY, _FILL, _ADDW, _MAXABS, _CHKSUM :: helpers

#ifdef SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ
CBOP
C     !ROUTINE: SEAICE_EVP_SUBCYCLES_FPADJ_AD
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ_AD(
     I     uIce, vIce,
     I     seaice_sigma1, seaice_sigma2, seaice_sigma12,
     I     e11, e22, e12, zetaZ, deltaC,
     I     stressDivergenceX, stressDivergenceY,
     I     DWATN,
     I     FORCEX, FORCEY, press, zeta, eta,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA, areaW, areaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     U     aduIce, advIce,
     U     adsigma1, adsigma2, adsigma12,
     U     ade11, ade22, ade12, adzetaZ, addeltaC,
     U     adstrDivX, adstrDivY,
     U     adDWATN,
     U     adFORCEX, adFORCEY, adpress, adzeta, adeta,
     U     adpress0,
     U     adseaiceMassC, adseaiceMassU, adseaiceMassV,
     U     adFORCEX0, adFORCEY0, aduIceNm1, advIceNm1,
     U     aduVel, advVel, adHEFF, adAREA,
     U     adareaW, adareaS,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ_AD
C     | o The fixed-point adjoint of the mEVP subcycle loop
C     |   SEAICE_EVP_SUBCYCLES: the hand-written adjoint of
C     |   SEAICE_EVP_SUBCYCLES_FPADJ. The AD tool's adjoint of
C     |   that routine calls it (Tapenade: the _B routine in
C     |   seaice_evp_subcycles_fpadj_tap.F; TAF: the ADNAME
C     |   routine in seaice_evp_subcycles_fpadj_mad.F).
C     *==========================================================*
C     | The forward arguments hold their values at the entry of
C     | the forward call (x_0 and the inputs p). The adjoint of
C     | the carry (aduIce .. adeta) holds the adjoint of the
C     | loop's output on entry and of its input on exit; the
C     | adjoints of the inputs (adpress0 .. adareaS) are
C     | incremented. tensileStrFac, denom1/2, evpAlphaC/Z and
C     | evpBetaU/V are constant fields here (no adjoint). With
C     | SEAICE_ALLOW_BOTTOMDRAG, CbotC is zero (SEAICE_CHECK
C     | refuses SEAICEbasalDragK2 > 0), so it is not an argument.
C     | The settings are the namelist's (SEAICE_fpOmega,
C     | SEAICE_fpExchAdj, SEAICE_fpExactReplay); the work is done
C     | in SEAICE_EVP_FPADJ_CORE.
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
C     uIce .. myThid :: the arguments of SEAICE_EVP_SUBCYCLES_FPADJ
C                       (values at the entry of the forward call)
C     aduIce .. adeta :: adjoint of the carry (overwritten)
C     adpress0 .. adareaS :: adjoint of the inputs (incremented)
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
      INTEGER kSrf
      INTEGER nEVPstep
      _RL COSWAT
      _RS SINWAT
      _RL ecc2, recip_ecc2, deltaMinSq
      _RL evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT
      _RL EVPcFac
      LOGICAL useAdaptiveEVP
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
      _RL adFORCEX      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEY      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adpress       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adzeta        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adeta         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adpress0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
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
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP )
#if (defined SEAICE_ALLOW_SIDEDRAG || defined SEAICE_ITD \
  || defined SEAICE_ALLOW_CLIPZETA)
      STOP 'S/R SEAICE_EVP_SUBCYCLES_FPADJ_AD: CPP options not covered'
#else
C     === Local variables ===
C     nC             :: number of carry fields (the loop's state:
C                       1 uIce, 2 vIce, 3-5 seaice_sigma1/2/12,
C                       6-8 e11/e22/e12, 9 zetaZ, 10 deltaC,
C                       11-12 stressDivergenceX/Y, 13 DWATN, 14 CbotC,
C                       15-16 FORCEX/Y, 17-19 press, zeta, eta)
C     xIn            :: the carry at the entry of the forward call x_0
C     lam            :: the adjoint of the carry
      INTEGER nC
      PARAMETER ( nC = 19 )
      _RL xIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
      _RL lam (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
C     tmpCb          :: the adjoint of CbotC (not used)
C     adTsf          :: the adjoint of tensileStrFac (not used)
      _RL tmpCb(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adTsf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      SAVE xIn, lam, tmpCb, adTsf
#define XIF(n) xIn(1-OLx,1-OLy,1,1,n)
#define LMF(n) lam(1-OLx,1-OLy,1,1,n)

C--   pack the carry and its adjoint
      CALL SEAICE_FPADJ_COPY(
     &     uIce,
     &     XIF(1),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     vIce,
     &     XIF(2),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     seaice_sigma1,
     &     XIF(3),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     seaice_sigma2,
     &     XIF(4),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     seaice_sigma12,
     &     XIF(5),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     e11,
     &     XIF(6),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     e22,
     &     XIF(7),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     e12,
     &     XIF(8),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     zetaZ,
     &     XIF(9),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     deltaC,
     &     XIF(10),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     stressDivergenceX,
     &     XIF(11),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     stressDivergenceY,
     &     XIF(12),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     DWATN,
     &     XIF(13),
     &     myThid )
C     CbotC is zero (see above)
      CALL SEAICE_FPADJ_FILL(
     &     0. _d 0,
     &     XIF(14),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     FORCEX,
     &     XIF(15),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     FORCEY,
     &     XIF(16),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     press,
     &     XIF(17),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     zeta,
     &     XIF(18),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     eta,
     &     XIF(19),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     aduIce,
     &     LMF(1),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     advIce,
     &     LMF(2),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adsigma1,
     &     LMF(3),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adsigma2,
     &     LMF(4),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adsigma12,
     &     LMF(5),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     ade11,
     &     LMF(6),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     ade22,
     &     LMF(7),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     ade12,
     &     LMF(8),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adzetaZ,
     &     LMF(9),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     addeltaC,
     &     LMF(10),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adstrDivX,
     &     LMF(11),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adstrDivY,
     &     LMF(12),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adDWATN,
     &     LMF(13),
     &     myThid )
      CALL SEAICE_FPADJ_FILL(
     &     0. _d 0,
     &     LMF(14),
     &     myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, adTsf, myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adFORCEX,
     &     LMF(15),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adFORCEY,
     &     LMF(16),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adpress,
     &     LMF(17),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adzeta,
     &     LMF(18),
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     adeta,
     &     LMF(19),
     &     myThid )

      CALL SEAICE_EVP_FPADJ_CORE(
     I     xIn,
     U     lam,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA, areaW, areaS,
     U     adpress0, adTsf,
     U     adseaiceMassC, adseaiceMassU, adseaiceMassV,
     U     adFORCEX0, adFORCEY0, aduIceNm1, advIceNm1,
     U     aduVel, advVel, adHEFF, adAREA, adareaW, adareaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     I     SEAICE_fpOmega, nEVPstep, SEAICE_fpExchAdj,
     I     SEAICE_fpExactReplay, .TRUE.,
     I     myTime, myIter, myThid )

C--   unpack the adjoint of the carry
      CALL SEAICE_FPADJ_COPY(
     &     LMF(1),
     &     aduIce,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(2),
     &     advIce,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(3),
     &     adsigma1,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(4),
     &     adsigma2,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(5),
     &     adsigma12,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(6),
     &     ade11,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(7),
     &     ade22,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(8),
     &     ade12,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(9),
     &     adzetaZ,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(10),
     &     addeltaC,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(11),
     &     adstrDivX,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(12),
     &     adstrDivY,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(13),
     &     adDWATN,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(14),
     &     tmpCb,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(15),
     &     adFORCEX,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(16),
     &     adFORCEY,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(17),
     &     adpress,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(18),
     &     adzeta,
     &     myThid )
      CALL SEAICE_FPADJ_COPY(
     &     LMF(19),
     &     adeta,
     &     myThid )
#undef XIF
#undef LMF
#endif /* SIDEDRAG, ITD, CLIPZETA */
#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_EVP_FPADJ_CORE
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_FPADJ_CORE(
     I     xIn,
     U     lam,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA, areaW, areaS,
     U     adpress0, adtensileStrFac,
     U     adseaiceMassC, adseaiceMassU, adseaiceMassV,
     U     adFORCEX0, adFORCEY0, aduIceNm1, advIceNm1,
     U     aduVel, advVel, adHEFF, adAREA, adareaW, adareaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     I     fpOmega, nIter, exchAdj, exactReplay, doPrint,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_FPADJ_CORE
C     | o The work of the fixed-point adjoint of the subcycle loop.
C     |   x = the carry of the loop (xIn = x_0 at the entry of
C     |   the forward call), p = its inputs, F = one modified
C     |   subcycle F_mod at x_N (SEAICE_EVP_FP_SUBCYCLE):
C     |   1. x_N: the loop SEAICE_EVP_SUBCYCLES is rerun from x_0
C     |      and p on work copies (the AD tools do not hand the
C     |      forward's output to a hand-written adjoint)
C     |   2. the frozen state at x_N, once (SEAICE_EVP_FP_FROZEN)
C     |   3. nIter times:
C     |        L = lam; L(uIce,vIce) = E^T L(uIce,vIce) (the
C     |        adjoint of the subcycle's final EXCH_UV_XY_RL,
C     |        SEAICE_FP_EXCH_UV_XY_RL_AD with exchAdj);
C     |        L = F_x^T L, acc = acc + F_p^T L (F_mod_b,
C     |        SEAICE_EVP_FP_SUBCYCLE_AD: overwrites L, increments
C     |        acc); lam = (1-omega) lam + omega L (omega = 1:
C     |        lam = L)
C     |   4. the adjoint of the inputs += omega acc; the anchor
C     |      uIceNm1/vIceNm1 is one of the inputs
C     | exactReplay (a check): the exact adjoint of the loop.
C     | Step 3 runs k = N..1 with F_mod_b at x_(k-1), recomputed
C     | from x_0 (k-1 subcycles; no storage, N(N-1)/2 subcycles
C     | in all), with the frozen state of x_(k-1) and the exact
C     | settings: theta = theta_rigid = 1, no raise, live
C     | replacement pressure, no near-rigid and no thin points,
C     | omega = 1. It must equal the AD tool's own adjoint of the
C     | loop (SEAICEuseEVPfixedPointAdj = .FALSE.) to roundoff.
C     | doPrint: one lambda-monitor line (max|lam| in, out and
C     | their ratio); with debugLevel >= debLevB also the
C     | checksum of the recomputed x_N.
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
      INTEGER nC
      PARAMETER ( nC = 19 )
      _RL xIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
      _RL lam (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
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
      INTEGER kSrf
      INTEGER nEVPstep
      _RL COSWAT
      _RS SINWAT
      _RL ecc2, recip_ecc2, deltaMinSq
      _RL evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT
      _RL EVPcFac
      LOGICAL useAdaptiveEVP
      _RL fpOmega
      INTEGER nIter, exchAdj
      LOGICAL exactReplay, doPrint
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP )
#if (defined SEAICE_ALLOW_SIDEDRAG || defined SEAICE_ITD \
  || defined SEAICE_ALLOW_CLIPZETA)
      STOP 'S/R SEAICE_EVP_FPADJ_CORE: CPP options not covered'
#else
C     === Local variables ===
C     xW             :: work copy of the carry (x_N, or x_(k-1))
C     lW             :: L, the transposed subcycle's work adjoint
C     acc*           :: sum of F_p^T L over the iterations
C     uFr .. rFz     :: the frozen state
C     om             :: omega (1 with exactReplay)
C     lamIn, lamOut  :: max |lam| in and out (lambda monitor)
      INTEGER i, j, k, bi, bj, n, it, kX, nAvg
      _RL xW  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
      _RL lW  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
      _RL acP0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acTS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acMC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acMU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acMV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acFX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acFY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acUN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acVN(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acHF(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL acUV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL acVV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL acAW(1:sNx,1:sNy,nSx,nSy)
      _RL acAS(1:sNx,1:sNy,nSx,nSy)
      _RL uFr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vFr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL dFr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL rFz(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      SAVE xW, lW
      SAVE acP0, acTS, acMC, acMU, acMV, acFX, acFY, acUN, acVN
      SAVE acHF, acAR, acUV, acVV, acAW, acAS
      SAVE uFr, vFr, dFr, rFz
      _RL om, lamIn, lamOut, ratio, chk(3)
      _RL thetaU, thetaRigidU, rigidDeltaU, thinAreaU, thinMassU
      LOGICAL raiseU, prLiveU
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#define XWF(n) xW(1-OLx,1-OLy,1,1,n)
#define LWF(n) lW(1-OLx,1-OLy,1,1,n)
#define LMF(n) lam(1-OLx,1-OLy,1,1,n)

C--   settings
      IF ( exchAdj.NE.1 .AND. exchAdj.NE.2 ) THEN
       WRITE(msgBuf,'(A,I12)') 'SEAICE_EVP_FPADJ_CORE: exchAdj = ',
     &      exchAdj
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R SEAICE_EVP_FPADJ_CORE'
      ENDIF
      IF ( exactReplay .AND. nIter.NE.nEVPstep ) THEN
       WRITE(msgBuf,'(A,2I8)') 'SEAICE_EVP_FPADJ_CORE: exactReplay '//
     &      'needs nIter = nEVPstep: ', nIter, nEVPstep
       CALL PRINT_ERROR( msgBuf, myThid )
       STOP 'ABNORMAL END: S/R SEAICE_EVP_FPADJ_CORE'
      ENDIF
      IF ( useAdaptiveEVP ) THEN
       STOP 'S/R SEAICE_EVP_FPADJ_CORE: adaptive EVP not covered'
      ENDIF
      IF ( exactReplay ) THEN
       om          = 1. _d 0
       thetaU      = 1. _d 0
       thetaRigidU = 1. _d 0
       rigidDeltaU = 0. _d 0
       thinAreaU   = -1. _d 0
       thinMassU   = -1. _d 0
       raiseU      = .FALSE.
       prLiveU     = .TRUE.
      ELSE
       om          = fpOmega
       thetaU      = SEAICE_fpTheta
       thetaRigidU = SEAICE_fpThetaRigid
       rigidDeltaU = SEAICE_fpRigidDelta
       thinAreaU   = SEAICE_fpThinArea
       thinMassU   = SEAICE_fpThinMass
       raiseU      = SEAICE_fpRaise
       prLiveU     = SEAICE_fpPrLive
      ENDIF

C--   1. x_N from x_0 on work copies (not needed with exactReplay)
      IF ( .NOT.exactReplay ) THEN
       CALL SEAICE_EVP_FPADJ_RERUN(
     I      xIn,
     O      xW,
     I      denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I      press0, tensileStrFac,
     I      seaiceMassC, seaiceMassU, seaiceMassV,
     I      FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I      uVel, vVel, HEFF, AREA, areaW, areaS,
     I      kSrf, nEVPstep, COSWAT, SINWAT,
     I      ecc2, recip_ecc2, deltaMinSq,
     I      evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I      EVPcFac, useAdaptiveEVP,
     I      myTime, myIter, myThid )
       IF ( doPrint .AND. debugLevel.GE.debLevB ) THEN
        CALL SEAICE_FPADJ_CHKSUM(
     &      XWF(1),
     &      XWF(2),
     &      XWF(3),
     &      chk,
     &      myThid )
        WRITE(msgBuf,'(A,I10,A,1P3E23.15)')
     &       'SEAICE_EVP_FPADJ: myIter=', myIter,
     &       ' recomputed x_N checksum uIce vIce sigma1:', chk
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT, myThid )
       ENDIF
      ENDIF

C--   2. the frozen state at x_N (once; with exactReplay at each
C     x_(k-1) in step 3)
      IF ( .NOT.exactReplay ) THEN
       CALL SEAICE_EVP_FP_FROZEN(
     I      XWF(1),
     I      XWF(2),
     I      press0, tensileStrFac, seaiceMassC,
     I      recip_ecc2, raiseU, SEAICE_fpKappa0,
     O      uFr, vFr, dFr, rFz,
     I      nEVPstep, myTime, myIter, myThid )
      ENDIF

C--   3. the transposed iteration
      CALL SEAICE_FPADJ_MAXABS( lam, nC, lamIn, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acP0, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acTS, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acMC, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acMU, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acMV, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acFX, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acFY, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acUN, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acVN, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acHF, myThid )
      CALL SEAICE_FPADJ_FILL( 0. _d 0, acAR, myThid )
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           acUV(i,j,k,bi,bj) = 0. _d 0
           acVV(i,j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
          acAW(i,j,bi,bj) = 0. _d 0
          acAS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO it=1,nIter
       IF ( exactReplay ) THEN
C     the subcycle k = N-it+1 at its own x_(k-1), recomputed from x_0
        kX = nIter - it + 1
        CALL SEAICE_EVP_FPADJ_RERUN(
     I       xIn,
     O       xW,
     I       denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I       press0, tensileStrFac,
     I       seaiceMassC, seaiceMassU, seaiceMassV,
     I       FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I       uVel, vVel, HEFF, AREA, areaW, areaS,
     I       kSrf, kX-1, COSWAT, SINWAT,
     I       ecc2, recip_ecc2, deltaMinSq,
     I       evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I       EVPcFac, useAdaptiveEVP,
     I       myTime, myIter, myThid )
        CALL SEAICE_EVP_FP_FROZEN(
     I       XWF(1),
     I       XWF(2),
     I       press0, tensileStrFac, seaiceMassC,
     I       recip_ecc2, raiseU, SEAICE_fpKappa0,
     O       uFr, vFr, dFr, rFz,
     I       kX, myTime, myIter, myThid )
       ELSE
        kX = nEVPstep
       ENDIF
C     L = lam
       DO n=1,nC
        CALL SEAICE_FPADJ_COPY(
     &     LMF(n),
     &     LWF(n),
     &     myThid )
       ENDDO
C     the adjoint of the subcycle's final exchange
C     (seaice_evp.F:1271)
       CALL SEAICE_FP_EXCH_UV_XY_RL_AD(
     &     LWF(1),
     &     LWF(2),
     &     .TRUE.,
     &      exchAdj, myThid )
C     L = F_x^T L, acc = acc + F_p^T L
       CALL SEAICE_EVP_FP_SUBCYCLE_AD(
     I      XWF(1),
     I      XWF(2),
     I      XWF(3),
     I      XWF(4),
     I      XWF(5),
     I      XWF(6),
     I      XWF(7),
     I      XWF(8),
     I      XWF(9),
     I      XWF(10),
     I      XWF(11),
     I      XWF(12),
     I      XWF(13),
     I      XWF(14),
     I      XWF(15),
     I      XWF(16),
     I      XWF(17),
     I      XWF(18),
     I      XWF(19),
     U      LWF(1),
     U      LWF(2),
     U      LWF(3),
     U      LWF(4),
     U      LWF(5),
     U      LWF(6),
     U      LWF(7),
     U      LWF(8),
     U      LWF(9),
     U      LWF(10),
     U      LWF(11),
     U      LWF(12),
     U      LWF(13),
     U      LWF(14),
     U      LWF(15),
     U      LWF(16),
     U      LWF(17),
     U      LWF(18),
     U      LWF(19),
     I      denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I      press0, tensileStrFac,
     I      seaiceMassC, seaiceMassU, seaiceMassV,
     I      FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I      uVel, vVel, HEFF, AREA,
     I      areaW, areaS,
     U      acP0, acTS,
     U      acMC, acMU, acMV,
     U      acFX, acFY, acUN, acVN,
     U      acUV, acVV, acHF, acAR,
     U      acAW, acAS,
     I      uFr, vFr, dFr, rFz,
     I      thetaU, thetaRigidU, rigidDeltaU,
     I      thinAreaU, thinMassU, raiseU, prLiveU,
     I      kSrf, COSWAT, SINWAT, recip_ecc2,
     I      evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I      kX, myTime, myIter, myThid )
C     lam = (1-omega) lam + omega L
       IF ( om.EQ.1. _d 0 ) THEN
        DO n=1,nC
         CALL SEAICE_FPADJ_COPY(
     &     LWF(n),
     &     LMF(n),
     &     myThid )
        ENDDO
       ELSE
C     the loop writes press, zeta, eta (17-19) only with the
C     diagnostics (seaice_evp.F:792-803); without, their
C     adjoint bypasses the loop unchanged (F_mod_b passes it through,
C     L = lam there), so it is not relaxed (an average of two equal
C     values can differ from them by one ulp)
        nAvg = nC
        IF ( .NOT.useDiagnostics ) nAvg = 16
        DO n=1,nAvg
         DO bj=myByLo(myThid),myByHi(myThid)
          DO bi=myBxLo(myThid),myBxHi(myThid)
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             lam(i,j,bi,bj,n) = (1. _d 0 - om)*lam(i,j,bi,bj,n)
     &                        + om*lW(i,j,bi,bj,n)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ENDDO

C--   4. the adjoint of the inputs += omega acc
      CALL SEAICE_FPADJ_ADDW( om, acP0, adpress0, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acTS, adtensileStrFac, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acMC, adseaiceMassC, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acMU, adseaiceMassU, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acMV, adseaiceMassV, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acFX, adFORCEX0, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acFY, adFORCEY0, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acUN, aduIceNm1, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acVN, advIceNm1, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acHF, adHEFF, myThid )
      CALL SEAICE_FPADJ_ADDW( om, acAR, adAREA, myThid )
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           aduVel(i,j,k,bi,bj) = aduVel(i,j,k,bi,bj)
     &                         + om*acUV(i,j,k,bi,bj)
           advVel(i,j,k,bi,bj) = advVel(i,j,k,bi,bj)
     &                         + om*acVV(i,j,k,bi,bj)
          ENDDO
         ENDDO
        ENDDO
        DO j=1,sNy
         DO i=1,sNx
          adareaW(i,j,bi,bj) = adareaW(i,j,bi,bj) + om*acAW(i,j,bi,bj)
          adareaS(i,j,bi,bj) = adareaS(i,j,bi,bj) + om*acAS(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   lambda monitor
      IF ( doPrint ) THEN
       CALL SEAICE_FPADJ_MAXABS( lam, nC, lamOut, myThid )
       ratio = 0. _d 0
       IF ( lamIn.GT.0. _d 0 ) ratio = lamOut/lamIn
       WRITE(msgBuf,'(A,I10,A,I2,L2,A,I4,A,1P3E12.4)')
     &      'SEAICE_EVP_FPADJ: myIter=', myIter, ' exchAdj,replay=',
     &      exchAdj, exactReplay, ' nIter=', nIter,
     &      ' max|lam| in, out, out/in:', lamIn, lamOut, ratio
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid )
      ENDIF
#undef XWF
#undef LWF
#undef LMF
#endif /* SIDEDRAG, ITD, CLIPZETA */
#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SEAICE_EVP_FPADJ_RERUN
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_FPADJ_RERUN(
     I     xIn,
     O     xW,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA, areaW, areaS,
     I     kSrf, nSub, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_FPADJ_RERUN
C     | o xW = the carry after nSub subcycles of the loop
C     |   SEAICE_EVP_SUBCYCLES from the carry xIn and the inputs;
C     |   the inputs are not changed (work copies of the fields
C     |   that SEAICE_EVP_SUBCYCLES may update), and no EVP
C     |   residual is printed. nSub = 0: xW = xIn.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     xIn, xW        :: the carry (the 19 fields of
C                       SEAICE_EVP_SUBCYCLES_FPADJ_AD), in and out
C     nSub           :: number of subcycles
C     the other arguments are those of SEAICE_EVP_SUBCYCLES
      INTEGER nC
      PARAMETER ( nC = 19 )
      _RL xIn (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
      _RL xW  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nC)
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
      INTEGER kSrf
      INTEGER nSub
      _RL COSWAT
      _RS SINWAT
      _RL ecc2, recip_ecc2, deltaMinSq
      _RL evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT
      _RL EVPcFac
      LOGICAL useAdaptiveEVP
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP )
#if !(defined SEAICE_ALLOW_SIDEDRAG || defined SEAICE_ITD \
  || defined SEAICE_ALLOW_CLIPZETA)
C     === Local variables ===
C     wD1 .. wBV     :: work copies of denom1 .. evpBetaV (arguments U
C                       of SEAICE_EVP_SUBCYCLES)
      INTEGER n, nS
      _RL wD1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wD2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wAC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wAZ(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wBU(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL wBV(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      SAVE wD1, wD2, wAC, wAZ, wBU, wBV
#ifdef ALLOW_SEAICE_EVP_RESIDUAL
      LOGICAL noResidual
#endif
#define XWF(n) xW(1-OLx,1-OLy,1,1,n)

      DO n=1,nC
       CALL SEAICE_FPADJ_COPY(
     &     xIn(1-OLx,1-OLy,1,1,n),
     &     XWF(n),
     &     myThid )
      ENDDO
      IF ( nSub.GE.1 ) THEN
       CALL SEAICE_FPADJ_COPY( denom1, wD1, myThid )
       CALL SEAICE_FPADJ_COPY( denom2, wD2, myThid )
       CALL SEAICE_FPADJ_COPY( evpAlphaC, wAC, myThid )
       CALL SEAICE_FPADJ_COPY( evpAlphaZ, wAZ, myThid )
       CALL SEAICE_FPADJ_COPY( evpBetaU, wBU, myThid )
       CALL SEAICE_FPADJ_COPY( evpBetaV, wBV, myThid )
       nS = nSub
#ifdef ALLOW_SEAICE_EVP_RESIDUAL
       noResidual = .FALSE.
#endif
       CALL SEAICE_EVP_SUBCYCLES(
     U      XWF(1),
     U      XWF(2),
     U      XWF(3),
     U      XWF(4),
     U      XWF(5),
     O      XWF(6),
     O      XWF(7),
     O      XWF(8),
     O      XWF(9),
     O      XWF(10),
     O      XWF(11),
     O      XWF(12),
     O      XWF(13),
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     O      XWF(14),
#endif
     O      XWF(15),
     O      XWF(16),
     O      XWF(17),
     O      XWF(18),
     O      XWF(19),
     U      wD1, wD2, wAC, wAZ, wBU, wBV,
     I      press0, tensileStrFac,
     I      seaiceMassC, seaiceMassU, seaiceMassV,
     I      FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I      uVel, vVel, HEFF, AREA,
     I      areaW, areaS,
     I      kSrf, nS, COSWAT, SINWAT,
     I      ecc2, recip_ecc2, deltaMinSq,
     I      evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I      EVPcFac, useAdaptiveEVP,
#ifdef ALLOW_SEAICE_EVP_RESIDUAL
     I      noResidual,
#endif
     I      myTime, myIter, myThid )
      ENDIF
#undef XWF
#endif /* not SIDEDRAG, ITD, CLIPZETA */
#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     helpers of the fixed-point adjoint (every grid point, halos included)

C     SEAICE_FPADJ_COPY: b = a
      SUBROUTINE SEAICE_FPADJ_COPY( a, b, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
      _RL a(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
      INTEGER i, j, bi, bj
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          b(i,j,bi,bj) = a(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END

C     SEAICE_FPADJ_FILL: a = c
      SUBROUTINE SEAICE_FPADJ_FILL( c, a, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
      _RL c
      _RL a(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
      INTEGER i, j, bi, bj
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          a(i,j,bi,bj) = c
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END

C     SEAICE_FPADJ_ADDW: b = b + w a
      SUBROUTINE SEAICE_FPADJ_ADDW( w, a, b, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
      _RL w
      _RL a(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      INTEGER myThid
      INTEGER i, j, bi, bj
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          b(i,j,bi,bj) = b(i,j,bi,bj) + w*a(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END

C     SEAICE_FPADJ_MAXABS: m = max over the nF fields and every grid point
C     of |a| (all processes)
      SUBROUTINE SEAICE_FPADJ_MAXABS( a, nF, m, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
      INTEGER nF
      _RL a(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,nF)
      _RL m
      INTEGER myThid
      INTEGER i, j, bi, bj, n
      m = 0. _d 0
      DO n=1,nF
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           m = MAX( m, ABS(a(i,j,bi,bj,n)) )
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _GLOBAL_MAX_RL( m, myThid )
      RETURN
      END

C     SEAICE_FPADJ_CHKSUM: c = sums over the interior of uIce, vIce and
C     sigma1 (per tile, then GLOBAL_SUM_TILE_RL: the same order for any
C     number of processes)
      SUBROUTINE SEAICE_FPADJ_CHKSUM( u, v, s, c, myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
      _RL u(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL v(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL s(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL c(3)
      INTEGER myThid
      INTEGER i, j, bi, bj
      _RL t1(nSx,nSy), t2(nSx,nSy), t3(nSx,nSy)
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        t1(bi,bj) = 0. _d 0
        t2(bi,bj) = 0. _d 0
        t3(bi,bj) = 0. _d 0
        DO j=1,sNy
         DO i=1,sNx
          t1(bi,bj) = t1(bi,bj) + u(i,j,bi,bj)
          t2(bi,bj) = t2(bi,bj) + v(i,j,bi,bj)
          t3(bi,bj) = t3(bi,bj) + s(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL GLOBAL_SUM_TILE_RL( t1, c(1), myThid )
      CALL GLOBAL_SUM_TILE_RL( t2, c(2), myThid )
      CALL GLOBAL_SUM_TILE_RL( t3, c(3), myThid )
      RETURN
      END

#endif /* SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ */
