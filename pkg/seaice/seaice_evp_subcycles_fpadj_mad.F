#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#include "AD_CONFIG.h"

C     TAF hook of the fixed-point adjoint of the mEVP subcycle loop.
C     pkg/seaice/seaice.flow names SEAICE_EVP_SUBCYCLES_FPADJ_MAD the
C     adjoint (ADNAME) of SEAICE_EVP_SUBCYCLES_FPADJ, with all 55
C     arguments DEPEND, the arguments 1-18, 25 and 27-39 ACTIVE, and
C     1-18 OUTPUT. The argument list below follows the convention of the
C     hand-written adjoints that TAF calls in MITgcm (CG2D_MAD,
C     ADSTREAMICE_CG_SOLVE): for each argument in order, its value if it
C     is DEPEND, then its adjoint if it is ACTIVE. The adjoint of an
C     argument that is INPUT and OUTPUT is overwritten (the adjoint of
C     the output on entry, of the input on exit), that of an INPUT-only
C     argument is incremented. SEAICE_EVP stores the arguments before
C     the call (CADJ STORE), so that their values are those at the entry
C     of the forward call.

#if ( defined ALLOW_AUTODIFF_TAMC && defined ALLOW_ADJOINT_RUN \
  && defined SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ )
CBOP
C     !ROUTINE: SEAICE_EVP_SUBCYCLES_FPADJ_MAD
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ_MAD(
     U     uIce, aduIce, vIce, advIce, seaice_sigma1, adseaice_sigma1,
     U     seaice_sigma2, adseaice_sigma2, seaice_sigma12,
     U     adseaice_sigma12, e11, ade11, e22, ade22, e12, ade12, zetaZ,
     U     adzetaZ, deltaC, addeltaC, stressDivergenceX,
     U     adstressDivergenceX, stressDivergenceY, adstressDivergenceY,
     U     DWATN, adDWATN, FORCEX, adFORCEX, FORCEY, adFORCEY, press,
     U     adpress, zeta, adzeta, eta, adeta,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, adpress0, tensileStrFac, seaiceMassC, adseaiceMassC,
     I     seaiceMassU, adseaiceMassU, seaiceMassV, adseaiceMassV,
     I     FORCEX0, adFORCEX0, FORCEY0, adFORCEY0, uIceNm1, aduIceNm1,
     I     vIceNm1, advIceNm1, uVel, aduVel, vVel, advVel, HEFF, adHEFF,
     I     AREA, adAREA, areaW, adareaW, areaS, adareaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ_MAD
C     | o TAF's adjoint of SEAICE_EVP_SUBCYCLES_FPADJ (ADNAME in
C     |   seaice.flow): the arguments hold the forward values at
C     |   the call's entry (DEPEND), the ad* arguments the
C     |   adjoints of the ACTIVE arguments (carry: overwritten,
C     |   inputs: incremented); calls the hand-written adjoint
C     |   SEAICE_EVP_SUBCYCLES_FPADJ_AD.
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
      _RL uIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL aduIce         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL advIce         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaice_sigma1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaice_sigma2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaice_sigma12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade11          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade22          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ade12          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetaZ         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adzetaZ        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL addeltaC       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adstressDivergenceX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adstressDivergenceY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DWATN         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adDWATN        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEX       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEY       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adpress        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zeta          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adzeta         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL eta           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adeta          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL denom1        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL denom2        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpAlphaC     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpAlphaZ     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpBetaU      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpBetaV      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adpress0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tensileStrFac (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaiceMassC  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassU   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaiceMassU  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassV   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adseaiceMassV  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEX0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adFORCEY0      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uIceNm1       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL aduIceNm1      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceNm1       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL advIceNm1      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uVel          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL aduVel         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVel          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL advVel         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HEFF          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adHEFF         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL AREA          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL adAREA         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL areaW         (1:sNx,1:sNy,nSx,nSy)
      _RL adareaW        (1:sNx,1:sNy,nSx,nSy)
      _RL areaS         (1:sNx,1:sNy,nSx,nSy)
      _RL adareaS        (1:sNx,1:sNy,nSx,nSy)
      INTEGER kSrf
      INTEGER nEVPstep
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
C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C     a different argument list of the generated call shifts the scalar
C     arguments: stop rather than run with wrong arguments
      IF ( nEVPstep.NE.SEAICEnEVPstarSteps
     &     .OR. ( kSrf.NE.1 .AND. kSrf.NE.Nr )
     &     .OR. myThid.LT.1 .OR. myThid.GT.MAX_NO_THREADS ) THEN
       WRITE(msgBuf,'(2A)') 'SEAICE_EVP_SUBCYCLES_FPADJ_MAD: the ',
     &      'generated call (seaice_evp_ad.f) does not match the'
       CALL PRINT_ERROR( msgBuf, 1 )
       WRITE(msgBuf,'(2A)') 'SEAICE_EVP_SUBCYCLES_FPADJ_MAD: ',
     &      'argument list of seaice_evp_subcycles_fpadj_mad.F'
       CALL PRINT_ERROR( msgBuf, 1 )
       STOP 'ABNORMAL END: S/R SEAICE_EVP_SUBCYCLES_FPADJ_MAD'
      ENDIF

      CALL SEAICE_EVP_SUBCYCLES_FPADJ_AD(
     I     uIce, vIce, seaice_sigma1, seaice_sigma2, seaice_sigma12,
     I     e11, e22, e12, zetaZ, deltaC, stressDivergenceX,
     I     stressDivergenceY, DWATN, FORCEX, FORCEY, press, zeta, eta,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac, seaiceMassC, seaiceMassU,
     I     seaiceMassV, FORCEX0, FORCEY0, uIceNm1, vIceNm1, uVel, vVel,
     I     HEFF, AREA, areaW, areaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     U     aduIce, advIce, adseaice_sigma1, adseaice_sigma2,
     U     adseaice_sigma12, ade11, ade22, ade12, adzetaZ, addeltaC,
     U     adstressDivergenceX, adstressDivergenceY, adDWATN, adFORCEX,
     U     adFORCEY, adpress, adzeta, adeta, adpress0, adseaiceMassC,
     U     adseaiceMassU, adseaiceMassV, adFORCEX0, adFORCEY0,
     U     aduIceNm1,
     U     advIceNm1, aduVel, advVel, adHEFF, adAREA, adareaW, adareaS,
     I     myTime, myIter, myThid )

      RETURN
      END
#endif /* ALLOW_AUTODIFF_TAMC, ALLOW_ADJOINT_RUN, FIXEDPOINT_ADJ */
