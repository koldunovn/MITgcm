#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#include "AD_CONFIG.h"

C     Tapenade hook of the fixed-point adjoint of the mEVP subcycle loop.
C     tools/TAP_support/flow_tap declares SEAICE_EVP_SUBCYCLES_FPADJ an
C     external routine (arguments 1-18 read and written, 19-55 read
C     only), so Tapenade calls SEAICE_EVP_SUBCYCLES_FPADJ_B in the
C     adjoint of SEAICE_EVP instead of differentiating the loop.
C     Tapenade passes each active argument followed by its adjoint and
C     each passive argument alone. The list below assumes that the
C     arguments 1-18 (the loop carry and outputs), press0, seaiceMassC/U/V,
C     FORCEX0/Y0, uIceNm1/vIceNm1, uVel/vVel, HEFF, AREA, areaW/S are
C     active at the call, and that tensileStrFac, denom1/2, evpAlphaC/Z,
C     evpBetaU/V and the scalars are passive. This holds when the sea-ice
C     dynamics and the ocean velocity depend on the control variables.
C     If Tapenade finds a different activity, the generated call in
C     seaice_evp_b.f has a different argument list. A Fortran 77 build
C     still links; the checks at the start of the routine stop the run.

#if ( defined ALLOW_TAPENADE && defined ALLOW_ADJOINT_RUN \
  && defined SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ )
CBOP
C     !ROUTINE: SEAICE_EVP_SUBCYCLES_FPADJ_B
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ_B(
     U     uIce, uIceb, vIce, vIceb, seaice_sigma1, seaice_sigma1b,
     U     seaice_sigma2, seaice_sigma2b, seaice_sigma12,
     U     seaice_sigma12b, e11, e11b, e22, e22b, e12, e12b, zetaZ,
     U     zetaZb, deltaC, deltaCb, stressDivergenceX,
     U     stressDivergenceXb, stressDivergenceY, stressDivergenceYb,
     U     DWATN, DWATNb, FORCEX, FORCEXb, FORCEY, FORCEYb, press,
     U     pressb, zeta, zetab, eta, etab,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, press0b, tensileStrFac, seaiceMassC, seaiceMassCb,
     I     seaiceMassU, seaiceMassUb, seaiceMassV, seaiceMassVb,
     I     FORCEX0, FORCEX0b, FORCEY0, FORCEY0b, uIceNm1, uIceNm1b,
     I     vIceNm1, vIceNm1b, uVel, uVelb, vVel, vVelb, HEFF, HEFFb,
     I     AREA, AREAb, areaW, areaWb, areaS, areaSb,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ_B
C     | o Tapenade's adjoint of SEAICE_EVP_SUBCYCLES_FPADJ: the
C     |   arguments hold the forward values at the call's entry,
C     |   the *b arguments the adjoints (carry: overwritten,
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
      _RL uIceb         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIce          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceb         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma1b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma2b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma12(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaice_sigma12b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e11b          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e22b          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL e12b          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetaZ         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetaZb        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL deltaCb       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceXb(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL stressDivergenceYb(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DWATN         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL DWATNb        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEXb       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEYb       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL pressb        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zeta          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL zetab         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL eta           (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL etab          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL denom1        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL denom2        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpAlphaC     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpAlphaZ     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpBetaU      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL evpBetaV      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL press0b       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tensileStrFac (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassCb  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassU   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassUb  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassV   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL seaiceMassVb  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEX0b      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY0       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL FORCEY0b      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uIceNm1       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uIceNm1b      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceNm1       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vIceNm1b      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uVel          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL uVelb         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVel          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL vVelb         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HEFF          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL HEFFb         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL AREA          (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL AREAb         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL areaW         (1:sNx,1:sNy,nSx,nSy)
      _RL areaWb        (1:sNx,1:sNy,nSx,nSy)
      _RL areaS         (1:sNx,1:sNy,nSx,nSy)
      _RL areaSb        (1:sNx,1:sNy,nSx,nSy)
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
       WRITE(msgBuf,'(2A)') 'SEAICE_EVP_SUBCYCLES_FPADJ_B: the ',
     &      'generated call (seaice_evp_b.f) does not match the'
       CALL PRINT_ERROR( msgBuf, 1 )
       WRITE(msgBuf,'(2A)') 'SEAICE_EVP_SUBCYCLES_FPADJ_B: ',
     &      'argument list of seaice_evp_subcycles_fpadj_tap.F'
       CALL PRINT_ERROR( msgBuf, 1 )
       STOP 'ABNORMAL END: S/R SEAICE_EVP_SUBCYCLES_FPADJ_B'
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
     U     uIceb, vIceb, seaice_sigma1b, seaice_sigma2b,
     U     seaice_sigma12b, e11b, e22b, e12b, zetaZb, deltaCb,
     U     stressDivergenceXb, stressDivergenceYb, DWATNb, FORCEXb,
     U     FORCEYb, pressb, zetab, etab, press0b, seaiceMassCb,
     U     seaiceMassUb, seaiceMassVb, FORCEX0b, FORCEY0b, uIceNm1b,
     U     vIceNm1b, uVelb, vVelb, HEFFb, AREAb, areaWb, areaSb,
     I     myTime, myIter, myThid )

      RETURN
      END
#endif /* ALLOW_TAPENADE, ALLOW_ADJOINT_RUN, FIXEDPOINT_ADJ */
