#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

#ifdef SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ
CBOP
C     !ROUTINE: SEAICE_EVP_SUBCYCLES_FPADJ
C     !INTERFACE:
      SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ(
     U     uIce, vIce,
     U     seaice_sigma1, seaice_sigma2, seaice_sigma12,
     U     e11, e22, e12, zetaZ, deltaC,
     U     stressDivergenceX, stressDivergenceY, DWATN,
     U     FORCEX, FORCEY, press, zeta, eta,
     I     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA, areaW, areaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
     I     myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SEAICE_EVP_SUBCYCLES_FPADJ
C     | o Wrapper of SEAICE_EVP_SUBCYCLES that SEAICE_EVP calls
C     |   instead of it when SEAICEuseEVPfixedPointAdj = .TRUE.
C     |   in adjoint builds; the forward is the same.
C     *==========================================================*
C     | This routine is the place where the hand-written
C     | fixed-point adjoint of the subcycle loop is attached:
C     | it is not differentiated by the AD tools (it is not in
C     | seaice_ad_diff.list). TAF finds its adjoint through the
C     | ADNAME directive in seaice.flow
C     | (SEAICE_EVP_SUBCYCLES_FPADJ_MAD), Tapenade through its
C     | entry in tools/TAP_support/flow_tap
C     | (SEAICE_EVP_SUBCYCLES_FPADJ_B). Both call
C     | SEAICE_EVP_SUBCYCLES_FPADJ_AD.
C     | The argument list is the same for every set of CPP
C     | options, because the flow directives of both tools are
C     | not run through CPP. The fields that the loop reads or
C     | writes only with some CPP options are not arguments:
C     | - SEAICE_ALLOW_SIDEDRAG, SEAICE_ITD, SEAICE_ALLOW_CLIPZETA
C     |   are not covered (SEAICE_CHECK stops; so does this
C     |   routine);
C     | - with SEAICE_ALLOW_BOTTOMDRAG, CbotC is zero, as
C     |   SEAICE_CHECK refuses SEAICEbasalDragK2 > 0;
C     | - with ALLOW_SEAICE_EVP_RESIDUAL, no residual is printed.
C     | denom1/2, evpAlphaC/Z and evpBetaU/V are input only: the
C     | loop changes them only with adaptive EVP, which
C     | SEAICE_CHECK refuses.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     The arguments keep the names of the SEAICE.h and DYNVARS.h fields
C     and of the local variables of SEAICE_EVP that are passed:
C     loop carry (read and written):
C     uIce, vIce     :: ice velocity
C     seaice_sigma1/2/12 :: EVP stresses
C     written (and read by the next subcycle):
C     e11, e22, e12  :: strain rates
C     zetaZ, deltaC  :: bulk viscosity at Z points, Delta at C points
C     stressDivergenceX/Y :: divergence of the stress tensor
C     DWATN          :: ocean drag coefficient
C     FORCEX, FORCEY :: forcing of the momentum equations
C     press, zeta, eta :: pressure, viscosities (diagnostics)
C     read:
C     denom1, denom2, evpAlphaC/Z, evpBetaU/V :: fields of SEAICE_EVP
C     press0, tensileStrFac, seaiceMassC/U/V, FORCEX0, FORCEY0,
C     uIceNm1, vIceNm1, HEFF, AREA :: SEAICE.h fields
C     uVel, vVel     :: ocean velocity (DYNVARS.h)
C     areaW, areaS   :: fractional area at velocity points
C     kSrf, nEVPstep, COSWAT, SINWAT, ecc2, recip_ecc2, deltaMinSq,
C     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT, EVPcFac,
C     useAdaptiveEVP :: scalars set in SEAICE_EVP
C     myTime         :: Simulation time
C     myIter         :: Simulation timestep number
C     myThid         :: my Thread Id. number
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
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if ( defined SEAICE_CGRID &&  defined SEAICE_ALLOW_EVP )
#if (defined SEAICE_ALLOW_SIDEDRAG || defined SEAICE_ITD \
  || defined SEAICE_ALLOW_CLIPZETA)
      STOP 'ABNORMAL END: S/R SEAICE_EVP_SUBCYCLES_FPADJ: CPP options'
#else /* not SIDEDRAG, ITD, CLIPZETA */
C     === Local variables ===
C     i, j, bi, bj   :: loop indices
C     CbotC          :: zero bottom drag coefficient
C     chk            :: checksum of the loop's output (debug print)
      INTEGER i, j, bi, bj
#ifdef SEAICE_ALLOW_BOTTOMDRAG
      _RL CbotC         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef ALLOW_SEAICE_EVP_RESIDUAL
      LOGICAL printResidual
#endif
      _RL chk(3)
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifdef SEAICE_ALLOW_BOTTOMDRAG
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          CbotC(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif
#ifdef ALLOW_SEAICE_EVP_RESIDUAL
      printResidual = .FALSE.
#endif

      CALL SEAICE_EVP_SUBCYCLES(
     U     uIce, vIce,
     U     seaice_sigma1, seaice_sigma2, seaice_sigma12,
     O     e11, e22, e12, zetaZ, deltaC,
     O     stressDivergenceX, stressDivergenceY,
     O     DWATN,
#ifdef SEAICE_ALLOW_BOTTOMDRAG
     O     CbotC,
#endif
     O     FORCEX, FORCEY, press, zeta, eta,
     U     denom1, denom2, evpAlphaC, evpAlphaZ, evpBetaU, evpBetaV,
     I     press0, tensileStrFac,
     I     seaiceMassC, seaiceMassU, seaiceMassV,
     I     FORCEX0, FORCEY0, uIceNm1, vIceNm1,
     I     uVel, vVel, HEFF, AREA,
     I     areaW, areaS,
     I     kSrf, nEVPstep, COSWAT, SINWAT,
     I     ecc2, recip_ecc2, deltaMinSq,
     I     evpStarFac, evpRevFac, recip_evpRevFac, recip_deltaT,
     I     EVPcFac, useAdaptiveEVP,
#ifdef ALLOW_SEAICE_EVP_RESIDUAL
     I     printResidual,
#endif
     I     myTime, myIter, myThid )

      IF ( debugLevel.GE.debLevB ) THEN
C     checksum of the loop's output x_N, to compare with the one of the
C     x_N that SEAICE_EVP_FPADJ_CORE recomputes
       CALL SEAICE_FPADJ_CHKSUM( uIce, vIce, seaice_sigma1, chk,
     &      myThid )
       WRITE(msgBuf,'(A,I10,A,1P3E23.15)')
     &      'SEAICE_EVP_FPADJ: myIter=', myIter,
     &      ' forward    x_N checksum uIce vIce sigma1:', chk
       CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid )
      ENDIF

#endif /* not SIDEDRAG, ITD, CLIPZETA */
#endif /* SEAICE_CGRID and SEAICE_ALLOW_EVP */

      RETURN
      END
#endif /* SEAICE_ALLOW_EVP_FIXEDPOINT_ADJ */
