! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE GPPWCVFE(YDGEOMETRY,KPROMA,KSTART,KPROF,KFLEV,PWCH,PQ,PRXP,PWCH_DP)

!**** *GPPWCVFE* - Computes half level PWC.

!     Purpose.
!     --------
!           Computes half level PWC.

!**   Interface.
!     ----------
!        *CALL* *GPPWCVFE(KPROMA,KSTART,KPROF,KFLEV,PWCH,PQ,PRXP,PWCH_DP)

!        Explicit arguments :
!        --------------------
!        KPROMA                     - HORIZ. DIMENSIONING      (INPUT)
!        KSTART                     - START OF WORK            (INPUT)
!        KPROF                      - DEPTH OF WORK            (INPUT)
!        KFLEV                      - NUMBER OF LEVELS         (INPUT)
!        PWCH (KPROMA,0:KFLEV)      - HALF LEVEL PWC           (OUTPUT)
!        PQ  (KPROMA,KFLEV)         - HUMIDITY                 (INPUT)
!        PRXP(KPROMA,0:KFLEV,NPPM)  - HALF,FULL AND LN HALF LEVEL
!                                     PRESSURES (SEE PPINIT)   (INPUT)
!        PWCH_DP (KPROMA,0:KFLEV)   - OPTIONAL DP HALF LEVEL PWC (OUTPUT)

!        Implicit arguments :    None.
!        --------------------

!     Method.
!     -------
!        See documentation.

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Michail Diamantakis, based on GPPWC(), extended to
!        handle Vertical Finite Element discretization

!     Modifications.
!     --------------
!   Original : May 2012
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : RG
USE PARDIM       , ONLY : JPNPPM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWCH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV,JPNPPM) 
REAL(KIND=JPRD)   ,OPTIONAL,INTENT(OUT) :: PWCH_DP(KPROMA,0:KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL, JLEV

REAL(KIND=JPRB) :: ZDELP, ZRGI
REAL(KIND=JPRD) :: ZRGI_DP
REAL(KIND=JPRB) :: ZWCIN(KPROMA,0:KFLEV+1)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPPWCVFE',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
!     ------------------------------------------------------------------

!*       1.    COMPUTES HALF LEVEL PWC
!              -----------------------

ZRGI=1.0_JPRB/RG
ZRGI_DP=1.0_JPRD/REAL(RG, KIND=JPRD)

IF (YDGEOMETRY%YRCVER%LVERTFE) THEN
  DO JLEV=1,KFLEV
    DO JL=KSTART,KPROF
      ZDELP=YDVAB%VDELA(JLEV) + YDVAB%VDELB(JLEV)*PRXP(JL,KFLEV,1)
      ZWCIN(JL,JLEV)=PQ(JL,JLEV)*ZDELP*YDVETA%VFE_RDETAH(JLEV)*ZRGI
    ENDDO
  ENDDO
  DO JL=KSTART,KPROF
    ZWCIN(JL,0)=0.0_JPRB
    ZWCIN(JL,KFLEV+1)=0.0_JPRB
  ENDDO
  IF (PRESENT(PWCH_DP)) THEN
    CALL VERDISINT(YDVFE,YDGEOMETRY%YRCVER,'ITOP','11',KPROMA,KSTART,KPROF,KFLEV,ZWCIN,PWCH,POUT_DP=PWCH_DP)
  ELSE
    CALL VERDISINT(YDVFE,YDGEOMETRY%YRCVER,'ITOP','11',KPROMA,KSTART,KPROF,KFLEV,ZWCIN,PWCH)
  ENDIF
ELSE
  DO JL=KSTART,KPROF
    PWCH(JL,0)=0.0_JPRB
    IF (PRESENT(PWCH_DP)) PWCH_DP(JL,0)=0.0_JPRD
  ENDDO
  DO JLEV=1,KFLEV
    DO JL=KSTART,KPROF
      ZDELP=PRXP(JL,JLEV,1)-PRXP(JL,JLEV-1,1)
      PWCH(JL,JLEV)=PWCH(JL,JLEV-1)+PQ(JL,JLEV)*ZDELP*ZRGI
      IF (PRESENT(PWCH_DP)) THEN
        PWCH_DP(JL,JLEV)=PWCH_DP(JL,JLEV-1)+REAL(PQ(JL,JLEV), KIND=JPRD)*REAL(ZDELP, KIND=JPRD)*ZRGI_DP
      ENDIF
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPPWCVFE',1,ZHOOK_HANDLE)
END SUBROUTINE GPPWCVFE
