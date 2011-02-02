MODULE atomdata

  TYPE OrbitInfo
     INTEGER :: nps, npp, npd ,npf, npg, norbit
     INTEGER, POINTER :: np(:),l(:)
     REAL(8), POINTER :: eig(:),occ(:),wfn(:,:)
  END TYPE OrbitInfo

  TYPE SCFInfo
     INTEGER :: iter
     REAL(8) :: delta,eone,ekin,estatic,ecoul,eexc,etot
  END TYPE SCFInfo

  TYPE FCinfo
     REAL(8), POINTER :: coreden(:),valeden(:)
     REAL(8) :: evale,zvale,zcore,corekin
     INTEGER :: norbit
     LOGICAL, POINTER :: iscore(:)
  END TYPE FCinfo

  LOGICAL :: scalarrelativistic
  LOGICAL :: finitenucleus
  LOGICAL :: gaussianshapefunction,besselshapefunction
  REAL(8) :: minlogderiv,maxlogderiv,gaussparam
  INTEGER :: coretailpoints,ivale,ivion,nlogderiv

CONTAINS

  SUBROUTINE InitOrbit1(Orbit,j)
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    INTEGER, INTENT(IN) :: j
    INTEGER :: ok
    ALLOCATE(Orbit%np(j),Orbit%l(j),&
&            Orbit%eig(j),Orbit%occ(j),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'Error in allocation of nl, l, occ...',ok
       STOP
    ENDIF
    Orbit%np=0;Orbit%l=0
    Orbit%eig=0.d0;Orbit%occ=0.d0
  END SUBROUTINE InitOrbit1

  SUBROUTINE InitOrbit2(Orbit,n,j)
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    INTEGER, INTENT(IN) :: n,j
    INTEGER :: ok
    ALLOCATE(Orbit%wfn(n,j),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'Error in allocation of wfn,...',ok
       STOP
    ENDIF
    Orbit%wfn=0.d0
  END SUBROUTINE InitOrbit2

  SUBROUTINE DestroyOrbit(Orbit)
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    IF (ASSOCIATED(Orbit%l)) DEALLOCATE(Orbit%l)
    IF (ASSOCIATED(Orbit%np)) DEALLOCATE(Orbit%np)
    IF (ASSOCIATED(Orbit%eig)) DEALLOCATE(Orbit%eig)
    IF (ASSOCIATED(Orbit%occ)) DEALLOCATE(Orbit%occ)
    IF (ASSOCIATED(Orbit%wfn)) DEALLOCATE(Orbit%wfn)
  END SUBROUTINE DestroyOrbit

  SUBROUTINE InitFC(FC,norbit,n)
    INTEGER, INTENT(IN) :: norbit,n
    TYPE (FCInfo), INTENT(INOUT) :: FC
    INTEGER :: ok
    FC%norbit=norbit
    ALLOCATE(FC%coreden(n),FC%valeden(n),FC%iscore(norbit),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'Error in allocation of coreden, valeden,...',ok
       STOP
    ENDIF
    FC%coreden=0.d0;FC%valeden=0.d0;FC%iscore=.true.
  END SUBROUTINE InitFC

  SUBROUTINE DestroyFC(FC)
    TYPE (FCInfo), INTENT(INOUT) :: FC
    IF (ASSOCIATED(FC%coreden)) DEALLOCATE(FC%coreden)
    IF (ASSOCIATED(FC%valeden)) DEALLOCATE(FC%valeden)
    IF (ASSOCIATED(FC%iscore)) DEALLOCATE(FC%iscore)
  END SUBROUTINE DestroyFC

END MODULE atomdata
