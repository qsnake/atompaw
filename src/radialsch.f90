MODULE radialsch
  USE GlobalMath
  USE gridmod
  USE atomdata
  USE calcpotential

  IMPLICIT NONE

CONTAINS
!*******************************************************************
!  FUNCTION wfninit(nz,l,v0,v0p,energy,r)
!*******************************************************************
  FUNCTION wfninit(nz,l,v0,v0p,energy,r)
   ! returns the solution of the Schroedinger equation near r=0
   !  using power series expansion
    REAL(8) :: wfninit
    REAL(8), INTENT(IN) :: nz
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: v0,v0p,energy,r

    REAL(8) :: c1,c2,c3

       c1=-nz/(l+1.d0)
       c2=((v0-energy)-2*nz*c1)/(4*l+6.d0)
       c3=(v0p+(v0-energy)*c1-2*nz*c2)/(6*l+12.d0)

        wfninit=(r**(l+1))*(1+r*(c1+r*(c2+r*c3)))
  End function wfninit
    
  
  !**********************************************************************
  !      subroutine unboundsch(Grid,Pot,nr,l,energy,wfn,nodes)
  !  pgm to solve radial schroedinger equation for unbound states
  !    at energy 'energy' and at angular momentum l
  !
  !    with potential rv/r, given in uniform mesh of n points
  !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
  !                               p((n+1)*h)=0
  !  nz=nuclear charge
  !
  !  uses Noumerov algorithm
  !
  !  For l=0,1 corrections are needed to approximate wfn(r=0)
  !     These depend upon:
  !         e0 (current guess of energy eigenvalue)
  !         l,nz
  !         v(0) == v0 electronic potential at r=0
  !         v'(0) == v0p derivative of electronic potential at r=0
  !
  ! also returns node == number of nodes for calculated state
  !************************************************************************
  SUBROUTINE unboundsch(Grid,Pot,nr,l,energy,wfn,nodes)
    TYPE(GridInfo), INTENT(IN)  :: Grid
    TYPE(PotentialInfo), INTENT(IN)  :: Pot
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(IN) :: energy
    REAL(8), INTENT(INOUT) :: wfn(:)
    INTEGER, INTENT(INOUT) :: nodes

    INTEGER :: n,nz,i,j,k,ierr
    REAL(8) :: zeroval,scale

    n=Grid%n
    IF (nr > n) THEN
       WRITE(6,*) 'Error in unboundsch -- nr > n', nr,n
       STOP
    ENDIF

    ! initialize wfn
    wfn=0
    wfn(2)=wfninit(Pot%nz,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
    zeroval=0
    if (l==0) zeroval=-2*Pot%nz
    if (l==1) zeroval=2

    call forward_numerov(Grid,l,nr,energy,Pot%rv,zeroval,wfn,nodes)
    !
    ! normalize to unity within integration range
    !
    scale=1.d0/overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
    scale=SIGN(SQRT(scale),wfn(nr-2))
    wfn(1:nr)=wfn(1:nr)*scale  

  END SUBROUTINE unboundsch

!******************************************************************
!  SUBROUTINE boundsch(Grid,Pot,Orbit,l,start,nroot,emin,ierr)
!******************************************************************
  SUBROUTINE boundsch(Grid,Pot,Orbit,l,start,nroot,emin,ierr)
    !  pgm to solve radial schroedinger equation for nroot bound state
    !    energies and wavefunctions for angular momentum l
    !    with potential rv/r, given in uniform mesh of n points
    !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
    !                               p((n+1)*h)=0
    !  nz=nuclear charge
    !  emin=is estimate of lowest eigenvalue; used if nz=0
    !     otherwise, set to the value of -(nz/(l+1))**2
    !
    !  It is assumed that the wavefunction has np-l-1 nodes, where
    !    np is the principle quantum number-- np=1,2,..nroot
    !
    !  uses Noumerov algorithm
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
    !
    !  Corrections are also needed for r>n*h, depending on:
    !         e0 (current guess of energy eigenvalue
    !         the extrapolated value of rv == r * v
    !
    ! ierr=an nroot digit number indicating status of each root
    !   a digit of 1 indicates success in converging root
    !              2 indicates near success in converging root
    !              9 indicates that root not found
    !
    ! first check how many roots expected =  ntroot (returned as argument)
    !
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
    INTEGER, INTENT(IN) :: l,start,nroot
    INTEGER, INTENT(INOUT) :: ierr
    REAL(8), INTENT(INOUT) :: emin

    REAL(8), PARAMETER :: convre=1.d-10,vlrg=1.d30
    INTEGER, PARAMETER :: niter=1000

    REAL(8), POINTER :: rv(:),eig(:),wfn(:,:)
    REAL(8), ALLOCATABLE :: p1(:),p2(:),dd(:)
    INTEGER :: nz,n
    REAL(8) :: h,v0,v0p
    REAL(8) :: err,convrez,energy,zeroval
    REAL(8) :: scale,emax,best,rout,ppp
    REAL(8) :: arg,r,r2,veff,pppp1,rin,dele,x,rvp1,pnp1,bnp1
    INTEGER :: iter,i,j,k,node,match,mxroot,ntroot,ir,iroot
    INTEGER :: least,many,ifac
    LOGICAL :: ok

    n=Grid%n
    h=Grid%h
    wfn=>Orbit%wfn
    eig=>Orbit%eig
    ALLOCATE(p1(n),p2(n),dd(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) ' Error in boundsch allocation ',i,n
       STOP
    ENDIF
    nz=Pot%nz
    v0=Pot%v0
    v0p=Pot%v0p
    rv=>Pot%rv
    err=n*nz*(h**4)
    convrez=convre
    IF (nz.GT.0) convrez=convre*nz
    !     write(6,*) 'expected error = ',err
    ierr=0

    WRITE(6,*) 'z , l = ',nz,l
    ! check how many roots expected by integration outward at
    !   energy = 0
    energy = 0
    !
    !  start outward integration
    !    correct behavior near r=0
    ! initialize wfn
    p1=0
    p1(2)=wfninit(Pot%nz,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
    zeroval=0
    if (l==0) zeroval=-2*Pot%nz
    if (l==1) zeroval=2

    call forward_numerov(Grid,l,n,energy,Pot%rv,zeroval,p1,node)

    WRITE(6,*) ' nodes at e=0  ', node
    mxroot=node+1
    ntroot=node
    IF (mxroot.LT.nroot) THEN
       WRITE(6,*)'error in boundsch - for l = ',l
       WRITE(6,*) nroot,' states requested but only',mxroot,' possible'
       DO ir=mxroot+1,nroot
          ierr=ierr+9*(10**(ir-1))
       ENDDO
    ENDIF
    mxroot=min0(mxroot,nroot)
    !
    IF (nz.EQ.0) energy=-ABS(emin)
    IF (nz.NE.0) energy=-(nz/(l+1.d0))**2
    emin=energy-err
    emax=0.d0

    DO iroot=1,mxroot
       best=1.d10; dele=1.d10
       energy=emin+err
       IF (energy.LT.emin) energy=emin
       IF (energy.GT.emax) energy=emax
       ok=.FALSE.
       !write(6,*) 'iter,iroot,energy',iter,iroot,energy
       !write(6,*) 'emin,max',emin,emax
       BigIter: DO iter=1,niter
          !write(6,*) 'In iter with energy', iter,energy,niter,l,iroot
          !  start inward integration
          !  start integration at n
          ! find classical turning point
          call ClassicalTurningPoint(Grid,Pot,l,energy,match)
          match=max(5,match)
          match=min(n-15,match)
          x=0.5d0*(rv(n)/Grid%r(n)+rv(n-1)/Grid%r(n-1))+l*(l+1)/(Grid%r(n)**2)
          ppp=SQRT(ABS(x-energy))
          p2=0
          p2(n)=1
          p2(n-1)=exp(-ppp*(Grid%r(n-1)-Grid%r(n)))
          call backward_numerov(Grid,l,match,energy,rv,p2)
          match=match+6
          call derivative(Grid,p2,dd,match-5,match+5)
          rin=dd(match)/p2(match)
              ! write(6,*) ' match point = ',match,rin,p2(match)
          !  start outward integration
          !    correct behavior near r=0
          ! initialize p1
          p1=0
          p1(2)=wfninit(Pot%nz,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
          zeroval=0
          if (l==0) zeroval=-2*Pot%nz
          if (l==1) zeroval=2

          call forward_numerov(Grid,l,match+6,energy,Pot%rv,zeroval,p1,node)
             
          call derivative(Grid,p1,dd,match-5,match+5)
          rout=dd(match)/p1(match)
             !write(6,*) 'node,match,rin,rout',node,(iroot-1),match,rin,rout
          ! check whether node = (iroot-1)
          !   not enough nodes -- raise energy
          IF (node.LT.iroot-1) THEN
             emin=MAX(emin,energy)-err
             energy=emax-(emax-energy)*ranx()
             ifac=9
             !   too many nodes -- lower energy
          ELSEIF (node.GT.iroot-1) THEN
             IF (energy.LE.emin) THEN
                ierr=ierr+9*(10**(iroot-1))
                WRITE(6,*) 'boundsch error -- emin too high',l,nz,emin,energy
                STOP
             ENDIF
             emax=MIN(emax,energy+err)
             energy=emin+(energy-emin)*ranx()
             !   correct number of nodes -- estimate correction
          ELSEIF (node.EQ.iroot-1) THEN
             DO j=1,match
                p1(j)=p1(j)/p1(match)
                         !write(6,*) 'j,p1',j,p1(j)
             ENDDO
             DO j=match,n
                p1(j)=p2(j)/p2(match)
                         !write(6,*) 'j,p2',j,p1(j)
             ENDDO
             scale=1.d0/overlap(Grid,p1,p1)
             dele=(rout-rin)*scale
                  !write(6,*) 'energy,dele,scale',energy,dele,scale
             x=ABS(dele)
             IF (x.LT.best) THEN
                scale=SQRT(scale)
                p1(1:n)=p1(1:n)*scale 
                k=start+iroot-1
                wfn(1:n,k)=p1(1:n) 
                eig(k)=energy
                !write(6,*) 'root',l,iroot,eig(k),emin,emax
                best=x
             ENDIF
             IF (ABS(dele).LE.convrez) THEN
                write(6,*) 'iter with dele' , iter,dele
                ok=.TRUE.
                !  eigenvalue found
                ierr=ierr+10**(iroot-1)
                IF (iroot+1.LE.mxroot) THEN
                   emin=energy+err
                   emax=0
                   energy=(emin+emax)/2
                   IF (energy.LT.emin) energy=emin
                   IF (energy.GT.emax) energy=emax
                   best=1.d10
                ENDIF
                EXIT BigIter
             ENDIF
             IF (ABS(dele).GT.convrez) THEN
                !write(6,*) 'iter with dele' , iter,dele
                energy=energy+dele
                ! if energy is out of range, pick random energy in correct range
                IF (emin-energy.GT.convrez.OR.energy-emax.GT.convrez)         &
&                    energy=emin+(emax-emin)*ranx()
                ifac=2
                !write(6,*) 'continuing with iter dele', iter,dele
             ENDIF
          ENDIF
       ENDDO BigIter !iter
       IF (.NOT.ok) THEN
          ierr=ierr+ifac*(10**(iroot-1))
          WRITE(6,*) 'no convergence in boundsch',iroot,l,dele,energy
          WRITE(6,*) ' best guess of eig, dele = ',eig(start+iroot-1),best
          IF (iroot.LT.mxroot) THEN
             DO ir=iroot+1,mxroot
                ierr=ierr+9*(10**(ir-1))
             ENDDO
          ENDIF
        ! reset wfn with hydrogenic form
        k=start+iroot-1; j=iroot+l+1
        wfn(:,k)=0
        ppp=(j)*sqrt(abs(eig(start+iroot-1)))
        do i=2,n
           wfn(i,k)=hwfn(ppp,j,l,Grid%r(i))
        enddo
       ENDIF
    ENDDO !iroot

    DEALLOCATE(p1,p2,dd)
             write(6,*) 'returning from boundsch -- ierr=',ierr
  END SUBROUTINE Boundsch

END MODULE radialsch
