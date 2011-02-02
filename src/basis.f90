MODULE basis
  USE atomdata
  USE aeatom
  USE calcpotential
  USE Globalmath
  USE gridmod
  USE pseudo
  USE radialsch
  USE radialsr
  USE anderson_realmix

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setbasis(Grid,Pot,Orbit,FC,PAW,ifinput)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    TYPE(FCInfo), INTENT(IN) :: FC
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER, INTENT(IN) :: ifinput


    INTEGER :: n,irc,norbit,nbase,l,lmax,mxbase
    INTEGER :: i,j,k,io,ok,nbl,nr,nodes,ib,loop,niter,iter
    REAL(8) :: h,rc,q00,energy,rat,delta,thisconv,qeff,tq
    REAL(8) :: ecoul,etxc,eexc
    CHARACTER(1) :: answer
    REAL(8), POINTER  :: r(:)
    CHARACTER(132) :: inputline

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    irc=PAW%irc
    nr=irc+20
    rc=PAW%rc
    lmax=PAW%lmax

    norbit=Orbit%norbit
    nbase=PAW%nbase
    mxbase=nbase+5*max(1,PAW%lmax)
    !
    nbase=0
    WRITE(6,*) '  basis functions:'
    WRITE(6,*)' No.   n     l         energy         occ   '
    DO l=0,lmax
       nbl=0
       DO io=1,norbit    ! cycle through all configuration
          IF (Orbit%l(io).EQ.l.AND.(.NOT.FC%iscore(io))) THEN
             nbl=nbl+1
             nbase=nbase+1
             PAW%np(nbase)=Orbit%np(io)
             PAW%l(nbase)=l
             PAW%eig(nbase)=Orbit%eig(io)
             PAW%occ(nbase)=Orbit%occ(io)
             PAW%phi(:,nbase)=Orbit%wfn(:,io)
             WRITE(6,'(3i6,1p,2e15.6)') nbase,PAW%np(nbase),l,             &
&                 PAW%eig(nbase),PAW%occ(nbase)
          ENDIF
       ENDDO
       generalizedloop: DO
          WRITE(6,*) 'For l = ',l,' there are currently ',nbl,&
&              'basis functions'
          WRITE(6,*) 'enter y to add additional functions or n to ' &
&              ,'go to next l'
          READ(5,'(a)') inputline
          WRITE(ifinput,'(a)') TRIM(inputline)
          READ(inputline,*) answer
          IF (answer.NE.'y') EXIT generalizedloop
          WRITE(6,*) 'enter energy for generalized function'
          READ(5,'(a)') inputline
          WRITE(ifinput,'(a)') TRIM(inputline)
          READ(inputline,*) energy
          IF(energy.LT.0.d0) THEN
             WRITE(6,*) 'energy is negative',energy,' -- WARNING WARNING !!!'
          ENDIF
          nbase=nbase+1
          IF (nbase > mxbase ) THEN
             WRITE(6,*) 'Error in  setbasis -- too many functions ', nbase,mxbase
             STOP
          ENDIF
          PAW%l(nbase)=l
          PAW%np(nbase)=999
          PAW%eig(nbase)=energy
          PAW%occ(nbase)=0.d0
          PAW%phi(1:n,nbase)=0.d0
          if (scalarrelativistic) then
             CALL unboundsr(Grid,Pot,n,l,energy,PAW%phi(:,nbase),nodes)
          else
             CALL unboundsch(Grid,Pot,n,l,energy,PAW%phi(:,nbase),nodes)
          endif
          rat=MAX(ABS(PAW%phi(irc,nbase)),ABS(PAW%phi(irc+1,nbase)))
          rat=DSIGN(rat,PAW%phi(irc,nbase))
          PAW%phi(1:n,nbase)=PAW%phi(1:n,nbase)/rat
          WRITE(6,'(3i6,1p,2e15.6)') nbase,PAW%np(nbase),l,             &
&              PAW%eig(nbase),PAW%occ(nbase)
          nbl=nbl+1
       ENDDO generalizedloop
       !
    ENDDO   ! end lmax loop

    WRITE(6,*) 'completed phi basis with ',nbase,' functions '
    PAW%nbase=nbase     ! reset nbase

  END SUBROUTINE setbasis

  !**************************************************************************
  !  Program to generate atomic basis functions
  !    Version using Bloechl's form of projector and orthogonalization procedure
  !     At the end of this subroutine, the basis functions and projectors are
  !     orthogonalized with a Gram-Schmidt like scheme
  !**************************************************************************
  SUBROUTINE makebasis_bloechl(Grid,AEPot,PAW,ifinput,option)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: AEPOT
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER,INTENT(IN) :: ifinput,option

    INTEGER :: n,irc,nr
    INTEGER :: i,j,k,io,ok,lmax
    REAL(8) :: h,rc
    REAL(8), ALLOCATABLE :: denout(:),tmp(:),VNC(:)
    TYPE(PotentialInfo), TARGET:: PS
    REAL(8), POINTER  :: r(:)

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    irc=PAW%irc
    nr=irc+20
    rc=PAW%rc
    lmax=PAW%lmax
    PS%nz=0

    call InitPotential(Grid,PS)
    ALLOCATE(tmp(n),VNC(n),stat=ok)
    IF (ok /= 0 ) THEN
       WRITE(6,*) 'Error in allocating  denout in makebasis ',n,ok
       STOP
    ENDIF

    ! find basis functions
    PS%rv=PAW%rveff
    call zeropot(Grid,PS%rv,PS%v0,PS%v0p)
    !write(6,*) 'VNC  v0 ', PS%v0,PS%v0p,PS%rv(5)
    CALL formprojectors(Grid,AEPot,PS,PAW,ifinput,option)

    DEALLOCATE(tmp,VNC,PS%den,PS%rv)

  END SUBROUTINE makebasis_bloechl


  !**************************************************************************
  !  Program to generate atomic basis functions
  !
  !   1) Pseudization of partial waves:
  !        - simple polynom scheme                                            [optps=1]
  !                   r^(l+1).Sum[Ci.r^2i]  0<=i<=4
  !   OR   - ultrasoft polynom scheme                                         [optps=2]
  !                   r^(l+1).{Sum[Ci.r^2i]+Sum[Cj.r^2j]}  0<=i<=3
  !                           3<j adjusted using Fourier filtering
  !   OR   - RRKJ scheme with 2 Bessel functions (PHYS REV B 41,1227 (1990))  [optps=3]
  !
  !   2) Build and orthogonalization of projectors
  !        - Vanderbilt generation method (PHYS REV B 41,7892 (1990))  [optorth=0]
  !   OR   - Gram-Schmidt like sheme                                   [optorth=1]
  !**************************************************************************
  SUBROUTINE makebasis_custom(Grid,Pot,PAW,ifinput,optps,optorth,pdeg,qcut)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER, INTENT(IN) :: ifinput,optps,optorth,pdeg
    REAL(8), INTENT(IN) :: qcut

    INTEGER :: i,j,k,l,io,jo,ok,lmax,nbase,n,irc,irc_vloc,nr,np,thisrc
    INTEGER :: icount,jcount,istart,ifinish,ibase,jbase
    REAL(8) :: choice,rc,xx,yy,gg,g,gp,gpp,al(2),ql(2)
    INTEGER, ALLOCATABLE :: omap(:)
    REAL(8), ALLOCATABLE :: VNC(:),Ci(:),aa(:,:),ai(:,:)
    TYPE(PotentialInfo), TARGET:: PS
    REAL(8), POINTER  :: r(:)
    CHARACTER(132) :: inputline

    if (optps<1.or.optps>3) stop 'bug: error calling makebasis_custom routine'
    if (optorth<0.or.optorth>1) stop 'bug: error calling makebasis_custom routine'

    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    irc_vloc=PAW%irc_vloc
    nbase=PAW%nbase
    lmax=PAW%lmax

    np=5;if (optps==2) np=pdeg+1
    if (optps==1.or.optps==2) allocate(Ci(np))

  ! Set screened local pseudopotential
    allocate(VNC(n),stat=i)
    if (i/=0) stop 'allocation error in makebasis_vanderbilt'
    VNC(2:n)=PAW%rveff(2:n)/r(2:n)
    call extrapolate(Grid,VNC)

    write(6,*) 'For each of the following basis functions enter rc'

  ! Loop on basis elements
    do io=1,nbase
       l=PAW%l(io)

     ! Read matching radius
       write(6,'(3i5,1pe15.7)') io,PAW%np(io),PAW%l(io),PAW%eig(io)
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       read(inputline,*) rc
       thisrc=FindGridIndex(Grid,rc)
       thisrc=MIN(thisrc,irc)       ! make sure rc<total rc
       rc=r(thisrc)
       write(6,*) 'rc for this wfn', rc
       if (thisrc<3.or.thisrc>irc.or. &
&          (optps==1.and.thisrc>n-3).or. &
&          (optps==2.and.thisrc>n-6)) then
          write(6,*) 'rc out of range', thisrc,n,irc
          stop
       endif

     ! Find partial wave pseudization
       if (optps==1) then
        call pspolyn(PAW%phi(:,io),Ci,r,l,np,thisrc,n)
       else if (optps==2) then
        call psuspolyn(PAW%phi(:,io),Ci,r,l,np,thisrc,n,qcut)
       else if (optps==3) then
        call psbes(PAW%phi(:,io),al,ql,Grid,l,thisrc,n)
       endif

     ! Compute pseudized partial wave and unnormalized projector
       PAW%tphi(:,io)=PAW%phi(:,io)
       PAW%tp(:,io)=0.d0
       if (optps==1.or.optps==2) then
        do i=1,thisrc-1
         xx=r(i)*r(i)
         PAW%tphi(i,io)=Ci(1)+Ci(2)*xx
         gp=2.d0*Ci(2)
         gpp=2.d0*Ci(2)
         do j=3,np
          PAW%tphi(i,io)=PAW%tphi(i,io)+Ci(j)*xx**(j-1)
          gp=gp+dble(2*j-2)*Ci( j)*xx**(j-2)
          gpp=gpp+dble((2*j-2)*(2*j-3))*Ci(j)*xx**(j-2)
         enddo
         PAW%tphi(i,io)=PAW%tphi(i,io)*r(i)**(l+1)
         PAW%tp(i,io)=(dble(2*(l+1))*gp+gpp)*r(i)**(l+1)+(PAW%eig(io)-VNC(i))*PAW%tphi(i,io)
        enddo
       else if (optps==3) then
        PAW%tphi(1,io)=0.d0
        do i=2,thisrc-1
         xx=ql(1)*r(i)
         call jbessel(g,gp,gpp,l,2,xx)
         PAW%tphi(i,io)=al(1)*g*r(i)
         gg=al(1)*(2.d0*ql(1)*gp+ql(1)*xx*gpp)
         xx=ql(2)*r(i)
         call jbessel(g,gp,gpp,l,2,xx)
         PAW%tphi(i,io)=PAW%tphi(i,io)+al(2)*g*r(i)
         gg=gg+al(2)*(2.d0*ql(2)*gp+ql(2)*xx*gpp)
         PAW%tp(i,io)=(PAW%eig(io)-VNC(i)-dble(l*(l+1))/(r(i)**2))*PAW%tphi(i,io)+gg
        enddo
       endif

       if (thisrc<irc_vloc) then
         do i=thisrc,irc_vloc-1
             gpp=Gsecondderiv(Grid,i,PAW%tphi(:,io))
             PAW%tp(i,io)=(PAW%eig(io)-VNC(i)-dble(l*(l+1))/(r(i)**2))*PAW%tphi(i,io) &
&                    +gpp
         enddo
       endif

     ! If Gram-Schmidt orthogonalization, form normalized projector
       if (optorth==1) then
        xx=overlap(Grid,PAW%tp(:,io),PAW%tphi(:,io),1,max(thisrc,irc_vloc))
        PAW%tp(:,io)=PAW%tp(:,io)/xx
       endif

    enddo  !nbase

    deallocate(VNC);if (optps==1.or.optps==2) deallocate(Ci)

  ! Form orthogonalized projector functions
    do io=1,nbase
       PAW%ophi(:,io)=PAW%phi(:,io)
       PAW%otphi(:,io)=PAW%tphi(:,io)
       PAW%Kop(1,io)=0
       PAW%Kop(2:n,io)=(PAW%eig(io)-Pot%rv(2:n)/Grid%r(2:n))*PAW%phi(2:n,io)
       if (optorth==1) PAW%otp(:,io)=PAW%tp(:,io)
    enddo

  ! First option : VANDERBILT SCHEME
    if (optorth==0) then
     do l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       write(6,*) 'For l = ', l,icount,' basis functions'
       if (icount==0) cycle
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0
       do io=1,nbase
        if (PAW%l(io)==l) then
          icount=icount+1
          omap(icount)=io
        endif
       enddo
       do i=1,icount
         io=omap(i)
         do j=1,icount
           jo=omap(j)
           aa(i,j)=overlap(Grid,PAW%otphi(:,io),PAW%tp(:,jo),1,irc)
         enddo
       enddo
       ai=aa;call minverse(ai,icount,icount,icount)

       do i=1,icount
         io=omap(i)
         PAW%ck(io)=ai(i,i)
         PAW%otp(:,io)=0
         do j=1,icount
           jo=omap(j)
           PAW%otp(:,io)=PAW%otp(:,io)+PAW%tp(:,jo)*ai(j,i)
         enddo
       enddo
       deallocate(aa,ai,omap)
     enddo

  ! Second option : GRAM-SCHMIDT SCHEME
    else if (optorth==1) then
     DO l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       if (icount==0) cycle
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       icount=0
       DO io=1,nbase
         IF (PAW%l(io)==l) THEN
           IF  (icount==0) istart=io
           IF  (icount>=0) ifinish=io
           icount=icount+1;omap(icount)=io
         ENDIF
       ENDDO
       DO ibase=istart,ifinish
         DO jbase=istart,ibase
           IF (jbase.LT.ibase) THEN
             xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
             yy=overlap(Grid,PAW%otphi(:,jbase),PAW%otp(:,ibase),1,irc)
             PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)-PAW%ophi(1:n,jbase)*xx
             PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)-PAW%Kop(1:n,jbase)*xx
             PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)-PAW%otphi(1:n,jbase)*xx
             PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)-PAW%otp(1:n,jbase)*yy
             aa(ibase-istart+1,jbase-istart+1)=xx
             aa(jbase-istart+1,ibase-istart+1)=xx
           ELSE IF (jbase.EQ.ibase) THEN
             xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
             choice=1.d0/SQRT(ABS(xx))
             PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)*DSIGN(choice,xx)
             PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)*choice
             PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)*choice
             PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)*choice
             aa(ibase-istart+1,ibase-istart+1)=xx
           ENDIF
         ENDDO
       ENDDO
       ai=aa;call minverse(aa,icount,icount,icount)
       do i=1,icount
        io=omap(i);PAW%ck(io)=ai(i,i)
       enddo
      deallocate(aa,ai,omap)
     ENDDO
    endif

  END SUBROUTINE makebasis_custom


!  subroutine not presently used
!  SUBROUTINE SCbasis(Grid,FC,nz,PAW)
!    TYPE(GridInfo), INTENT(IN) :: Grid
!    TYPE(FCInfo), INTENT(IN) :: FC
!    INTEGER, INTENT(IN) :: nz
!    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
!
!    INTEGER :: iv,niter,loop,i,n,nbase,io,index
!    REAL(8) :: q00,tq,h,conv,ecoul,etxc,eexc,delta,v1,v2,v3,v4
!    REAL(8), ALLOCATABLE :: denout(:),dum(:)
!    REAL(8), PARAMETER  :: conv1=1.d13,conv2=2.d13,conv3=3.d13,conv4=4.d13
!    TYPE(PotentialInfo) :: PS
!
!    h=Grid%h
!    n=Grid%n
!    nbase=PAW%nbase
!
!    v1=conv1;v2=conv2;v3=conv3;v4=conv4
!
!    DO
!       WRITE(6,*) 'Input index of basis function for setting vloc'
!         READ(5,*) index
!         IF (index>=1.AND.index<=nbase) EXIT
!      ENDDO
!
!      WRITE(PAW%Vloc_description,&
!&          '("Fixed shape Vloc with amplitude set for basis index",i4)')&
!&          index
!      WRITE(6,*) PAW%Vloc_description
!
!      ALLOCATE(PS%rv(n),PS%den(n),denout(n),dum(n),stat=i)
!      IF (i/=0) THEN
!         WRITE(6,*) 'Allocation error in SCbasis ', n,i
!         STOP
!      ENDIF
!
!      PS%nz=0
!      PAW%vloc=0
!      DO io=1,nbase
!         PAW%tphi(:,io)=PAW%phi(:,io)
!      ENDDO
!
!      PAW%den=FC%coreden; PAW%tden=PAW%tcore
!      DO io=1,nbase
!         PAW%den(:)=PAW%den(:)+PAW%occ(io)*(PAW%phi(:,io))**2
!         PAW%tden(:)=PAW%tden(:)+PAW%occ(io)*(PAW%tphi(:,io))**2
!      ENDDO
!
!      DO iv=1,2
!
!         conv=nbase*1.d-7
!         niter=1000
!         loop=0
!         delta=1.d10
!         denout=PAW%tden
!
!         tphiloop: DO
!            loop=loop+1
!            IF(loop > niter) THEN
!               WRITE(6,*) 'SCbasis pgm terminating in tphiloop ', loop
!               STOP
!            ENDIF
!
!
!            tq=integrator(Grid,PAW%tden)
!            WRITE(6,*) 'tq = ',tq
!            dum=PAW%den-PAW%tden
!            q00=-nz + integrator(Grid,dum,1,PAW%irc)
!            WRITE(6,*) 'q00 = ',q00
!
!            PS%q=tq
!            PS%den=PAW%tden
!            CALL potential(Grid,PS,ecoul,etxc,eexc)
!            PAW%rveff(:)=PS%rv(:)+PAW%vloc(:)*Grid%r(:)+q00*PAW%hatpot(:)
!
!            IF (loop>=4) THEN
!               IF ((delta.LT.1.e-11).OR.(.NOT.(v4.LE.v3.AND.v3.LE.v2 &
!&                   .AND.v2.LE.v1).AND.v4.LE.conv)) THEN
!                  !
!                  !  converged result
!                  !
!                  WRITE(6,*) '  SCprojectors converged in',loop,' iterations'
!
!                  EXIT tphiloop
!               ENDIF
!            ENDIF
!
!
!            CALL makebasis(Grid,PAW)
!
!            PAW%tden=PAW%tcore
!            DO io=1,nbase
!               PAW%tden(:)=PAW%tden(:)+PAW%occ(io)*(PAW%tphi(:,io))**2
!            ENDDO
!
!            dum=denout-PAW%tden
!            delta=SQRT(DOT_PRODUCT(dum(1:n),dum(1:n)))
!            CALL shift4(v1,v2,v3,v4,delta)
!            WRITE(6,*) 'loop = ', loop,' delta = ', delta
!            denout=PAW%tden
!
!         ENDDO tphiloop
!
!         IF (iv==1) THEN
!            WRITE(6,*) 'Completed first tphiloop', PAW%ck(1:nbase)
!            WRITE(6,*) 'Resetting vloc with vlocfac == ', -PAW%ck(index)
!            PAW%vloc=-PAW%ck(index)*PAW%projshape
!         ENDIF
!
!      ENDDO   !iv loop
!      DEALLOCATE(PS%den,PS%rv,denout,dum)
!    END SUBROUTINE SCbasis
!
    SUBROUTINE FindVlocfromVeff(Grid,FC,nz,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(FCInfo), INTENT(IN) :: FC
      REAL(8), INTENT(IN) :: nz
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      TYPE(PotentialInfo), TARGET:: PS
      REAL(8), POINTER  :: r(:)
      REAL(8) :: h,qeff,tq,rat,q00,ecoul,etxc,eexc
      INTEGER :: n,i,irc,io,nbase
      REAL(8), allocatable :: d(:),v(:),dd(:),vv(:) ! temporary arrays for abinit

      n=Grid%n
      allocate(d(n),dd(n),v(n),vv(n),STAT=i)
          if (i /= 0) then
             write(6,*) 'Error in allocating abinit arrays',n
             stop
          endif
      h=Grid%h
      r=>Grid%r
      irc=max(PAW%irc,PAW%irc_shap,PAW%irc_vloc,PAW%irc_core)
      d=FC%coreden
      dd=PAW%tcore
      PAW%den=0.d0
      PAW%tden=0.d0
      nbase=PAW%nbase
      DO io=1,nbase
         IF (PAW%occ(io) > 1.d-8)  THEN
            DO i=1,n
               PAW%den(i)= PAW%den(i)+PAW%occ(io)*(PAW%phi(i,io)**2)
               PAW%tden(i)= PAW%tden(i)+PAW%occ(io)*(PAW%tphi(i,io)**2)
            ENDDO
         ENDIF
      ENDDO

      d=d+PAW%den
      dd=dd+PAW%tden

      ALLOCATE(PS%den(n),PS%rv(n),stat=i)
      IF (i /= 0) THEN
         WRITE(6,*) 'Error in FindVlocfromVeff -- ',n,i
         STOP
      ENDIF

      PS%den=d-dd

      qeff=-nz+integrator(Grid,PS%den,1,irc)
      WRITE(6,*) 'qeff = ',qeff
      tq=integrator(Grid,dd)
      WRITE(6,*) 'tq  = ',tq
      q00=qeff
      WRITE(6,*) 'q00 = ',q00
      !  generate effective potential
      !
      PS%den=dd; PS%q=tq; PS%nz=0
      CALL potential(Grid,PS,ecoul,etxc,eexc)
      !
      rat=0
      DO i=2,n
         PAW%vloc(i)=(PAW%rveff(i)-PS%rv(i)-q00*PAW%hatpot(i))/r(i)
         IF (i>=irc) rat=rat+ABS(PAW%vloc(i))
      ENDDO
      call extrapolate(Grid,PAW%vloc)
      PAW%vloc(irc:n)=0
      WRITE(6,*) 'Error in vloc -- ', rat

   !  Construct ionic local potential for abinit from screened
   !     pseudopotential
   !     in addition to ionic unscreening include hathat density
   !        in exchange-correlation functional
   !     also generate version without hathat density in vxc
      d=0;dd=0
      DO io=1,nbase
         IF (PAW%occ(io) > 1.d-8)  THEN
            DO i=1,n
               d(i)= d(i)+PAW%occ(io)*(PAW%phi(i,io)**2)
               dd(i)= dd(i)+PAW%occ(io)*(PAW%tphi(i,io)**2)
            ENDDO
         ENDIF
      ENDDO

      d=d-dd
      tq=integrator(Grid,d,1,irc)
      write(6,*) ' abinit tq = ', tq
      d=dd+tq*PAW%hatden
      write(6,*) ' check valence ', FC%zvale,integrator(Grid,dd)
      !   tq=FC%zvale
      CALL poisson(Grid,q00,d,v,rat)   ! Coul(\tilde(n) + \hat(n))
      write(6,*) 'Check poisson (should equal ',FC%zvale,') ',q00
      d=d+PAW%tcore
      CALL exch(Grid,d,vv,etxc,eexc)    ! with \hat(n)
      d=vv                      ! store vxc(\tilde(n)+\tilde(n_c)+\hat(n))
      dd=dd+PAW%tcore
      CALL exch(Grid,dd,vv,etxc,eexc)   !store vxc(\tilde(n)+\tilde(n_c))

      do i=2,irc-1
         PAW%abinitvloc(i)=(PAW%rveff(i)-v(i)-d(i))/r(i)  ! with nhat
         PAW%abinitnohat(i)=(PAW%rveff(i)-v(i)-vv(i))/r(i)  ! without nhat
      enddo
      call extrapolate(Grid,PAW%abinitvloc)
      call extrapolate(Grid,PAW%abinitnohat)

      dd=FC%coreden-PAW%tcore
      qeff=-nz+integrator(Grid,dd,1,PAW%irc_core)
      dd=PAW%tcore+qeff*PAW%hatden
      CALL poisson(Grid,q00,dd,vv,rat)
      write(6,*) 'Check poisson (should equal ',-FC%zvale,') ',q00

      do i=irc,n
         PAW%abinitvloc(i)=vv(i)/r(i)  ! still in Rydberg units
         PAW%abinitnohat(i)=vv(i)/r(i)  ! still in Rydberg units
      enddo

      open(123,file='compare.abinit', form='formatted')
       do i=2,n
          write(123,'(1p,5e16.7)') r(i),PAW%rveff(i)/r(i),&
&                 PAW%vloc(i),PAW%abinitvloc(i),PAW%abinitnohat(i)
       enddo
      close(123)

      DEALLOCATE(PS%rv,PS%den,d,dd,v,vv)
    END SUBROUTINE FindVlocfromVeff

    !*************************************************************************
    !  program to generate projector functions for Blochl's paw formalism
    !    starting with smooth functions
    !    for every basis function phi, choose smooth function with
    !     form for r<rc: (r**(l+1))*sum(n)*(cn*(r**(2*n)))
    !       where, rc==r(iiirc)
    !  phi(i,io), eval(io), tocc(io)  original basis function, energy, occupancy
    !  tphi(i,io) smooth basis function corresponding to phi(i,io)
    !  tp(i,io)=constant*normalized Harmonic Oscillator function
    !  tp == 0 for r>r(iiirc)
    !  functions defined to be identically zero for r>r(iiirc))
    !*************************************************************************
    SUBROUTINE formprojectors(Grid,AEPot,PS,PAW,ifinput,option)
      TYPE(GridInfo),  INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: AEPot,PS
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER,INTENT(IN) :: ifinput,option

      INTEGER :: nbase,lmax,l,io,irc,wantednodes,nb,n,icount
      INTEGER ::  istart,ifinish,ibase,jbase
      REAL(8) :: h,xx,yy,choice,rc
      INTEGER,allocatable :: irc_io(:)
      CHARACTER(132) :: inputline

      n=Grid%n
      h=Grid%h
      lmax=PAW%lmax
      nbase=PAW%nbase
      irc=PAW%irc
      !
      !   form projector functions for each l
      !

      if (option==1) then
       allocate(irc_io(nbase))
       write(6,*) 'For each of the following basis functions enter rc'
       do io=1,nbase
        write(6,'(3i5,1pe15.7)') io,PAW%np(io),PAW%l(io),PAW%eig(io)
        READ(5,'(a)') inputline
        WRITE(ifinput,'(a)') TRIM(inputline)
        read(inputline,*) rc
        irc_io(io)=FindGridIndex(Grid,rc)
        rc=Grid%r(irc_io(io))
        write(6,*) 'rc for this wfn', rc
        if(irc_io(io)>PAW%irc) then
         write(6,*) 'rc out of range', irc_io(io),n,PAW%irc
         stop
        endif
       enddo
      endif

      DO l=0,lmax
         icount=0
         DO io=1,nbase
            IF (PAW%l(io)==l) THEN
               IF  (icount==0) istart=io
               IF  (icount >= 0) ifinish=io
               icount=icount+1
               wantednodes=icount-1
               ! form unorthonormalized projector functions tp
               WRITE(6,*) '******* projector for l = ',l
               if (option==1) then
                CALL  bsolv(Grid,PS,PAW,io,wantednodes,irc_io(io))
               else
                CALL  bsolv(Grid,PS,PAW,io,wantednodes)
               endif
               PAW%ophi(:,io)=PAW%phi(:,io)
               PAW%otphi(:,io)=PAW%tphi(:,io)
               PAW%otp(:,io)=PAW%tp(:,io)
               PAW%Kop(1,io)=0
               PAW%Kop(2:n,io)=(PAW%eig(io)-AEPot%rv(2:n)/Grid%r(2:n))&
&                       *PAW%phi(2:n,io)
            ENDIF
         ENDDO
         !write(6,*) 'orthnormalization'
         !write(6,*) 'start orthogonalization',istart,ifinish
         DO ibase=istart,ifinish
            DO jbase=istart,ibase
               IF (jbase.LT.ibase) THEN
                  xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
                  yy=overlap(Grid,PAW%otphi(:,jbase),PAW%otp(:,ibase),1,irc)
                  !write(6,*) 'before',jbase,ibase,xx,yy
                  PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)-PAW%ophi(1:n,jbase)*xx
                  PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)-PAW%Kop(1:n,jbase)*xx
                  PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)-PAW%otphi(1:n,jbase)*xx
                  PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)-PAW%otp(1:n,jbase)*yy
               ELSE IF (jbase.EQ.ibase) THEN
                  xx=overlap(Grid,PAW%otp(:,ibase),PAW%otphi(:,ibase),1,irc)
                !write(6,*) 'before',jbase,ibase,xx
                  choice=1.d0/SQRT(ABS(xx))
                  PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)*DSIGN(choice,xx)
                  PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)*choice
                  PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)*choice
                  PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)*choice
               ENDIF
            ENDDO
         ENDDO

      ENDDO

      if (option==1) deallocate(irc_io)

    END SUBROUTINE formprojectors
    !*************************************************************************
    !  on input tphi=phi
    !  on output tphi recalculated for r<nrc*h
    !  l = angular momentum quantum number
    !  nodes = the number of desired nodes in tp and tphi
    !  shape function
    !*************************************************************************
    SUBROUTINE bsolv(Grid,Pot,PAW,io,wantednodes,irc_io)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER, INTENT(IN) :: io,wantednodes
      INTEGER, OPTIONAL :: irc_io

      REAL(8), ALLOCATABLE :: f(:),chi(:),fakerv(:)
      REAL(8), ALLOCATABLE :: tphi(:),tp(:)
      REAL(8), POINTER :: r(:),rv(:)
      REAL(8) :: h,en,v0,v0p,ol,h2,eh2,veh2,xnorm,rc_io
      INTEGER ::  l,n,num,nrc,i,j,k,ii,jj,iter,irc,nodes,ok
      REAL(8) :: cpmin,cpmax,zeroval
      REAL(8) :: cp,cph2,del,deriv,tderiv
      REAL(8), PARAMETER :: small=1.d-10,step=1.d0
      INTEGER, PARAMETER :: mxiter=1000

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      rv=>Pot%rv
      if (present(irc_io)) then
       irc=irc_io
      else
       irc=PAW%irc
      endif
      ALLOCATE(fakerv(n),f(n),chi(n),tphi(n),tp(n),stat=ok)
      IF (ok /= 0) THEN
         WRITE(6,*) 'Error in bsolv allocation ', n,ok
         STOP
      ENDIF
      tphi=PAW%phi(:,io)
      en=PAW%eig(io)
      l=PAW%l(io)
      deriv=(Gfirstderiv(Grid,irc,tphi))/tphi(irc) ! log derivative
      tderiv=0    ! initial log derive
      cp=0         ! initial projector constant
      del=1000     ! initial error
      cpmin=-100
      cpmax=100
      if (present(irc_io)) then
       rc_io=r(irc_io)
       f(1)=1.d0;f(irc_io:n)=0.d0
       DO i=2,irc_io-1
        f(i)=(SIN(pi*r(i)/rc_io)/(pi*r(i)/rc_io))**2
       ENDDO
      else
       DO i=1,n
          f(i)=PAW%projshape(i)
       ENDDO
      endif

      WRITE(6,*) 'in bsolv -- l, en, n',l,en,wantednodes

      iter=0

      DO

         iter=iter+1

         WRITE(6,*) 'bsolv iter cp',iter,cp

         fakerv=rv-cp*f*Grid%r
         call zeropot(Grid,fakerv,v0,v0p)
         ! initialize chi
         chi=0
         chi(2)=wfninit(0.d0,l,v0,v0p,en,Grid%r(2))
         zeroval=0
         if (l==1) zeroval=2

         call forward_numerov(Grid,l,irc+5,en,fakerv,zeroval,chi,nodes)

         WRITE(6,'("iter nodes cp cpmin cpmax",2i5,1p,3e15.7)')&
&             iter,nodes,cp,cpmin,cpmax
         IF (nodes.EQ.wantednodes) THEN
            tderiv=(Gfirstderiv(Grid,irc,chi))/chi(irc) ! log derivative
            tp=0
            tp(1:irc)=chi(1:irc)/chi(irc)
            chi(1:irc)=f(1:irc)*(tp(1:irc))**2
            xnorm=integrator(Grid,chi(1:irc),1,irc)
            del=(tderiv-deriv)/xnorm
            WRITE(6,*) 'iter nodes del', iter,nodes,del

            IF (ABS(del).LT.small) EXIT

            IF (iter.GE.mxiter) THEN
               WRITE(6,*)' terminating projector',iter
               STOP
            ENDIF

            IF (ABS(del).GT.step) del=DSIGN(step,del)
            cp=cp+del
            IF (cp.GT.cpmax) cp=cpmax-ranx()*step
            IF (cp.LT.cpmin) cp=cpmin+ranx()*step

         ELSE IF(nodes.GT.wantednodes) THEN
            cpmax=cp
            cp=cpmax-ranx()*step
         ELSE IF(nodes.LT.wantednodes) THEN
            cpmin=cp
            cp=cpmin+ranx()*step
         ENDIF

      ENDDO

      tphi(1:irc-1)=tphi(irc)*tp(1:irc-1)
      tphi(irc:n)=PAW%phi(irc:n,io)
      tp(1:irc)=f(1:irc)*tphi(1:irc)
      chi(1:irc)=tphi(1:irc)*tp(1:irc)
      xnorm=integrator(Grid,chi(1:irc),1,irc)
      WRITE(6,*) 'normalization for projector l,n=',l,nodes,xnorm
      tp(1:irc)=tp(1:irc)/xnorm
      tp(irc+1:n)=0

      PAW%tphi(:,io)=tphi
      PAW%tp(:,io)=tp
      PAW%ck(io)=cp

      WRITE(6,*) 'completed bsolv',io,cp

      DEALLOCATE(fakerv,f,chi,tphi,tp)
    END SUBROUTINE bsolv

    SUBROUTINE trunk(Grid,f,rstart,rend)
      TYPE (GridInfo), INTENT(IN) :: Grid
      REAL(8),  INTENT(INOUT) :: f(:)
      REAL(8),  INTENT(IN) :: rstart,rend

      INTEGER :: i,n
      REAL(8) :: r,delta,arg

      delta=rend-rstart

      DO i=1,Grid%n
         r=Grid%r(i)
         IF (r>rstart .AND. r<=rend) THEN
            arg=pi*(r-rstart)/delta
            f(i)=f(i)*(SIN(arg)/arg)**2
         ENDIF
         IF (r>rend) f(i)=0
      ENDDO
    END SUBROUTINE trunk

    !***********************************************************************
    !  program to calculate Fourier transform of paw product tp*tphi
    !   in order to get an idea of the convergence
    !   and output them to files tprod.l
    !***********************************************************************
    SUBROUTINE ftprod(Grid,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,tphij
      INTEGER :: i,ib,l,iq,n,irc
      CHARACTER(4) flnm
      !
      WRITE(6,*) 'calculating Fourier transforms of tp*tphi products  ',&
&          'For bound states only '

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      irc=PAW%irc
      ALLOCATE(q(nq),dum(n),dum1(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in allocating space in ftprod',n,nq,i
         STOP
      ENDIF

      dq=qmax/nq
      DO ib=1,PAW%nbase
         IF (PAW%eig(ib)< 0.d0) THEN
            l=PAW%l(ib)
            CALL mkname(ib,flnm)
            OPEN(55,file='tprod.'//TRIM(flnm),form='formatted')
            DO iq=1,nq
               q(iq)=iq*dq
               DO i=1,n
                  dum(i)=q(iq)*r(i)
               ENDDO
               CALL sphbes(l,n,dum)
               DO i=1,n
                  dum1(i)=dum(i)*r(i)*PAW%otphi(i,ib)
                  IF (i.LE.irc) dum(i)=dum(i)*r(i)*PAW%otp(i,ib)
               ENDDO
               tphij=integrator(Grid,dum1(1:n))*&
&                    integrator(Grid,dum(1:irc),1,irc)

               WRITE(55,'(1p,2e16.7)') q(iq),tphij
            ENDDO
            CLOSE(55)
         ENDIF  !(e<=0)
      ENDDO  ! ib
      DEALLOCATE(q,dum,dum1)
    END SUBROUTINE ftprod
    !***********************************************************************
    !  program to calculate Fourier transform of hatpot functions
    !   in order to get an idea of the convergence
    !   and output them to files hatpot.l
    !***********************************************************************
    SUBROUTINE fthatpot(Grid,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:),dum2(:),dum3(:),arg(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,fthatp,fthatd,fthattd
      INTEGER :: i,ib,l,iq,n,irc,ll
      CHARACTER(4) flnm
      !
      WRITE(6,*) 'calculating Fourier transforms of hatpot for  each l '

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      irc=PAW%irc
      ALLOCATE(q(nq),dum(n),dum1(n),dum2(n),dum3(n),arg(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in allocating space in fthatpot',n,nq,i
         STOP
      ENDIF

      dq=qmax/nq
      ll=2*MAXVAL(PAW%l(:))
      DO l=0,ll
         CALL mkname(l,flnm)
         OPEN(55,file='fthatpot.'//TRIM(flnm),form='formatted')
         CALL hatL(Grid,PAW,l,dum2)
         DO iq=1,nq
            q(iq)=iq*dq
            DO i=1,n
               dum(i)=q(iq)*r(i)
            ENDDO
            CALL sphbes(l,n,dum)
            DO i=1,n
               dum3(i)=dum(i)*dum2(i)
            ENDDO
            fthatd=integrator(Grid,dum3(1:n))
            fthatp=8*PI*fthatd/q(iq)**2

            fthattd=0
            IF(l==0) THEN
               arg(1:n)=dum(1:n)*PAW%tden(1:n)
               fthattd=integrator(Grid,arg)
            ENDIF

            WRITE(55,'(1p,5e16.7)') q(iq),fthatp,fthatd,fthatp*fthatd,fthatp*fthattd
         ENDDO
         CLOSE(55)
      ENDDO  ! l
      DEALLOCATE(q,dum,dum1,dum2,dum3,arg)
    END SUBROUTINE fthatpot

    !***********************************************************************
    !  program to calculate Fourier transform of tphi*q^2
    !   in order to get an idea of the convergence of kinetic energy
    !   and output them to files ftkin.l
    !***********************************************************************
    SUBROUTINE ftkin(Grid,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,kin
      INTEGER :: i,ib,l,iq,n,irc
      CHARACTER(4) flnm
      !
      WRITE(6,*) 'calculating Fourier transforms of tphi  ',&
&          'For bound states only '

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      ALLOCATE(q(nq),dum(n),dum1(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in allocating space in ftkin',n,nq,i
         STOP
      ENDIF

      dq=qmax/nq
      DO ib=1,PAW%nbase
         IF (PAW%eig(ib)<=0.d0) THEN
            l=PAW%l(ib)
            CALL mkname(ib,flnm)
            OPEN(55,file='ftkin.'//TRIM(flnm),form='formatted')
            DO iq=1,nq
               q(iq)=iq*dq
               DO i=1,n
                  dum(i)=q(iq)*r(i)
               ENDDO
               CALL sphbes(l,n,dum)
               DO i=1,n
                  dum1(i)=dum(i)*r(i)*PAW%tphi(i,ib)
               ENDDO
               kin=integrator(Grid,dum1(1:n))*(q(iq))
               kin=kin*kin

               WRITE(55,'(1p,2e16.7)') q(iq),kin
            ENDDO
            CLOSE(55)
         ENDIF  !(e<=0)
      ENDDO  ! ib
      DEALLOCATE(q,dum,dum1)
    END SUBROUTINE ftkin

    !***********************************************************************
    !  program to calculate Fourier transform of vloc and  tden
    !   in order to get an idea of their convergence
    !   and output them to files ftvloc
    !***********************************************************************
    SUBROUTINE ftvloc(Grid,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,vloc,tden
      INTEGER :: i,ib,l,iq,n,irc
      CHARACTER(4) flnm
      !
      WRITE(6,*) 'calculating Fourier transforms of vloc and tden'

      n=Grid%n
      h=Grid%h
      irc=PAW%irc
      r=>Grid%r
      ALLOCATE(q(nq),dum(n),dum1(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in allocating space in ftvloc',n,nq,i
         STOP
      ENDIF

      OPEN(55,file='ftvloc',form='formatted')
      dq=qmax/nq
      l=0
      DO iq=1,nq
         q(iq)=iq*dq
         DO i=1,n
            dum(i)=q(iq)*r(i)
         ENDDO
         CALL sphbes(l,n,dum)
         DO i=1,n
            dum1(i)=dum(i)*PAW%vloc(i)*(r(i)**2)
            dum(i)=dum(i)*PAW%tden(i)
         ENDDO
         vloc=integrator(Grid,dum1(1:n),1,irc)
         tden=integrator(Grid,dum(1:n))

         WRITE(55,'(1p,4e16.7)') q(iq),vloc,tden,vloc*tden
      ENDDO
      CLOSE(55)

      DEALLOCATE(q,dum,dum1)
    END SUBROUTINE ftvloc

    !*******************************************************************
    !  function to calculated <wfn|O|wfn> for smooth paw wavefunction
    !*******************************************************************
    FUNCTION sepnorm(Grid,PAW,nr,l,wfn)
      REAL(8) :: sepnorm
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: nr,l
      REAL(8), INTENT(IN) :: wfn(:)

      INTEGER :: n,ib,ic,nbase,irc
      REAL(8) :: h
      REAL(8), ALLOCATABLE :: b(:)

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc
      ALLOCATE(b(nbase),stat=ib)
      IF (ib/=0) THEN
         WRITE(6,*) 'Error in sepnorm allocation', nbase,ib
         STOP
      ENDIF

      IF (nr<irc) THEN
         WRITE(6,*) 'Error in sepnorm -- nr < irc'
         STOP
      ENDIF
      sepnorm=overlap(Grid,wfn,wfn,1,nr)
      b=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) b(ib)=overlap(Grid,wfn,PAW%otp(:,ib),1,irc)
      ENDDO
      DO ib=1,nbase
         DO ic=1,nbase
            sepnorm=sepnorm+b(ib)*PAW%oij(ib,ic)*b(ic)
         ENDDO
      ENDDO

      DEALLOCATE(b)
    END FUNCTION sepnorm


    !***************************************************************************
    !  pgm to solve separable radial schroedinger equation
    !    at energy 'energy' and at angular momentum l
    !
    !    with smooth potential rveff/r, given in uniform mesh of n points
    !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
    !                               p((n+1)*h)=0
    !
    !  uses Noumerov algorithm
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz==0
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
    !
    ! also returns node == number of nodes for calculated state
    !
    !  proj == projector functions
    !  hij and qij == hamiltonianian and overlap matrix elements
    !***************************************************************************
    SUBROUTINE unboundsep(Grid,Pot,PAW,nr,l,energy,wfn,node)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l,nr
      REAL(8), INTENT(IN) :: energy
      REAL(8), INTENT(INOUT) :: wfn(:)
      INTEGER, INTENT(INOUT) :: node

      INTEGER :: n,ia,ib,ic,nbase,icount,jcount,lcount,irc
      REAL(8) :: summ,h,scale,zeroval
      REAL(8), ALLOCATABLE :: y(:,:),b(:),a(:,:)

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc

      IF (nr<irc) THEN
         WRITE(6,*) 'error in unboundsep -- nr < irc', nr,irc
         STOP
      ENDIF
      !
      ! initialize wfn
      wfn=0
      wfn(2)=wfninit(0.d0,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
      zeroval=0
      if (l==1) zeroval=2

      call forward_numerov(Grid,l,nr,energy,Pot%rv,zeroval,wfn,node)

      ALLOCATE(y(nr,nbase),b(nbase),stat=ib)
      IF (ib/=0) THEN
         WRITE(6,*) 'Error in unboundsep allocation',  nr,nbase,ib
         STOP
      ENDIF

      lcount=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) THEN
            lcount=lcount+1
            CALL inhomogeneous_numerov(Grid,l,nr,energy,&
&              PAW%rveff,PAW%otp(:,ib),y(:,lcount))
         ENDIF
      ENDDO
      !
      IF(lcount>0) THEN
         ALLOCATE(a(lcount,lcount),stat=ib)
         IF (ib/=0) THEN
            WRITE(6,*) 'Error in unboundsep allocation',  lcount,ib
            STOP
         ENDIF

         icount=0
         DO ib=1,nbase
            summ=0.d0
            IF (l==PAW%l(ib)) THEN
               icount=icount+1
               DO ic=1,nbase
                  IF (l==PAW%l(ic)) THEN
                     summ=summ+(PAW%dij(ib,ic)-energy*PAW%oij(ib,ic))*&
&                         overlap(Grid,PAW%otp(:,ic),wfn,1,irc)
                  ENDIF
               ENDDO
               b(icount)=-summ
            ENDIF
         ENDDO
         !
         icount=0
         DO ia=1,nbase
            IF (l==PAW%l(ia)) THEN
               icount=icount+1
               jcount=0
               DO ib=1,nbase
                  IF (l==PAW%l(ib)) THEN
                     jcount=jcount+1
                     summ=0.d0
                     IF (ia.EQ.ib) summ=1.d0
                     DO ic=1,nbase
                        IF (l==PAW%l(ic)) THEN
                           summ=summ+(PAW%dij(ia,ic)-energy*PAW%oij(ia,ic))*&
&                          overlap(Grid,PAW%otp(:,ic),y(:,jcount),1,irc)
                        ENDIF
                     ENDDO
                     a(icount,jcount)=summ
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         !
         CALL linsol(a,b,lcount,lcount,lcount,nbase)

         icount=0
         DO ib=1,nbase
            IF(l==PAw%l(ib)) THEN
               icount=icount+1
               wfn(1:nr)=wfn(1:nr)+b(icount)*y(1:nr,icount)
            ENDIF
         ENDDO
         DEALLOCATE(a)
      ENDIF
      !
      ! normalize to unity within integration range
      !
      scale=1.d0/sepnorm(Grid,PAW,nr,l,wfn)
      IF (scale.LE.0.d0) THEN
         WRITE(6,*) 'warning -- negative norm for l=',l
         scale=-scale
         IF (scale.EQ.0.d0) scale=1.d0
      ENDIF
      scale=DSIGN(SQRT(scale),wfn(nr-2))
      wfn(1:nr)=wfn(1:nr)*scale
      DEALLOCATE(b,y)
    END SUBROUTINE unboundsep

    !***************************************************************************
    !  pgm to solve separable radial schroedinger equation
    !    for bound state near energy 'energy' and at angular momentum l
    !
    !    with smooth potential rveff/r, given in uniform mesh of n points
    !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
    !                               p((n+1)*h)=0
    !
    !  uses Noumerov algorithm
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz==0
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
    !
    ! also returns node == number of nodes for calculated state
    !
    !  proj == projector functions
    !  hij and qij == hamiltonianian and overlap matrix elements
    !***************************************************************************
    SUBROUTINE boundsep(Grid,Pot,PAW,l,node,energy,emin,emax,wfn)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l,node
      REAL(8), INTENT(INOUT) :: energy,emin,emax
      REAL(8), INTENT(INOUT) :: wfn(:)

      INTEGER, PARAMETER :: niter=100
      INTEGER :: n,i,j,ia,ib,ic,nbase,icount,jcount,lcount
      INTEGER :: match,mmatch,irc,node1,iter
      REAL(8), PARAMETER :: convre=1.d-10,err=1.d-9
      REAL(8) :: summ,h,scale,best,energy0,dele,x,rin,rout,zeroval,ppp
      REAL(8), ALLOCATABLE :: p1(:),u(:),y(:,:),b(:),a(:,:)

      n=Grid%n; h=Grid%h; nbase=PAW%nbase;  irc=PAW%irc
      !
      ALLOCATE(p1(n),u(n),y(n,nbase),b(nbase),stat=ib)
      IF (ib/=0) THEN
         WRITE(6,*) 'Error in boundsep allocation',  n,nbase,ib
         STOP
      ENDIF
      !WRITE(6,*) 'in boundsep with ', l,node,energy,emin,emax
      lcount=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) THEN
            lcount=lcount+1
         ENDIF
      ENDDO
      !
      IF(lcount>0) THEN
         ALLOCATE(a(lcount,lcount),stat=ib)
         IF (ib/=0) THEN
            WRITE(6,*) 'Error in boundsep allocation',  lcount,ib
            STOP
         ENDIF
      ENDIF

      IF (emax.GT.0.d0) emax=0.d0
      best=1.d10
      IF (energy.LT.emin) energy=emin+err
      IF (energy.GT.emax) energy=emax
      energy0=energy

      iter=0


      DO
         match=MIN(irc+10,n-10)
         x=0.5d0*(Pot%rv(n)/Grid%r(n)+Pot%rv(n-1)/Grid%r(n-1))+&
&               l*(l+1)/(Grid%r(n)**2)
         ppp=SQRT(ABS(x-energy))
         p1=0
         p1(n)=1
         p1(n-1)=exp(-ppp*(Grid%r(n-1)-Grid%r(n)))

         !write(6,*) 'before backward', n,p1(n-1),p1(n)
         !write(6,*) 'before backward', x,ppp,exp(-ppp*(Grid%r(n-1)-Grid%r(n)))
         !write(6,*) 'x,energy', x,energy,ABS(x-energy),SQRT(ABS(x-energy))
         !call flush(6)

         CALL backward_numerov(Grid,l,match-5,energy,Pot%rv,p1)
         rin=Gfirstderiv(Grid,match,p1)/p1(match)
         mmatch=match+1
         !WRITE(6,*) 'match, rin ' ,match,rin
         !
         !  perform outward integration until match point -- it is assumed
         !   that projector functions proj are zero for r>r(match)
         !
         ! initialize u
         u=0
         u(2)=wfninit(0.d0,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
         zeroval=0
         if (l==1) zeroval=2

         call forward_numerov(Grid,l,mmatch+5,energy&
                  ,Pot%rv,zeroval,u,node1)
         !
         lcount=0
         DO ib=1,nbase
            IF (l==PAW%l(ib)) THEN
               lcount=lcount+1
               CALL inhomogeneous_numerov(Grid,l,mmatch+5,energy,&
&                  Pot%rv,PAW%otp(:,ib),y(:,lcount))
            ENDIF
         ENDDO
         !
         IF(lcount>0) THEN

            icount=0
            DO ib=1,nbase
               summ=0.d0
               IF (l==PAW%l(ib)) THEN
                  icount=icount+1
                  DO ic=1,nbase
                     IF (l==PAW%l(ic)) THEN
                        summ=summ+(PAW%dij(ib,ic)-energy*PAW%oij(ib,ic))*&
&                            overlap(Grid,PAW%otp(:,ic),u,1,irc)
                     ENDIF
                  ENDDO
                  b(icount)=-summ
               ENDIF
            ENDDO
            !
            icount=0
            DO ia=1,nbase
               IF (l==PAW%l(ia)) THEN
                  icount=icount+1
                  jcount=0
                  DO ib=1,nbase
                     IF (l==PAW%l(ib)) THEN
                        jcount=jcount+1
                        summ=0.d0
                        IF (ia.EQ.ib) summ=1.d0
                        DO ic=1,nbase
                           IF (l==PAW%l(ic)) THEN
                              summ=summ+(PAW%dij(ia,ic)-energy*PAW%oij(ia,ic))*&
&                             overlap(Grid,PAW%otp(:,ic),y(:,jcount),1,irc)
                           ENDIF
                        ENDDO
                        a(icount,jcount)=summ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            !
            CALL linsol(a,b,lcount,lcount,lcount,nbase)

            wfn=0
            wfn(1:mmatch+5)=u(1:mmatch+5)
            icount=0
            DO ib=1,nbase
               IF(l==PAw%l(ib)) THEN
                  icount=icount+1
                  wfn(1:mmatch+5)=wfn(1:mmatch+5)+b(icount)*y(1:mmatch+5,icount)
               ENDIF
            ENDDO
         ENDIF

         rout=Gfirstderiv(Grid,match,wfn)/wfn(match)
         !WRITE(6,'("node,match,rin,rout",2i8,1p,2e15.7)') node1,match,rin,rout
         !   -- estimate correction
         node1=0
         wfn(:)=wfn(:)/wfn(match)
         DO j=3,match
               IF (wfn(j)*wfn(j-1).LT.0.d0) node1=node1+1
         ENDDO
         !WRITE(6,*) 'actual number of nodes', node1

!This test is obsolete: pseudo-WFs do not have to be orthogonal
!          IF (node1<node) THEN
!             ! too few nodes -- raise energy
!             emin=MAX(energy+err,emin)
!             energy=emax-(emax-energy)*ranx()
!             !WRITE(6,*) 'too few nodes -- energy raised', energy,emin,emax
!          ELSEIF (node1>node) THEN
!             ! too many nodes -- lower energy
!             emax=MIN(energy-err,emax)
!             energy=emin+(energy-emin)*ranx()
!             !WRITE(6,*) 'too many nodes -- energy lowered', energy,emin,emax
!             !do i=1,mmatch
!             !write(200+iter,'(1p,7e15.7)')Grid%r(i),wfn(i)
!             !enddo
!          ELSEIF (node1==node) THEN
            DO j=match,n
               wfn(j)=p1(j)/p1(match)
            ENDDO
            !  normalization
            scale=1.d0/sepnorm(Grid,PAW,n,l,wfn)
            dele=(rout-rin)*scale
            !WRITE(6,*) 'dele,scale',dele,scale
            scale=SQRT(scale)
            wfn=scale*wfn

            !do i=1,n
            ! write(100+iter,'(1p,2E15.7)') Grid%r(i),wfn(i)
            !enddo

            x=ABS(dele)
            IF (x.LT.best) THEN
               energy0=energy
               best=x
            ENDIF
            IF (ABS(dele).LE.convre) EXIT
            energy=energy+dele
            !WRITE(6,*) 'next energy' , energy,dele
            ! if energy is out of range, pick random energy in correct range
            IF (emin-energy.GT.convre.OR.energy-emax.GT.convre) THEN
               energy=emin+(emax-emin)*ranx()+err
            !   WRITE(6,*) 'energy out of range -- reranged --', energy
            ENDIF
!          ENDIF
         iter=iter+1
         !WRITE(6,*) 'Energy for next iteration ', iter,energy
         IF (iter.GT.niter) THEN
            WRITE(6,*) 'no convergence in boundsep',l,dele,energy
            WRITE(6,*) ' best guess of eig, dele = ',energy0,best
            STOP
         ENDIF
      ENDDO
      !
      ! normalize to unity within integration range
      !
      CALL filter(n,wfn,1.d-11)
      scale=1.d0/sepnorm(Grid,PAW,n,l,wfn)
      IF (scale.LE.0.d0) THEN
         WRITE(6,*) 'warning -- negative norm for l=',l
         scale=-scale
         IF (scale.EQ.0.d0) scale=1.d0
      ENDIF
      scale=DSIGN(SQRT(scale),wfn(n-2))
      wfn(1:n)=wfn(1:n)*scale
      !WRITE(6,*) 'exiting boundsep with energy ', l,energy
      DEALLOCATE(a,b,y,u,p1)
    END SUBROUTINE boundsep

    !*********************************************************
    !  program to to transform smooth wavefunction to all-electron
    !     wavefunction within Blochl's paw formalism
    !   otp == projector function
    !   odphi == ophi - otphi
    !*********************************************************
    SUBROUTINE PStoAE(Grid,PAW,nr,l,tpsi,psi)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: nr,l
      REAL(8), INTENT(IN) :: tpsi(:)
      REAL(8), INTENT(INOUT) :: psi(:)

      INTEGER :: n,ib,nbase,irc
      REAL(8) :: h,scale

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc

      IF (nr<irc) THEN
         WRITE(6,*) 'error in PStoAE  -- nr < irc', nr,irc
         STOP
      ENDIF

      psi(1:nr)=tpsi(1:nr)
      DO ib=1,nbase
         IF (l==PAW%l(ib)) psi(1:nr)=psi(1:nr)+&
&             (PAW%ophi(1:nr,ib)-PAW%otphi(1:nr,ib))*&
&             overlap(Grid,PAW%otp(:,ib),tpsi,1,irc)
      ENDDO
    END SUBROUTINE PStoAE

    !************************************************************************
    !  program to calculate logerivatives of paw wavefunctions
    !   and to compare them with all electron wavefunctions
    !  optionally, wavefunctions are written to a file
    !************************************************************************
    SUBROUTINE logderiv(Grid,Pot,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      TYPE(PotentialInfo) :: PS
      INTEGER :: n,l,ke,ie,nbase,ib,ic,nr,i,nodes,mbase,irc
      REAL(8), PARAMETER :: e0=-0.05d0 ! Where to print WFn in file
      REAL(8) :: de,h,x,dwdr,dcwdr,scale,energy
      REAL(8), ALLOCATABLE :: psi(:),tpsi(:),ttpsi(:)
      CHARACTER(4)  :: flnm

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc;  nr=irc+10

      ALLOCATE(psi(nr),tpsi(nr),ttpsi(nr),PS%rv(n),stat=ie)
      IF (ie/=0) THEN
         WRITE(6,*) 'Error in logderiv allocation',n,ie
         STOP
      ENDIF

     ! PAW%dij and PAW%oij already loaded --
      !!!!  load PAW%dij and PAW%oij

      !!!! PAW%dij=0; PAW%oij=0

      !!!! DO ib=1,nbase
         !!!! DO ic=ib,nbase
            !!!! IF (PAW%l(ib)==PAW%l(ic)) THEN
               !!!! CALL dqij(Grid,PAW,ib,ic,PAW%oij(ib,ic))
               !!!! !CALL dtij(Grid,PAW,ib,ic,PAW%dij(ib,ic))
               !!!! CALL altdtij(Grid,Pot,PAW,ib,ic,PAW%dij(ib,ic))
               !!!! CALL avij(Grid,Pot,PAW,ib,ic,x)
               !!!! !WRITE(6,*) 'dij',ib,ic,PAW%dij(ib,ic),x
               !!!! PAW%dij(ib,ic)=PAW%dij(ib,ic)+x
               !!!! IF  (ic>ib) THEN
                  !!!! PAW%oij(ic,ib)=PAW%oij(ib,ic)
                  !!!! PAW%dij(ic,ib)=PAW%dij(ib,ic)
               !!!! ENDIF
            !!!! ENDIF
         !!!! ENDDO
      !!!! ENDDO

      ! load  PS
      PS%rv=PAW%rveff ; PS%nz=0
      call zeropot(Grid,PS%rv,PS%v0,PS%v0p)
      !
      !   calculate logderivatives at irc
      !
      WRITE(6,*) 'calculating log derivatives at irc',Grid%r(irc)
      !
      de=(maxlogderiv-minlogderiv)/dble(nlogderiv-1)
      ke=1+anint((e0-minlogderiv)/de)

      mbase=nbase
      DO l=0,PAW%lmax+1
         CALL mkname(l,flnm)
         OPEN(56,file='logderiv.'//TRIM(flnm),form='formatted')

         DO ie=1,nlogderiv
            energy=minlogderiv+de*(ie-1)
            psi=0;tpsi=0;ttpsi=0
            if (scalarrelativistic) then
               CALL unboundsr(Grid,Pot,nr,l,energy,psi,nodes)
            else
               CALL unboundsch(Grid,Pot,nr,l,energy,psi,nodes)
            endif
            CALL unboundsep(Grid,PS,PAW,nr,l,energy,tpsi,nodes)
            CALL PStoAE(Grid,PAW,nr,l,tpsi,ttpsi)
            !
            dwdr=Gfirstderiv(Grid,irc,psi)/psi(irc)
            dcwdr=Gfirstderiv(Grid,irc,ttpsi)/ttpsi(irc)

            WRITE(56,'(1p,5e12.4)') energy,dwdr,dcwdr
            IF (ie.EQ.ke) THEN
               mbase=mbase+1
               CALL mkname(mbase,flnm)
               OPEN(57,file='wfn'//TRIM(flnm),form='formatted')
               WRITE(57,*) '# l=',l,'energy=',energy
               !
               ! form converted wavefunction and rescale exact wavefunction
               !
               scale=ttpsi(irc)/psi(irc)
               DO i=1,nr
                  WRITE(57, '(1p,5e12.4)') Grid%r(i),tpsi(i),ttpsi(i),psi(i)*scale
               ENDDO
               CLOSE(57)
            ENDIF
         ENDDO !ie
         CLOSE(56)

      ENDDO !l

      DEALLOCATE(psi,tpsi,PS%rv)
    END SUBROUTINE logderiv

    !**************************************************************************
    !**************************************************************************
    SUBROUTINE SCFPAW(Grid,nz,PAW,FC,AEOrbit,PSOrbit,Etotal,newconfig)
      TYPE (GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: nz
      TYPE (PseudoInfo), INTENT(INOUT) :: PAW
      TYPE (FCinfo), INTENT(INOUT) :: FC
      TYPE (OrbitInfo), INTENT(IN) :: AEOrbit
      TYPE (OrbitInfo), INTENT(INOUT) :: PSOrbit
      REAL(8), INTENT(OUT) :: Etotal
      LOGICAL, INTENT(IN), OPTIONAL :: newconfig

      TYPE (Anderson_Context) , POINTER :: PSd,PSw

      !  program to calculate the self-consistent density functional
      !    atom for a given configuration within the PAW formalism


      INTEGER, PARAMETER :: mxloop=1000
      INTEGER :: n,io,nbase,ip,l,nfix,i,j,k,mnbase,norbit
      INTEGER :: loop,irc,ib,ic,id,ie,node,ns,nb
      INTEGER , ALLOCATABLE :: mmap(:),mp(:)
      REAL(8), PARAMETER :: rimix=0.5,worst=1.d-5
      REAL(8), PARAMETER  :: conv1=1.d13,conv2=2.d13,conv3=3.d13,conv4=4.d13
      REAL(8) :: occ,xocc,rmix,cnvrg,ecoul,etxc,eexc,Q00,dEdQ,h
      REAL(8) :: energy, emin,emax,x,delta,w1,w2,w3,w4,dcore,tcore
      REAL(8), ALLOCATABLE :: v1(:),v2(:),d1(:),d2(:),wij(:,:),Psi(:,:),Eig(:)
      REAL(8), ALLOCATABLE :: den(:),dendiff(:),dij(:,:),wijold(:),wijdiff(:)
      REAL(8), ALLOCATABLE :: vtcore(:)
      TYPE(PotentialInfo) :: PSPot

      w1=conv1;w2=conv2;w3=conv3;w4=conv4
      nbase=PAW%nbase;  irc=PAW%irc;  n=Grid%n; mnbase=(nbase*(nbase+1))/2
      h=Grid%h
      !WRITE(6,*) ' in SCFPAW ', nbase,irc,n,mnbase
      !CALL flush(6)

      if (present(newconfig)) then
         norbit=0
         WRITE(6,*) 'Current occupancies:'
         WRITE(6,*) ' n l   occupancy        energy    '
         DO io=1,AEOrbit%norbit
            IF (.NOT.FC%iscore(io)) THEN
               WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') AEOrbit%np(io),&
&                   AEOrbit%l(io),AEOrbit%occ(io),&
&                   AEOrbit%eig(io)
               norbit=norbit+1
            ENDIF
         ENDDO
         PSOrbit%norbit=norbit
         ALLOCATE(PSOrbit%np(norbit),PSOrbit%l(norbit),PSOrbit%occ(norbit),&
&             PSOrbit%eig(norbit),PSOrbit%wfn(n,norbit),stat=ib)
         IF (ib /= 0 ) THEN
            WRITE(6,*) 'Allocation error in SCFPAW -- norbit,n',norbit,n,ib
            STOP
         ENDIF
         norbit=0
         DO io=1,AEOrbit%norbit
            IF (.NOT.FC%iscore(io)) THEN
               norbit=norbit+1
               PSOrbit%np(norbit)=AEOrbit%np(io)
               PSOrbit%l(norbit)=AEOrbit%l(io)
               PSOrbit%occ(norbit)=AEOrbit%occ(io)
               PSOrbit%eig(norbit)=AEOrbit%eig(io)
            ENDIF
         ENDDO

         ALLOCATE(d1(n),d2(n),v1(n),v2(n),PSPot%rv(n),PSPot%den(n),&
&             den(n),dendiff(n),vtcore(n),&
&             wij(nbase,nbase),dij(nbase,nbase),wijold(mnbase),wijdiff(mnbase),&
&             Psi(n,norbit),Eig(norbit),mmap(norbit),mp(norbit),stat=ib)
         IF (ib /= 0 ) THEN
            WRITE(6,*) 'Allocation error in SCFPAW -- nbase,n',nbase,n,ib
            STOP
         ENDIF


         WRITE(6,*) 'enter np  l   occ    for all revisions'
         WRITE(6,*) ' enter 0 0 0 to end'

         DO
            READ(5,*) ip,l,xocc
            IF (ip.LE.0) EXIT
            nfix=-100
            DO io=1,norbit
               IF (ip==PSOrbit%np(io).AND.l==PSOrbit%l(io)) THEN
                  nfix=io
                  EXIT
               ENDIF
            ENDDO
            IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
               WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
&                   ip,l,xocc,nfix,norbit
               STOP
            ENDIF
            PSOrbit%occ(nfix)=xocc
         ENDDO

         PSPot%nz=0
         PSPot%den=0
         FC%zvale=0
         WRITE(6,*) 'New configuration:'
         DO io=1,norbit
            WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') PSOrbit%np(io),&
&                PSOrbit%l(io),PSOrbit%occ(io), PSOrbit%eig(io)
            FC%zvale=FC%zvale+PSOrbit%occ(io)
            DO ib=1,PAW%nbase
               IF (PSOrbit%l(io)==PAW%l(ib).AND. &
&                   PSOrbit%np(io)==PAW%np(ib)) THEN
                  PSOrbit%wfn(:,io)=PAW%tphi(:,ib)
               ENDIF
            ENDDO
         ENDDO
    else
         norbit=0
         DO io=1,AEOrbit%norbit
            IF (.NOT.FC%iscore(io)) THEN
               norbit=norbit+1
            ENDIF
         ENDDO
         PSOrbit%norbit=norbit
         ALLOCATE(PSOrbit%np(norbit),PSOrbit%l(norbit),PSOrbit%occ(norbit),&
&             PSOrbit%eig(norbit),PSOrbit%wfn(n,norbit),stat=ib)
         IF (ib /= 0 ) THEN
            WRITE(6,*) 'Allocation error in SCFPAW -- norbit,n',norbit,n,ib
            STOP
         ENDIF
         norbit=0
         DO io=1,AEOrbit%norbit
            IF (.NOT.FC%iscore(io)) THEN
               norbit=norbit+1
               PSOrbit%np(norbit)=AEOrbit%np(io)
               PSOrbit%l(norbit)=AEOrbit%l(io)
               PSOrbit%occ(norbit)=AEOrbit%occ(io)
               PSOrbit%eig(norbit)=AEOrbit%eig(io)
            ENDIF
         ENDDO

         ALLOCATE(d1(n),d2(n),v1(n),v2(n),PSPot%rv(n),PSPot%den(n),&
&             den(n),dendiff(n),vtcore(n),&
&             wij(nbase,nbase),dij(nbase,nbase),wijold(mnbase),wijdiff(mnbase),&
&             Psi(n,norbit),Eig(norbit),mmap(norbit),mp(norbit),stat=ib)
         IF (ib /= 0 ) THEN
            WRITE(6,*) 'Allocation error in SCFPAW -- nbase,n',nbase,n,ib
            STOP
         ENDIF

         PSPot%nz=0
         PSPot%den=0
         FC%zvale=0
         WRITE(6,*) 'Recalculating configuration:'
         DO io=1,norbit
            WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') PSOrbit%np(io),&
&                PSOrbit%l(io),PSOrbit%occ(io), PSOrbit%eig(io)
            FC%zvale=FC%zvale+PSOrbit%occ(io)
            DO ib=1,PAW%nbase
               IF (PSOrbit%l(io)==PAW%l(ib).AND. &
&                   PSOrbit%np(io)==PAW%np(ib)) THEN
                  PSOrbit%wfn(:,io)=PAW%tphi(:,ib)
               ENDIF
            ENDDO
         ENDDO
    endif

      WRITE(6,*) ' Total charge, Valence charge = ', FC%zvale+FC%zcore, FC%zvale

      !  initialize Anderson Mix arrays
      rmix=rimix
      CALL InitAnderson(PSd,6,5,n,rmix,1.d5)
      CALL InitAnderson(PSw,6,5,mnbase,rmix,1.d5)
      cnvrg=worst
      WRITE(6,*) 'Density convergence parameter set to at least',cnvrg
      den=0; dendiff=0; wijold=0; wijdiff=0

      d1=FC%coreden-PAW%tcore
      dcore=integrator(Grid,d1,1,PAW%irc_core)
      WRITE(6,*) 'Delta core den = ', dcore
      tcore=integrator(Grid,PAW%tcore)
      WRITE(6,*) 'tcore total = ', tcore
      CALL poisson(Grid,tcore,PAW%tcore,vtcore,ecoul)
      DO loop = 1,mxloop
         delta=0
         IF (loop > 1) THEN
            den=PSPot%den
            i=0
            do ib=1,PAW%nbase
               do ic=ib,PAW%nbase
                  i=i+1
                  wijold(i)=wij(ib,ic)
               enddo
            enddo
         ENDIF

         PSPot%den=0;   wij=0
         Etotal=0
         WRITE(6,*) '  results for loop = ',loop
         WRITE(6,*) ' n  l     occupancy       energy'
         DO io=1,norbit
            PSPot%den=PSPot%den+PSOrbit%occ(io)*(PSOrbit%wfn(:,io))**2
            CALL calcwij(Grid,PAW,PSOrbit%l(io),PSOrbit%occ(io),&
&                PSOrbit%wfn(:,io),wij)
            CALL kinetic(Grid,PSOrbit%l(io),PSOrbit%wfn(:,io),x)
            Etotal=Etotal+PSOrbit%occ(io)*x
            WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&              PSOrbit%np(io),PSOrbit%l(io),&
&              PSOrbit%occ(io),PSOrbit%eig(io)
         ENDDO
            i=0
            do ib=1,PAW%nbase
               do ic=ib,PAW%nbase
                  i=i+1
                  wijdiff(i)=wij(ib,ic)-wijold(i)
               enddo
            enddo

         !  open(88,file='wijtest',form='formatted')
         !    Do ib=1,nbase
         !       Do ic=1,nbase
         !          write(88,'(2i5,1pe15.7)') ib,ic,wij(ib,ic)
         !       Enddo
         !    Enddo
         !  close(88)

         Q00=-nz+dcore
         DO ib=1,nbase
            DO ic=1,nbase
               Q00=Q00+wij(ib,ic)*PAW%Oij(ib,ic)
            ENDDO
         ENDDO

         WRITE(6,*) 'loop = ',loop,'  Q00 = ',Q00

         x=integrator(Grid,PSPot%den)
         !WRITE(6,*) '  x, density scaling ', x,-(Q00+tcore)/x
         !x=-(Q00+tcore)/x
         WRITE(6,*) '  x, density scaling ', &
&              x,(FC%zvale-nz+dcore-Q00)/x
         x=(FC%zvale-nz+dcore-Q00)/x
         PSPot%den=PSPot%den*x
         dendiff=PSPot%den-den

         IF (loop > 1 ) THEN
            delta=SUM(ABS(dendiff))
            WRITE(6,*) ' tden delta ',delta
            delta=delta+SUM(ABS(wijdiff))
            WRITE(6,*) ' second delta ' , delta
            CALL shift4(w1,w2,w3,w4,delta)
            CALL Anderson_Mix(PSd,den,dendiff)
            do i=1,n
               if (den(i) < 1.d-11) den(i)=0
            enddo
            PSPot%den=den
            CALL Anderson_Mix(PSw,wijold,wijdiff)
            CALL filter(mnbase,wijold,1.d-11)
            i=0
            DO ib=1,nbase
               DO ic=ib,nbase
                  i=i+1
                  wij(ib,ic)=wijold(i)
                  wij(ic,ib)=wijold(i)
               ENDDO
            ENDDO
            !Rescale densities
            Q00=-nz+dcore
            DO ib=1,nbase
               DO ic=1,nbase
                  Q00=Q00+wij(ib,ic)*PAW%Oij(ib,ic)
               ENDDO
            ENDDO
            WRITE(6,*) ' Q00 after Andersonmix = ',Q00
            x=integrator(Grid,PSPot%den)
            !WRITE(6,*) '  x, density scaling ', x,-(Q00+tcore)/x
            !x=-(Q00+tcore)/x
            WRITE(6,*) '  x, density scaling ', &
&               x,(FC%zvale-nz+dcore-Q00)/x
            x=(FC%zvale-nz+dcore-Q00)/x
            PSPot%den=PSPot%den*x
         ENDIF

         PSPot%q=integrator(Grid,PSPot%den)
         WRITE(6,*) 'tq = ', PSPot%q

         CALL poisson(Grid,PSPot%q,PSPot%den,PSPot%rv,ecoul)
         v1(2:n)=PSPot%rv(2:n)/Grid%r(2:n)
         call extrapolate(Grid,v1)
         dEdQ=overlap(Grid,v1,PAW%hatden,1,PAW%irc_shap)
         write(6,*) 'dEdQ old', dEdQ
         d1(2:n)=PAW%hatpot(2:n)/Grid%r(2:n)
         call extrapolate(Grid,d1)
         dEdQ=overlap(Grid,d1,PSPot%den)
         write(6,*) 'dEdQ new', dEdQ
         DO ib=1,nbase
            DO ic=1,nbase
               dEdQ=dEdQ-wij(ib,ic)*PAW%vhatij(ib,ic)
            ENDDO
         ENDDO

         WRITE(6,*) ' dEdQ = ', dEdQ
         d1=PAW%tcore+PSPot%den
         x=integrator(Grid,d1)
         write(6,*) 'tcore+tden total ', x
         CALL poisson(Grid,x,d1,v1,ecoul)
         CALL exch(Grid,d1,v2,etxc,eexc)
         PSPot%rv(:)=v1+v2+Grid%r(:)*PAW%vloc(:)+Q00*PAW%hatpot(:)
         !do i=1,n
         !write(122,'(1p,2e15.7)') Grid%r(i),PSPot%rv(i)
         !enddo
         !stop
         call zeropot(Grid,PSPot%rv,PSPot%v0,PSPot%v0p)

         ! if((loop/10)*10==loop) then
         !  do i=1,n
         !  write(100+(loop/10),'(1p12e15.7)') Grid%r(i),PSPot%rv(i),rv(i),rvdiff(i)
         !  enddo
         ! endif

         v1(2:n)=PAW%vloc(2:n)+(vtcore(2:n)+Q00*PAW%hatpot(2:n))/Grid%r(2:n)
         call extrapolate(Grid,v1)
         Etotal=Etotal+overlap(Grid,PSPot%den,v1)
         CALL poisson(Grid,PSPot%q,PSPot%den,v1,ecoul)
         Etotal=Etotal+ecoul+eexc
         PAW%dij=0

         DO ib=1,nbase
            DO ic=1,nbase
               PAW%dij(ib,ic)=PAW%tvij(ib,ic)-Q00*PAW%vhatij(ib,ic) &
&                   +dEdQ*PAW%v0ij(ib,ic)
               ! write(6,'("dij part",2i3,1p,10e15.7)') ib,ic,PAW%dij(ib,ic),PAW%tvij(ib,ic),PAW%vhatij(ib,ic),PAW%v0ij(ib,ic)
               ! call flush(6)
               Etotal=Etotal+wij(ib,ic)*(PAW%tvij(ib,ic)-Q00*PAW%vhatij(ib,ic))
               DO id=1,nbase
                  DO ie=1,nbase
                     PAW%dij(ib,ic)=PAW%dij(ib,ic)+&
&                         wij(id,ie)*PAW%vhijkl(ib,ic,id,ie)
                     !       write(6,'("dij part",4i3,1p,10e15.7)') ib,ic,id,ie,PAW%dij(ib,ic),PAW%vhijkl(ib,ic,id,ie)
                     !       call flush(6)
                     Etotal=Etotal+0.5d0*wij(ib,ic)*wij(id,ie)*PAW%vhijkl(ib,ic,id,ie)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         !  exchange-correlation part
         d1=FC%coreden; d2=PAW%tcore
         DO ib=1,nbase
            DO ic=1,nbase
               d1=d1+wij(ib,ic)*PAW%ophi(:,ib)*PAW%ophi(:,ic)
               d2=d2+wij(ib,ic)*PAW%otphi(:,ib)*PAW%otphi(:,ic)
            ENDDO
         ENDDO

         CALL exch(Grid,d1,v1,etxc,eexc,irc)
         Etotal=Etotal+eexc
         CALL exch(Grid,d2,v2,etxc,eexc,irc)
         Etotal=Etotal-eexc
         !CALL exch(n,h,FC%coreden,d2,x,eexc)
         !Etotal=Etotal-eexc
         !CALL exch(n,h,PAW%tcore,d2,x,eexc)
         !Etotal=Etotal+eexc


         DO ib=1,nbase
            DO ic=1,nbase
               IF (PAW%l(ib)==PAW%l(ic)) THEN
                  d1=v1(:)*PAW%ophi(:,ib)*PAW%ophi(:,ic) &
&                      -v2(:)*PAW%otphi(:,ib)*PAW%otphi(:,ic)
                  d1(2:n)=d1(2:n)/Grid%r(2:n)
                  call extrapolate(Grid,d1)
                  PAW%dij(ib,ic)=PAW%dij(ib,ic)+integrator(Grid,d1,1,irc)
                  !  write(6,*)' after vxc',ib,ic,PAW%dij(ib,ic)
               ENDIF
            ENDDO
         ENDDO
         WRITE(6,'(" loop, Etotal, delta",i8,1p,2e15.7)') loop,Etotal,delta


         DO l=0,PAW%lmax
            ns=0; Psi=0; Eig=0;mmap=0;mp=0;nb=0
            nb=0
            DO io=1,PAW%nbase
               IF (l == PAW%l(io)) THEN
                  nb=nb+1
                  mp(nb)=io
               ENDIF
            ENDDO
            DO io=1,norbit
               IF (l == PSOrbit%l(io)) THEN
                  ns=ns+1
                  mmap(ns)=io
                  Psi(:,ns)=PSOrbit%wfn(:,io)
                  CALL calcHPsi(Grid,PSPot,PAW,l,nb,mp,Psi(:,ns),d1)
                  Eig(ns)=overlap(Grid,Psi(:,ns),d1)
                  !write(6,*) 'Eig',ns,Eig(ns)
                  !Eig(ns)=PSOrbit%eig(io)
               ENDIF
            ENDDO
            !CALL BlockDavidson(Grid,PSPot,PAW,l,ns,Psi,Eig)
            emin=-500.d0
            DO io=1,ns
               If (io==1) then
                  do i=1,AEOrbit%norbit
                     if (AEOrbit%l(i)==l.AND.&
&                       FC%iscore(i)) &
&                      emin=AEOrbit%eig(i)+0.00001
                  enddo
               Endif
               node=io-1; emax=0
               !write(6,*) 'calling boundsep', node,emax,Eig(io)
               CALL boundsep(Grid,PSPot,PAW,l,node,&
&                   Eig(io),emin,emax,Psi(:,io))
               !PSOrbit%eig(mmap(io))=Eig(io)
               !PSOrbit%wfn(:,mmap(io))=Psi(:,io)
               emin=Eig(io)+0.00001
            ENDDO
            CALL gramschmidt(Grid,PAW,nb,mp,ns,Psi(:,1:ns))
            PSOrbit%eig(mmap(1:ns))=Eig(1:ns)
            PSOrbit%wfn(:,mmap(1:ns))=Psi(:,1:ns)
         ENDDO

         IF (loop>=4) THEN
            IF (.NOT.(w4.LE.w3.AND.w3.LE.w2 &
&                .AND.w2.LE.w1).AND.w4.LE.cnvrg) THEN
               !
               !  converged result
               !
               WRITE(6,*) ' SCFPAW converged in',loop,&
&               ' iterations -- delta = ',delta
               EXIT
            ENDIF
         ENDIF

      ENDDO


      CALL FreeAnderson(PSd)
      CALL FreeAnderson(PSw)
      DEALLOCATE(d1,d2,v1,v2,PSPot%rv,PSPot%den,den,dendiff,wij,&
&          wijold,wijdiff,Psi,Eig,mmap,mp,vtcore)
    END SUBROUTINE SCFPAW

    !*** NOTE: Block-Davidson routines are not currently working
    !***************************************************************************
    ! Routines needed for Block-Davidson solver
    !***************************************************************************

    SUBROUTINE gramschmidt(Grid,PAW,nb,mp,many,wfn)
      TYPE (GridInfo), INTENT(IN) :: Grid
      TYPE (PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: nb,mp(:)
      INTEGER, INTENT(INOUT) :: many
      REAL(8), INTENT(INOUT) :: wfn(:,:)

      REAL(8) :: smallest
      REAL(8), ALLOCATABLE :: b(:),bp(:),res(:)
      REAL(8) :: h,x
      INTEGER :: i,j,k,l,n,irc

      smallest=1.e-11
      n=Grid%n;  h=Grid%h;  irc=PAW%irc

      ALLOCATE(b(nb),bp(nb),res(many),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in graschmidt allocation', nb,i
         STOP
      ENDIF

      DO i=1,many
         IF (i>1) THEN
            DO j=1,i-1
               IF (res(j) >= smallest) THEN
                  DO k=1,nb
                     b(k)=overlap(Grid,wfn(:,i),PAW%otp(:,mp(k)),1,irc)
                  ENDDO
                  DO k=1,nb
                     bp(k)=overlap(Grid,wfn(:,j),PAW%otp(:,mp(k)),1,irc)
                  ENDDO
                  x=overlap(Grid,wfn(:,i),wfn(:,j))
                  DO k=1,nb
                     DO l=1,nb
                        x=x+b(k)*bp(l)*PAW%Oij(mp(k),mp(l))
                     ENDDO
                  ENDDO
                  wfn(:,i)=wfn(:,i)-x*wfn(:,j)
               ENDIF
            ENDDO
         ENDIF
         DO k=1,nb
            b(k)=overlap(Grid,wfn(:,i),PAW%otp(:,mp(k)),1,irc)
         ENDDO
         res(i)=overlap(Grid,wfn(:,i),wfn(:,i))
         DO k=1,nb
            DO l=1,nb
               res(i)=res(i)+b(k)*b(l)*PAW%Oij(mp(k),mp(l))
            ENDDO
         ENDDO
         !write(6,*) 'gramschmidt norm',i,res(i)

         IF (res(i) >= smallest) THEN
            x=SQRT(res(i))
            wfn(:,i)=wfn(:,i)/x
         ENDIF
   ENDDO

   i=0
   DO j=1,many
      IF (res(j) >= smallest) THEN
         i=i+1
         IF (j/=i) wfn(:,j)=wfn(:,i)
      ENDIF
   ENDDO
   !WRITE(6,*) 'Returning from gramschmidt with ',i,many
   many=i

   DEALLOCATE(b,bp,res)
 END SUBROUTINE gramschmidt


    !*** NOTE: Block-Davidson routines are not currently working
 SUBROUTINE calcHPsi(Grid,Pot,PAW,l,nb,mp,Psi,HPsi)
   TYPE(GridInfo), INTENT(IN) :: Grid
   TYPE(PotentialInfo), INTENT(IN) :: Pot
   TYPE(PseudoInfo), INTENT(INOUT) :: PAW
   INTEGER, INTENT(in) :: l,nb,mp(:)
   REAL(8), INTENT(IN) :: Psi(:)
   REAL(8), INTENT(INOUT) :: HPsi(:)

   INTEGER :: n,irc,ib,ic
   REAL(8) :: h
   REAL(8), ALLOCATABLE :: b(:)

   ALLOCATE(b(nb),stat=ib)
   IF (ib /= 0) THEN
      WRITE(6,*) 'Allocation Error in calcHpsi --',ib,nb
      STOP
   ENDIF

   n=GRID%n; h=GRID%h;  irc=PAW%irc
   HPsi=0
   CALL laplacian(Grid,l,Psi,HPsi)
   HPsi(2:n)=(Pot%rv(2:n)*Psi(2:n)-Hpsi(2:n))/Grid%r(2:n)
   DO ib=1,nb
      b(ib)=overlap(Grid,Psi,PAW%otp(:,mp(ib)),1,irc)
   ENDDO

   DO ib=1,nb
      DO ic=1,nb
         HPsi=HPsi+PAW%otp(:,mp(ib))*PAW%Dij(mp(ib),mp(ic))*b(ic)
      ENDDO
   ENDDO

   DEALLOCATE(b)
 END SUBROUTINE calcHPsi

    !*** NOTE: Block-Davidson routines are not currently working
 SUBROUTINE calcOPsi(Grid,PAW,nb,mp,Psi,OPsi)
   TYPE(GridInfo), INTENT(IN) :: Grid
   TYPE(PseudoInfo), INTENT(INOUT) :: PAW
   INTEGER, INTENT(in) :: nb,mp(:)
   REAL(8), INTENT(IN) :: Psi(:)
   REAL(8), INTENT(INOUT) :: OPsi(:)

   INTEGER :: n,irc,ib,ic
   REAL(8) :: h
   REAL(8), ALLOCATABLE :: b(:)

   ALLOCATE(b(nb),stat=ib)
   IF (ib /= 0) THEN
      WRITE(6,*) 'Allocation Error in calcOpsi --',ib,nb
      STOP
   ENDIF

   n=GRID%n; h=GRID%h;  irc=PAW%irc
   OPsi=Psi
   DO ib=1,nb
      b(ib)=overlap(Grid,Psi,PAW%otp(:,mp(ib)),1,irc)
   ENDDO

   DO ib=1,nb
      DO ic=1,nb
         OPsi=OPsi+PAW%otp(:,mp(ib))*PAW%Oij(mp(ib),mp(ic))*b(ic)
      ENDDO
   ENDDO

   DEALLOCATE(b)
 END SUBROUTINE calcOPsi


!******************************************************************************
!
!  Diagonalizer - Diagonalizes Hbase & OBase matrices
!
!  VecSize      - Number of basis vectors used to construct Hbase, OBase
!  ArraySize    - Dimension of Hbase, OBase
!  NewSize      - Number of EigenValues and EigenVectors returned
!  Eigen        - List of EigenValues
!  Vec          - MAtrix containing the new eigenvectors
!  DoOrthog     - If true, assume Obase is identity matrix
!  Hbase, Obase - Assumed to be Hermitian
!
!   Based on Alan Tackett's routine
!******************************************************************************

    !*** NOTE: Block-Davidson routines are not currently working
Subroutine Diagonalizer(VecSize, ArraySize, NewSize, Hbase, Obase, &
&                       Eigen,  Vec)
  Integer,          Intent(IN)  :: VecSize
  Integer,          Intent(IN)  :: ArraySize
  Integer,          Intent(OUT) :: NewSize
  REAL(8),          Intent(INOUT) :: Hbase(:,:)
  REAL(8),          Intent(INOUT) :: Obase(:,:)
  REAL(8),          Intent(OUT) :: Eigen(:)
  REAL(8),          Intent(OUT) :: Vec(:,:)

  Integer ::  i, j, k, LWork,  LSize
  Integer ::  Info
  REAL(8), allocatable  :: Omat(:,:),Hmat(:,:),VecR(:,:),Work(:),JUNK(:,:)
  Real(8) ,allocatable    ::  Lambda(:)
  Real(8)     :: tol,val

  Allocate(Hmat(ArraySize, ArraySize), VecR(ArraySize, ArraySize), &
&    Omat(ArraySize, ArraySize), Work(4*ArraySize), &
&    JUNK(ArraySize,ArraySize),  Lambda(ArraySize))
  tol=1.d-8

  LWork = 4*ArraySize
  NewSize=VecSize

  Hmat = Hbase;  Omat = Obase;  VecR = 0

  Hmat = Hbase;  Omat = Obase

  Info =13; NewSize=VecSize

        Call DSYEV('V', 'U', VecSize, Omat(1,1), ArraySize, Lambda, &
&         Work, LWork, Info)

        write(6,*) ' completed Omat diagonalization with Info=',Info

        j=0 ; VecR=0
        Do i=1,VecSize
           write(6,*) 'i lambda',i,Lambda(i)
           If (Lambda(i) > tol) then
              j=j+1
              VecR(1:VecSize,j)=Omat(1:VecSize,i)/SQRT(Lambda(i))
           EndIf
        Enddo


        If (j > 0) then
          NewSize=j
        Else
          Write(6,*) 'O matrix is singular', Lambda
          Stop
        Endif

        Write(6,*) 'NewSize = ',NewSize
        Omat=0
        do k=1,NewSize
           do i=1,VecSize
              do j=1,VecSize
                 Omat(i,k)=Omat(i,k)+(Hmat(j,i))*VecR(j,k)
              Enddo
           Enddo
         Enddo
        Hmat=0
        do k=1,NewSize
           do i=1,NewSize
              Do j=1,VecSize
                 Hmat(i,k)=Hmat(i,k)+(VecR(j,i))*Omat(j,k)
              Enddo
           Enddo
        Enddo

        Call DSYEV('V', 'U', NewSize, Hmat(1,1), ArraySize, Eigen, &
&            Work, LWork, Info)

        write(6,*) ' completed Hmat diagonalization with Info=',Info

        Omat=Hmat
        Hmat=0
        JUNK(1:NewSize,1:VecSize)=TRANSPOSE(VecR(1:VecSize,1:NewSize))
        do k=1,NewSize
           do i=1,VecSize
              Do j=1,NewSize
                 Hmat(i,k)=Hmat(i,k)+JUNK(j,i)*Omat(j,k)
              Enddo
           Enddo
        Enddo

        write(6,*) ' completed Hmat diagonalization with Info=',Info
        if (info /= 0) then
            write(6,*) 'Stopping due to diagonalizer error'
            stop
        endif


  Vec = Hmat

  DeAllocate(Hmat,Omat,VecR,Work,JUNK,Lambda)

End Subroutine


    !*** NOTE: Block-Davidson routines are not currently working
 SUBROUTINE BlockDavidson(Grid,Pot,PAW,l,ns,Psi,Eig)
   TYPE(GridInfo), INTENT(IN) :: Grid
   TYPE(PotentialInfo), INTENT(IN) :: Pot
   TYPE(PseudoInfo), INTENT(INOUT) :: PAW
   INTEGER, INTENT(in) :: l,ns
   REAL(8), INTENT(INOUT) :: Psi(:,:),Eig(:)

   INTEGER, ALLOCATABLE :: mp(:)
   INTEGER :: nb,na,last,start,finish,nbase
   INTEGER :: i,j,k,ib,ic,in,iter,n,irc,lwork
   INTEGER, PARAMETER :: niter=1000,repeat=1
   REAL(8), PARAMETER :: small=1.e-8
   REAL(8), PARAMETER  :: conv1=4.d13,conv2=3.d13,conv3=2.d13,conv4=1.d13
   REAL(8), ALLOCATABLE :: A(:,:),O(:,:),Vec(:,:),fn(:,:),en(:),w(:),work(:)
   REAL(8) :: delta,h,v1,v2,v3,v4

   v1=conv1;v2=conv2;v3=conv3;v4=conv4
   IF (ALLOCATED(mp)) DEALLOCATE(mp)
   nbase=PAW%nbase

   nb=0
   DO i=1,nbase
      IF (l==PAW%l(i)) THEN
         nb=nb+1
      ENDIF
   ENDDO

   IF (nb > 0) THEN
      ALLOCATE(mp(nb),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in BlockDavidson-makemp',nb,i
         STOP
      ENDIF
      j=0
      DO i=1,nbase
         IF (l==PAW%l(i)) THEN
            j=j+1
            mp(j)=i
         ENDIF
      ENDDO
   ELSE
      WRITE(6,*) 'Error in BlockDavidson-makemap',l,PAW%l
      STOP
   ENDIF

   n=Grid%n; h=Grid%h;  irc=PAW%irc
   na=ns*(repeat+1); lwork=na*na+3*na
   ALLOCATE(A(na,na),O(na,na),Vec(na,na),fn(n,na),en(na),w(n),work(lwork),stat=i)
   IF (i /=0) THEN
      WRITE(6,*) 'Allocation error in BlockDavidson',na,i
      STOP
   ENDIF
   A=0;O=0;fn=0;en=0

   in=ns
   CALL gramschmidt(Grid,PAW,nb,mp,in,Psi)
   IF (in<ns) THEN
      WRITE(6,*) 'Error in initial states ', in,ns
      STOP
   ENDIF

   delta=1.e10
   DO iter=1,niter
      IF ((ABS(delta)<small).OR.  &
&            (.NOT.(v4.LE.v3.AND.v3.LE.v2 &
&                .AND.v2.LE.v1).AND.delta.LE.100*small)) THEN
         WRITE(6,*) 'returning from BlockDavidson with',&
&             ' iter, delta = ', iter,delta
         DEALLOCATE(mp,A,O,Vec,fn,en,w,work)
         RETURN
      ELSE
         start=1; finish=ns; last=finish
         fn=0; en=0
         DO i=1,ns
            fn(:,i)=Psi(:,i)
            en(i)=Eig(i)
         ENDDO

         w=0
         DO k=1,repeat
            last=finish
            DO i=start,finish
               CALL calcHPsi(Grid,Pot,PAW,l,nb,mp,fn(:,i),w)
               en(i)=overlap(Grid,fn(:,i),w)
               last=last+1
               fn(:,last)=w
               CALL calcOPsi(Grid,PAW,nb,mp,fn(:,i),w)
               en(i)=en(i)/overlap(Grid,fn(:,i),w)
               fn(:,last)=fn(:,last)-en(i)*w
            ENDDO
            start=finish+1; finish=last
         ENDDO
         in=last
         !WRITE(6,*) 'calling gramschmidt with in = ',in
         !CALL gramschmidt(Grid,PAW,nb,mp,in,fn)
         !WRITE(6,*) 'returned gramschmidt with in = ',in
         A=0 ; O=0
         DO i=1,in
            CALL calcHPsi(Grid,Pot,PAW,l,nb,mp,fn(:,i),w)
            DO j=1,i
               A(j,i)=overlap(Grid,fn(:,j),w)
               A(i,j)=A(j,i)
               write(6,*) 'A(ij)', i,j,A(i,j)
            ENDDO
            CALL calcOPsi(Grid,PAW,nb,mp,fn(:,i),w)
            DO j=1,i
               O(j,i)=overlap(Grid,fn(:,j),w)
               O(i,j)=O(j,i)
               write(6,*) 'O(ij)', i,j,O(i,j)
            ENDDO
         ENDDO

         CALL Diagonalizer(in,na,i,A,O,w,Vec)

         write(6,*) 'Finished Diagonalizer with i = ', i
         if (i < ns) then
            write(6,*) 'Too few eigenvalues',i,ns
            stop
         endif
         in=i

         delta=0
         DO i=1,ns
            delta=delta+ABS(w(i)-Eig(i))
            WRITE(6,'("EIGEN UPDATE", i5, 1p,4e15.7)')i,w(i),Eig(i),w(i)-Eig(i)
            Eig(i)=w(i)
         ENDDO
         WRITE(6,*) 'loop, delta ', iter,delta
         CALL shift4(v1,v2,v3,v4,delta)
         Psi=0
         DO i=1,ns
            DO j=1,in
               Psi(:,i)=Psi(:,i)+fn(:,j)*Vec(j,i)
            ENDDO
         ENDDO
      ENDIF
   ENDDO

   WRITE(6,*) ' BlockDavidson did not converge '
   STOP

 END SUBROUTINE BlockDavidson

!******************************************************************************
! Pseudization routine: PSPOLYN
!  Pseudize a function with a polynom
!                   tfunc(r)=r^(l+1).Sum[Ci.r^2i]  0<=i<=np-1  if r<=rc
!                   tfunc(r)=func(r)                           if r>rc
!  Ci coefficients are returned
!******************************************************************************

 SUBROUTINE pspolyn(func,Ci,r,l,np,irc,n)

   INTEGER,INTENT(IN) :: irc,l,n,np
   REAL(8),INTENT(IN) :: func(n),r(n)
   REAL(8),INTENT(OUT) :: Ci(np)

   INTEGER :: i,j,np2
   REAL(8) :: rc,xx,y,scale
   REAL(8),ALLOCATABLE :: A(:,:),X(:)

   if (irc<3.or.irc>irc.or.irc>n-3) stop 'pspolyn: rc out of range'
   if (np<1) stop 'pspolyn: p out of range'

   allocate(A(np,np),X(np),stat=i)
   if (i/=0) stop 'allocation error in pspolyn'

   rc=r(irc);np2=np/2
   scale=(rc/(rc-r(irc-1)))**2 ! Scale to limit rounding error in linsol
   do i=1,np
    xx=r(i+irc-np2-1)/rc
    y=xx*xx
    A(i,1)=scale
    do j=2,np
     A(i,j)=A(i,j-1)*y
    enddo
    X(i)=scale*func(i+irc-np2-1)/(xx**(l+1))
   enddo

   call linsol(A,X,np,np,np,np)
   write(6,*) 'Completed linsol with coefficients'
   write(6,'(1p,10e15.7)') (X(i),i=1,np)

   do i=1,np
    Ci(i)=X(i)/rc**(l+2*i-1)
   enddo

   deallocate(A,X)

 END SUBROUTINE pspolyn

!******************************************************************************
! Pseudization routine: PSUSPOLYN
!  Pseudize a function with a polynom
!                   tfunc(r)=r^(l+1).Sum[Ci.r^2i]  0<=i<=np-1  if r<=rc
!                   tfunc(r)=func(r)                           if r>rc
!      For i>3, Ci coefficients are computed so that to minimize
!      Fourier coefficients of pseudized function for q>qcut
!  Ci coefficients are returned
!******************************************************************************

 SUBROUTINE psuspolyn(func,Ci,r,l,np,irc,n,qcut)

   INTEGER,INTENT(IN) :: irc,l,n,np
   REAL(8),INTENT(IN) :: qcut
   REAL(8),INTENT(IN) :: func(n),r(n)
   REAL(8),INTENT(OUT) :: Ci(np)

   INTEGER,PARAMETER :: nq=2001
   REAL(8),PARAMETER :: lfact(0:5)=(/1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/)
   INTEGER :: i,j,k,ip,iq,jq,il,ix,jx
   REAL(8) :: qh,qrc2,rc,xx,scale,yy(6,6),zz(6)
   REAL(8),ALLOCATABLE :: A(:,:),X(:),Q(:,:),qq(:),ff(:),gg(:)

   if (irc<3.or.irc>irc.or.irc>n-6) stop 'psuspolyn: rc out of range'
   if (np<4) stop 'psuspolyn: polynomial degree for pseudization was too small'

   allocate(A(np+4,np+4),X(np+4),Q(nq,np+4),qq(nq),ff(nq),gg(nq),stat=i)
   if (i/=0) stop 'allocation error in psuspolyn'

   il=l+1;rc=r(irc)

   qh=qcut/dfloat(nq-1)
   do iq=1,nq
    qq(iq)=0.5d0*qcut+qh*dble(iq-1)
   enddo

   Q(:,:)=0.d0
   do ip=1,np
    do iq=1,nq
     qrc2=0.5d0*(qq(iq)*rc)**2
     xx=1.d0;ix=0
     do while (abs(xx)>1.d-20)
      ix=ix+1
      xx=xx*qrc2/dble(ix)/dble(2*(ix+l)+1)
     enddo
     xx=0.d0
     do jx=ix,1,-1
      xx=qrc2/dble(jx)/dble(2*(jx+l)+1)*(1.d0/dble(2*(jx+ip+l)+1)-xx)
     enddo
     Q(iq,ip)=4.d0*pi*rc**dble(2*ip+il)*(qq(iq)*rc)**dble(l)/lfact(l) &
&                 *(1.0/dble(2*(l+ip)+1)-xx)
    enddo
   enddo

   A(:,:)=0.d0
   do iq=1,np
    do ip=1,np
     do jq=1,nq
      ff(jq)=qq(jq)**4*Q(jq,iq)*Q(jq,ip)
     enddo
     A(iq,ip)=overint(nq,qh,ff)
    enddo
   enddo
   do iq=1,np
    ix=2*(iq-1)+il
    A(np+1,iq)=rc**dble(ix)
    A(np+2,iq)=dble(ix)*rc**(ix-1)
    A(np+3,iq)=dble(ix*(ix-1))*rc**(ix-2)
    A(np+4,iq)=dble(ix*(ix-1)*(ix-2))*rc**(ix-3)
    A(iq,np+1)=A(np+1,iq)
    A(iq,np+2)=A(np+2,iq)
    A(iq,np+3)=A(np+3,iq)
    A(iq,np+4)=A(np+4,iq)
   enddo

   yy(:,:)=0.d0;zz(:)=0.d0
   do i=1,6
    do j=1,6
     if (i==1.and.j==1) then
      yy(i,j)=12.d0
     else
      do k=1,12
       yy(i,j)=yy(i,j)+(r(irc+k-6)-rc)**(i+j-2)
      enddo
     endif
    enddo
   enddo
   do k=1,12
    zz(1)=zz(1)+func(irc+k-6)
   end do
   do i=2,6
    do k=1,12
     zz(i)=zz(i)+func(irc+k-6)*(r(irc+k-6)-rc)**(i-1)
    end do
   end do
   scale=1/(rc-r(irc-1))**3;yy=yy*scale;zz=zz*scale ! Scale to limit rounding error in linsol
   call linsol(yy,zz,6,6,6,6)
   zz(3)=2.d0*zz(3);zz(4)=6.d0*zz(4)
   X(np+1:np+4)=zz(1:4)

   do iq=1,nq
    gg(iq)=4.d0*pi*intjl(rc,qq(iq),zz,l)
   enddo
   do ip=1,np
    do iq=1,nq
     ff(iq)=qq(iq)**4*gg(iq)*Q(iq,ip)
    enddo
    X(ip)=-overint(nq,qh,ff)
   enddo

   scale=(rc/(rc-r(irc-1)))**2;A=A*scale;X=X*scale ! Scale to limit rounding error in linsol
   call linsol(A,X,np+4,np+4,np+4,np+4)
   write(6,*) 'Completed linsol with coefficients'
   write(6,'(1p,10e15.7)') (X(i),i=1,np+4)

   Ci(1:np)=X(1:np)

   deallocate(A,X,Q,qq,ff,gg)

 END SUBROUTINE psuspolyn

!******************************************************************************
! Pseudization routine: PSBES
!  Pseudize a function with a sum of 2 Bessel functions
!                (following PHYS REV B 41,1227 (1990))
!                   tfunc(r)=[al(1)*jl(ql(1)*r)+al(2)*jl(ql(2)*r)]*r if r<=rc
!                   tfunc(r)=func(r)                                 if r>rc
!  al and ql coefficients are returned
!******************************************************************************

 SUBROUTINE psbes(func,al,ql,Grid,l,irc,n)

   INTEGER,INTENT(IN) :: irc,l,n
   REAL(8),INTENT(IN) :: func(n)
   REAL(8),INTENT(OUT) :: al(2),ql(2)
   TYPE(GridInfo),INTENT(IN) :: Grid

   INTEGER :: i
   REAL(8) :: alpha,beta,det,qr,jbes,jbesp,jbespp,rc
   REAL(8) :: amat(2,2),bb(2)

   rc=Grid%r(irc)

   beta=1.D0
   alpha=1.D0-Gfirstderiv(Grid,irc,func)*rc/func(irc)
   call solvbes(ql,alpha,beta,l,2)
   ql(1:2)=ql(1:2)/rc

   do i=1,2
    qr=ql(i)*rc
    call jbessel(jbes,jbesp,jbespp,l,2,qr)
    jbespp=2.d0*ql(i)*jbesp+jbespp*ql(i)*ql(i)*rc
    jbesp=jbes+jbesp*ql(i)*rc
    jbes=jbes*rc
    amat(1,i)=jbes
    amat(2,i)=jbespp
   enddo

   bb(1)=func(irc)
   bb(2)=Gsecondderiv(Grid,irc,func)

   det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
   al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
   al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det

 END SUBROUTINE psbes

END MODULE basis
