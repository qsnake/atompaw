MODULE pseudo
  USE GlobalMath
  USE atomdata
  USE aeatom
  USE calcpotential
  USE radialsch
  USE radialsr
  USE anderson_realmix

  IMPLICIT NONE

  TYPE  Pseudoinfo
     INTEGER  :: lmax,irc,irc_shap,irc_vloc,irc_core
     CHARACTER(132) :: Vloc_description
     CHARACTER(132) :: Proj_description
     CHARACTER(132) :: Comp_description
     REAL(8) :: rc,rc_shap,rc_vloc,rc_core,energyoflmax
     REAL(8), POINTER :: vloc(:),rveff(:),abinitvloc(:),abinitnohat(:)
     REAL(8), POINTER :: projshape(:),hatshape(:),hatden(:),hatpot(:)
     REAL(8), POINTER :: den(:),tden(:),tcore(:)
     INTEGER :: nbase
     INTEGER, POINTER :: np(:),l(:)
     CHARACTER(8), POINTER :: label(:)
     REAL(8), POINTER :: phi(:,:),tphi(:,:),tp(:,:) ! before orthog
     REAL(8), POINTER :: ophi(:,:),otphi(:,:),otp(:,:) ! after orthog
     REAL(8), POINTER :: Kop(:,:)    ! for storing K|phi>
     REAL(8), POINTER :: eig(:),occ(:),ck(:),vrc(:)
     REAL(8), POINTER :: oij(:,:),dij(:,:),multipole(:,:,:)
     !** L=0 matrix elements for atomic SC calculations
     REAL(8), POINTER :: tvij(:,:),vhatij(:,:),kin(:,:)
     REAL(8), POINTER :: v0ij(:,:),vhijkl(:,:,:,:)
  END  TYPE Pseudoinfo


CONTAINS

  !***************************************************************
  ! SUBROUTINE troullier(lmax,Grid,Pot)
  !  Creates  screened norm-conserving pseudopotential following
  !    approach of N. Troullier and J. L. Martins, PRB 43, 1993 (1991)
  !    Uses p(r)=a0+f(r); f(r)=SUMm(Coef(m)*r^(2*m), where
  !          m=1,2..6
  !    Psi(r) = r^(l+1)*exp(p(r))
  !***************************************************************
  SUBROUTINE Troullier(Grid,Pot,PAW,l,e)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
    INTEGER,INTENT(IN) :: l
    REAL(8),INTENT(IN) :: e

    REAL(8), ALLOCATABLE :: VNC(:)
    REAL(8) :: A0,A,B,B0,C,C0,D,F,S
    REAL(8) :: Coef(6),Coef0,Coef0old
    REAL(8) :: h,rc,delta,x,pp,dpp,ddpp,dddpp,ddddpp
    REAL(8) :: gam,bet
    INTEGER :: i,j,k,n,iter,nr,nodes,irc,ok,m,wavetype
    INTEGER, PARAMETER :: niter=5000
    REAL(8), PARAMETER :: small=1.0d-9
    REAL(8), ALLOCATABLE ::  wfn(:),p(:),dum(:)
    REAL(8), POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    ALLOCATE(VNC(n),wfn(nr),p(nr),dum(nr),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in troullier  -- in allocating wfn,p', nr,ok
       STOP
    ENDIF

    !write(6,*) ' Troullier ', n,nr,irc
    !call flush(6)
    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot,nr,l,e,wfn,nodes)
    endif

    IF (wfn(irc)<0) wfn=-wfn
    dum(1:irc)=(wfn(1:irc)**2)
    S=integrator(Grid,dum(1:irc),1,irc)
    A0=LOG(wfn(irc)/(rc**(l+1)))
    B0=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
    C0=rc*(rv(irc)-rc*e)-B0*(B0+2*l+2)
    D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))-2*B0*C0-2*(l+1)*(C0-B0)
    F=rc*(2*rv(irc)-rc*(2*Gfirstderiv(Grid,irc,rv) &
&        -rc*Gsecondderiv(Grid,irc,rv)))+&
&        4*(l+1)*(C0-B0)-2*(l+1)*D-2*C0**2-2*B0*D

    WRITE(6,*) 'In troullier -- matching parameters',S,A0,B0,C0,D,F

    delta=1.d10
    iter=0
    Coef0=0

    DO WHILE(delta>small.AND.iter<=niter)
       iter=iter+1
       A=A0-Coef0
       B=B0
       C=C0
       CALL EvaluateTp(l,A,B,C,D,F,coef)

       dum=0
       DO  i=1,irc
          x=(r(i)/rc)**2
          p(i)=x*(Coef(1)+x*(Coef(2)+x*(Coef(3)+&
&              x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
          dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
       ENDDO
       Coef0old=Coef0

       x=integrator(Grid,dum(1:irc),1,irc)
       Coef0=(LOG(S/x))/2

       delta=ABS(Coef0-Coef0old)
       WRITE(6,'(" VNC: iter Coef0 delta",i5,1p,2e15.7)') iter,Coef0,delta
    ENDDO

    WRITE(6,*) '  VNC converged in ', iter,'  iterations'
    WRITE(6,*) '  Coefficients  -- ', Coef0,Coef(1:6)
    !
    ! Now  calculate VNC
    OPEN(88,file='NC',form='formatted')
    !
    VNC=0
    DO  i=2,nr
       x=(r(i)/rc)**2
       p(i)=Coef0+x*(Coef(1)+x*(Coef(2)+&
&           x*(Coef(3)+x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
       dpp=2*r(i)/(rc**2)*(Coef(1)+x*(2*Coef(2)+x*(3*Coef(3)+&
            x*(4*Coef(4)+x*(5*Coef(5)+x*6*Coef(6))))))
       ddpp=(1/(rc**2))*(2*Coef(1)+x*(12*Coef(2)+x*(30*Coef(3)+&
&           x*(56*Coef(4)+x*(90*Coef(5)+x*132*Coef(6))))))
       dddpp=(r(i)/rc**4)*(24*Coef(2)+x*(120*Coef(3)+x*(336*Coef(4)+&
&           x*(720*Coef(5)+x*1320*Coef(6)))))
       ddddpp=(1/(rc**4)*(24*Coef(2)+x*(360*Coef(3)+x*(1680*Coef(4)+&
&           x*(5040*Coef(5)+x*11880*Coef(6))))))
       IF (i==irc) THEN
          WRITE(6,*) 'check  dp ', dpp,  B0/rc
          WRITE(6,*) 'check ddp ', ddpp, C0/rc**2
          WRITE(6,*) 'check dddp', dddpp, D/rc**3
          WRITE(6,*) 'check ddddp', ddddpp, F/rc**4
       ENDIF
       VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
       dum(i)=(r(i)**(l+1))*EXP(p(i))
       WRITE(88,'(1p,5e15.7)') r(i),wfn(i),dum(i),VNC(i)*r(i),rv(i)
    ENDDO
    CLOSE(88)
    x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
    WRITE(6,*) 'check norm ',x,S

    VNC(irc:n)=rv(irc:n)/r(irc:n)
    PAW%rveff(1:n)=VNC(1:n)*r(1:n)

    DEALLOCATE(VNC,wfn,p,dum)
  END SUBROUTINE troullier

  !***************************************************************
  ! SUBROUTINE kerker(lmax,Grid,Pot)
  !  Creates  screened norm-conserving pseudopotential following
  !    approach of G. P. Kerker, J. Phys. C. 13,L189-L194 (1980)
  !    Uses p(r)=a0+f(r); f(r)=SUMi(Coef(i)*r^m(i)), where m(i)
  !          are input powers
  !    Psi(r) = r^(l+1)*exp(p(r)) if PStype = EXPF
  !    Psi(r) = r^(l+1)*(p(r))    if PStype = POLY
  !***************************************************************
  SUBROUTINE kerker(Grid,Pot,PAW)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW

    REAL(8), ALLOCATABLE :: VNC(:)
    REAL(8) :: A0,A,B,C,D,S,Coef(4),Coef0,Coef0old
    REAL(8) :: h,e,rc,delta,x,pp,dpp,ddpp,dddpp
    REAL(8) :: gam,bet
    INTEGER :: i,j,k,n,iter,nr,nodes,irc,l,ok,m(4),wavetype
    INTEGER, PARAMETER :: EXPF=1, POLY=2
    INTEGER, PARAMETER :: niter=5000
    REAL(8), PARAMETER :: small=1.0d-12
    CHARACTER(10) :: vtype
    REAL(8), ALLOCATABLE ::  wfn(:),p(:),dum(:)
    REAL(8), POINTER :: r(:),rv(:)

    DO
       WRITE(6,*) 'Input "EXPF" or "POLY"  pseudowave form'
       READ(5,*) vtype
       IF (TRIM(vtype)=="EXPF".OR.TRIM(vtype)=="expf") THEN
          wavetype=EXPF
          EXIT
       ELSE IF (TRIM(vtype)=="POLY".OR.TRIM(vtype)=="poly") THEN
          wavetype=POLY
          EXIT
       ENDIF
    ENDDO

    DO
       WRITE(6,*) 'Input angular momentum l and energy e to set VNC '
       READ(5,*) l,e
       IF (l >= 0 .AND. l < 10) EXIT
    ENDDO

    m=0
    DO
       WRITE(6,*) 'Input the 4 powers for the polynomial f(r)'
       READ(5,*) m(1),m(2),m(3),m(4)
       IF (m(1)>0.AND.m(2)>0.AND.m(3)>0.AND.m(4)>0) EXIT
    ENDDO

    IF (wavetype==EXPF) THEN
       WRITE(PAW%Vloc_description,&
&           '("Norm-conserving Exp Vloc; l = ",i1,"; powers = ",4i3,"; e = ",1pe12.3)')&
&           l, m(1),m(2),m(3),m(4),e
       WRITE(6,*) PAW%Vloc_description
    ENDIF
    IF (wavetype==POLY) THEN
       WRITE(PAW%Vloc_description,&
&           '("Norm-conserving Poly Vloc; l = ",i1,"; powers = ",4i3,"; e = ",1pe12.3)')&
&           l, m(1),m(2),m(3),m(4),e
       WRITE(6,*) PAW%Vloc_description
    ENDIF

    n=Grid%n
    n=Grid%n
    h=Grid%h
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc
    ALLOCATE(VNC(n),wfn(nr),p(nr),dum(nr),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in kerker  -- in allocating wfn,p', nr,ok
       STOP
    ENDIF

    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot,nr,l,e,wfn,nodes)
    endif



    IF (wfn(irc)<0) wfn=-wfn
    dum(1:irc)=(wfn(1:irc)**2)
    S=integrator(Grid,dum(1:irc),1,irc)
    IF (wavetype==EXPF) THEN
       A0=LOG(wfn(irc)/(rc**(l+1)))
       B=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
       C=rc*(rv(irc)-rc*e)-B*(B+2*l+2)
       D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))-2*B*C-2*(l+1)*(C-B)
    ENDIF

    IF (wavetype==POLY) THEN
       A0=(wfn(irc)/(rc**(l+1)))
       B=(rc*Gfirstderiv(Grid,irc,wfn))/(rc**(l+1))-(l+1)*A0
       C=rc*(rv(irc)-rc*e)*A0-2*(l+1)*B
       D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))*A0+2*(l+1)*(B-C)+&
&           rc*(rv(irc)-rc*e)*B
    ENDIF


    WRITE(6,*) 'In kerker -- matching parameters',S,A0,B,C,D

    delta=1.d10
    iter=0
    Coef0=0

    DO WHILE(delta>small.AND.iter<=niter)
       iter=iter+1
       A=A0-Coef0
       CALL EvaluateP(m,A,B,C,D,Coef)

       dum=0
       DO  i=1,irc
          x=(r(i)/rc)
          p(i)=(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+(x**m(3))*Coef(3)+(x**m(4))*Coef(4)
          IF (wavetype==EXPF)dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
          IF (wavetype==POLY)dum(i)=(wfn(i))**2-((r(i)**(l+1))*(p(i)))**2
       ENDDO
       Coef0old=Coef0
       IF (wavetype==EXPF) THEN
          x=integrator(Grid,dum(1:irc),1,irc)
          Coef0=(LOG(S/x))/2
       ENDIF
       IF (wavetype==POLY) THEN
          gam=(2*l+3)*integrator(Grid,dum(1:irc),1,irc)/(rc**(2*l+3))
          bet=(2*l+3)*(Coef(1)/(2*l+3+m(1))+Coef(2)/(2*l+3+m(2))+&
&              Coef(3)/(2*l+3+m(3))+Coef(4)/(2*l+3+m(4)))
          WRITE(6,'("VNC: iter -- bet,gam = ",i5,1p,4e15.7)') iter,bet,gam
          x=bet**2+gam
          Coef0old=Coef0
          IF (x<0.d0) THEN
             WRITE(6,*) 'Warning in Kerker subroutine x = ',x
               Coef0=Coef0+0.1*A0
            ELSE
               Coef0=SQRT(x)-bet
            ENDIF
         ENDIF
         delta=ABS(Coef0-Coef0old)
         WRITE(6,*) '  VNC: iter  Coef0  delta', iter,Coef0,delta
      ENDDO

      WRITE(6,*) '  VNC converged in ', iter,'  iterations'
      WRITE(6,*) '  Coefficients  -- ', Coef0,Coef(1:4)
      !
      ! Now  calculate VNC
      OPEN(88,file='NC',form='formatted')
      !
      VNC=0
      DO  i=1,nr
         x=(r(i)/rc)
         p(i)=Coef0+(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+&
&             (x**m(3))*Coef(3)+(x**m(4))*Coef(4)
         dpp=(m(1)*(x**(m(1)-1))*Coef(1)+m(2)*(x**(m(2)-1))*Coef(2)+&
&             m(3)*(x**(m(3)-1))*Coef(3)+m(4)*(x**(m(4)-1))*Coef(4))/rc
         ddpp=(m(1)*(m(1)-1)*(x**(m(1)-2))*Coef(1)+&
&             m(2)*(m(2)-1)*(x**(m(2)-2))*Coef(2)+&
&             m(3)*(m(3)-1)*(x**(m(3)-2))*Coef(3)+&
&             m(4)*(m(4)-1)*(x**(m(4)-2))*Coef(4))/(rc**2)
         dddpp=(m(1)*(m(1)-1)*(m(1)-2)*(x**(m(1)-3))*Coef(1)+&
&             m(2)*(m(2)-1)*(m(2)-2)*(x**(m(2)-3))*Coef(2)+&
&             m(3)*(m(3)-1)*(m(3)-2)*(x**(m(3)-3))*Coef(3)+&
&             m(4)*(m(4)-1)*(m(4)-2)*(x**(m(4)-3))*Coef(4))/(rc**3)
         IF (i==irc) THEN
            WRITE(6,*) 'check  dp ', dpp,  B/rc
            WRITE(6,*) 'check ddp ', ddpp, C/rc**2
            WRITE(6,*) 'check dddp', dddpp,  D/rc**3
         ENDIF
         IF (wavetype==EXPF) THEN
            VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
            dum(i)=(r(i)**(l+1))*EXP(p(i))
         ENDIF
         IF (wavetype==POLY) THEN
            VNC(i)=e+(ddpp+2*(l+1)*dpp/r(i))/p(i)
            dum(i)=(r(i)**(l+1))*(p(i))
         ENDIF
         WRITE(88,'(1p,5e15.7)') r(i),wfn(i),dum(i),VNC(i)*r(i),rv(i)
      ENDDO
      CLOSE(88)
      x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
      WRITE(6,*) 'check norm ',x,S

      VNC(irc:n)=rv(irc:n)/r(irc:n)
      PAW%rveff(1:n)=VNC(1:n)*r(1:n)

      DEALLOCATE(VNC,wfn,p,dum)
    END SUBROUTINE kerker

  !***************************************************************
  ! SUBROUTINE nonncps(lmax,Grid,Pot)
  !  Creates  screened pseudopotential by inverting Schroedinger
  !    equation from a pseudized radial wave function of the form:
  !        Psi(r) = r**(l+1) * exp (a + b*r**2 + c*r**4 + d*r**6)
  !  No norm-conserving condition is imposed on Psi
  !***************************************************************
  SUBROUTINE nonncps(Grid,Pot,PAW,l,e)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
    INTEGER,INTENT(IN) :: l
    REAL(8),INTENT(IN) :: e

    INTEGER :: i,irc,n,nr,ok,nodes,i1,i2,i3,i4
    REAL(8) :: rc,x,y1,y2,y3,p0,p1,p2,p3,sgn
    REAL(8) :: b(4),c(4),d(4),amat(4,4)
    REAL(8),ALLOCATABLE ::  VNC(:),wfn(:)
    REAL(8),POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

   !Polynomial definitions
    p0(x,y1,y2,y3)=(x-y1)*(x-y2)*(x-y3)
    p1(x,y1,y2,y3)=(x-y2)*(x-y3)+(x-y1)*(x-y3)+(x-y1)*(x-y2)
    p2(x,y1,y2,y3)=2.0d0*((x-y1)+(x-y2)+(x-y3))
    p3(x,y1,y2,y3)=6.0d0

    n=Grid%n
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    ALLOCATE(VNC(n),wfn(nr),stat=ok)
    IF (ok/=0) stop 'Error in uspseudo -- allocating arrays'

    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot,nr,l,e,wfn,nodes)
    endif
    IF (wfn(irc)<0) wfn=-wfn

    DO i=2,nr
      wfn(i)=wfn(i)/r(i)**dble(l+1)
    ENDDO

    i1=irc-1;i2=i1+1;i3=i2+1;i4=i3+1
    c(1)=wfn(i1)/p0(r(i1),r(i2),r(i3),r(i4))
    c(2)=wfn(i2)/p0(r(i2),r(i3),r(i4),r(i1))
    c(3)=wfn(i3)/p0(r(i3),r(i4),r(i1),r(i2))
    c(4)=wfn(i4)/p0(r(i4),r(i1),r(i2),r(i3))
    d(1)=c(1)*p0(rc,r(i2),r(i3),r(i4)) + c(2)*p0(rc,r(i3),r(i4),r(i1)) + &
    &    c(3)*p0(rc,r(i4),r(i1),r(i2)) + c(4)*p0(rc,r(i1),r(i2),r(i3))
    d(2)=c(1)*p1(rc,r(i2),r(i3),r(i4)) + c(2)*p1(rc,r(i3),r(i4),r(i1)) + &
     &   c(3)*p1(rc,r(i4),r(i1),r(i2)) + c(4)*p1(rc,r(i1),r(i2),r(i3))
    d(3)=c(1)*p2(rc,r(i2),r(i3),r(i4)) + c(2)*p2(rc,r(i3),r(i4),r(i1)) + &
     &   c(3)*p2(rc,r(i4),r(i1),r(i2)) + c(4)*p2(rc,r(i1),r(i2),r(i3))
    d(4)=c(1)*p3(rc,r(i2),r(i3),r(i4)) + c(2)*p3(rc,r(i3),r(i4),r(i1)) + &
     &   c(3)*p3(rc,r(i4),r(i1),r(i2)) + c(4)*p3(rc,r(i1),r(i2),r(i3))

    sgn=d(1)/abs(d(1));d(1:4)=d(1:4)*sgn
    b(1)=log(d(1));b(2:4)=d(2:4)
    amat(1,1)= 1.0d0
    amat(2:4,1)= 0.0d0
    amat(1,2)= rc**2
    amat(2,2)= 2.0d0*d(1)*rc
    amat(3,2)= 2.0d0*d(1)   +2.0d0*d(2)*rc
    amat(4,2)=               4.0d0*d(2)   +2.0d0*d(3)*rc
    amat(1,3)= rc**4
    amat(2,3)=  4.0d0*d(1)*rc**3
    amat(3,3)= 12.0d0*d(1)*rc**2+ 4.0d0*d(2)*rc**3
    amat(4,3)= 24.0d0*d(1)*rc   +24.0d0*d(2)*rc**2+4.0d0*d(3)*rc**3
    amat(1,4)= rc**6
    amat(2,4)=   6.0d0*d(1)*rc**5
    amat(3,4)=  30.0d0*d(1)*rc**4+ 6.0d0*d(2)*rc**5
    amat(4,4)= 120.0d0*d(1)*rc**3+60.0d0*d(2)*rc**4+6.0d0*d(3)*rc**5

    CALL linsol(amat,b,4,4,4,4)
    write(6,*) 'Completed linsol with coefficients'
    write(6,'(1p,10e15.7)') (b(i),i=1,4)

    PAW%rveff(1)=0.d0
    DO i=2,irc-1
     c(1)=2.0d0*b(2)*r(i)+ 4.0d0*b(3)*r(i)**3+ 6.0d0*b(4)*r(i)**5
     c(2)=2.0d0*b(2)     +12.0d0*b(3)*r(i)**2+30.0d0*b(4)*r(i)**4
     PAW%rveff(i)=r(i)*(e+dble(2*l+2)*c(1)/r(i)+c(1)**2+c(2))
    ENDDO
    PAW%rveff(irc:n)=rv(irc:n)

    DEALLOCATE(VNC,wfn)

  END SUBROUTINE nonncps


  !***************************************************************
  ! SUBROUTINE besselps(lmax,Grid,Pot)
  !  Creates screened pseudopotential by simply pseudizing the
  !    AE potential with a l=0 spherical Bessel function:
  !                                     Vps(r) = a.sin(qr)/r
  !***************************************************************
  SUBROUTINE besselps(Grid,Pot,PAW)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW

    INTEGER :: i,irc,l,n
    REAL(8) :: e,rc,alpha,beta,vv,vvp,AA,QQ,xx(1)
    REAL(8),ALLOCATABLE ::  VNC(:),wfn(:)
    REAL(8),POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

    n=Grid%n
    r=>Grid%r
    rv=>Pot%rv
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    vv=rv(irc);vvp=Gfirstderiv(Grid,irc,rv)

    alpha=1.D0-rc*vvp/vv;beta=1.D0
    call solvbes(xx,alpha,beta,0,1);QQ=xx(1)
    AA=vv/sin(QQ);QQ=QQ/rc

    PAW%rveff(1)=0.d0
    PAW%rveff(irc+1:n)=rv(irc+1:n)
    do i=2,irc
     PAW%rveff(i)=AA*sin(QQ*r(i))
    enddo

  END SUBROUTINE besselps


    !***************************************************************
    ! SUBROUTINE EvaluateP
    !   Inverts 4x4 matrix used  by kerker subroutine
    !***************************************************************
    SUBROUTINE EvaluateP(m,A,B,C,D,coef)
      INTEGER, INTENT(IN) :: m(4)
      REAL(8), INTENT(IN) :: A,B,C,D
      REAL(8), INTENT(OUT) ::  coef(4)

      REAL(8) :: t(4,4)
      INTEGER :: i,n

      t=0
      Coef(1)=A; Coef(2)=B;  Coef(3)=C;    Coef(4)=D
      t(1,1:4)=1
      t(2,1:4)=m(1:4)
      DO i=1,4
         t(3,i)=m(i)*(m(i)-1)
      ENDDO
      DO i=1,4
         t(4,i)=m(i)*(m(i)-1)*(m(i)-2)
      ENDDO
      n=4
      CALL linsol(t,Coef,n,4,4,4)
    END SUBROUTINE EvaluateP

    !***************************************************************
    ! SUBROUTINE EvaluateTp
    !   Inverts 5x5 matrix used  by troullier subroutine
    !***************************************************************
    SUBROUTINE EvaluateTp(l,A,B,C,D,F,coef)
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(IN) :: A,B,C,D,F
      REAL(8), INTENT(OUT) ::  coef(6)

      REAL(8) :: t(6,6),coef10,old
      REAL(8), PARAMETER :: small=1.e-10
      INTEGER :: i,n,iter
      INTEGER, PARAMETER :: niter=1000

      old=-1.e30; Coef10=-1; iter=-1
      DO WHILE (iter < niter .AND. ABS(old-coef10)> small)
         iter=iter+1
         t=0
         Coef(1)=A-Coef10; Coef(2)=B-2*Coef10;  Coef(3)=C-2*Coef10;
         Coef(4)=D;    Coef(5)=F
         Coef(6)=-Coef10**2
         DO i=1,6
            t(1,i)=1
            t(2,i)=2*i
            t(3,i)=2*i*(2*i-1)
            t(4,i)=2*i*(2*i-1)*(2*i-2)
            t(5,i)=2*i*(2*i-1)*(2*i-2)*(2*i-3)
         ENDDO
         t(6,1)=2*Coef10;  t(6,2)=2*l+5

         n=6
         CALL linsol(t,Coef,n,6,6,6)

         old=Coef10; Coef10=Coef10+Coef(1)
         WRITE(6,'("EvaluateTp: iter",i5,1p,2e15.7)') iter,Coef(1),Coef10
         WRITE(6,'("Coef: ",1p,6e15.7)')Coef10,(Coef(i),i=2,6)
         Coef(1)=Coef10
      ENDDO

      IF (iter >= niter) THEN
         WRITE(6,*) 'Error in EvaluateTP -- no convergence'
         STOP
      ENDIF
    END SUBROUTINE EvaluateTp

    SUBROUTINE checkghosts(Grid,Orbit,FC,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(OrbitInfo), INTENT(IN) :: Orbit
      TYPE(FCInfo), INTENT(IN) :: FC
      TYPE(PseudoInfo), INTENT(in) :: PAW

      INTEGER :: l,nr,nodes,i,io
      REAL(8) :: energy,h
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: wfn(:),VNC(:)
      TYPE(PotentialInfo) :: Pot

      h=Grid%h
      r=>Grid%r
      nr=min(PAW%irc+5,Grid%n)
      ALLOCATE(VNC(nr),wfn(nr),POT%rv(nr),stat=i)
      IF (i /= 0) THEN
         WRITE(6,*) 'Error in checkghosts allocation',nr,i
         STOP
      ENDIF

      POT%rv(1:nr)=PAW%rveff(1:nr)
      call zeropot(Grid,POT%rv,POT%v0,POT%v0p)

      DO l=0,PAW%lmax
         DO io=1,Orbit%norbit
            IF((.NOT.FC%iscore(io)).AND.(Orbit%l(io)==l)) THEN
               energy=Orbit%eig(io)
               WRITE(6,*) 'Check  for ghosts with  l', l,energy
                  CALL unboundsch(Grid,Pot,nr,l,energy,wfn,nodes)
               !DO i=1,nr
               !   WRITE(l+17,'(1p,2e15.7)') Grid%r(i),wfn(i)
               !ENDDO
               WRITE(6,*) '    Found # nodes = ', nodes
               EXIT
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(VNC,wfn,POT%rv)
    END SUBROUTINE checkghosts

    SUBROUTINE InitPseudopot(Grid,PAW,Orbit,FC)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      TYPE(OrbitInfo), INTENT(IN) :: Orbit
      TYPE(FCInfo), INTENT(IN) :: FC

      INTEGER :: io,l,n,mxbase,nbase,ok

!     Compute initial size of basis
      n=Grid%n
      nbase=0
      DO l=0,PAW%lmax
         DO io=1,Orbit%norbit    ! cycle through all configurations
            IF (Orbit%l(io).EQ.l.AND.(.NOT.FC%iscore(io))) THEN
               nbase=nbase+1
            ENDIF
         ENDDO
      ENDDO
      mxbase=nbase+5*max(1,PAW%lmax)   !estimate excess
      PAW%nbase=nbase
      WRITE(6,*) 'Found ', nbase,' valence basis functions '
      WRITE(6,*) 'Allocating for ', mxbase, ' total basis functions'

      ALLOCATE(PAW%projshape(n),PAW%hatden(n),PAW%hatpot(n),&
&          PAW%hatshape(n),PAW%vloc(n),PAW%rveff(n),PAW%abinitvloc(n),&
&          PAW%abinitnohat(n),PAW%den(n),PAW%tden(n),PAW%tcore(n),stat=ok)
      IF (ok /=  0) WRITE(6,*) 'Allocation error 1 in initpseudopot',n,ok
      PAW%projshape=0.d0;PAW%hatden=0.d0;PAW%hatpot=0.d0
      PAW%hatshape=0.d0;PAW%vloc=0.d0;PAW%rveff=0.d0
      PAW%abinitvloc=0.d0;PAW%abinitnohat=0.d0
      PAW%den=0.d0;PAW%tden=0.d0;PAW%tcore=0.d0

      ALLOCATE(PAW%tvij(mxbase,mxbase),PAW%vhatij(mxbase,mxbase),&
&          PAW%v0ij(mxbase,mxbase),PAW%kin(mxbase,mxbase),&
&          PAW%vhijkl(mxbase,mxbase,mxbase,mxbase),&
&          PAW%multipole(mxbase,mxbase,2*PAW%lmax+1),stat=ok)
      IF (ok /= 0) WRITE(6,*) 'Allocation error 2 in initpseudopot',nbase,ok
      PAW%tvij=0.d0;PAW%vhatij=0.d0;PAW%v0ij=0.d0
      PAW%kin=0.d0;PAW%vhijkl=0.d0;PAW%multipole=0.d0

      ALLOCATE(PAW%phi(n,mxbase),PAW%tphi(n,mxbase),PAW%tp(n,mxbase),&
&          PAW%ophi(n,mxbase),PAW%otphi(n,mxbase),PAW%otp(n,mxbase),&
&          PAW%np(mxbase),PAW%l(mxbase),PAW%eig(mxbase),PAW%occ(mxbase),&
&          PAW%ck(mxbase),PAW%vrc(mxbase),PAW%Kop(n,mxbase),stat=ok)
      IF (ok /= 0) WRITE(6,*) 'Allocation error 3 in initpseudopot',n,mxbase,ok
      PAW%phi=0.d0;PAW%tphi=0.d0;PAW%tp=0.d0
      PAW%ophi=0.d0;PAW%otphi=0.d0;PAW%otp=0.d0
      PAW%eig=0.d0;PAW%occ=0.d0;PAW%vrc=0.d0;PAW%ck=0.d0;PAW%Kop=0.d0
      PAW%np=0;PAW%l=0

      ALLOCATE(PAW%oij(mxbase,mxbase),PAW%dij(mxbase,mxbase),stat=ok)
      IF (ok/=0) WRITE(6,*) 'Allocation error 4 in initpseudopot',mxbase,ok
      PAW%oij=0.d0; PAW%dij=0.d0

    END SUBROUTINE InitPseudopot

    SUBROUTINE DestroyPseudopot(PAW)
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      IF (ASSOCIATED(PAW%projshape)) DEALLOCATE(PAW%projshape)
      IF (ASSOCIATED(PAW%hatden)) DEALLOCATE(PAW%hatden)
      IF (ASSOCIATED(PAW%hatpot)) DEALLOCATE(PAW%hatpot)
      IF (ASSOCIATED(PAW%hatshape)) DEALLOCATE(PAW%hatshape)
      IF (ASSOCIATED(PAW%projshape)) DEALLOCATE(PAW%projshape)
      IF (ASSOCIATED(PAW%vloc)) DEALLOCATE(PAW%vloc)
      IF (ASSOCIATED(PAW%rveff)) DEALLOCATE(PAW%rveff)
      IF (ASSOCIATED(PAW%abinitvloc)) DEALLOCATE(PAW%abinitvloc)
      IF (ASSOCIATED(PAW%abinitnohat)) DEALLOCATE(PAW%abinitnohat)
      IF (ASSOCIATED(PAW%den)) DEALLOCATE(PAW%den)
      IF (ASSOCIATED(PAW%tden)) DEALLOCATE(PAW%tden)
      IF (ASSOCIATED(PAW%tcore)) DEALLOCATE(PAW%tcore)
      IF (ASSOCIATED(PAW%phi)) DEALLOCATE(PAW%phi)
      IF (ASSOCIATED(PAW%tphi)) DEALLOCATE(PAW%tphi)
      IF (ASSOCIATED(PAW%tp)) DEALLOCATE(PAW%tp)
      IF (ASSOCIATED(PAW%ophi)) DEALLOCATE(PAW%ophi)
      IF (ASSOCIATED(PAW%otphi)) DEALLOCATE(PAW%otphi)
      IF (ASSOCIATED(PAW%otp)) DEALLOCATE(PAW%otp)
      IF (ASSOCIATED(PAW%np)) DEALLOCATE(PAW%np)
      IF (ASSOCIATED(PAW%l)) DEALLOCATE(PAW%l)
      IF (ASSOCIATED(PAW%eig)) DEALLOCATE(PAW%eig)
      IF (ASSOCIATED(PAW%occ)) DEALLOCATE(PAW%occ)
      IF (ASSOCIATED(PAW%ck)) DEALLOCATE(PAW%ck)
      IF (ASSOCIATED(PAW%vrc)) DEALLOCATE(PAW%vrc)
      IF (ASSOCIATED(PAW%Kop)) DEALLOCATE(PAW%Kop)
      IF (ASSOCIATED(PAW%oij)) DEALLOCATE(PAW%oij)
      IF (ASSOCIATED(PAW%dij)) DEALLOCATE(PAW%dij)
      IF (ASSOCIATED(PAW%multipole)) DEALLOCATE(PAW%multipole)
      IF (ASSOCIATED(PAW%v0ij)) DEALLOCATE(PAW%v0ij)
      IF (ASSOCIATED(PAW%tvij)) DEALLOCATE(PAW%tvij)
      IF (ASSOCIATED(PAW%vhatij)) DEALLOCATE(PAW%vhatij)
      IF (ASSOCIATED(PAW%vhijkl)) DEALLOCATE(PAW%vhijkl)
      IF (ASSOCIATED(PAW%kin)) DEALLOCATE(PAW%kin)
    END SUBROUTINE DestroyPseudopot

    SUBROUTINE sethat(Grid,PAW,gaussparam,besselopt)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER,INTENT(IN), OPTIONAL :: besselopt
      REAL(8),INTENT(IN), OPTIONAL :: gaussparam

      INTEGER :: n,irc,irc_shap,i
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,con,rc,rc_shap,selfen,d,dd,jbes1,jbes2,qr
      REAL(8) :: al(2),ql(2)

      n=Grid%n
      h=Grid%h
      irc=PAW%irc
      rc=PAW%rc
      irc_shap=PAW%irc_shap
      rc_shap=PAW%rc_shap
      r=>Grid%r

      PAW%hatden=0
      PAW%projshape=0
      PAW%hatshape=0
      PAW%projshape(1)=1
      PAW%hatshape(1)=1
      DO i=2,irc-1
       PAW%projshape(i)=(SIN(pi*r(i)/rc)/(pi*r(i)/rc))**2
      ENDDO
      if(present(gaussparam)) then
       d=rc_shap/SQRT(LOG(1.d0/gaussparam))
       DO i=2,irc
        PAW%hatshape(i)=EXP(-(r(i)/d)**2)
       ENDDO
       PAW%irc_shap=PAW%irc
       PAW%rc_shap=PAW%rc
      else if(present(besselopt)) then
       call shapebes(al,ql,0,rc_shap)
       DO i=1,irc_shap-1
        qr=ql(1)*r(i);CALL jbessel(jbes1,d,dd,0,0,qr)
        qr=ql(2)*r(i);CALL jbessel(jbes2,d,dd,0,0,qr)
        PAW%hatshape(i)=al(1)*jbes1+al(2)*jbes2
       ENDDO
      else
       DO i=2,irc_shap-1
        PAW%hatshape(i)=(SIN(pi*r(i)/rc_shap)/(pi*r(i)/rc_shap))**2
       ENDDO
      endif
      PAW%hatden(1:irc)=PAW%hatshape(1:irc)*(r(1:irc)**2)

      !  normalize
      if (.not.besselshapefunction) then
       con=integrator(Grid,PAW%hatden,1,PAW%irc_shap)
       WRITE(6,*) ' check hatden normalization', con
       PAW%hatden=PAW%hatden/con
      endif

      CALL poisson(Grid,con,PAW%hatden,PAW%hatpot,selfen)
      WRITE(6,*) 'Self energy for L=0 hat density  ', selfen

    END SUBROUTINE sethat

    SUBROUTINE coretailselfenergy(Grid,PAW,ctctse,cthatse)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      REAL(8), INTENT(OUT) :: ctctse,cthatse

      INTEGER :: i,irc,n
      REAL(8) :: rc,h,x,y,z
      REAL(8), allocatable :: d1(:),d2(:)

      n=Grid%n
      h=Grid%h
      irc=PAW%irc

      allocate(d1(n),d2(n),stat=i)
          if (i /= 0) then
             write(6,*) 'coretailselfenergy: allocation error -- ', n,i
             stop
          endif

      x=integrator(Grid,PAW%tcore)
      write(6,*) 'tcore charge ' , x
      CALL poisson(Grid,x,PAW%tcore,d1,ctctse)
      d2(2:n)=PAW%hatden(2:n)*d1(2:n)/Grid%r(2:n)
      d2(1)=0
      cthatse=integrator(Grid,d2(1:irc),1,PAW%irc_shap)
      write(6,*) 'ctctse,cthatse = ', ctctse,cthatse

      deallocate(d1,d2)

    END SUBROUTINE coretailselfenergy


    SUBROUTINE setcoretail(Grid,coreden,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: coreden(:)
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      REAL(8) :: rc,h,x,y,z,u0,u2,u4
      REAL(8), allocatable :: d1(:),d2(:)
      INTEGER :: i,j,k,n,irc

      n=Grid%n
      h=Grid%h
      irc=PAW%irc_core
      rc=PAW%rc_core

      allocate(d1(n),d2(n),stat=i)
          if (i /= 0) then
             write(6,*) 'setcoretail: allocation error -- ', n,i
             stop
          endif
      CALL derivative(Grid,coreden,d1)
      CALL derivative(Grid,d1,d2)

      x=coreden(irc)
      y=d1(irc)*rc
      z=d2(irc)*(rc*rc)
      write(6,*) 'setcoretail: x,y,z = ', x,y,z

      u0=3*x - 9*y/8 + z/8
      u2=-3*x + 7*y/4 - z/4
      u4=x - 5*y/8 + z/8

      write(6,*) 'setcoretail: u0,u2,u4 = ', u0,u2,u4

      PAW%tcore=coreden

      do i=1,irc
         x=(Grid%r(i)/rc)**2
         PAW%tcore(i)= x*(u0+x*(u2+x*u4))
      enddo

      deallocate(d1,d2)

    END SUBROUTINE setcoretail

    !**************************************************************
    ! subroutine hatpotL
    !   Calculates potential associated with L component
    !    of unit hat density
    !**************************************************************
    SUBROUTINE hatpotL(Grid,PAW,l,vhat)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(OUT) :: vhat(:)

      INTEGER :: n,irc,i
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: den(:),a(:)
      REAL(8) :: h,con
      REAL(8) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)

      n=Grid%n
      h=Grid%h
      r=>Grid%r

      irc=PAW%irc

      ALLOCATE(den(n),a(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in hatpotL allocation',n,i
         STOP
      ENDIF


      if (besselshapefunction) then
       call shapebes(al,ql,l,PAW%rc_shap)
       DO i=1,PAW%irc_shap
        qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
        qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
        den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
       ENDDO
       if (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
      else
       DO i=1,n
        den(i)=(r(i)**l)*PAW%hatden(i)
       ENDDO
       a(1:n)=den(1:n)*(r(1:n)**l)
       con=integrator(Grid,a,1,PAW%irc_shap)
       den=den/con
      endif

       a(1:n)=den(1:n)*(r(1:n)**l)
       write(6,*) 'hatpotL check l ', l,integrator(Grid,a)
      vhat=0

      CALL apoisson(Grid,l,n,den,vhat(1:n))

      ! apoisson returns vhat*r
      !DO i=1,n
      !   WRITE (78+l,'(i5,1p,5e15.7)') i,Grid%r(i),den(i),vhat(i)
      !ENDDO
      vhat(2:n)=vhat(2:n)/r(2:n)
      call extrapolate(Grid,vhat)

      DEALLOCATE(den,a)
    END SUBROUTINE hatpotL

    !**************************************************************
    ! subroutine hatL
    !   Calculates density associated with L component
    !    normalized to unity
    !**************************************************************
    SUBROUTINE hatL(Grid,PAW,l,dhat)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(OUT) :: dhat(:)

      INTEGER :: n,irc,i
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: den(:),a(:)
      REAL(8) :: h,con
      REAL(8) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)

      n=Grid%n
      h=Grid%h
      r=>Grid%r

      irc=PAW%irc

      ALLOCATE(den(n),a(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in hatL allocation',irc,i
         STOP
      ENDIF

      if (besselshapefunction) then
       call shapebes(al,ql,l,PAW%rc_shap)
       DO i=1,PAW%irc_shap
        qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
        qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
        den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
       ENDDO
       if (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
      else
       DO i=1,n
        den(i)=(r(i)**l)*PAW%hatden(i)
       ENDDO
       a(1:n)=den(1:n)*(r(1:n)**l)
       con=integrator(Grid,a,1,PAW%irc_shap)
       den=den/con
      endif

      dhat=0

      dhat(1:n)=den(1:n)

      DEALLOCATE(den,a)
    END SUBROUTINE hatL

    SUBROUTINE selfhatpot(Grid,PAW,l,eself)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(OUT) :: eself


      INTEGER :: n,irc,i
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: den(:),a(:)
      REAL(8) :: h,con
      REAL(8) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      irc=PAW%irc

      ALLOCATE(den(n),a(n),stat=i)
      IF (i/=0) THEN
         WRITE(6,*) 'Error in selfhatpot allocation',irc,i
         STOP
      ENDIF

      if (besselshapefunction) then
       call shapebes(al,ql,l,PAW%rc_shap)
       DO i=1,PAW%irc_shap
        qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
        qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
        den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
       ENDDO
       if (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
      else
       DO i=1,n
        den(i)=(r(i)**l)*PAW%hatden(i)
       ENDDO
       a(1:n)=den(1:n)*(r(1:n)**l)
       con=integrator(Grid,a,1,PAW%irc_shap)
       den=den/con
      endif

      a=0

      CALL apoisson(Grid,l,n,den,a)

      ! apoisson returns a*r
      a(2:n)=a(2:n)/r(2:n)
      a(1)=0

      eself=0.5d0*overlap(Grid,a,den)

      DEALLOCATE(den,a)

    END SUBROUTINE selfhatpot

    !***********************************************************************88
    ! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
    ! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
    !   for r > rc, f1=t1, f2=t2
    ! on output: qqqq is difference overlap matrix element
    !   qqqq=<f1|f2>-<t1|t2>
    !***********************************************************************88
    SUBROUTINE dqij(Grid,PAW,ib,ic,qqqq)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: ib,ic
      REAL(8), INTENT(OUT) :: qqqq

      INTEGER :: n,i,ok,irc
      REAL(8) :: h
      REAL(8), ALLOCATABLE :: dum(:)

      qqqq=0
      IF (PAW%l(ib)/=PAW%l(ic)) RETURN
      n=Grid%n; h=Grid%h;  irc=PAW%irc
      ALLOCATE(dum(n),stat=ok)
      IF (ok /=0) THEN
         WRITE(6,*) 'Error in dqij allocation', n,ok
         STOP
      ENDIF
      DO i=1,n
         dum(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)-PAW%otphi(i,ib)*PAW%otphi(i,ic)
      ENDDO
      qqqq=integrator(Grid,dum,1,irc)

      DEALLOCATE(dum)
    END SUBROUTINE dqij

   !***********************************************************************
   ! SUBROUTINE dtij
   ! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
   ! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
   !   for r > rc, f1=t1, f2=t2
   ! on output: tij is difference kinetic energy matrix element in Rydberg units
   !   tij =<f1|T|f2>-<t1|T|t2>
   !************************************************************************
    SUBROUTINE dtij(Grid,PAW,ib,ic,tij)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: ib,ic
      REAL(8), INTENT(OUT) :: tij

      INTEGER :: n,i,ok,l,irc
      REAL(8) :: angm
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: dum(:),del1(:),del2(:),tdel1(:),tdel2(:)

      tij=0
      IF (PAW%l(ib)/=PAW%l(ic)) RETURN
      n=Grid%n;  r=>Grid%r;  l=PAW%l(ib);  irc=PAW%irc
      ALLOCATE(dum(n),del1(n),tdel1(n),del2(n),tdel2(n),stat=ok)
      IF (ok /=0) THEN
         WRITE(6,*) 'Error in dtij allocation', n,ok
         STOP
      ENDIF
      CALL derivative(Grid,PAW%ophi(:,ib),del1)
      CALL derivative(Grid,PAW%ophi(:,ic),del2)
      CALL derivative(Grid,PAW%otphi(:,ib),tdel1)
      CALL derivative(Grid,PAW%otphi(:,ic),tdel2)
      dum=0 ;   angm=l*(l+1)
      DO i=1,irc
         dum(i)=del1(i)*del2(i)-tdel1(i)*tdel2(i)
      ENDDO
         del1=0;del2=0;tdel1=0;tdel2=0
         del1(2:irc)=PAW%ophi(2:irc,ib)/Grid%r(2:irc)
         del2(2:irc)=PAW%ophi(2:irc,ic)/Grid%r(2:irc)
         tdel1(2:irc)=PAW%otphi(2:irc,ib)/Grid%r(2:irc)
         tdel2(2:irc)=PAW%otphi(2:irc,ic)/Grid%r(2:irc)
      DO i=1,irc
         dum(i)=dum(i)+angm*(del1(i)*del2(i)-tdel1(i)*tdel2(i))
      ENDDO
      tij=integrator(Grid,dum,1,irc)

      DEALLOCATE(dum,del1,del2,tdel1,tdel2)
    END SUBROUTINE dtij

   !***********************************************************************
   ! SUBROUTINE altdtij
   ! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
   ! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
   !   for r > rc, f1=t1, f2=t2
   ! on output: tij is difference kinetic energy matrix element in Rydberg units
   !   tij =<f1|T|f2>-<t1|T|t2>
   !************************************************************************
    SUBROUTINE altdtij(Grid,PAW,ib,ic,tij)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: ib,ic
      REAL(8), INTENT(OUT) :: tij

      INTEGER :: n,i,ok,l,irc
      REAL(8) :: angm
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: dum(:),tdel1(:),tdel2(:)

      tij=0
      IF (PAW%l(ib)/=PAW%l(ic)) RETURN
      n=Grid%n;  r=>Grid%r;  l=PAW%l(ib);  irc=PAW%irc
      ALLOCATE(dum(n),tdel1(n),tdel2(n),stat=ok)
      IF (ok /=0) THEN
         WRITE(6,*) 'Error in dtij allocation', n,ok
         STOP
      ENDIF
      dum=0
      do i=2,irc
        !dum(i)=(PAW%eig(ic)-AEPot%rv(i)/Grid%r(i))*PAW%ophi(i,ib)*PAW%ophi(i,ic)
        dum(i)=PAW%ophi(i,ib)*PAW%Kop(i,ic)
      enddo
      CALL derivative(Grid,PAW%otphi(:,ic),tdel1)
      CALL derivative(Grid,tdel1,tdel2)
      angm=l*(l+1)
      DO i=2,irc
         dum(i)=dum(i)+PAW%otphi(i,ib)*(tdel2(i)-&
&                   angm*PAW%otphi(i,ic)/(Grid%r(i)**2))
      ENDDO
      tij=integrator(Grid,dum,1,irc)

      DEALLOCATE(dum,tdel1,tdel2)
    END SUBROUTINE altdtij

    SUBROUTINE dvij(Grid,PAW,FC,nz,ib,ic,vij)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      TYPE(FCInfo), INTENT(IN) :: FC
      REAL(8), INTENT(IN) :: nz
      INTEGER, INTENT(IN) :: ib,ic
      REAL(8), INTENT(OUT) :: vij

      INTEGER :: n,i,ok,irc
      REAL(8) :: h,en,q,qt
      REAL(8), ALLOCATABLE :: dum(:),d1(:)
      REAL(8), POINTER :: r(:)

      vij=0
      IF (PAW%l(ib)/=PAW%l(ic)) RETURN
      n=Grid%n; h=Grid%h;  r=>Grid%r;  irc=PAW%irc
      ALLOCATE(dum(n),d1(n),stat=ok)
      IF (ok /=0) THEN
         WRITE(6,*) 'Error in dvij allocation', n,ok
         STOP
      ENDIF

      q=integrator(Grid,FC%coreden)
      WRITE(6,*) 'core electrons ',q,FC%zcore
      CALL poisson(Grid,q,FC%coreden,dum,en)
      dum=dum-2*nz
      qt=integrator(Grid,PAW%tcore)
      WRITE(6,*) 'coretail electrons ',qt
      CALL poisson(Grid,qt,PAW%tcore,d1,en)
      dum(1)=0
      DO i=2,irc
         dum(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)*dum(i)/r(i)-&
&             PAW%otphi(i,ib)*PAW%otphi(i,ic)*(PAW%vloc(i)+d1(i)/r(i))
      ENDDO
      vij=integrator(Grid,dum,1,irc)

      DEALLOCATE(dum)
    END SUBROUTINE dvij

    !****************************************************************
    ! SUBROUTINE avij -- potential part of Dij coefficients for
    !   estimating logderiv's
    !****************************************************************
    SUBROUTINE avij(Grid,Pot,PAW,ib,ic,vij)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: ib,ic
      REAL(8), INTENT(OUT) :: vij

      INTEGER :: n,i,ok,irc
      REAL(8) :: h,en,q
      REAL(8), ALLOCATABLE :: dum(:)
      REAL(8), POINTER :: r(:)

      vij=0
      IF (PAW%l(ib)/=PAW%l(ic)) RETURN
      n=Grid%n; h=Grid%h;  r=>Grid%r;  irc=PAW%irc
      ALLOCATE(dum(n),stat=ok)
      IF (ok /=0) THEN
         WRITE(6,*) 'Error in avij allocation', n,ok
         STOP
      ENDIF

      dum=0
      DO i=2,n
         dum(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)*Pot%rv(i)/r(i)-&
&             PAW%otphi(i,ib)*PAW%otphi(i,ic)*PAW%rveff(i)/r(i)
      ENDDO
      vij=integrator(Grid,dum,1,irc)

      DEALLOCATE(dum)
    END SUBROUTINE avij

    !********************************************************
    ! SUBROUTINE calcwij
    !  subroutine to accumulate the wij coefficience for an input
    !     smooth wavefunction twfn and occupancy and l
    !********************************************************
    SUBROUTINE calcwij(Grid,PAW,l,occ,twfn,wij)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(IN) :: twfn(:),occ
      REAL(8), INTENT(INOUT) :: wij(:,:)

      INTEGER :: n,i,ok,irc,ib,ic,nbase
      REAL(8) :: h
      REAL(8), ALLOCATABLE :: bm(:)

      n=Grid%n; h=Grid%h ;  irc=PAW%irc
      nbase=PAW%nbase

      ALLOCATE(bm(nbase))

      bm=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) bm(ib)=overlap(Grid,PAW%otp(:,ib),twfn,1,irc)
         !IF (l==PAW%l(ib)) write(6,*) 'accum wij',l,ib,bm(ib)
      ENDDO
      DO ib=1,nbase
         DO ic=1,nbase
            wij(ib,ic)=wij(ib,ic) + occ*bm(ib)*bm(ic)
         ENDDO
      ENDDO
      DEALLOCATE(bm)
    END SUBROUTINE calcwij

  END MODULE pseudo
