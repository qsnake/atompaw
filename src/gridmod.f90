MODULE gridmod
  Use globalmath

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: lineargrid=1  ! r(i)=h*(i-1)
  INTEGER, PARAMETER, PRIVATE :: loggrid=2  ! r(i)=r0*(exp(h*(i-1))-1)


  TYPE GridInfo
     INTEGER :: TYPE
     INTEGER :: n
     INTEGER :: ishift
     REAL(8) :: h
     REAL(8), POINTER :: r(:)
     REAL(8), POINTER :: drdu(:)    ! for loggrid -- dr/du
     REAL(8), POINTER :: pref(:)    ! for loggrid -- r0*exp(u/2)
     REAL(8), POINTER :: rr02(:)    ! for loggrid -- (r+r0)**2
  END TYPE GridInfo

CONTAINS

  !**********************************************************************
  ! function usingloggrid(Grid)
  !**********************************************************************
  FUNCTION usingloggrid(Grid)
       logical :: usingloggrid
       Type(GridInfo), INTENT(IN) :: Grid

       usingloggrid=.false.
       if (Grid%type==loggrid) usingloggrid=.true.
   END FUNCTION

  !**********************************************************************
  FUNCTION overint(n,h,f1,icorr)
    !   function to calculate the integral of one vectors f1
    !      using simpsons rule assuming a regular grid with
    !      spacing of h and n total points
    !
    !      icorr: optional parameter: used only when n is even
    !             if icorr<0,  a trapezoidal correction is applied
    !                          at the start of interval
    !             if icorr>=0, a trapezoidal correction is applied
    !                          at the end of interval
    !             default (if missing) is icorr=0
    REAL(8) :: overint
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: h,f1(:)
    INTEGER, OPTIONAL :: icorr


    REAL(8),PARAMETER :: tol=1.D-14
    INTEGER :: i,j,istart,m

    overint=0

    !Eliminate zeros at end of interval
    i=n;do while(abs(f1(i))<machine_zero.and.i>2);i=i-1;enddo
    m=min(i+1,n)

    IF (m<=1) THEN
       RETURN
    ELSEIF (m==2) THEN
       overint=(f1(1)+f1(2))*(h/2)   ! Trapezoidal rule
       RETURN
    ENDIF

    istart=1
    if (present(icorr)) then
     if (icorr<0.and.mod(m,2)==0) istart=2
    endif
    overint=f1(istart)+4*f1(istart+1)+f1(istart+2)
    j=((m-istart)/2)*2+istart
    IF (j>=istart+4) THEN
       DO i=istart+4,j,2
          overint=overint+f1(i-2)+4*f1(i-1)+f1(i)
       ENDDO
    ENDIF
    overint=overint*(h/3)
    IF (m>j) overint=overint+(f1(j)+f1(m))*(h/2)
    IF (istart==2) overint=overint+(f1(1)+f1(2))*(h/2)
    RETURN
  END FUNCTION overint

  !*******************************************************************
  ! function integrator(Grid,arg)
  !*******************************************************************
  FUNCTION integrator(Grid,arg,str,fin)
    REAL(8) :: integrator
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: arg(:)
    INTEGER, INTENT(IN), OPTIONAL :: str,fin

    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: i,n,i1,i2

    n=Grid%n
    i1=1;i2=n

    IF (PRESENT(str).AND.PRESENT(fin)) THEN
       i1=str; i2=fin; n=i2-i1+1
    ENDIF

    SELECT CASE(Grid%type)
    CASE default
       WRITE(6,*) 'Error in integrator -- grid ', Grid%type
       STOP
    CASE(lineargrid)
       integrator=overint(n,Grid%h,arg(i1:i2))
    CASE(loggrid)
       ALLOCATE(dum(i1:i2),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in integrator -- allocation ', Grid%n,i
          STOP
       ENDIF
       dum(i1:i2)=arg(i1:i2)*Grid%drdu(i1:i2)
       integrator=overint(n,Grid%h,dum(i1:i2),-1)
       DEALLOCATE(dum)
    END SELECT

  END FUNCTION integrator

  !*****************************************************************
  FUNCTION overlap(Grid,f1,f2,str,fin)
    !   function to calculate the overlap between two vectors f1 and f2
    REAL(8) :: overlap
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: f1(:),f2(:)
    INTEGER, INTENT(IN), OPTIONAL :: str,fin

    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: i,n,i1,i2

    n=Grid%n
    i1=1;i2=n
    IF (PRESENT(str).AND.PRESENT(fin)) THEN
       i1=str; i2=fin; n=i2-i1+1
    ENDIF


    ALLOCATE(dum(i1:i2),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in overlap allocation ', n,i
       STOP
    ENDIF
    dum(1:n)=f1(i1:i2)*f2(i1:i2)
    overlap=integrator(Grid,dum(1:n),1,n)
    DEALLOCATE(dum)
  END FUNCTION overlap
  !*****************************************************

  SUBROUTINE nderiv(h,y,z,ndim,ierr)
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(INOUT) :: ierr
    REAL(8) , INTENT(IN) :: h,y(:)
    REAL(8) , INTENT(INOUT) :: z(:)
    !     subroutine ddet5(h,y,z,ndim)
    !      ssp routine modified by nawh 6/8/76

    REAL(8) :: hh,yy,a,b,c
    INTEGER :: i

    ierr=-1
    IF (ndim.LT.5) RETURN
    !        prepare differentiation loop
    hh=.08333333333333333d0/h
    yy=y(ndim-4)
    b=hh*(-25.d0*y(1)+48.d0*y(2)-36.d0*y(3)+16.d0*y(4)-3.d0*y(5))
    c=hh*(-3.d0*y(1)-10.d0*y(2)+18.d0*y(3)-6.d0*y(4)+y(5))

    !
    !        start differentiation loop
    DO  i=5,ndim
       a=b
       b=c
       c=hh*(y(i-4)-y(i)+8.d0*(y(i-1)-y(i-3)))
       z(i-4)=a
    ENDDO
    !        end of differentiation loop
    !
    !        normal exit
    a=hh*(-yy+6.d0*y(ndim-3)-18.d0*y(ndim-2)+10.d0*y(ndim-1)          &
&        +3.d0*y(ndim))
    z(ndim)=hh*(3.d0*yy-16.d0*y(ndim-3)+36.d0*y(ndim-2)               &
&        -48.d0*y(ndim-1)+25.d0*y(ndim))
    z(ndim-1)=a
    z(ndim-2)=c
    z(ndim-3)=b
    !
    ierr=0
    RETURN
  END SUBROUTINE nderiv

  !**********************************************************************
  ! subroutine derivative(Grid,f,dfdr)
  !*********************************************************************
  SUBROUTINE derivative(Grid,f,dfdr,begin,bend)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: f(:)
    REAL(8), INTENT(OUT) :: dfdr(:)
    INTEGER, OPTIONAL, INTENT(IN) :: begin,bend

    INTEGER :: i,n,i1,i2

    i1=1;i2=Grid%n;n=i2-i1+1
    IF (PRESENT(begin).OR.PRESENT(bend)) THEN
       IF (begin>=1.AND.bend<= Grid%n) THEN
          i1=begin;i2=bend;n=i2-i1+1
       ELSE
          WRITE(6,*) 'Error in derivative', begin,bend,Grid%n
          STOP
       ENDIF
    ENDIF


    SELECT CASE(Grid%type)
    CASE default
       WRITE(6,*) 'Error in derivative -- grid ', Grid%type
       STOP
    CASE(lineargrid)
       CALL nderiv(Grid%h,f(i1:i2),dfdr(i1:i2),n,i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in derivative -nderiv problem', i
          STOP
       ENDIF
    CASE(loggrid)
       CALL nderiv(Grid%h,f(i1:i2),dfdr(i1:i2),n,i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in derivative -nderiv problem', i
          STOP
       ENDIF
       dfdr(i1:i2)=dfdr(i1:i2)/Grid%drdu(i1:i2)
    END SELECT

  END SUBROUTINE derivative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! simplederiv(Grid,f,dfdr,begin,bend)
!   low order formula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE simplederiv(Grid,f,dfdr,begin,bend)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: f(:)
    REAL(8), INTENT(OUT) :: dfdr(:)
    INTEGER, OPTIONAL, INTENT(IN) :: begin,bend

    INTEGER :: i,n,i1,i2
    REAL(8) :: HH

    i1=1;i2=Grid%n;n=i2-i1+1
    IF (PRESENT(begin).OR.PRESENT(bend)) THEN
       IF (begin>=1.AND.bend<= Grid%n) THEN
          i1=begin;i2=bend;n=i2-i1+1
          if (n<3) then
            WRITE(6,*) 'Error in simplederiv -- n too small',n,i1,i2
            stop
          endif
       ELSE
          WRITE(6,*) 'Error in simplederive', begin,bend,Grid%n
          STOP
       ENDIF
    ENDIF

    HH=0.5d0/Grid%h
    dfdr(i1)=HH*(-3*f(i1)+4*f(i1+1)-f(i1+2))
    do i=i1,i2-2
       dfdr(i+1)=HH*(f(i+2)-f(i))
    enddo
    dfdr(i2)=HH*(3*f(i2)-4*f(i2-1)+f(i2-2))

    if (Grid%type==loggrid) then
       dfdr(i1:i2)=dfdr(i1:i2)/Grid%drdu(i1:i2)
    endif

  end subroutine simplederiv

  SUBROUTINE laplacian(Grid,l,wfn,del,fin)
    TYPE(GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: wfn(:)
    REAL(8), INTENT(INOUT) :: del(:)
    INTEGER, INTENT(IN),OPTIONAL :: fin

    INTEGER :: OK,i,lfac,n
    REAL(8), ALLOCATABLE :: dum1(:),dum2(:)

    n=Grid%n
    if (present(fin)) n=fin
    ALLOCATE(dum1(n),dum2(n),stat=ok)
    IF (ok /= 0) THEN
       WRITE(6,*) 'Error in laplace allocation',n
       STOP
    ENDIF
    lfac=l*(l+1)
    del=0
    DO i=2,n
       dum2(i)=wfn(i)/Grid%r(i)
       del(i)=-lfac*dum2(i)
    ENDDO
    call extrapolate(Grid,dum2)
    call derivative(Grid,dum2,dum1,1,n)
    del(:)=del(:)+2*Grid%r(:)*dum1(:)
    call derivative(Grid,dum1,dum2,1,n)
    del(:)=del(:)+(Grid%r(:)**2)*dum2(:)

    DEALLOCATE(dum1,dum2)
  END SUBROUTINE laplacian

  !******************************************************
  ! SUBROUTINE poisson(Grid,q,den,rv,ecoul,v00)
  !*****************************************************
  SUBROUTINE poisson(Grid,q,den,rv,ecoul,v00)
    !  use Numerov algorithm to solve poisson equation
    !  den(n) is electron density * (4*pi*r**2)
    !  rv(n) is returned as electrostatic potential * r
    !  ecoul is the coulomb interaction energy

    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN):: den(:)
    REAL(8), INTENT(INOUT) :: rv(:),ecoul,q
    REAL(8), optional, INTENT(OUT) :: v00

    REAL(8), ALLOCATABLE :: a(:),b(:)
    REAL(8) :: sd,sl,h,h2
    INTEGER :: i,n

    n=Grid%n
    h=Grid%h

    rv=0
    q=integrator(Grid,den)

    ALLOCATE(a(n),b(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in allocating arrays in poisson',i,n
       STOP
    ENDIF

    IF (Grid%type==lineargrid) THEN

       sd=2
       sl=-1
       a(1)=0
       DO i=2,n
          a(i)=h*den(i)/(6*(i-1))
       ENDDO
       rv(1)=0
       rv(2)=10*a(2)+a(3)
       DO i=3,n-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       rv(n)=10*a(n)+a(n-1)+2*q

    ELSEIF (Grid%type==loggrid) THEN

       sd=2+10*h*h/48
       sl=-1+h*h/48
       a(1)=0
       h2=h*h
       DO i=2,n
          a(i)=h2*Grid%rr02(i)*den(i)/(6*Grid%r(i))/Grid%pref(i)
       ENDDO
       rv(1)=0
       rv(2)=10*a(2)+a(3)
       DO i=3,n-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       !   last term is boundary value at point n+1
       rv(n)=10*a(n)+a(n-1)-2*q*sl/(Grid%pref(n)*EXP(h/2))

    ENDIF

    CALL conthomas(n,sl,sd,rv)
    IF (Grid%type==loggrid) rv=rv*Grid%pref

    !
    !  calculate ecoul
    !
    DO i=2,n
       a(i)=den(i)*rv(i)/Grid%r(i)
    ENDDO
    a(1)=0
    ecoul=0.5d0*integrator(Grid,a)
    WRITE(6,*) ' from poisson: ecoul = ',ecoul

    if (present(v00)) then
       a=0
       a(2:n)=den(2:n)/Grid%r(2:n)
       v00=2*integrator(Grid,a)
    endif
    DEALLOCATE(a,b)
    !
  END SUBROUTINE poisson

!*****************************************************************
!Alternative form of poisson solver written by Marc Torrent 6/9/06
!  works well for loggrid
!******************************************************************
  SUBROUTINE poisson_marc(Grid,q,den,rv,ecoul)
    !  use Numerov algorithm to solve poisson equation
    !  den(n) is electron density * (4*pi*r**2)
    !  rv(n) is returned as electrostatic potential * r
    !  ecoul is the coulomb interation energy

    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN):: den(:)
    REAL(8), INTENT(INOUT) :: rv(:),ecoul,q

    REAL(8), ALLOCATABLE :: aa(:),bb(:),cc(:),dd(:)
    REAL(8) :: sd,sl,h,h2
    INTEGER :: i,n,ir,jr

    n=Grid%n
    h=Grid%h

    rv=0.d0
    q=integrator(Grid,den)

    ALLOCATE(aa(n),bb(n),cc(n),dd(n),stat=i)
    IF (i/=0) THEN
      WRITE(6,*) 'Error in allocating arrays in poisson ',i,n
      STOP
    ENDIF

    do jr=n,2,-1
     ir=n-jr+1
     aa(ir)=den(jr)*Grid%drdu(jr)
     bb(ir)=den(jr)*Grid%drdu(jr)/Grid%r(jr)
    end do

    cc=0.d0
    cc(5)=aa(n-4);cc(4)=aa(n-3);cc(3)=aa(n-2);cc(2)=aa(n-1)
    call extrapolate(Grid,cc);aa(n)=cc(1)
    cc(5)=bb(n-4);cc(4)=bb(n-3);cc(3)=bb(n-2);cc(2)=bb(n-1)
    call extrapolate(Grid,cc);bb(n)=cc(1)

    cc(1)=0.d0;dd(1)=0.d0
    do ir=3,n,2
     cc(ir)  =cc(ir-2)+h/3.d0*(aa(ir-2)+4.d0*aa(ir-1)+aa(ir))
     cc(ir-1)=cc(ir-2)+h/3.d0*(1.25d0*aa(ir-2)+2.0d0*aa(ir-1)-0.25d0*aa(ir))
     dd(ir)  =dd(ir-2)+h/3.d0*(bb(ir-2)+4.d0*bb(ir-1)+bb(ir))
     dd(ir-1)=dd(ir-2)+h/3.d0*(1.25d0*bb(ir-2)+2.d0*bb(ir-1)-0.25d0*bb(ir))
    end do
    if (mod(n,2)==0) then
     cc(n)=cc(n-2)+h/3.d0*(aa(n-2)+4.d0*aa(n-1)+aa(n))
     dd(n)=dd(n-2)+h/3.d0*(bb(n-2)+4.d0*bb(n-1)+bb(n))
    end if

    rv(1)=0.d0
    do ir=2,n
     jr=n-ir+1
     rv(ir)=2.d0*(dd(jr)*Grid%r(ir)+(cc(n)-cc(jr)))
    end do

    !  calculate ecoul
    aa(1)=0.d0
    do i=2,n
     aa(i)=den(i)*rv(i)/Grid%r(i)
    end do
    ecoul=0.5d0*integrator(Grid,aa)
!   WRITE(6,*) 'ecoul = ',ecoul

    deallocate(aa,bb,cc,dd)

 end subroutine poisson_marc

  !********************************************************************
  !  use Numerov algorithm to solve poisson equation
  !  for angularly dependent charge distribution of angular momentum l
  !  den(n) is electron density * (4*pi*r**2) appropriate for l
  !  rv(n) is returned as electrostatic potential * r
  !  a(n), b(n), and c(n) are work arrays
  !********************************************************************

  SUBROUTINE apoisson(Grid,l,irc,den,rv)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,irc
    REAL(8), INTENT(IN) :: den(:)
    REAL(8), INTENT(OUT) :: rv(:)

    INTEGER :: i,j
    REAL(8) :: angm,r,q,h,h2
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:)


    ALLOCATE(a(irc),b(irc),c(irc),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in allocating arrays in apoisson ',i,irc
       STOP
    ENDIF

    b=den(1:irc)*(Grid%r(1:irc))**L
    q=integrator(Grid,b,1,irc)/(2*l+1)
    h=Grid%h
    WRITE(6,*) 'check l multipole',l,q
    IF (Grid%type==lineargrid) THEN
       a=0;b=0; c=0; rv=0
       DO i=2,irc
          a(i)=0.2d0*h*den(i)/((i-1))
       ENDDO
       DO i=2,irc-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       !  set up Tridiagonal equations
       a=0; angm=l*(l+1)
       DO i=2,irc-1
          if (i>2) a(i)=-1.2d0+0.1d0*angm/((i-2)**2)
          c(i)=-1.2d0+0.1d0*angm/((i)**2)
          b(i)=2.4d0+angm/((i-1)**2)
       ENDDO
       a(2)=0
       IF (l==1) b(2)=b(2)+0.2d0
       rv(irc)=2*q/(Grid%r(irc)**l)
       rv(irc-1)=rv(irc-1)-c(irc-1)*rv(irc)
       c(irc-1)=0;rv(1)=0
       CALL thomas(irc-2,a(2:irc-1),b(2:irc-1),c(2:irc-1),rv(2:irc-1))
    ELSEIF (Grid%type==loggrid) THEN
       a=0;b=0; c=0; rv=0; h2=h*h
       DO i=2,irc
          a(i)=0.2d0*h2*Grid%rr02(i)*den(i)/(Grid%pref(i)*Grid%r(i))
       ENDDO
       DO i=2,irc-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       !  set up Tridiagonal equations
       a=0; angm=l*(l+1)
       DO i=2,irc-1
          if (i>2) a(i)=-1.2d0+&
&             0.1d0*h2*(0.25d0+angm*Grid%rr02(i-1)/Grid%r(i-1)**2)
          c(i)=-1.2d0+0.1d0*h2*(0.25d0+angm*Grid%rr02(i+1)/Grid%r(i+1)**2)
          b(i)=2.4d0+h2*(0.25d0+angm*Grid%rr02(i)/Grid%r(i)**2)
       ENDDO
       a(2)=0
       IF (l==1) b(2)=b(2)+0.2d0*Grid%rr02(1)/(Grid%r(2)**2)
       rv(irc)=2*q/(Grid%r(irc)**l)/Grid%pref(irc)
       rv(irc-1)=rv(irc-1)-c(irc-1)*rv(irc)
       c(irc-1)=0;rv(1)=0
       CALL thomas(irc-2,a(2:irc-1),b(2:irc-1),c(2:irc-1),rv(2:irc-1))
       rv(1:irc)=Grid%pref(1:irc)*rv(1:irc)
    ENDIF
    !
    DEALLOCATE(a,b,c)
  END SUBROUTINE apoisson

!******************************************************************
! pgm to determine r=0 form of potential assuming that
!   nuclear contribution (-2*nz) is not yet included
!******************************************************************

  SUBROUTINE zeropot(Grid,rv,v0,v0p)
    ! extrapolate potential to value at r=0
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(INOUT) :: rv(:)    ! Note: rv(1) corresponds to r=0
    REAL(8), INTENT(OUT) :: v0,v0p

    REAL(8) :: tmp(15),tmp1(15)
    INTEGER :: i,n

    tmp(2:15)=rv(2:15)/Grid%r(2:15)
    call extrapolate(Grid,tmp(1:15))
    v0=tmp(1)

    CALL derivative(Grid,tmp(1:15),tmp1(1:15),2,15)

    call extrapolate(Grid,tmp1(1:15))
    v0p=tmp1(1)
  END SUBROUTINE zeropot

  SUBROUTINE extrapolate(Grid,v)
   ! extrapolate array v to r=0 at v(1)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(INOUT) :: v(:)  ! assume v(2),v(3)...  given

    v(1)=v(4)+3.d0*(v(2)-v(3))

  END SUBROUTINE extrapolate

  !*************************************************************
  ! subroutine forward_numerov(Grid,l,many,energy,rv,zeroval,wfn,nodes)
  !*************************************************************
  SUBROUTINE forward_numerov(Grid,l,many,energy,rv,zeroval,wfn,nodes)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,many
    REAL(8), INTENT(IN) :: energy,zeroval,rv(:)
    ! zeroval == lim r-->0 A(r)*P(r)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(1) and wfn(2) given
    ! on ouput wfn(i) given for
    !      i<=many
    INTEGER, INTENT(OUT) :: nodes

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    ALLOCATE(a(many),b(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error forward_numerov ', many,i
       STOP
    ENDIF

    p(1)=wfn(1)
    p(2)=wfn(2)
    xx=zeroval
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(1:2)=Grid%pref(1:2)*wfn(1:2)
       xx=Grid%pref(1)*Grid%rr02(1)*xx
       a=0.25d0+Grid%rr02(1:many)*a
    ENDIF

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a
    p(3)=(b(2)*p(2)+0.1d0*h2*xx)/a(3)

    nodes=0
    DO i=4,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2))/a(i)
       IF (p(i)*p(i-1) < 0.d0) nodes=nodes+1
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(1:i)=scale*p(1:i)
       ENDIF
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE forward_numerov

  !*************************************************************
  ! subroutine shifted_forward_numerov(Grid,many,istart,ww,wfn,nodes)
  !     designed for use with scalar relativistic case in which
  !       l and energy information is stored in ww
  !*************************************************************
  SUBROUTINE shifted_forward_numerov(Grid,many,istart,ww,wfn,nodes)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: many,istart
    REAL(8), INTENT(IN) :: ww(:)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(1:5) given
    ! on ouput wfn(i) given for
    !      i<=many
    INTEGER, INTENT(OUT) :: nodes

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    if (istart>many) then
         write(6,*) 'shifted_forward_numerov:  error istart many',istart,many
         stop
    endif
    ALLOCATE(a(many),b(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error forward_numerov ', many,i
       STOP
    ENDIF

    p(1:istart)=wfn(1:istart)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=ww(i)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(1:istart)=Grid%pref(1:istart)*wfn(1:istart)
       a=0.25d0+Grid%rr02(1:many)*a
    ENDIF

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a

    nodes=0
    DO i=istart+1,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2))/a(i)
       IF (p(i)*p(i-1) < 0.d0) nodes=nodes+1
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(1:i)=scale*p(1:i)
       ENDIF
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE shifted_forward_numerov
   !**********************************************************************
   !  pgm to integrate outward the radial schroedinger inhomogeneous equation
   !    at energy 'energy' and at angular momentum l
   !    with potential smooth rv/r
   !        proj/r == inhomogeneous function
   !     It is assumed that for r~0, proj~~(r**(l+1))*(c0+r**2*c2+...)
   !       and wfn~C*r**(l+3)*polynomial(r) for r~0;
   !
   !  uses Noumerov algorithm
   !*************************************************************************
  !*************************************************************
  ! subroutine inhomogeneous_numerov(Grid,l,many,energy,rv,proj,wfn)
  !*************************************************************
  SUBROUTINE inhomogeneous_numerov(Grid,l,many,energy,rv,proj,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,many
    REAL(8), INTENT(IN) :: energy,rv(:),proj(:)
    REAL(8), INTENT(OUT) :: wfn(:)
             ! initial values of wfn determined from proj(r=0)
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    INTEGER :: i,j,k,n

    ALLOCATE(a(many),b(many),c(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error forward_numerov ', many,i
       STOP
    ENDIF

    a=0
    a(2:3)=proj(2:3)/(Grid%r(2:3)**(l+1))
    call extrapolate(Grid,a)
    wfn=0
    wfn(2)=-a(1)*(Grid%r(2)**(l+3))/(4*l+6.d0)

    p(1)=wfn(1)
    p(2)=wfn(2)
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    b(1:many)=0.1d0*h2*proj(1:many)
    IF (Grid%type==loggrid) THEN
       p(1:2)=Grid%pref(1:2)*wfn(1:2)
       a=0.25d0+Grid%rr02(1:many)*a
       b=Grid%rr02(1:many)*b/Grid%pref(1:many)
    ENDIF

    c=0
    do i=2,many-1
      c(i)=10*b(i)+b(i-1)+b(i+1)
    enddo

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a
    p(3)=(b(2)*p(2)-c(2))/a(3)

    DO i=4,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2)-c(i-1))/a(i)
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,c,p)

  END SUBROUTINE inhomogeneous_numerov

  !*********************************************************
  ! subroutine backward_numerov(Grid,l,match,energy,rv,wfn)
  !*********************************************************
  SUBROUTINE backward_numerov(Grid,l,match,energy,rv,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,match
    REAL(8), INTENT(IN) :: energy,rv(:)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(n-1) and wfn(n) given
    ! on ouput wfn(i) given for
    !      start<=i<=n

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    n=Grid%n
    ALLOCATE(a(n),b(n),p(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error backward_numerov ', n,i
       STOP
    ENDIF

    p(n)=wfn(n)
    p(n-1)=wfn(n-1)
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=match,n
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(n-1:n)=Grid%pref(n-1:n)*p(n-1:n)
       a(match:n)=0.25d0+Grid%rr02(match:n)*a(match:n)
    ENDIF

    b(match:n)=2.4d0+h2*a(match:n)
    a(match:n)=1.2d0-0.1d0*h2*a(match:n)

    DO i=n-2,match,-1
       p(i)=(b(i+1)*p(i+1)-a(i+2)*p(i+2))/a(i)
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(i:n)=scale*p(i:n)
       ENDIF
    ENDDO

    wfn(match:n)=p(match:n)

    IF (Grid%type==loggrid) THEN
       wfn(match:n)=wfn(match:n)*Grid%pref(match:n)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE backward_numerov

! Subroutine from David Vanderbilt's USPS code, modified by Marc
!     Torrent and Francois Jollet, further modified by NAWH
!===========================================================================
!      subroutine cfdsol(zz,yy,jj1,jj2,mesh)
!===========================================================================

!     routine for solving coupled first order differential equations
!
!      d yy(x,1)
!      ---------   =  zz(x,1,1) * yy(x,1) + zz(x,1,2) * yy(2,1)
!         dx
!
!      d yy(x,2)
!      ---------   =  zz(x,2,1) * yy(x,1) + zz(x,2,2) * yy(2,1)
!         dx
!
!
!     using fifth order predictor corrector algorithm
!
!     routine integrates from jj1 to jj2 and can cope with both cases
!     jj1 < jj2 and jj1 > jj2.  first five starting values of yy must
!     be provided by the calling program.

   subroutine cfdsol(Grid,zz,yy,jj1,jj2)
      Type(gridinfo), INTENT(IN) :: Grid
      real(8), INTENT(IN):: zz(:,:,:)
      real(8), INTENT(INOUT):: yy(:,:)
      integer, INTENT(IN)  :: jj1,jj2


      real(8):: fa(0:5),fb(0:5),abp(1:5),amc(0:4)
      integer :: isgn,i,j,ip,mesh
      real(8):: arp,brp
      real(8), allocatable :: tmpz(:,:,:)
      real(8), parameter :: verylarge=1.d30
      real(8) :: scale

      mesh=size(yy(2,:))
      !write (6,*) ' in cdfdol with mesh jj1,j22', mesh, jj1,jj2
      if (size(zz(2,2,:))/=mesh) then
         write(6,*) 'cfdsol error - incompatible arrays', mesh,size(zz(2,2,:))
         stop
      endif
      isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
      if ( isgn .eq. + 1 ) then
        if ( jj1 .le. 5 .or. jj2 .gt. mesh ) then
          write(6,10) isgn,jj1,jj2,mesh
          call exit(1)
        endif
      elseif ( isgn .eq. - 1 ) then
        if ( jj1 .ge. ( mesh - 4 ) .or. jj2 .lt. 1 ) then
          write(6,10) isgn,jj1,jj2,mesh
          call exit(1)
        endif
      else
        write(6,10) isgn,jj1,jj2,mesh
      endif

  10  format(' ***error in subroutine difsol',/,&
     &' isgn =',i2,' jj1 =',i5,' jj2 =',i5,' mesh =',i5,&
     &' are not allowed')

      allocate(tmpz(2,2,mesh))
      tmpz=zz

         do i=1,2
            do j=1,2
               tmpz(i,j,:)=tmpz(i,j,:)*Grid%h
               if (Grid%TYPE==loggrid) tmpz(i,j,1:mesh)=tmpz(i,j,1:mesh)*Grid%drdu(1:mesh)
            enddo
         enddo

      abp(1) = 1901.d0 / 720.d0
      abp(2) = -1387.d0 / 360.d0
      abp(3) = 109.d0 / 30.d0
      abp(4) = -637.d0 / 360.d0
      abp(5) = 251.d0 / 720.d0
      amc(0) = 251.d0 / 720.d0
      amc(1) = 323.d0 / 360.d0
      amc(2) = -11.d0 / 30.d0
      amc(3) = 53.d0 / 360.d0
      amc(4) = -19.d0 / 720.d0

      do j = 1,5
        ip = jj1 - isgn * j
        fa(j) = tmpz(1,1,ip) * yy(1,ip) + tmpz(1,2,ip) * yy(2,ip)
        fb(j) = tmpz(2,1,ip) * yy(1,ip) + tmpz(2,2,ip) * yy(2,ip)
      enddo

      do j = jj1,jj2,isgn
        arp = yy(1,j-isgn)
        brp = yy(2,j-isgn)
        if (abs(arp)>verylarge.or.brp>verylarge) then
           scale=1.d0/(abs(arp)+abs(brp))
           arp=arp*scale
           brp=brp*scale
           fa(:)=fa(:)*scale; fb(:)=fb(:)*scale
           yy=yy*scale
        endif
        do  i = 1,5
          arp = arp + dble(isgn) * abp(i) * fa(i)
          brp = brp + dble(isgn) * abp(i) * fb(i)
        enddo

        fa(0) = tmpz(1,1,j) * arp + tmpz(1,2,j) * brp
        fb(0) = tmpz(2,1,j) * arp + tmpz(2,2,j) * brp
        yy(1,j) = yy(1,j-isgn)
        yy(2,j) = yy(2,j-isgn)
        do  i = 0,4,1
          yy(1,j) = yy(1,j) + dble(isgn) * amc(i) * fa(i)
          yy(2,j) = yy(2,j) + dble(isgn) * amc(i) * fb(i)
        enddo

        do i = 5,2,-1
          fa(i) = fa(i-1)
          fb(i) = fb(i-1)
        enddo
        fa(1) = tmpz(1,1,j) * yy(1,j) + tmpz(1,2,j) * yy(2,j)
        fb(1) = tmpz(2,1,j) * yy(1,j) + tmpz(2,2,j) * yy(2,j)
      enddo

     deallocate(tmpz)
 end subroutine cfdsol


  !*********************************************************
  ! subroutine mod_backward_numerov(Grid,match,ww,wfn)
  !    version modified for scalar relativistic case when
  !     l and energy terms represented in array ww
  !*********************************************************
  SUBROUTINE mod_backward_numerov(Grid,match,ww,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: match
    REAL(8), INTENT(IN) :: ww(:)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(n-1) and wfn(n) given
    ! on ouput wfn(i) given for
    !      start<=i<=n

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    n=Grid%n
    ALLOCATE(a(n),b(n),p(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error backward_numerov ', n,i
       STOP
    ENDIF

    p(n)=wfn(n)
    p(n-1)=wfn(n-1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=match,n
       a(i)=ww(i)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(n-1:n)=Grid%pref(n-1:n)*p(n-1:n)
       a(match:n)=0.25d0+Grid%rr02(match:n)*a(match:n)
    ENDIF

    b(match:n)=2.4d0+h2*a(match:n)
    a(match:n)=1.2d0-0.1d0*h2*a(match:n)

    DO i=n-2,match,-1
       p(i)=(b(i+1)*p(i+1)-a(i+2)*p(i+2))/a(i)
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(i:n)=scale*p(i:n)
       ENDIF
    ENDDO

    wfn(match:n)=p(match:n)

    IF (Grid%type==loggrid) THEN
       wfn(match:n)=wfn(match:n)*Grid%pref(match:n)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE mod_backward_numerov
  !******************************************************************
  !  subroutine kinetic(Grid,wfn,l,ekin)
  !       calculates expectation value of kinetic energy for wfn
  !        with orbital angular momentum l
  !        wfn == r*radialwfn in Schroedinger Equation
  !        assumes wfn=(constant)*r^(l+1) at small r
  !*****************************************************************
  SUBROUTINE kinetic(Grid,l,wfn,ekin)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn(:)
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(OUT) :: ekin

    REAL(8), ALLOCATABLE :: dfdr(:),arg(:)
    INTEGER :: i,n

    n=Grid%n

    ALLOCATE(dfdr(n),arg(n),stat=i)

    CALL derivative(Grid,wfn,dfdr)

    arg=0
    DO i=2,n
       arg(i)=wfn(i)/Grid%r(i)
    ENDDO

    DO i=1,n
       arg(i)=(dfdr(i))**2+(l*(l+1))*(arg(i))**2
    ENDDO

    ekin=integrator(Grid,arg)

    DEALLOCATE(dfdr,arg)
  END SUBROUTINE kinetic

  !******************************************************************
  !  subroutine altkinetic(Grid,wfn,energy,rv,ekin)
  !       calculates expectation value of kinetic energy for wfn
  !        with orbital wfn by integrating
  !          int(wfn**2 * (energy-rv/r), r=0..rmax)
  !*****************************************************************
  SUBROUTINE altkinetic(Grid,wfn,energy,rv,ekin)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn(:),rv(:),energy
    REAL(8), INTENT(OUT) :: ekin

    REAL(8), ALLOCATABLE :: arg(:)
    INTEGER :: i,n

    n=Grid%n

    ALLOCATE(arg(n),stat=i)

    arg=0
    DO i=2,n
       arg(i)=(wfn(i)**2)*(energy-rv(i)/Grid%r(i))
    ENDDO

    ekin=integrator(Grid,arg)

    DEALLOCATE(arg)
  END SUBROUTINE altkinetic

  !****************************************************************
  ! function FindGridIndex(Grid,rpoint)
  !****************************************************************
  FUNCTION FindGridIndex(Grid,rpoint)
    INTEGER :: FindGridIndex
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: rpoint

    REAL(8) :: r0

    FindGridIndex=0
    IF (Grid%type==lineargrid) THEN
       FindGridIndex=rpoint/Grid%h+1
       IF (Grid%h*(FindGridIndex-1)<rpoint-1.d-10) FindGridIndex=FindGridIndex+1
    ELSEIF (Grid%type==loggrid) THEN
       r0=Grid%drdu(1)
       FindGridIndex=LOG(rpoint/r0+1)/Grid%h+1
       IF (r0*EXP(Grid%h*(FindGridIndex-1))<rpoint-1.d-10) &
&           FindGridIndex=FindGridIndex+1
    ENDIF
  END FUNCTION FindGridIndex

  !****************
  ! finite difference second derivative
  !*****************
  FUNCTION secondderiv(index,f,h)
    REAL(8) :: secondderiv
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: f(:),h

    IF (index<3.OR.index+2>SIZE(f)) THEN
       WRITE(6,*) 'Error in secondderiv', index, SIZE(f)
       STOP
    ENDIF

    secondderiv=0

    secondderiv=-(f(index-2)+f(index+2))/12 + &
&        4*(f(index-1)+f(index+1))/3 - 5*f(index)/2
    secondderiv=secondderiv/(h*h)
  END FUNCTION secondderiv



  !****************
  ! finite difference first derivative
  !*****************
  FUNCTION firstderiv(index,f,h)
    REAL(8) :: firstderiv
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: f(:),h

    IF (index<3.OR.index+2>SIZE(f)) THEN
       WRITE(6,*) 'Error in firstderiv', index, SIZE(f)
       STOP
    ENDIF

    firstderiv=0

    firstderiv=(f(index-2)-f(index+2))/12-2*(f(index-1)-f(index+1))/3
    firstderiv=firstderiv/h
  END FUNCTION firstderiv


  !*****************
  ! Second derivative for general grid
  !*****************

  FUNCTION Gsecondderiv(Grid,index,g)
    REAL(8) :: Gsecondderiv
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: g(:)

    REAL(8), ALLOCATABLE :: f(:)
    INTEGER :: i,n

    IF (index<3.OR.index+2>SIZE(g)) THEN
       WRITE(6,*) 'Error in secondderiv', index, SIZE(f)
       STOP
    ENDIF

    Gsecondderiv=0

    IF (Grid%type==lineargrid) THEN
       Gsecondderiv=secondderiv(index,g,Grid%h)
    ELSEIF  (Grid%type==loggrid) THEN
       n=Grid%n
       ALLOCATE(f(index-2:index+2),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in secondderiv ', i,n
          STOP
       ENDIF
       f(index-2:index+2)=g(index-2:index+2)/Grid%pref(index-2:index+2)
       Gsecondderiv=secondderiv(3,f(index-2:index+2),Grid%h)
       Gsecondderiv=Grid%pref(index)*(Gsecondderiv-&
&           0.25d0*f(index))/Grid%rr02(index)
       DEALLOCATE(f)
    ENDIF

  END FUNCTION Gsecondderiv


  !*****************
  ! First  derivative for general grid
  !*****************

  FUNCTION Gfirstderiv(Grid,index,g)
    REAL(8) :: Gfirstderiv
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: g(:)

    REAL(8), ALLOCATABLE :: f(:)
    INTEGER :: i,n

    IF (index<3.OR.index+2>SIZE(g)) THEN
       WRITE(6,*) 'Error in firstderiv', index, SIZE(f)
       STOP
    ENDIF

    Gfirstderiv=0

    IF (Grid%type==lineargrid) THEN
       Gfirstderiv=firstderiv(index,g,Grid%h)
    ELSEIF  (Grid%type==loggrid) THEN
       n=Grid%n
       ALLOCATE(f(index-2:index+2),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in firstderiv ', i,n
          STOP
       ENDIF
       f(index-2:index+2)=g(index-2:index+2)/Grid%pref(index-2:index+2)
       Gfirstderiv=firstderiv(3,f(index-2:index+2),Grid%h)
       Gfirstderiv=Grid%pref(index)*(Gfirstderiv+&
&           0.5d0*f(index))/SQRT(Grid%rr02(index))
       DEALLOCATE(f)
    ENDIF

  END FUNCTION Gfirstderiv


  !*****************************************************************
  ! subroutine reportgrid(Grid,unit)
  !*****************************************************************
  SUBROUTINE reportgrid(Grid,unit)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: unit

    IF (Grid%type==lineargrid) THEN
       WRITE(unit,*) ' Radial integration grid is linear '
       WRITE(unit,*)  ' h = ', Grid%h,'   n = ', Grid%n
    ELSEIF (Grid%type==loggrid) THEN
       WRITE(unit,*) ' Radial integration grid is logarithmic '
       WRITE(unit,*)  'r0 = ', Grid%drdu(1),' h = ', Grid%h,'   n = ', Grid%n
    ENDIF
  END SUBROUTINE reportgrid


  !******************************************************************
  ! function gridindex(Grid,r)
  !******************************************************************
  FUNCTION gridindex(Grid,r)
    INTEGER :: gridindex
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: r

    gridindex=0
    if (Grid%type==lineargrid) then
       gridindex=r/Grid%h +0.1d0 +1
    elseif (Grid%type==loggrid) then
       gridindex=LOG(1.d0+r/Grid%drdu(1))/Grid%h +0.1d0 +1
    endif

  END function gridindex

  !*********************************************************************
  !  subroutine findh(Z,range,n,hval)
  !    find hval for fixed number of input grid points n in loggrid case
  !    assumes form r(i)=(h/Z)*(exp(h*(i-1))-1)
  !*********************************************************************
   SUBROUTINE findh(Z,range,n,hval)
     REAL(8), INTENT(IN) :: Z
     INTEGER, INTENT(IN) :: n
     REAL(8), INTENT(IN) :: range
     REAL(8), INTENT(INOUT) :: hval

     REAL(8) :: h0,dh,f,df
     INTEGER :: i,j,k
     INTEGER, parameter :: iter=1000
     REAL(8), parameter :: eps=1.e-15
     LOGICAL :: success

     h0=hval
     success=.false.
     do i=1,iter
        f=LOG(Z*range/h0+1.d0)/h0
        df=-f/h0-(Z*range/h0**3)/(Z*range/h0+1.d0)
        dh=(n-1-f)/df
        if (ABS(dh)< eps) then
           success=.true.
           exit
        endif
        if (h0+dh<0.d0) then
           h0=h0/2
        else
           h0=h0+dh
        endif
      enddo

      if (.not.success) then
        write(6,*) 'Warning in findh -- dh > eps ', dh,h0
      endif
      hval=h0

  end subroutine findh



  !******************************************************************
  ! subroutine initgrid(Grid,h,range,r0)
  !******************************************************************
  SUBROUTINE initgrid(Grid,h,range,r0)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    REAL(8), INTENT(IN) :: range
    REAL(8), INTENT(IN) :: h
    REAL(8), OPTIONAL, INTENT(IN) :: r0

    INTEGER :: i,n

    IF (PRESENT(r0)) THEN

       Grid%type=loggrid
       n=LOG(range/r0+1)/h+1
       Grid%ishift=5
       IF (r0*(EXP(h*(n-1))-1)<range-1.d-5) n=n+1
       Grid%h=h
       Grid%n=n
       WRITE(6,*) 'InitGrid: -- logarithmic ',n, h,range,r0
       ALLOCATE(Grid%r(n),Grid%drdu(n),Grid%pref(n),Grid%rr02(n),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Allocation error in initgrid ', n,i
          STOP
       ENDIF
       DO i=1,n
          Grid%r(i)=r0*(EXP(Grid%h*(i-1))-1)
          Grid%drdu(i)=r0*EXP(Grid%h*(i-1))
          Grid%pref(i)=r0*EXP(Grid%h*(i-1)/2.d0)
          Grid%rr02(i)=(Grid%r(i)+r0)**2
       ENDDO
    ELSE
       Grid%type=lineargrid
       n=range/h+1
       Grid%ishift=25
       IF (h*(n-1)<range-1.d-5) n=n+1
       Grid%n=n
       Grid%h=h
       WRITE(6,*) 'InitGrid: -- linear  ', n,h,range
       ALLOCATE(Grid%r(n),Grid%drdu(n),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Allocation error in initgrid ', n,i
          STOP
       ENDIF
       DO i=1,n
          Grid%r(i)=(Grid%h*(i-1))
          Grid%drdu(i)=1.d0
       ENDDO
       NULLIFY(Grid%pref)
       NULLIFY(Grid%rr02)
    ENDIF

  END SUBROUTINE initgrid

  !******************************************************************
  ! subroutine destroygrid(Grid)
  !******************************************************************
  SUBROUTINE DestroyGrid(Grid)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    IF (ASSOCIATED(Grid%r)) DEALLOCATE(Grid%r)
    IF (ASSOCIATED(Grid%drdu)) DEALLOCATE(Grid%drdu)
    IF (ASSOCIATED(Grid%pref)) DEALLOCATE(Grid%pref)
    IF (ASSOCIATED(Grid%rr02)) DEALLOCATE(Grid%rr02)
  END SUBROUTINE DestroyGrid

  !******************************************************************
  ! subroutine nullifygrid(Grid)
  !******************************************************************
  SUBROUTINE NullifyGrid(Grid)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    NULLIFY(Grid%r)
    NULLIFY(Grid%drdu)
    NULLIFY(Grid%pref)
    NULLIFY(Grid%rr02)
  END SUBROUTINE NullifyGrid

END MODULE gridmod
