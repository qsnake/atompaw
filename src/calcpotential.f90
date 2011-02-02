MODULE calcpotential
  USE gridmod
  USE excor
  USE atomdata, only : finitenucleus
  USE globalmath, only : pi,derf

  IMPLICIT NONE

  TYPE PotentialInfo
     REAL(8) :: nz       !  nz is nuclear charge
     CHARACTER(2) :: sym
     REAL(8) :: q,v0,v0p  !  q is total electron charge
                          !  v0,v0p are potential value and deriv at r=0
     REAL(8) , POINTER :: den(:),rv(:)
                       !  den(n) is electron density * (4*pi*r**2)
                       !  rv(n) is  veff * r
  END TYPE PotentialInfo




CONTAINS
  SUBROUTINE potential(Grid,Pot,ecoul,etxc,eexc)
    !  calculate veff for density functional theory from electron density
    !  ecoul and etxc are coulomb and exchange-correlation contributions
    !    to the total energy
    !  eexc is the total exchange energy (int(den*exc))
    !  v0=v(0)
    !  v0p=dv/dr(0)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
    REAL(8), INTENT(INOUT) :: ecoul,etxc,eexc

    REAL(8), POINTER :: den(:),rv(:)
    REAL(8), allocatable :: rvxc(:)

    INTEGER :: n
    REAL(8) :: h,q,v0,v0p,nz
    REAL(8) :: r,RR
    INTEGER :: i,j,k
    !  INTEGER :: fcount=0

    n=Grid%n
    h=Grid%h
    den=>Pot%den
    rv=>Pot%rv
    nz=Pot%nz
    v0=0;v0p=0


    !
    !  Coulomb contribution
    !
    CALL poisson(Grid,Pot%q,den,rv,ecoul,v0)
    !
    !  Exchange-correlation contribution  exc-vxc
    !
    allocate(rvxc(n),stat=i)
       if (i /= 0) then
         write(6,*) 'Error in potential -- allocation of rvxc', n,i
         stop
       endif
    Pot%v0=v0
    CALL exch(Grid,den,rvxc,etxc,eexc,n,v0,v0p)
    !
    !fcount=fcount+1
    !do i=1,n
    !   write(400+fcount,'(1p,7e15.7)') Grid%r(i),den(i),rv(i),rvxc(i)
    !enddo
    rv(1:n)=rv(1:n)+rvxc(1:n)
    !
    !  calculate v0 and v0p
    !
    Pot%v0=Pot%v0+v0
    Pot%v0p=v0p
    !call zeropot(Grid,rv,Pot%v0,Pot%v0p)
    !
    !  Add in Nuclear contribution
    !
    If (.not.finitenucleus) then
       DO i=1,n
          rv(i)=rv(i)-2*nz
       ENDDO
    Else
       RR=nz
       RR=2.9d-5*(RR**0.3333333333333333333333333d0)
       DO i=1,n
          rv(i)=rv(i)-2*nz*derf(Grid%r(i)/RR)
       ENDDO
       Pot%v0=Pot%v0-4*nz/(sqrt(pi)*RR)
    Endif
    WRITE(6,*) 'v0, v0p = ', Pot%v0,Pot%v0p
    DEALLOCATE(rvxc)
  END SUBROUTINE potential

  subroutine ClassicalTurningPoint(Grid,Pot,l,energy,turningpoint)
    Type(GridInfo), INTENT(IN) :: Grid
    Type(PotentialInfo), INTENT(IN) :: Pot
    Integer, INTENT(IN) :: l
    Real(8), INTENT(IN) :: energy
    Integer, INTENT(OUT) :: turningpoint

    integer :: i,n
    Real(8), allocatable :: v(:)

    n=Grid%n
    allocate(v(n), stat=i)
       if (i /= 0) then
          write(6,*) 'Allocation error in ClassicalTurningPoint ', i,n
          stop
       endif

    v=0
    v(2:n)=Pot%rv(2:n)/Grid%r(2:n)+l*(l+1)/(Grid%r(2:n)**2)

    turningpoint=n
    do i=n,2,-1
       if (v(i)<energy) exit
    enddo
    turningpoint=i

    !write(6,*) 'Found turning point at ', turningpoint, Grid%r(turningpoint)

    deallocate(v)

  End subroutine ClassicalTurningPoint


  Subroutine InitPotential(Grid,Pot)
    Type(GridInfo), INTENT(IN) :: Grid
    Type(PotentialInfo), INTENT(OUT) :: Pot

    Integer :: i

    Allocate(Pot%den(Grid%n),Pot%rv(Grid%n), stat=i)
         if ( i /= 0) then
            write(6,*) 'InitPotential allocation error ', i,Grid%n
            stop
         endif
    Pot%den=0.d0;Pot%rv=0.d0
  End subroutine InitPotential

  Subroutine DestroyPotential(Pot)
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
    IF (ASSOCIATED(Pot%rv)) DEALLOCATE(Pot%rv)
    IF (ASSOCIATED(Pot%den)) DEALLOCATE(Pot%den)
  End subroutine DestroyPotential

END MODULE calcpotential
