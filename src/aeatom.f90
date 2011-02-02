MODULE AEatom
  USE GlobalMath
  USE excor
  USE radialsch
  USE radialsr
  USE calcpotential
  USE anderson_realmix
  USE atomdata

  IMPLICIT NONE
  REAL(8), PARAMETER, PRIVATE ::small0=1.d-15,rimix=0.5d0,worst=1.d-5
  REAL(8), PARAMETER, PRIVATE ::linrange=50.d0,linh=0.0025d0,mxgridlin=20001
  REAL(8), PARAMETER, PRIVATE ::logrange=80.d0,logh=0.020d0,mxgridlog=2001
  REAL(8), PARAMETER, PRIVATE ::logder_min=-5.d0,logder_max=4.95d0,logder_pts=200
  INTEGER, PARAMETER, PRIVATE :: niter=1000,mxloop=500
  REAL(8), PARAMETER, PRIVATE :: conv1=4.d13,conv2=3.d13,conv3=2.d13,conv4=1.d13
  INTEGER, PRIVATE :: nps,npp,npd,npf,npg,norbit,n
  REAL(8), PRIVATE :: h, electrons,nz

CONTAINS

  SUBROUTINE iSCFatom(AEGrid,AEPot,AEOrbit,AESCF,ifinput)
    TYPE (GridInfo), INTENT(INOUT) :: AEGrid
    TYPE (PotentialInfo), INTENT(INOUT) :: AEPot
    TYPE (OrbitInfo), INTENT(INOUT) :: AEOrbit
    TYPE (SCFInfo), INTENT(INOUT) :: AESCF
    INTEGER, INTENT(IN),OPTIONAL :: ifinput

    !  program to calculate the self-consistent density functional
    !    atom ground state for atom with atomic number nz
    !    for self-consistent potential rv


    REAL(8) :: xocc,qf,small,zeff,hval,gridmatch,gridrange,logdmin,logdmax
    REAL(8) :: qcal, rescale,nzeff
    INTEGER :: icount,i,j,k,it,start,np,ierr,gridpoints,logdpoints
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir
    INTEGER :: ilin,ilog,inrl,iscl,ipnt,ifin,iend,ilgd
    INTEGER :: igrid,irelat,ilogder
    INTEGER, ALLOCATABLE :: nl(:,:)
    CHARACTER(128) :: exchangecorrelationandgridline,gridkey,relkey
    CHARACTER(132) :: inputline,inputword

    scalarrelativistic=.false.; finitenucleus=.false.
    WRITE(6,*) 'enter atomic symbol and atomic number'
    If (present(ifinput)) then
      READ(5,'(a)') inputline
      WRITE(ifinput,'(a)') TRIM(inputline)
      READ(inputline,*) AEPot%sym,nz
    else
      READ(5,*) AEPot%sym,nz
    endif
    AEPot%nz=nz

    WRITE(6,*) 'exchange-correlation type -- LDA-PW(default) or GGA-PBE or libXC keyword (if linked)'
    WRITE(6,*) 'further optionally (space) "nonrelativistic/scalarrelativistic" keyword'
    WRITE(6,*) 'further optionally (space) "point-nucleus/finite-nucleus" keyword'
    WRITE(6,*) 'optionally (space) "loggrid/lineargrid" keyword if appropriate'
    WRITE(6,*) '    further optionally n (number of grid points)'
    WRITE(6,*) '                       r_max (max. grid radius)'
    WRITE(6,*) '                       r_match (exact value of r(n))'
    WRITE(6,*) 'optionally (space) "logderivrange" keyword'
    WRITE(6,*) '    further optionally emin (minimum energy for log. deriv. plot)'
    WRITE(6,*) '                       emax (maximum energy for log. deriv. plot)'
    WRITE(6,*) '                       ne   (#  of energies for log. deriv. plot)'
    READ(5,'(a)') exchangecorrelationandgridline
    if(present(ifinput)) WRITE(ifinput,'(a)')&
&          TRIM(exchangecorrelationandgridline)
    call Uppercase(exchangecorrelationandgridline)
    exchangecorrelationandgridline=trim(exchangecorrelationandgridline)
    READ(unit=exchangecorrelationandgridline,fmt=*) exctype
    ilin=INDEX(exchangecorrelationandgridline,'LINEARGRID')
    ilog=INDEX(exchangecorrelationandgridline,'LOGGRID')
    inrl=INDEX(exchangecorrelationandgridline,'NONRELATIVISTIC')
    iscl=INDEX(exchangecorrelationandgridline,'SCALARRELATIVISTIC')
    ipnt=INDEX(exchangecorrelationandgridline,'POINT-NUCLEUS')
    ifin=INDEX(exchangecorrelationandgridline,'FINITE-NUCLEUS')
    ilgd=INDEX(exchangecorrelationandgridline,'LOGDERIVRANGE')
    if (iscl>0.and.inrl==0) scalarrelativistic=.true.
    if (ifin>0.and.ipnt==0) finitenucleus=.true.
    gridkey='LINEAR';gridpoints=mxgridlin;gridrange=linrange;gridmatch=linrange
    if (ilog>0.and.ilin==0) then
     gridkey='LOGGRID';gridpoints=mxgridlog;gridrange=logrange;gridmatch=logrange
    endif
    minlogderiv=logder_min;maxlogderiv=logder_max;nlogderiv=logder_pts

    igrid=max(ilin,ilog)
    irelat=max(inrl,iscl)
    ilogder=ilgd

    if (igrid>0) then
     iend=128
     if (irelat >igrid.and.irelat-1 <iend) iend=irelat -1
     if (ilogder>igrid.and.ilogder-1<iend) iend=ilogder-1
     inputline=""
     if (ilog>0.and.iend>igrid+7) inputline=trim(exchangecorrelationandgridline(igrid+7:iend))
     if (ilin>0.and.iend>igrid+10)inputline=trim(exchangecorrelationandgridline(igrid+10:iend))
     if (inputline/="") then
      call extractword(1,inputline,inputword);inputword=trim(inputword)
      if (inputword/="") then
       read(inputword,*) gridpoints
       call extractword(2,inputline,inputword);inputword=trim(inputword)
       if (inputword/="") then
        read(inputword,*) gridrange
        gridmatch=gridrange
        call extractword(3,inputline,inputword);inputword=trim(inputword)
        if (inputword/="") read(inputword,*) gridmatch
       endif
      endif
     endif
    endif

    if (ilogder>0) then
     iend=128
     if (igrid >ilogder.and.igrid-1 <iend) iend=igrid -1
     if (irelat>ilogder.and.irelat-1<iend) iend=irelat-1
     inputline=""
     if (iend>ilogder+13) inputline=trim(exchangecorrelationandgridline(ilogder+13:iend))
     if (inputline/="") then
      call extractword(1,inputline,inputword);inputword=trim(inputword)
      if (inputword/="") then
       read(inputword,*) minlogderiv
       call extractword(2,inputline,inputword);inputword=trim(inputword)
       if (inputword/="") then
        read(inputword,*) maxlogderiv
        call extractword(3,inputline,inputword);inputword=trim(inputword)
        if (inputword/="") read(inputword,*) nlogderiv
       endif
      endif
     endif
    endif

    WRITE(6,*) 'Calculation for atomic number = ',nz
    WRITE(6,*) 'enter maximum principal quantum numbers for s,p,d,f,g'
    if(present(ifinput)) then
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       READ(inputline,*) nps,npp,npd,npf,npg
    else
       READ(5,*) nps,npp,npd,npf,npg
    endif
    IF(nps<0)nps=0
    IF(npp<0)npp=0
    IF(npd<0)npd=0
    IF(npf<0)npf=0
    IF(npg<0)npg=0
    WRITE(6,'(5i4)') nps,npp,npd,npf,npg
    i=MAX(nps,npp,npd,npf,npg)
    j=nps+npp+npd+npf+npg
    CALL InitOrbit1(AEOrbit,j)
    ALLOCATE(nl(i,j));nl=0
    icount=0
    IF (nps.GT.0) THEN
       DO is=1,nps
          icount=icount+1
          nl(is,1)=icount
          AEOrbit%occ(icount)=2.d0
          AEOrbit%np(icount)=is
          AEOrbit%l(icount)=0
       ENDDO
    ENDIF
    IF (npp.GT.1) THEN
       DO ip=2,npp
          icount=icount+1
          nl(ip,2)=icount
          AEOrbit%occ(icount)=6.d0
          AEOrbit%np(icount)=ip
          AEOrbit%l(icount)=1
       ENDDO
    ENDIF
    IF (npd.GT.2) THEN
       DO id=3,npd
          icount=icount+1
          nl(id,3)=icount
          AEOrbit%occ(icount)=10.d0
          AEOrbit%np(icount)=id
          AEOrbit%l(icount)=2
       ENDDO
    ENDIF
    IF (npf.GT.3) THEN
       DO jf=4,npf
          icount=icount+1
          nl(jf,4)=icount
          AEOrbit%occ(icount)=14.d0
          AEOrbit%np(icount)=jf
          AEOrbit%l(icount)=3
       ENDDO
    ENDIF
    IF(npg.GT.4) THEN
       DO ig=5,npg
          icount=icount+1
          nl(ig,5)=icount
          AEOrbit%occ(icount)=18.d0
          AEOrbit%np(icount)=ig
          AEOrbit%l(icount)=4
       ENDDO
    ENDIF
    norbit=icount
    AEOrbit%nps=nps
    AEOrbit%npp=npp
    AEOrbit%npd=npd
    AEOrbit%npf=npf
    AEOrbit%npg=npg
    AEOrbit%norbit=norbit
    WRITE(6,*) norbit, ' orbitals will be calculated'
    !
    WRITE(6,*)' Below are listed the default occupations '
    WRITE(6,"(' n  l     occupancy')")
    DO io=1,norbit
       WRITE(6,'(i2,1x,i2,4x,1pe15.7)') &
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
    ENDDO
    !
    WRITE(6,*)' enter np l occ for all occupations for all revisions'
    WRITE(6,*)'  enter 0 0 0. to end'

    DO
      if(present(ifinput)) then
         READ(5,'(a)') inputline
         WRITE(ifinput,'(a)') TRIM(inputline)
         READ(inputline,*) ip,l,xocc
      else
         READ(5,*) ip,l,xocc
      endif
       IF (ip.LE.0) EXIT
       nfix=nl(ip,l+1)
       IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
          WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
&              ip,l,xocc,nfix,norbit
          STOP
       ENDIF
       AEOrbit%occ(nfix)=xocc
    ENDDO

    !
    WRITE(6,*) ' Corrected occupations are: '
    WRITE(6,"(' n  l     occupancy')")
    electrons=0.d0
    DO io=1,norbit
       WRITE(6,'(i2,1x,i2,4x,1pe15.7)')&
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
       electrons=electrons+AEOrbit%occ(io)
    ENDDO
    AEPot%q=electrons
    qf=nz-electrons
    WRITE(6,*)
    WRITE(6,*) 'nuclear charge    = ' , nz
    WRITE(6,*) 'electronic charge = ', electrons
    WRITE(6,*) 'net charge        = ', qf
    !
    if (trim(gridkey)=='LOGGRID') then
       hval=logh;call findh(nz,gridmatch,gridpoints,hval)
       Call initgrid(AEGrid,hval,gridrange,hval/nz)
    else
       hval=gridmatch/dble(gridpoints-1)
       Call initgrid(AEGrid,hval,gridrange)
    endif

    CALL InitOrbit2(AEOrbit,AEGrid%n,j)
    call InitPotential(AEGrid,AEPot)
    if (scalarrelativistic) Call Allocate_Scalar_Relativistic(AEGrid)

    !
    !  calculate initial charge density from hydrogen-like functions
    !    also initial energies
    !
    AEPot%den(1:n)=0.d0
    zeff=nz
    n=AEGrid%n
    DO io=1,norbit
       np=AEOrbit%np(io)
       l=AEOrbit%l(io)
       xocc=AEOrbit%occ(io)
       zeff=zeff-0.5d0*xocc
       zeff=MAX(zeff,1.d0)
       nzeff=zeff+0.1d0    ! no longer used
       AEOrbit%eig(io)=-(zeff/(np))**2
       WRITE(6,*) io,np,l,xocc,AEOrbit%eig(io)
       DO ir=1,n
          AEOrbit%wfn(ir,io)=hwfn(zeff,np,l,AEGrid%r(ir))
          AEPot%den(ir)=AEPot%den(ir)+xocc*(AEOrbit%wfn(ir,io)**2)
       ENDDO
       zeff=zeff-0.5d0*xocc
    ENDDO
    !
    !  check charge
    !
    qcal=integrator(AEGrid,AEPOT%den)
    qf=qcal
    WRITE(6,*) 'qcal electrons = ',qcal, electrons
    !  rescale density
    rescale=electrons/qcal
    AEPot%den(1:n)=AEPot%den(1:n)*rescale
    !
    !  choose form of exchange-correlation potential
    CALL initexch
    !
    DEALLOCATE(nl)
    CALL SCFloop(AEGrid,AEPot,AEOrbit,AESCF)
  END SUBROUTINE iSCFatom

  SUBROUTINE SCFloop(AEGrid,AEPot,AEOrbit,AESCF)
    TYPE (GridInfo), INTENT(INOUT) :: AEGrid
    TYPE (PotentialInfo), INTENT(INOUT) :: AEPot
    TYPE (OrbitInfo), INTENT(INOUT) :: AEOrbit
    TYPE (SCFInfo), INTENT(INOUT) :: AESCF

    !  program to perform self-consistency loop

    TYPE (Anderson_Context) , POINTER :: AErho

    REAL(8) :: rmix,xocc,qf,small,zeff,delta,v1,v2,v3,v4,x,y,nzeff
    REAL(8) :: qcal, rescale,cnvrg,emin,ecoul,eexc,etxc,eone,etot,ekin
    INTEGER :: icount,i,j,k,it,start,np,ierr,nroot
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,loop,jierr
    REAL(8), ALLOCATABLE :: denout(:)
    INTEGER :: fcount=0

    v1=conv1;v2=conv2;v3=conv3;v4=conv4
    n=AEGrid%n; h=AEGrid%h
    rmix=rimix
    cnvrg=worst
    small=small0

    WRITE(6,*) 'Density convergence parameter set to at least',cnvrg
    !
    ALLOCATE(denout(n),STAT=k)
    IF (k /= 0) THEN
       WRITE(6,*) 'Error in denout allocation ', n,k
       STOP
    ENDIF

    !    start iteration loop
    !
    DO   loop=1,mxloop
       !  calculate  potential * r == rv for given density
       !
       CALL potential(AEGrid,AEPot,ecoul,etxc,eexc)
       denout=0
       denout(2:n)=-2*nz*AEPot%den(2:n)/AEGrid%r(2:n)
       AESCF%estatic=integrator(AEGrid,denout)+ecoul
       !
       !  solve for bound states of Schroedinger equation
       !
       icount=0
       qf=nz-electrons
       jierr=0
       it=0
       !  s states :
       IF (nps.GT.0) THEN
          it=it+1
          emin=-nz*nz
          l=0
          nroot=nps
          start=1
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  p states :
       IF (npp.GT.1) THEN
          it=it+1
          emin=-nz*nz/4.d0
          l=1
          nroot=npp-1
          start=start+nps
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  d states :
       IF (npd.GT.2) THEN
          it=it+1
          emin=-nz*nz/9.d0
          l=2
          nroot=npd-2
          start=start+npp-1
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  f states :
       IF (npf.GT.3) THEN
          it=it+1
          emin=-nz*nz/16.d0
          l=3
          nroot=npf-3
          start=start+npd-2
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  g states :
       IF (npg.GT.4) THEN
          it=it+1
          emin=-nz*nz/25.d0
          l=4
          nroot=npg-4
          start=start+npf-3
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !
       !  calculate new density
       !
       denout(1:n)=0.d0
       WRITE(6,*) '  results for loop = ',loop
       WRITE(6,*) ' n  l     occupancy       energy'
       DO io=1,norbit
          WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&              AEOrbit%np(io),AEOrbit%l(io),&
&              AEOrbit%occ(io),AEOrbit%eig(io)
          IF (AEOrbit%occ(io).GT.small) THEN
             DO i=1,n
                denout(i)=denout(i)+AEOrbit%occ(io)*(AEOrbit%wfn(i,io)**2)
             ENDDO
          ENDIF
       ENDDO
       qcal=integrator(AEGrid,denout)
       WRITE(6,*) 'qcal electrons = ',qcal, electrons
       !  rescale density
       rescale=electrons/qcal
       denout(1:n)=denout(1:n)*rescale
        !fcount=fcount+1
        !do i=1,n
        ! write(fcount+100,'(1p,2e15.7)') AEGrid%r(i),denout(i)
        !enddo
       !
       denout=denout-AEPot%den
       delta=SUM(ABS(denout))
       CALL shift4(v1,v2,v3,v4,delta)
       WRITE(6,*) 'density iter',loop,delta
       IF (loop.EQ.1) CALL InitAnderson(AErho,6,5,n,rmix,1.d5)
       CALL Anderson_Mix(AErho,AEPot%den(1:n),denout(1:n))
       !  correct and recale density
       DO i=1,n
          IF (AEPot%den(i).LT.0.d0) AEPot%den(i)=0.d0
       ENDDO
       qcal=integrator(AEGrid,AEPot%den)
       WRITE(6,*) 'qcal electrons = ',qcal, electrons
       !  rescale density
       rescale=electrons/qcal
       AEPot%den(1:n)=AEPot%den(1:n)*rescale
       !
       WRITE(6,*) '  results for loop ,delta = ',loop,delta
       WRITE(6,*) ' n  l     occupancy       energy'
       eone=0.d0
       DO io=1,norbit
          WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&              AEOrbit%np(io),AEOrbit%l(io),&
&              AEOrbit%occ(io),AEOrbit%eig(io)
          eone=eone+AEOrbit%occ(io)*AEOrbit%eig(io)
       ENDDO
       WRITE(6,*)
       WRITE(6,*) ' Total energies'
       WRITE(6,*) '    One-electron contribution:  ',eone
       WRITE(6,*) '    Coulomb contribution     :  ',ecoul
       WRITE(6,*) '    Exch-correl contribution :  ',eexc
       etot=eone-ecoul+etxc
       WRITE(6,*) '    Total                    :  ',etot
       !
       IF (loop>=4) THEN
          IF (.NOT.(v4.LE.v3.AND.v3.LE.v2 &
&              .AND.v2.LE.v1).AND.v4.LE.cnvrg) THEN
             !
             !  converged result
             !
             WRITE(6,*) '  dfatom converged in',loop,' iterations'
             AESCF%iter=loop
             WRITE(6,*) '     for nz = ',nz
             WRITE(6,*) '    delta(density)  = ', delta
             AESCF%delta=delta
             WRITE(6,*) '  results for loop = ',loop
             WRITE(6,*) ' n  l     occupancy       energy'
             eone=0.d0
             ekin=0.0
             DO io=1,norbit
                WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&                    AEOrbit%np(io),AEOrbit%l(io),&
&                    AEOrbit%occ(io),AEOrbit%eig(io)
                eone=eone+AEOrbit%occ(io)*AEOrbit%eig(io)
                !CALL kinetic(AEGrid,AEOrbit%l(io),&
                !&  AEOrbit%wfn(:,io),x)
                !CALL altkinetic(AEGrid,AEOrbit%wfn(:,io),AEOrbit%eig(io)&
                !&       ,AEPot%rv,y)
                !WRITE(6,*) 'Kinetic compare',io,x,y
                CALL altkinetic(AEGrid,AEOrbit%wfn(:,io),AEOrbit%eig(io)&
&                      ,AEPot%rv,x)
                ekin=ekin+AEOrbit%occ(io)*x
             ENDDO
             WRITE(6,*)
             WRITE(6,*) ' Total energies'
             WRITE(6,*) '    One-electron contribution:  ',eone
             AESCF%eone=eone
             WRITE(6,*) '    Kinetic energy contribution:',ekin
             AESCF%ekin=ekin
             WRITE(6,*) '    Coulomb contribution     :  ',ecoul
             AESCF%ecoul=ecoul
             WRITE(6,*) '    Electrostatic contribution: ',AESCF%estatic
             WRITE(6,*) '    Exch-correl contribution :  ',eexc
             AESCF%eexc=eexc
             etot=eone-ecoul+etxc
             AESCF%etot=etot
             WRITE(6,*) '    Total                    :  ',etot
             WRITE(6,*) '    Total (alt form)         :  ',&
&                   AESCF%ekin+AESCF%estatic+AESCF%eexc
             CALL FreeAnderson(AErho)
             DEALLOCATE(denout)
             RETURN
          ENDIF
       ENDIF
    ENDDO   ! mxloop

    WRITE(6,*)'calculation terminating without density convergence'
    WRITE(6,*) 'delta, cnvrg =', delta,cnvrg
    WRITE(6,*) 'loop,mxloop = ',loop,mxloop
    CALL FreeAnderson(AErho)
    DEALLOCATE(denout)
    STOP
  END SUBROUTINE SCFloop

  SUBROUTINE cSCFatom(AEGrid,AEPot,AEOrbit,AESCF)
    TYPE (GridInfo), INTENT(INOUT) :: AEGrid
    TYPE (PotentialInfo), INTENT(INOUT) :: AEPot
    TYPE (OrbitInfo), INTENT(INOUT) :: AEOrbit
    TYPE (SCFInfo), INTENT(INOUT) :: AESCF

    !  program to calculate the self-consistent density functional
    !    atom ground state for atom with atomic number nz
    !    for self-consistent potential rv
    !    version for changing configuration after initial
    !    SCF run

    REAL(8) :: xocc,qf,small,zeff
    REAL(8) :: qcal, rescale
    INTEGER :: icount,i,j,k,it,start,np,ierr
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,nzeff
    INTEGER, ALLOCATABLE :: nl(:,:)

    WRITE(6,*) norbit, ' orbitals will be calculated'
    !
    WRITE(6,*)' Below are listed the current occupations '
    WRITE(6,"(' n  l     occupancy')")
    DO io=1,norbit
       WRITE(6,'(i2,1x,i2,4x,1pe15.7)') &
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
    ENDDO
    !
    WRITE(6,*)' enter np l occ for all occupations for all revisions'
    WRITE(6,*)'  enter 0 0 0. to end'

    DO
       READ(5,*) ip,l,xocc
       IF (ip.LE.0) EXIT
       nfix=-100
       DO io=1,norbit
          IF (ip==AEOrbit%np(io).AND.l==AEOrbit%l(io)) THEN
             nfix=io
             EXIT
          ENDIF
       ENDDO
       IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
          WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
&              ip,l,xocc,nfix,norbit
          STOP
       ENDIF
       AEOrbit%occ(nfix)=xocc
    ENDDO

    !
    WRITE(6,*) ' Corrected occupations are: '
    WRITE(6,"(' n  l     occupancy')")
    electrons=0.d0
    DO io=1,norbit
       WRITE(6,'(i2,1x,i2,4x,1pe15.7)')&
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
       electrons=electrons+AEOrbit%occ(io)
    ENDDO
    AEPot%q=electrons
    qf=nz-electrons
    WRITE(6,*)
    WRITE(6,*) 'nuclear charge    = ' , nz
    WRITE(6,*) 'electronic charge = ', electrons
    WRITE(6,*) 'net charge        = ', qf
    !
    !
    !  calculate initial charge density from stored wavefunctions
    !    also initial energies
    !
    AEPot%den(1:n)=0.d0
    DO io=1,norbit
       xocc=AEOrbit%occ(io)
       DO ir=1,n
          AEPot%den(ir)=AEPot%den(ir)+xocc*(AEOrbit%wfn(ir,io)**2)
       ENDDO
    ENDDO
    !
    !  check charge
    !
    qcal=integrator(AEGrid,AEPOT%den)
    qf=qcal
    WRITE(6,*) 'qcal electrons = ',qcal, electrons
    !  rescale density
    rescale=electrons/qcal
    AEPot%den(1:n)=AEPot%den(1:n)*rescale
    !
    CALL SCFloop(AEGrid,AEPot,AEOrbit,AESCF)
  END SUBROUTINE cSCFatom

 SUBROUTINE ChooseValence(Grid,Orbit,FC,ifinput)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Orbitinfo), INTENT(IN) :: Orbit
    TYPE(FCinfo), INTENT(INOUT) :: FC
    INTEGER, INTENT(IN), OPTIONAL :: ifinput

    CHARACTER(1) :: answer
    INTEGER :: io,n,ok,norbit
    CHARACTER(132) :: inputline

    norbit=Orbit%norbit;n=Grid%n
    call InitFC(FC,norbit,n)

    WRITE(6,*) 'for each state enter c for core or v for valence'
    FC%zcore=0; FC%zvale=0
    DO io=1,FC%norbit
       WRITE(6,'(3i5,1p,2e15.7)') io,Orbit%np(io),Orbit%l(io),&
&           Orbit%occ(io),Orbit%eig(io)
       DO
         if (present(ifinput)) then
          READ(5,'(a)') inputline
          WRITE(ifinput,'(a)') TRIM(inputline)
          READ(inputline,*) answer
         else
          READ(5,*) answer
         endif
          IF (answer.NE.'c'.AND.answer.NE.'C'.AND.answer.NE.'v'&
&              .AND.answer.NE.'V')  THEN
             WRITE(6,*) 'Please input c or v'
          ELSE
             EXIT
          ENDIF
       ENDDO
       IF (answer.EQ.'c'.OR.answer.EQ.'C') FC%iscore(io)=.true.
       IF (answer.EQ.'v'.OR.answer.EQ.'V') FC%iscore(io)=.false.
       IF (FC%iscore(io)) then
           FC%zcore=FC%zcore+Orbit%occ(io)
           FC%coreden=FC%coreden+Orbit%occ(io)*(Orbit%wfn(:,io))**2
       ENDIF
       IF (.NOT.FC%iscore(io)) then
           FC%zvale=FC%zvale+Orbit%occ(io)
           FC%valeden=FC%valeden+Orbit%occ(io)*(Orbit%wfn(:,io))**2
       ENDIF
    ENDDO
  END SUBROUTINE ChooseValence

  !***********************************************************************
  ! Subroutine FCenergy --
  !   Calculates valence energy of atom given fixed FC%coreden
  !   reference: U. von Barth and C. D. Gelatt, Phys. Rev. B 21, 2222(1980)
  !***********************************************************************
  SUBROUTINE FCenergy(Grid,Pot,Orbit,SCF,FC)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Orbitinfo), INTENT(IN) :: Orbit
    TYPE(SCFInfo), INTENT(INOUT) :: SCF
    TYPE(FCinfo), INTENT(INOUT) :: FC

    INTEGER :: n,i,io,ok
    REAL(8) :: h,nz
    REAL(8) :: qtot,qchk,qcore,qvale,y
    REAL(8) :: etotal,ec,ev,ecv,tc,tv,ekin,cext,vext,tmp
    REAL(8) :: ccoul,vcoul,cvcoul,chkcvcoul1,chkcvcoul2
    REAL(8) :: cexc,vexc,texc
    REAL(8), ALLOCATABLE :: den(:),dum1(:),dum2(:),dum3(:)

    n=Grid%n
    h=Grid%h
    nz=Pot%nz
    ALLOCATE(den(n),dum1(n),dum2(n),dum3(n), stat=ok)
    IF(ok/=0) THEN
       WRITE(6,*) ' Error in FCenergy allocation ' , n
       STOP
    ENDIF
    FC%valeden(1:n)=0.d0
    DO io=1,FC%norbit
       IF (.NOT.FC%iscore(io)) THEN
          DO i=1,n
             FC%valeden(i)=FC%valeden(i)+Orbit%occ(io)*(Orbit%wfn(i,io)**2)
          ENDDO
       ENDIF
    ENDDO
    DO i=1,n
       den(i)=FC%coreden(i)+FC%valeden(i)
    ENDDO
    !
    etotal=0
    ec=0
    ev=0
    ecv=0
    ! kinetic energy
    tc=0
    tv=0
    DO io=1,FC%norbit
       IF (FC%iscore(io)) THEN
          !CALL kinetic(Grid,Orbit%l(io),Orbit%wfn(:,io),ekin)
          !CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io)&
          !&             ,Pot%rv,y)
          !      WRITE(6,*) 'Kinetic compare',io,ekin,y
          CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io)&
&                      ,Pot%rv,ekin)
          tc=tc+Orbit%occ(io)*ekin
       ELSE IF (.NOT.FC%iscore(io)) THEN
          !CALL kinetic(Grid,Orbit%l(io),Orbit%wfn(:,io),ekin)
          !CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io)&
          !&             ,Pot%rv,y)
          !      WRITE(6,*) 'Kinetic compare',io,ekin,y
          CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io)&
&                      ,Pot%rv,ekin)
          tv=tv+Orbit%occ(io)*ekin
       ELSE
          WRITE(6,*) 'Error in Kinetic energy loop of FCenergy',io
          STOP
       ENDIF
    ENDDO
    WRITE(6,*) 'core and valence kinetic energies', tc,tv
    FC%corekin=tc
    ! external potential interaction
    dum1=0;dum2=0
    DO i=2,n
       dum1(i)=FC%coreden(i)/Grid%r(i)
       dum2(i)=FC%valeden(i)/Grid%r(i)
    ENDDO
    cext=-2*nz*integrator(Grid,dum1)
    vext=-2*nz*integrator(Grid,dum2)
    WRITE(6,*) 'core and valence external potential energies',cext,vext
    !
    qcore=integrator(Grid,FC%coreden)
    qvale=integrator(Grid,FC%valeden)
    qtot=qcore+qvale
    qchk=integrator(Grid,den)
    WRITE(6,*) 'qcore,qvale,qtot,qchk=', qcore,qvale,qtot,qchk
    WRITE(6,*) 'zcore,zvale = ', FC%zcore,FC%zvale
    CALL poisson(Grid,qcore,FC%coreden,dum3,ccoul)
    chkcvcoul1=overlap(Grid,dum3,dum2)

    CALL poisson(Grid,qvale,FC%valeden,dum3,vcoul)
    chkcvcoul2=overlap(Grid,dum3,dum1)
    CALL poisson(Grid,qtot,den,dum3,cvcoul)
    CALL exch(Grid,FC%coreden,dum3,tmp,cexc)
    CALL exch(Grid,FC%valeden,dum3,tmp,vexc)
    CALL exch(Grid,den,dum3,tmp,texc)
    cvcoul=cvcoul-vcoul-ccoul
    WRITE(6,*) 'core , valence and interaction coulomb energies'
    WRITE(6,*) ccoul,vcoul,cvcoul,chkcvcoul1,chkcvcoul2
    WRITE(6,*) 'core , valence and interaction exc energies'
    WRITE(6,*) cexc,vexc,texc
    ec=tc+ccoul+cext+cexc
    ev=tv+vcoul+vext+vexc
    ecv=cvcoul+texc-cexc-vexc
    etotal=ec+ev+ecv
    SCF%etot=etotal
    !FC%evale=tv+vcoul+vext+cvcoul+texc-cexc
    FC%evale=tv+vcoul+vext+cvcoul+texc
    WRITE(6,*) 'ec = ',ec
    WRITE(6,*) 'ev = ',ev
    WRITE(6,*) 'ecv = ', ecv
    WRITE(6,*) 'etotal = ',etotal
    WRITE(6,*) 'evale = ',FC%evale

    DEALLOCATE(den,dum1,dum2,dum3)
  END SUBROUTINE FCenergy

  SUBROUTINE FCselfenergy(Grid,Orbit,FC,selfenergy)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Orbitinfo), INTENT(IN) :: Orbit
    TYPE(FCinfo), INTENT(IN) :: FC
    REAL(8), INTENT(OUT) :: selfenergy

    INTEGER :: n,i,io,norbit
    REAL(8) :: h,x,y
    REAL(8), allocatable :: dum1(:),dum2(:)

    n=Grid%n
    h=Grid%h
    ALLOCATE(dum1(n),dum2(n),stat=i)
       if(i/=0) then
          write(6,*) 'Allocation error in FCselfenergy', i,n
          stop
       endif

    selfenergy=0
    DO io=1,FC%norbit
       IF (.NOT.FC%iscore(io).and.Orbit%occ(io)>1.d-6) THEN
          dum1=0
          DO i=1,n
             dum1(i)=dum1(i)+(Orbit%wfn(i,io)**2)
          ENDDO
          x=integrator(Grid,dum1)
          dum1=dum1/x
          x=1.d0
          call poisson(Grid,x,dum1,dum2,y)
          write(6,*) 'Self energy contribution ', io,Orbit%occ(io),y
          selfenergy=selfenergy+Orbit%occ(io)*y
       ENDIF
    ENDDO

    write(6,*) 'Total valence self energy contribution ', selfenergy
    deallocate(dum1,dum2)
  END SUBROUTINE FCselfenergy

  SUBROUTINE FCSCFatom(AEGrid,AEPot,AEOrbit,FCOrbit,AESCF,FC)
    TYPE (GridInfo), INTENT(INOUT) :: AEGrid
    TYPE (PotentialInfo), INTENT(INOUT) :: AEPot
    TYPE (OrbitInfo), INTENT(INOUT) :: AEOrbit
    TYPE (OrbitInfo), INTENT(INOUT) :: FCOrbit
    TYPE (SCFInfo), INTENT(INOUT) :: AESCF
    TYPE (FCinfo), INTENT(INOUT) :: FC

    !  program to calculate the self-consistent density functional
    !    atom ground state for atom with atomic number nz
    !    for self-consistent potential rv
    !    version for changing configuration after initial
    !    SCF run and for fixing "frozen" core

    REAL(8) :: xocc,qf,small,zeff
    REAL(8) :: qcal, rescale,nzeff
    INTEGER :: icount,i,j,k,it,start,np,ierr,n
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir
    INTEGER, save :: init=0

    if (init==0) then
       call ChooseValence(AEGrid,AEOrbit,FC)
       init=1
    endif

    write(6,*) 'completed ChooseValence -- norbit = ', FC%norbit,AEOrbit%norbit
    norbit=AEOrbit%norbit ; n=AEGrid%n
    WRITE(6,*) 'Frozen core calculation '
    ! setup FCOrbit data structure
    FCOrbit%norbit=AEOrbit%norbit
    FCOrbit%nps=AEOrbit%nps
    FCOrbit%npp=AEOrbit%npp
    FCOrbit%npd=AEOrbit%npd
    FCOrbit%npf=AEOrbit%npf
    FCOrbit%npg=AEOrbit%npg
    Allocate(FCOrbit%np(norbit),FCOrbit%l(norbit),FCOrbit%eig(norbit),&
&       FCOrbit%occ(norbit),FCOrbit%wfn(n,norbit),stat=i)
       if (i /= 0) then
          write(6,*) ' Allocation error in FCSCFatom ',i,norbit,n
          stop
       endif
   FCOrbit%np(1:norbit)=AEOrbit%np(1:norbit)
   FCOrbit%l(1:norbit)=AEOrbit%l(1:norbit)
   FCOrbit%eig(1:norbit)=AEOrbit%eig(1:norbit)
   FCOrbit%occ(1:norbit)=AEOrbit%occ(1:norbit)
   FCOrbit%wfn(:,1:norbit)=AEOrbit%wfn(:,1:norbit)
    !
    WRITE(6,*)' Below are listed the core states'
    WRITE(6,"(' n  l     occupancy')")
    DO io=1,norbit
       if (FC%iscore(io)) WRITE(6,'(i2,1x,i2,4x,1pe15.7)') &
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
    ENDDO

    WRITE(6,*)' Below are listed the current valence occupations '
    WRITE(6,"(' n  l     occupancy')")
    DO io=1,norbit
       if (.NOT.FC%iscore(io)) WRITE(6,'(i2,1x,i2,4x,1pe15.7)') &
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
    ENDDO
    !
    WRITE(6,*)' enter np l occ for all occupations for all revisions'
    WRITE(6,*)'  enter 0 0 0. to end'

    DO
       READ(5,*) ip,l,xocc
       IF (ip.LE.0) EXIT
       nfix=-100
       DO io=1,norbit
          IF (ip==AEOrbit%np(io).AND.l==AEOrbit%l(io)&
&             .AND.(.NOT.FC%iscore(io))) THEN
             nfix=io
             EXIT
          ENDIF
       ENDDO
       IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
          WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
&              ip,l,xocc,nfix,norbit
          STOP
       ENDIF
       AEOrbit%occ(nfix)=xocc; FCOrbit%occ(nfix)=xocc
    ENDDO

    !
    WRITE(6,*) ' Corrected occupations are: '
    WRITE(6,"(' n  l     occupancy')")
    electrons=0.d0
    DO io=1,norbit
       if (.NOT.FC%iscore(io))then
         WRITE(6,'(i2,1x,i2,4x,1pe15.7)')&
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
         electrons=electrons+AEOrbit%occ(io)
       endif
    ENDDO
    AEPot%q=electrons+FC%zcore
    FC%zvale=electrons
    qf=nz-electrons-FC%zcore
    WRITE(6,*)
    WRITE(6,*) 'nuclear charge    = ' , nz
    WRITE(6,*) 'core charge       = ' , FC%zcore
    WRITE(6,*) 'electronic charge = ', electrons
    WRITE(6,*) 'net charge        = ', qf
    !
    small=small0

    !
    !  calculate initial charge density from stored wavefunctions
    !    also initial energies
    !
    FC%valeden=0
    DO io=1,norbit
       if (.NOT.FC%iscore(io)) then
          xocc=AEOrbit%occ(io)
          DO ir=1,n
             FC%valeden(ir)=FC%valeden(ir)+xocc*(AEOrbit%wfn(ir,io)**2)
          ENDDO
       endif
    ENDDO
    !
    !  check charge
    !
    qcal=integrator(AEGrid,FC%valeden)
    qf=qcal
    WRITE(6,*) 'qcal electrons = ',qcal, electrons
    !  rescale density
    rescale=electrons/qcal
    FC%valeden(1:n)=FC%valeden(1:n)*rescale
    !
    CALL FCSCFloop(AEGrid,AEPot,AEOrbit,AESCF,FC)
    ! reset FCOrbit datastructure
    Do io=1,norbit
       if (.NOT.FC%iscore(io)) then
          FCOrbit%occ(io)=AEOrbit%occ(io)
          FCOrbit%eig(io)=AEOrbit%eig(io)
          FCOrbit%wfn(:,io)=AEOrbit%wfn(:,io)
       endif
    ENDDO

    CALL FCenergy(AEGrid,AEPot,FCOrbit,AESCF,FC)
  END SUBROUTINE FCSCFatom

  SUBROUTINE FCSCFloop(AEGrid,AEPot,AEOrbit,AESCF,FC)
    TYPE (GridInfo), INTENT(INOUT) :: AEGrid
    TYPE (PotentialInfo), INTENT(INOUT) :: AEPot
    TYPE (OrbitInfo), INTENT(INOUT) :: AEOrbit
    TYPE (SCFInfo), INTENT(INOUT) :: AESCF
    TYPE (FCinfo), INTENT(INOUT) :: FC

    !  program to perform self-consistency loop -- frozencore case

    TYPE (Anderson_Context) , POINTER :: FCrho

    REAL(8) :: rmix,xocc,qf,small,zeff,delta,v1,v2,v3,v4,nzeff
    REAL(8) :: qcal, rescale,cnvrg,emin,ecoul,eexc,etxc,eone,etot
    INTEGER :: icount,i,j,k,it,start,np,ierr,nroot
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,loop,jierr
    REAL(8), ALLOCATABLE :: denout(:)

    v1=conv1;v2=conv2;v3=conv3;v4=conv4
    n=AEGrid%n; h=AEGrid%h
    rmix=rimix
    cnvrg=worst
    IF (small.LT.0.d0.OR.small.GT.1.d0) small=small0

    WRITE(6,*) 'Density convergence parameter set to at least',cnvrg
    !
    ALLOCATE(denout(n),STAT=k)
    IF (k /= 0) THEN
       WRITE(6,*) 'Error in denout allocation ', n,k
       STOP
    ENDIF

    !    start iteration loop
    !
    DO   loop=1,mxloop
       !  calculate  potential * r == rv for given density
       !
       AEPot%den=FC%coreden+FC%valeden
       AEPot%q=FC%zcore+FC%zvale
       CALL potential(AEGrid,AEPot,ecoul,etxc,eexc)
       !
       !  solve for bound states of Schroedinger equation
       !    all states are actually calculated, but only
       !      valence states are used
       !
       icount=0
       jierr=0
       it=0
       !  s states :
       IF (nps.GT.0) THEN
          it=it+1
          emin=-nz*nz
          l=0
          nroot=nps
          start=1
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  p states :
       IF (npp.GT.1) THEN
          it=it+1
          emin=-nz*nz/4.d0
          l=1
          nroot=npp-1
          start=start+nps
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  d states :
       IF (npd.GT.2) THEN
          it=it+1
          emin=-nz*nz/9.d0
          l=2
          nroot=npd-2
          start=start+npp-1
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  f states :
       IF (npf.GT.3) THEN
          it=it+1
          emin=-nz*nz/16.d0
          l=3
          nroot=npf-3
          start=start+npd-2
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !  g states :
       IF (npg.GT.4) THEN
          it=it+1
          emin=-nz*nz/25.d0
          l=4
          nroot=npg-4
          start=start+npf-3
          if (scalarrelativistic) then
            CALL boundsr(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          else
            CALL boundsch(AEGrid,AEPot,AEOrbit,l,start,nroot,emin,ierr)
          endif
       ENDIF
       !
       !  calculate new density
       !
       denout(1:n)=0.d0
       WRITE(6,*) '  results for loop = ',loop
       WRITE(6,*) ' n  l     occupancy       energy'
       DO io=1,norbit
          if (.NOT.FC%iscore(io)) then
               WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&             AEOrbit%np(io),AEOrbit%l(io),&
&              AEOrbit%occ(io),AEOrbit%eig(io)
             IF (AEOrbit%occ(io).GT.small) THEN
                DO i=1,n
                   denout(i)=denout(i)+AEOrbit%occ(io)*(AEOrbit%wfn(i,io)**2)
                ENDDO
             ENDIF
          endif
       ENDDO
       qcal=integrator(AEGrid,denout)
       WRITE(6,*) 'qcal valence electrons = ',qcal, FC%zvale
       !  rescale density
       rescale=FC%zvale/qcal
       denout(1:n)=denout(1:n)*rescale
       !
       denout=denout-FC%valeden
       delta=SUM(ABS(denout))
       CALL shift4(v1,v2,v3,v4,delta)
       WRITE(6,*) 'density iter',loop,delta
       IF (loop.EQ.1) CALL InitAnderson(FCrho,6,5,n,rmix,1.d5)
       CALL Anderson_Mix(FCrho,FC%valeden(1:n),denout(1:n))
       !  correct and recale density
       DO i=1,n
          IF (FC%valeden(i).LT.0.d0) FC%valeden(i)=0.d0
       ENDDO
       qcal=integrator(AEGrid,FC%valeden)
       WRITE(6,*) 'qcal valence electrons = ',qcal, FC%zvale
       !  rescale density
       rescale=FC%zvale/qcal
       FC%valeden(1:n)=FC%valeden(1:n)*rescale
       !
       WRITE(6,*) '  results for loop ,delta = ',loop,delta
       WRITE(6,*) ' n  l     occupancy       energy'
       DO io=1,norbit
          if (.NOT.FC%iscore(io)) WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&              AEOrbit%np(io),AEOrbit%l(io),&
&              AEOrbit%occ(io),AEOrbit%eig(io)
       ENDDO
       !
       IF (loop>=4) THEN
          IF (.NOT.(v4.LE.v3.AND.v3.LE.v2 &
&              .AND.v2.LE.v1).AND.v4.LE.cnvrg) THEN
             !
             !  converged result
             !
             WRITE(6,*) '  FCSCFatom converged in',loop,' iterations'
             AESCF%iter=loop
             WRITE(6,*) '     for nz = ',nz
             WRITE(6,*) '    delta(density)  = ', delta
             AESCF%delta=delta
             WRITE(6,*) '  results for loop = ',loop
             WRITE(6,*) ' n  l     occupancy       energy'
             DO io=1,norbit
                if(.NOT.FC%iscore(io)) WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
&                    AEOrbit%np(io),AEOrbit%l(io),&
&                    AEOrbit%occ(io),AEOrbit%eig(io)
             ENDDO
             CALL FreeAnderson(FCrho)
             DEALLOCATE(denout)
             RETURN
          ENDIF
       ENDIF
    ENDDO   ! mxloop

    WRITE(6,*)'calculation terminating without density convergence'
    WRITE(6,*) 'delta, cnvrg =', delta,cnvrg
    WRITE(6,*) 'loop,mxloop = ',loop,mxloop
    CALL FreeAnderson(FCrho)
    DEALLOCATE(denout)
    STOP
  END SUBROUTINE FCSCFloop

END MODULE AEatom
