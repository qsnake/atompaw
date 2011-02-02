PROGRAM atompaw
  !***************************************************************
  !  08-13-07 New options for the grid definition (r_max, r_match)
  !           Printing of pseudo valence density in atomicdata file
  !  12-21-06 Set version to be  2.0
  !  11-30-06 Many new options implemented by Marc Torrent -- see
  !      user guide
  !   9-26-06 Upgraded scalarrelativistic option by (1) using solver
  !      cfdsol obtained from David Vanderbilt's USPS code and
  !      adapted by Marc Torrent and Francois Jollet (2) introducing
  !      finite-nucleus option which can replace -2*Z/r potential
  !      with -2*Z*erf(r/RR)/r, where RR is a nuclear size parameter.
  !   7-07-06 added option to use Gaussian shape for hat density
  !      also simplified other options for projector generation
  !          keywords Bloechl for Peter Bloechl's scheme (VNCT)
  !                   Vanderbilt for David Vanderbilt' scheme (VNCTV)
  !      To specify Gaussian shape functions use
  !             Bloechl   Gaussian     1.d-4    or
  !             Vanderbilt   Gaussian  1.d-4
  !         where Gaussian = exp(-(r/d)^2) and d=rc*ln(1/1.d-4)
  !      To specify Sinc^2 shape function, line should be blank
  !          after projector keyword
  !           Bloechl (or VNCT)   or
  !           Vandervilt (or VNCTV)
  !   6-14-06 corrected xml portion to be consistent with FSATOM
  !     standard on website
  !       http://dcwww.camp.dtu.dk/campos//pawxml/pawxml.xhtml
  !   6-10-06 worked with Marc Torrent to validate Atompaw2abinit interface
  !   5-15-06 added option to use logarithmic grid for radial functions
  !   1-21-06 added ionic local potential output needed by
  !     abinit code to [atom].atomicdata file  -- function not
  !     added to xml output;   not currently clear that xml output
  !     is correctly implemented for current "standard" described in
  !     http://dcwww.camp.dtu.dk/campos//atomic_setup/atomic_setup.xhtml
  !  12-20-05 implemented option to construct projector and basis
  !     functions similar to David Vanderbilt's ultra-soft
  !     pseudopotentials (PRB 41, 7892 (1990)
  !   1-03-05 implemented PAW-XML output in conformance with
  !     http://www.fysik.dtu.dk/campos/atomic_setup/paw_setup.html
  !  12-31-04 introduced CoreTail density
  !   4-30-04 minor changes to simplify options
  !   2-29-04 major changes to code structure -- adding possibility
  !    of adjusting Vloc to norm conserving pseudopotential at given
  !      l
  !  4-19-00 pgm written by N. A. W. Holzwarth
  !     Calculates projector and basis functions needed by pwpaw pgm
  !       for electronic structure calculations using the PAW method
  !       of Blochl
  !     Modified verion of original genproj pgm
  !***************************************************************
  USE GlobalMath
  USE atomdata
  USE aeatom
  USE gridmod
  USE pseudo
  USE basis
  USE abinitinterface
  USE pwscfinterface
  USE excor
  USE libxc_mod
  USE pkginfo

  IMPLICIT NONE
  CHARACTER (len=4) :: flnm
  CHARACTER (len=20) :: nm
  CHARACTER (len=2) :: sym
  CHARACTER (len=1) :: syml
  REAL(8), POINTER :: r(:),den(:),rv(:),wfn(:,:)
  INTEGER, POINTER :: n,norbit,nps,npp,npd,npf,npg
  REAL(8), POINTER :: h
  !INTEGER, PARAMETER :: VNCT=1,VNCK=2,VSHAPE=3,VNCTV=4
  INTEGER, PARAMETER :: BLOECHL=1, VANDERBILT=2, CUSTOM=3
  INTEGER, PARAMETER :: BLOECHLPS=0, POLYNOM=1, POLYNOM2=2, RRKJ=3
  INTEGER, PARAMETER :: VANDERBILTORTHO=0, GRAMSCHMIDTORTHO=1
  INTEGER, PARAMETER :: MTROULLIER=1, ULTRASOFT=2, BESSEL=3
  CHARACTER(50) :: Projectortype,inputfilename
  CHARACTER(80) :: inputfileline
  CHARACTER(3) :: gridtype
  INTEGER :: Projectorindex,PSindex,Orthoindex,Vlocalindex
  INTEGER :: i,j,io,many,l,istart,irc,ishift,OK,nbase,ib,ic,icount
  INTEGER :: llmin, llmax, ll,lcount,jcount,id,ie,iop,lp,pdeg
  REAL(8) :: rc1,rc2,rc3,rc4,rc,x,rr,e,eself,cc,q00,qcut,evale,storeself,qcore,Etotal
  REAL(8) :: tildekin,tildepot,onehat,onehartree,ctexc,ctexc1,cexc1,texc,vexc1
  REAL(8) :: vtexc,vtexc1,oneenergy,ekin,vlocal,tq,fac,stuff,term,sqr4pi
  REAL(8) :: ctctse,cthatse,selfenergy
  REAL(8), ALLOCATABLE :: ttphi(:),soij(:),stij(:),svij(:)
  REAL(8), ALLOCATABLE :: shij(:),snij(:),dum(:),dum1(:),rh(:),rth(:)
  REAL(8), ALLOCATABLE :: shartree(:),sshartree(:,:,:),wf(:),twf(:)
  !REAL(8), ALLOCATABLE :: scself1(:),scself2(:,:,:)
  REAL(8), ALLOCATABLE :: bm(:),cm(:),dm(:,:)
  INTEGER, PARAMETER :: ifen=11,ifatompaw=12,ifout=13,ifxml=14,ifinput=15
  !INTEGER, PARAMETER :: ifself=16
  INTEGER :: lcao_points,lcao_i,iskip
  REAL(8) :: tphirange,hlcao
  REAL(8), PARAMETER :: coretailtol=1.d-7, gausstol=1.d-4
  LOGICAL :: even,fileexists,multi_rc
  CHARACTER(132) :: inputline,inputword
  TYPE (GridInfo), TARGET :: AEGrid
  TYPE (PotentialInfo), TARGET :: AEPot
  TYPE (OrbitInfo), TARGET :: AEOrbit
  TYPE (OrbitInfo), TARGET :: PSOrbit
  TYPE (SCFInfo), TARGET :: AESCF
  TYPE (FCInfo), TARGET :: FC
  TYPE (PseudoInfo), TARGET :: PAW

! First write out details on AtomPAW
  WRITE(6,*) atp_package,' v',atp_version
  WRITE(6,*) 'Compiled for ',atp_target
  WRITE(6,*)

  OPEN(ifinput,file='dummy',form='formatted')

  CALL Init_GlobalConstants
  CALL iSCFatom(AEGrid,AEPot,AEOrbit,AESCF,ifinput)
  CALL ChooseValence(AEGrid,AEOrbit,FC,ifinput)
  CALL FCenergy(AEGrid,AEPot,AEOrbit,AESCF,FC)
  CALL FCselfenergy(AEGrid,AEOrbit,FC,selfenergy)

  OPEN(ifen,file=TRIM(AEPot%sym),form='formatted')
  WRITE(ifen,'("Atom = ",a2,"  Z = ",f4.0)')  AEPot%sym, AEPot%nz
  SELECT CASE(TRIM(exctype))
  CASE default
   if (have_libxc) then
!   call libxc_print_func(6)
    call libxc_print_func(ifen)
   else
!   WRITE(6,*) 'Perdew-Wang correlation'
    WRITE(ifen,*) 'Perdew-Wang correlation'
   end if
  CASE('LDA-PW')
!    WRITE(6,*) 'Perdew-Wang correlation'
     WRITE(ifen,*) 'Perdew-Wang correlation'
  CASE('GGA-PBE')
!    WRITE(6,*) 'Perdew - Burke - Ernzerhof GGA'
     WRITE(ifen,*) 'Perdew - Burke - Ernzerhof GGA'
  END SELECT

  n=>AEGrid%n
  h=>AEGrid%h
  r=>AEGrid%r
  rv=>AEPot%rv
  ishift=AEGrid%ishift
  if (usingloggrid(AEGrid)) then
      gridtype="log"
      WRITE(ifen,'("Log grid -- n,r0,rmax = ",i5,1p,2e15.7)') &
&      AEGrid%n,AEGrid%drdu(1),AEGrid%r(n)
  else
      gridtype="lin"
      WRITE(ifen,'("Linear grid --  n,rmax = ",i5,1p,2e15.7)') &
&      AEGrid%n,AEGrid%r(n)
  endif
  if (scalarrelativistic) then
      if(.not.finitenucleus) then
          WRITE(ifen,*) 'Scalar relativistic calculation -- point nucleus'
      else
          WRITE(ifen,*) &
&            'Scalar relativistic calculation -- finite (Gaussian) nucleus'
      endif
  else
      WRITE(ifen,*) 'Non-relativistic calculation'
  endif

  WRITE(ifen,*) ' all-electron results '
  WRITE(ifen,*) ' core states (zcore) = ',FC%zcore
  DO io=1,FC%norbit
     IF (FC%iscore(io))THEN
        WRITE(ifen,'(3i5,1p,2e15.7)') io,AEOrbit%np(io),AEOrbit%l(io),&
&            AEOrbit%occ(io),AEOrbit%eig(io)
     ENDIF
  ENDDO
  WRITE(ifen,*) ' valence states (zvale) = ',FC%zvale
  DO io=1,FC%norbit
     IF (.NOT.FC%iscore(io)) THEN
        WRITE(ifen,'(3i5,1p,2e15.7)') io,AEOrbit%np(io),AEOrbit%l(io),&
&            AEOrbit%occ(io),AEOrbit%eig(io)
     ENDIF
  ENDDO
  WRITE(ifen,*) 'evale = ', FC%evale
  WRITE(ifen,*) 'selfenergy contribution = ', selfenergy

  WRITE(6,*) 'Enter maximum L for basis and projector functions'
  READ(5,'(a)') inputline
  WRITE(ifinput,'(a)') TRIM(inputline)
  READ(inputline,*) PAW%lmax

  PAW%rc=0;multi_rc=.false.
  PAW%rc_shap=0;PAW%rc_vloc=0;PAW%rc_core=0
  WRITE(6,*) 'enter  rc [and eventually: rc_shape, rc_vloc, rc_core]'
  READ(5,'(a)') inputline
  WRITE(ifinput,'(a)') TRIM(inputline)
  CALL extractword(1,inputline,inputword);inputword=trim(inputword)
  IF (inputword/="") READ(inputword,*) rc1
  rc2=rc1;rc3=rc1;rc4=rc1
  CALL extractword(2,inputline,inputword);inputword=trim(inputword)
  IF (inputword/="") THEN
   multi_rc=.true.
   READ(inputword,*) rc2
   CALL extractword(3,inputline,inputword);inputword=trim(inputword)
   IF (inputword/="") THEN
    READ(inputword,*) rc3
    CALL extractword(4,inputline,inputword);inputword=trim(inputword)
    IF (inputword/="") THEN
     READ(inputword,*) rc4
    ELSE
     WRITE(6,*) 'error -- rc(core) is missing '
     STOP
    ENDIF
   ELSE
    WRITE(6,*) 'error -- rc(Vloc) is missing '
    STOP
   ENDIF
  ENDIF
  IF (multi_rc) THEN
   IF (rc1.LE.0.d0.or.rc2.LE.0.d0.or.rc3.LE.0.d0.or.rc4.LE.0.d0) THEN
     WRITE(6,*) 'error -- one rc is too small !'
     STOP
   ENDIF
   IF (rc2.GT.rc1.or.rc3.GT.rc1.or.rc4.GT.rc1) THEN
     WRITE(6,*) 'error -- rc_shape, rc_vloc and rc_core must be <rc !'
     STOP
   ENDIF
  ELSEIF (rc1.LE.0.d0) THEN
     WRITE(6,*) 'error -- rc too small ',rc1
     STOP
  ENDIF
  PAW%irc     =FindGridIndex(AEGrid,rc1) ; PAW%rc     =r(PAW%irc)
  PAW%irc_shap=FindGridIndex(AEGrid,rc2) ; PAW%rc_shap=r(PAW%irc_shap)
  PAW%irc_vloc=FindGridIndex(AEGrid,rc3) ; PAW%rc_vloc=r(PAW%irc_vloc)
  PAW%irc_core=FindGridIndex(AEGrid,rc4) ; PAW%rc_core=r(PAW%irc_core)
  irc=PAW%irc;rc=PAW%rc
  if (irc>n-ishift) stop 'error -- rc is too big !'
  WRITE(6,*) ' adjusted rc ',rc, r(irc)
  WRITE(6,*) ' irc,rc = ',irc,rc
  if (multi_rc) then
   WRITE(6,*) ' adjusted rc_shape ',PAW%rc_shap
   WRITE(6,*) ' adjusted rc_vloc  ',PAW%rc_vloc
   WRITE(6,*) ' adjusted rc_core  ',PAW%rc_core
  endif
  WRITE(ifen,*) ' paw parameters: '
  WRITE(ifen,*) '      lmax = ',PAW%lmax
  WRITE(ifen,*) '        rc = ',rc
  WRITE(ifen,*) '       irc = ',irc
  if (multi_rc) then
   WRITE(ifen,*) '  rc_shape = ',PAW%rc_shap
   WRITE(ifen,*) '   rc_vloc = ',PAW%rc_vloc
   WRITE(ifen,*) '   rc_core = ',PAW%rc_core
  endif
  !
  CALL initpseudopot(AEGrid,PAW,AEOrbit,FC)
  CALL setbasis(AEGrid,AEPot,AEOrbit,FC,PAW,ifinput)

  ALLOCATE(ttphi(n),dum(n),dum1(n),stat=OK)
  IF(OK /= 0) THEN
     WRITE(6,*) 'Error in allocating arrays', n
     STOP
  ENDIF
  CALL setcoretail(AEGrid,FC%coreden,PAW)
  ! Find coretailpoints
     coretailpoints=irc+ishift
        do i=irc+ishift,n
           if(ABS(FC%zcore-integrator(AEGrid,FC%coreden,1,i))<coretailtol) then
             coretailpoints=i
             exit
           endif
        enddo
     write(6,*) 'coretailpoints = ',coretailpoints

!  simplified input for most useful options
!  WRITE(6,*) ' This code uses SINC2 shape for projector '&
!&      ,' generation and for hat densities'
!  CALL sethat(AEGrid,PAW)
!  WRITE(6,*) ' Enter type of Vloc -- ',&
!&  'VNCT (Norm conserving - Troullier-Martins)',&
!&  'VNCTV (Norm conserving - Troullier-Martins with Vanderbilt projectors)',&
!&  'VNCK (Norm conserving - Kerker)',&
!&  ' or VSHAPE (for fixed shape form)'
!  READ(5,*) Vlocaltype
!  Projectorindex=BLOECHL;Vlocalindex=MTROULLIER
!  IF (TRIM(Vlocaltype)=="VNCT") Projectorindex=BLOECHL
!  IF (TRIM(Vlocaltype)=="VNCTV") Projectorindex=VANDERBILT
!  IF (TRIM(Vlocaltype)=="VNCK") vlocalindex=VNCK
!  IF (TRIM(Vlocaltype)=="VSHAPE") vlocalindex=VSHAPE

   WRITE(6,*) 'Enter "Bloechl", "Vanderbilt", or "custom" keywords',&
&            ' for projector generation method.'
   WRITE(6,*) ' In case of "custom" choice, enter additional (optional) keywords:'
   WRITE(6,*) ' - for partial waves pseudization scheme:'
   WRITE(6,*) '                       "bloechlps", "polynom", "polynom2 p qcut" or "RRKJ"'
   WRITE(6,*) ' - for orthogonalization scheme: "GramSchmidtOrtho" or "VanderbiltOrtho"'
   WRITE(6,*) 'Compensation charge shape defaults set to "sinc^2";'
   WRITE(6,*) ' - Gaussian shape can be specified by adding "Gaussian" keyword',&
&             '   and tol (1.d-4, for ex).'
   WRITE(6,*) ' - Bessel shape can be specified by adding "Besselshape" keyword'
   READ(5,'(a)') inputfileline
   WRITE(ifinput,'(a)') TRIM(inputfileline)
   call Uppercase(inputfileline)
   inputfileline=TRIM(inputfileline)
   Projectorindex=BLOECHL;PSindex=BLOECHLPS;Orthoindex=GRAMSCHMIDTORTHO
   pdeg=4;qcut=10.d0
   read(unit=inputfileline,fmt=*) Projectortype
   if (TRIM(Projectortype)=='BLOECHL'.or.TRIM(Projectortype)=='VNCT') then
    Projectorindex=BLOECHL;PSindex=BLOECHLPS;Orthoindex=GRAMSCHMIDTORTHO
   else if (TRIM(Projectortype)=='VANDERBILT'.or.TRIM(Projectortype)=='VNCTV') then
    Projectorindex=VANDERBILT;PSindex=POLYNOM;Orthoindex=VANDERBILTORTHO
   else if (TRIM(Projectortype)=='CUSTOM') then
    Projectorindex=CUSTOM
    i=0;i=INDEX(inputfileline,'BLOECHLPS')
    if (i>0) then
     PSindex=BLOECHLPS;Orthoindex=GRAMSCHMIDTORTHO
    else
     i=INDEX(inputfileline,'POLYNOM2')
     if (i>0) then
      PSindex=POLYNOM2
      read(unit=inputfileline(i+8:80),fmt=*,err=111,end=111,iostat=i) pdeg,qcut
111   continue
     else
      i=INDEX(inputfileline,'POLYNOM')
      if (i>0) then
       PSindex=POLYNOM
      else
       i=INDEX(inputfileline,'RRKJ')
       if (i>0) PSindex=RRKJ
      endif
     endif
     i=INDEX(inputfileline,'GRAMSCHMIDTORTHO')
     if (i>0) Orthoindex=GRAMSCHMIDTORTHO
     i=INDEX(inputfileline,'VANDERBILTORTHO')
     if (i>0) Orthoindex=VANDERBILTORTHO
    endif
   endif
   write(PAW%Proj_description,'("Projector method:")')
   if (PSindex==BLOECHLPS) then
     if (Orthoindex==VANDERBILTORTHO) stop &
&     'Vanderbilt orthogonalization not compatible with Bloechl s projector scheme !'
    write(PAW%Proj_description,'(a," Bloechl")') trim(PAW%Proj_description)
   else
    if (Orthoindex==VANDERBILTORTHO) &
&    write(PAW%Proj_description,'(a," Vanderbilt (")') trim(PAW%Proj_description)
    if (PSindex==POLYNOM) then
     write(PAW%Proj_description,'(a,"polynomial pseudization")') trim(PAW%Proj_description)
    else if (PSindex==POLYNOM2) then
     write(PAW%Proj_description,'(a,"improved polynomial pseudization")') trim(PAW%Proj_description)
    else if (PSindex==RRKJ) then
     write(PAW%Proj_description,'(a,"RRKJ pseudization")') trim(PAW%Proj_description)
    endif
    if (Orthoindex==VANDERBILTORTHO) then
     write(PAW%Proj_description,'(a,")")') trim(PAW%Proj_description)
    else
     write(PAW%Proj_description,'(a," + Gram-Schmidt ortho.")') trim(PAW%Proj_description)
    endif
   endif

   gaussianshapefunction=.false.;besselshapefunction=.false.
   i=0;i=INDEX(inputfileline,'GAUSSIAN')
   if (i>0) then
      gaussianshapefunction=.true.
      gaussparam=gausstol
      read(unit=inputfileline(i+8:80),fmt=*) x
      if (x>0) gaussparam=x
      CALL sethat(AEGrid,PAW,gaussparam=gaussparam)    ! Gaussian shape function
      write(PAW%Comp_description,&
&      '("Gaussian compensation charge shape with gausstol = ",1pe12.4)')&
&         gaussparam
   else
    i=0;i=INDEX(inputfileline,'BESSELSHAPE')
    if (i>0) then
     besselshapefunction=.true.
     CALL sethat(AEGrid,PAW,besselopt=i)               ! Bessel shape function
     if (PAW%irc_shap/=PAW%irc) then
      write(PAW%Comp_description,&
&      '("Bessel compensation charge shape zeroed at ",1pe12.4)') PAW%rc_shap
     else
      write(PAW%Comp_description,&
&      '("Bessel compensation charge shape zeroed at rc")')
     endif
    else
     CALL sethat(AEGrid,PAW)                          ! sinc^2 shape function
     if (PAW%irc_shap/=PAW%irc) then
      write(PAW%Comp_description,&
&      '("Sinc^2 compensation charge shape zeroed at ",1pe12.4)') PAW%rc_shap
     else
      write(PAW%Comp_description,&
&      '("Sinc^2 compensation charge shape zeroed at rc")')
     endif
    endif
   endif

   WRITE(6,*) 'To generate the local pseudopotential, this code can use:'
   WRITE(6,*) '  1- a Troullier-Martins scheme for specified l value and energy'
   WRITE(6,*) '  2- a non norm-conserving pseudopotential scheme for specified l value and energy'
   WRITE(6,*) '  3- a simple pseudization scheme using a single spherical Bessel function'
   WRITE(6,*) 'For choice 1, enter (high) l value and energy e'
   WRITE(6,*) 'For choice 2, enter (high) l value, energy e and "ultrasoft"'
   WRITE(6,*) 'For choice 3, enter "bessel"'
   READ(5,'(a)') inputfileline
   WRITE(ifinput,'(a)') TRIM(inputfileline)
   call Uppercase(inputfileline)
   Vlocalindex=MTROULLIER
   i=0;i=INDEX(inputfileline,'BESSEL')
   if (i>0) then
    Vlocalindex=BESSEL
    WRITE(PAW%Vloc_description,'("Vloc: truncated form - Vps(r)=A.sin(qr)/r for r<rc")')
   else
    READ(inputfileline,*) l,e
    if (l<0.or.l>10) stop 'Error while reading Vloc parameters'
    i=0;i=INDEX(inputfileline,'ULTRASOFT')
    if (i>0) then
     Vlocalindex=ULTRASOFT
     WRITE(PAW%Vloc_description,&
&        '("Vloc: Non norm-conserving form with l= ",i1,";e= ",1pe12.4)')l,e
    else
     Vlocalindex=MTROULLIER
     WRITE(PAW%Vloc_description,&
&        '("Vloc: Norm-conserving Troullier-Martins form; l= ",i1,";e= ",1pe12.4)')l,e
    endif
   endif
   WRITE(6,*) PAW%Vloc_description

   IF (Vlocalindex==MTROULLIER) CALL troullier(AEGrid,AEPot,PAW,l,e)
   IF (Vlocalindex==ULTRASOFT) CALL nonncps(AEGrid,AEPot,PAW,l,e)
   IF (Vlocalindex==BESSEL) CALL besselps(AEGrid,AEPot,PAW)

   IF (Projectorindex==BLOECHL) THEN
    CALL makebasis_bloechl(AEGrid,AEPot,PAW,ifinput,0)
   ELSE IF (Projectorindex==CUSTOM.AND.PSindex==BLOECHLPS) THEN
    CALL makebasis_bloechl(AEGrid,AEPot,PAW,ifinput,1)
   ELSE IF (Projectorindex==VANDERBILT.OR.Projectorindex==CUSTOM) THEN
    CALL makebasis_custom(AEGrid,AEPot,PAW,ifinput,PSindex,Orthoindex,pdeg,qcut)
   ENDIF

   CALL FindVlocfromVeff(AEGrid,FC,AEPot%nz,PAW)

  !IF (vlocalindex==VNCK) THEN
  !   CALL kerker(AEGrid,AEPot,PAW)
  !   CALL makebasis(AEGrid,PAW)
  !   CALL FindVlocfromVeff(AEGrid,FC,AEPot%nz,PAW)
  !   PAW%Vloc_description=&
  !&        'Bloechl projectors with '//PAW%Vloc_description(1:250)
  !ENDIF
!
!  IF (vlocalindex==VSHAPE) THEN
!     CALL SCbasis(AEGrid,FC,AEPot%nz,PAW)
!     PAW%Vloc_description=&
!&           'Bloechl projectors with '//PAW%Vloc_description(1:250)
!  ENDIF

      WRITE(ifen,*) TRIM(PAW%Vloc_description)
      WRITE(ifen,*) TRIM(PAW%Proj_description)
      WRITE(ifen,*) TRIM(PAW%Comp_description)

  CALL checkghosts(AEGrid,AEOrbit,FC,PAW)

  OPEN(ifout,file='density', form='formatted')
  DO i=1,n
     IF (FC%coreden(i)<machine_zero) FC%coreden(i)=0
     IF (PAW%tcore(i)<machine_zero) PAW%tcore(i)=0
     IF (PAW%den(i)<machine_zero) PAW%den(i)=0
     IF (PAW%tden(i)<machine_zero) PAW%tden(i)=0
     WRITE(ifout,'(1p,1e15.7,1p,4e25.17)') r(i),FC%coreden(i),&
&         PAW%den(i),PAW%tden(i),PAW%tcore(i)
  ENDDO
  CLOSE(ifout)

  OPEN(ifout,file='potential', form='formatted')
  DO i=1,n
     IF (ABS(rv(i))<machine_zero) rv(i)=0
     IF (ABS(PAW%rveff(i))<machine_zero) PAW%rveff(i)=0
     WRITE(ifout,'(1p,1e15.7,1p,3e25.17)') r(i),AEPot%rv(i),PAW%rveff(i)
  ENDDO
  CLOSE(ifout)

  OPEN(ifout,file='vloc', form='formatted')
  DO i=1,irc+10
     WRITE(ifout,'(1p,1e15.7,1p,3e25.17)') r(i),PAW%vloc(i)
  ENDDO
  CLOSE(ifout)

  nbase=PAW%nbase
  WRITE(ifen,'(/"Number of basis functions ",i5)') nbase
  WRITE(ifen,*)'No.   n    l      Energy         Cp coeff         Occ'

  DO io=1,nbase
     WRITE(ifen,'(3i5,1p,3e15.7)') io,PAW%np(io),PAW%l(io),PAW%eig(io),&
&         PAW%ck(io),PAW%occ(io)
     CALL mkname(io,flnm)
     OPEN(ifout,file='wfn'//TRIM(flnm),form='formatted')
     WRITE(ifout,*) '# l=',PAW%l(io),'basis function with energy  ',PAW%eig(io)
       DO i=1,irc+50
          WRITE(ifout,'(1p,5e12.4)') r(i),PAW%ophi(i,io),PAW%otphi(i,io),PAW%otp(i,io)
       ENDDO
       CLOSE(ifout)
    ENDDO

!   Find radii ensuring:
!     - Abs(density)<10e-10 (for the pseudo-valence density)
!     - r>=10 bohr (for the ionic potential)
    sqr4pi=sqrt(4*pi)*1.d-10;ivion=gridindex(AEGrid,10.d0);ivale=ivion
    do while (ivale<AEGrid%n.and.abs(PAW%tden(ivale))>sqr4pi*AEGrid%r(ivale)**2)
      ivale=ivale+1
    end do

    !Begin  constructing output  for pwpaw pgm

    qcore=FC%zcore
    !OPEN(ifself,file=TRIM(AEpot%sym)//'.scself',form='formatted')
    OPEN(ifatompaw,file=TRIM(AEpot%sym)//'.atomicdata',form='formatted')
    WRITE(ifatompaw,'("  ATOMTYPE     ",a2)') AEpot%sym
    WRITE(ifatompaw,'("  ATOMXCTYPE     ",a)') trim(exctype)
    WRITE(ifatompaw,'("  ATOMIC_CHARGE    ",f5.0)') AEpot%nz
    WRITE(ifatompaw,'("  CORE_CHARGE    ",1pe20.13)') qcore
    WRITE(ifatompaw,'("  RC         ",1pe20.13)') PAW%rc
    if (gaussianshapefunction) then
     WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1pe20.13)') 'gaussian', &
&          PAW%rc_shap/SQRT(LOG(1.d0/gaussparam))
    else if (besselshapefunction) then
     if (multi_rc) then
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1pe20.13)') 'bessel',PAW%rc_shap
     else
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20)') 'bessel'
     endif
    else
     if (multi_rc) then
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1pe20.13)') 'sinc2',PAW%rc_shap
     else
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20)') 'sinc2'
     endif
    endif
    WRITE(ifatompaw,'("  BASIS_SIZE    ",i5)') PAW%nbase
    WRITE(ifatompaw,'("  ORBITALS      ")')
    WRITE(ifatompaw,'(20i4)') (PAW%l(ib),ib=1,nbase)
    WRITE(ifatompaw,'("  END     ")')
    WRITE(ifatompaw,'("  INITOCC       ")')
    WRITE(ifatompaw,'(8f10.6)') (PAW%occ(ib),ib=1,nbase)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("  MESH_SIZE    ",i10)') PAW%irc+ishift
    WRITE(ifatompaw,'("  MESH_STEP    ",1pe20.13)') AEGrid%h
    if (usingloggrid(AEGrid)) then
       WRITE(ifatompaw,'("  LOG_GRID    ",1pe20.13)') AEGrid%drdu(1)
    endif

    WRITE(ifatompaw,'("  CORETAIL_POINTS   ",i10)') coretailpoints

! xml part still in flux
    sqr4pi=sqrt(4*pi)
    OPEN(ifxml,file=TRIM(AEpot%sym)//'.'//TRIM(exctype)//'.xml',&
&          form='formatted')
    WRITE(ifxml,'("<?xml  version=""1.0""?>")')
    !WRITE(ifxml,'("<!DOCTYPE paw_setup SYSTEM ""paw_setup.dtd"">")')
    WRITE(ifxml,'("<paw_setup version=""0.5"">")')
    !note energy units are Hartrees
    i=qcore+0.2
    flnm=stripchar('"'//AEPot%sym//'"')
    !WRITE(ifxml,'("<atom symbol=",a4," Z=""",i3,""" core=""", i3,""" valence=""",i3"""/>")')&
    !&       TRIM(flnm),AEpot%nz,i,AEpot%nz-i
    WRITE(ifxml,'("<atom symbol=",a4," Z=""",i3,$)') TRIM(flnm),INT(AEpot%nz)
    WRITE(ifxml,'(""" core=""", i3,""" valence=""",i3"""/>")')i,INT(AEpot%nz)-i
    if (TRIM(exctype)=='LDA-PW') WRITE(ifxml,&
&         '("<xc_functional type=""LDA"" name=""PW""/>")')
    if (TRIM(exctype)=='GGA-PBE') WRITE(ifxml,&
&         '("<xc_functional type=""GGA"" name=""PBE""/>")')
    if (scalarrelativistic) then
    WRITE(ifxml,'("<generator type=""scalar-relativistic"" name=""atompaw"">")')
    else
    WRITE(ifxml,'("<generator type=""non-relativistic"" name=""atompaw"">")')
    endif
    WRITE(ifxml,'("</generator>)")')
    WRITE(ifxml,'("<!-- Contact info: email: natalie@wfu.edu")')
    WRITE(ifxml,'("               web: pwpaw.wfu.edu")')
    WRITE(ifxml,'(" Energy units=Hartree, length units=bohr")')
    WRITE(ifxml,'(" Note: consistent with 06-14-06 standard")')
    WRITE(ifxml,'(" As discussed at CECAM PAW Workshop     ")')
    Call PrintDate(ifxml, ' PAW functions generated on ')
    WRITE(ifxml,'("  ",a)') Trim(PAW%Vloc_description)
    WRITE(ifxml,'("  ",a)') Trim(PAW%Proj_description)
    WRITE(ifxml,'("  ",a)') Trim(PAW%Comp_description)
    WRITE(ifxml,'("  Program:  atompaw - input data follows: ")')
     rewind(ifinput)
         do
           read(ifinput,'(a)',iostat=i) inputfileline
           if (i/=0) exit
           write(ifxml,'("   ",a)') TRIM(inputfileline)
        enddo
    WRITE(ifxml,'("  Program:  atompaw - input end -->")')
!
    WRITE(ifxml,'("<ae_energy kinetic=""",1pe25.17,""" xc=""",1pe25.17,"""")')&
&           AESCF%ekin/2,AESCF%eexc/2
    WRITE(ifxml,'("  electrostatic=""",1pe25.17,""" total=""",1pe25.17,"""/>")') &
&           AESCF%estatic/2,AESCF%etot/2
    WRITE(ifxml,'("<core_energy kinetic=""",1pe25.17,"""/>")') FC%corekin/2
!
    WRITE(ifxml,'("<valence_states>")')
    do ib=1,nbase
       call mkname(ib,flnm)
       nm=stripchar('"'//AEPot%sym//flnm//'"')
       i=min(ABS(PAW%np(ib)),20)
       !WRITE(ifxml,&
       !& '("  <state n=""",i2,""" l=""",i1,""" f=""",1pe14.7,""" rc=""",1pe17.10,""" e=""",1pe14.7,""" id=",a6,"/>")')&
       !& i,PAW%l(ib),PAW%occ(ib),PAW%rc,PAW%eig(ib)/2,TRIM(nm)
       WRITE(ifxml,'("  <state n=""",i2,""" l=""",i1,""" f=""",1pe14.7,$)')&
&        i,PAW%l(ib),PAW%occ(ib)
       WRITE(ifxml,'(""" rc=""",1pe17.10,""" e=""",1pe14.7,""" id=",a6,"/>")')&
&        PAW%rc,PAW%eig(ib)/2,TRIM(nm)
    enddo
    WRITE(ifxml,'("</valence_states>")')
    if (usingloggrid(AEGrid)) then
     !WRITE(ifxml,'("<radial_grid eq=""r=a*(exp(d*i)-1)"" a=""",1pe17.10,""" d=""",1pe17.10,""" istart=""0"" iend=""",i5,""" id=""log1""/>")')&
     !&   AEGrid%drdu(1),AEGrid%h,coretailpoints-1
     WRITE(ifxml,'("<radial_grid eq=""r=a*(exp(d*i)-1)"" a=""",1pe17.10,$)')&
&        AEGrid%drdu(1)
     WRITE(ifxml,'(""" d=""",1pe17.10,""" istart=""0"" iend=""",i5,$)')&
&        AEGrid%h,coretailpoints-1
     WRITE(ifxml,'(""" id=""log1""/>")')
     !WRITE(ifxml,'("<radial_grid eq=""r=a*(exp(d*i)-1)"" a=""",1pe17.10,""" d=""",1pe17.10,""" istart=""0""  iend=""",i5,""" id=""log2""/>")')&
     !&   AEGrid%drdu(1),AEGrid%h,ivion-1
     WRITE(ifxml,'("<radial_grid eq=""r=a*(exp(d*i)-1)"" a=""",1pe17.10,$)')&
&        AEGrid%drdu(1)
     WRITE(ifxml,'(""" d=""",1pe17.10,""" istart=""0"" iend=""",i5,$)')&
&        AEGrid%h,ivion-1
     WRITE(ifxml,'(""" id=""log2""/>")')
     WRITE(ifxml,'("<radial_grid eq=""r=a*(exp(d*i)-1)"" a=""",1pe17.10,$)')&
&        AEGrid%drdu(1)
     WRITE(ifxml,'(""" d=""",1pe17.10,""" istart=""0"" iend=""",i5,$)')&
&        AEGrid%h,ivale-1
     WRITE(ifxml,'(""" id=""log3""/>")')
    else
     !WRITE(ifxml,'("<radial_grid eq=""r=d*i"" d=""",1pe17.10,""" istart=""0"" iend=""",i4,""" id=""lin1""/>")')&
     !&   AEGrid%h,coretailpoints-1
     WRITE(ifxml,'("<radial_grid eq=""r=d*i"" d=""",1pe17.10,$)') AEGrid%h
     WRITE(ifxml,'(""" istart=""0"" iend=""",i4,""" id=""lin1""/>")')&
&        coretailpoints-1
     !WRITE(ifxml,'("<radial_grid eq=""r=d*i"" d=""",1pe17.10,""" istart=""0"" iend=""",i5,""" id=""lin2""/>")')&
     !&   AEGrid%h,ivion-1
     WRITE(ifxml,'("<radial_grid eq=""r=d*i"" d=""",1pe17.10,$)') AEGrid%h
     WRITE(ifxml,'(""" istart=""0"" iend=""",i4,""" id=""lin2""/>")')&
&        ivion-1
     WRITE(ifxml,'("<radial_grid eq=""r=d*i"" d=""",1pe17.10,$)') AEGrid%h
     WRITE(ifxml,'(""" istart=""0"" iend=""",i4,""" id=""lin3""/>")')&
&        ivale-1
    endif
!    dum=0
!    dum(1:irc)=PAW%hatshape(1:irc)
!    WRITE(ifxml,'("<shape_function grid=""",a,i1,""">")') gridtype,1
    if (gaussianshapefunction) then
     WRITE(ifxml,'("<shape_function type=""gauss"" rc=""",1pe17.10"""/>")')&
&          PAW%rc_shap/SQRT(LOG(1.d0/gaussparam))
    else if (besselshapefunction) then
     WRITE(ifxml,'("<shape_function type=""bessel"" rc=""",1pe17.10"""/>")') PAW%rc_shap
    else
     WRITE(ifxml,'("<shape_function type=""sinc"" rc=""",1pe17.10"""/>")') PAW%rc_shap
    endif
!    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
!    WRITE(ifxml,'("</shape_function>")')
    dum=0
    dum(2:coretailpoints)=sqr4pi*FC%coreden(2:coretailpoints)&
&                /(4*pi*(AEGrid%r(2:coretailpoints))**2)
    call extrapolate(AEGrid,dum)
    WRITE(ifxml,'("<ae_core_density grid=""",a,i1""">")') gridtype,1
    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
    WRITE(ifxml,'("</ae_core_density>")')
    dum=0
    dum(2:coretailpoints)=sqr4pi*PAW%tcore(2:coretailpoints)&
&                /(4*pi*(AEGrid%r(2:coretailpoints))**2)
    call extrapolate(AEGrid,dum)
    WRITE(ifxml,'("<pseudo_core_density grid=""",a,i1,""">")') gridtype,1
    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
    WRITE(ifxml,'("</pseudo_core_density>")')
    dum=0
    dum(2:ivale)=sqr4pi*PAW%tden(2:ivale)&
&                /(4*pi*(AEGrid%r(2:ivale))**2)
    call extrapolate(AEGrid,dum)
    WRITE(ifxml,'("<pseudo_valence_density grid=""",a,i1,""">")') gridtype,3
    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,ivale)
    WRITE(ifxml,'("</pseudo_valence_density>")')
    dum=0
    dum(1:irc)=sqr4pi*PAW%vloc(1:irc)/2
    WRITE(ifxml,'("<zero_potential grid=""",a,i1,""">")') gridtype,1
    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
    WRITE(ifxml,'("</zero_potential>")')
    dum=0
    dum(1:ivion)=sqr4pi*PAW%abinitvloc(1:ivion)/2
    WRITE(ifxml,'("<kresse_joubert_local_ionic_potential grid=""",a,i1,""">")')&
&       gridtype,2
    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,ivion)
    WRITE(ifxml,'("</kresse_joubert_local_ionic_potential>")')
    dum=0
    dum(1:ivion)=sqr4pi*PAW%abinitnohat(1:ivion)/2
    WRITE(ifxml,'("<blochl_local_ionic_potential grid=""",a,i1,""">")')&
&       gridtype,2
    WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,ivion)
    WRITE(ifxml,'("</blochl_local_ionic_potential>")')

    Do ib=1,nbase
       call mkname(ib,flnm)
       nm=stripchar('"'//AEPot%sym//flnm//'"')
       dum=0
       dum(2:coretailpoints)=&
&           PAW%ophi(2:coretailpoints,ib)/AEGrid%r(2:coretailpoints)
       call extrapolate(AEGrid,dum)
       WRITE(ifxml,'("<ae_partial_wave state=",a6," grid=""",a,i1,""">")')&
&          TRIM(nm),gridtype,1
       WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
       WRITE(ifxml,'("</ae_partial_wave>")')
       dum=0
       dum(2:coretailpoints)=&
&            PAW%otphi(2:coretailpoints,ib)/AEGrid%r(2:coretailpoints)
       call extrapolate(AEGrid,dum)
       WRITE(ifxml,'("<pseudo_partial_wave state=",a6," grid=""",a,i1,""">")')&
&          TRIM(nm),gridtype,1
       WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
       WRITE(ifxml,'("</pseudo_partial_wave>")')
       dum=0
       dum(2:coretailpoints)=&
&            PAW%otp(2:coretailpoints,ib)/AEGrid%r(2:coretailpoints)
       call extrapolate(AEGrid,dum)
       WRITE(ifxml,'("<projector_function state=",a6," grid=""",a,i1,""">")') &
&            TRIM(nm),gridtype,1
       WRITE(ifxml,'(1p,3e25.17)') (dum(i),i=1,coretailpoints)
       WRITE(ifxml,'("</projector_function>")')
    Enddo



    ! new version of LCAO  -- assume range of 10 bohr sufficient
    !!!tphirange=10.d0       ! tphirange=8.d0
    !!!lcao_points=800       ! lcao_points=600
    !!!lcao_i=tphirange/AEGrid%h
    !!!iskip=MAX(1,lcao_i/lcao_points)
    !!!lcao_i=iskip*(lcao_i/iskip)+iskip
    !!!hlcao=AEGrid%h*iskip
    !!!lcao_points=lcao_i/iskip + 1
!   Find index for AEGrid%r(i)>10
    j=gridindex(AEGrid,10.d0)
    lcao_points=j
    hlcao=AEGrid%h
    WRITE(ifatompaw,'("  LCAO_SIZE  ",i10)') lcao_points
    WRITE(ifatompaw,'("  LCAO_STEP   ",1pe20.13)') hlcao

    WRITE(ifatompaw,'("   CORE_DENSITY   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (FC%coreden(i),i=1,PAW%irc+ishift)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   CORETAIL_DENSITY   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%tcore(i),i=1,coretailpoints)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   PSEUDO_VALENCE_DENSITY   ",3x,i8)') ivale
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%tden(i),i=1,ivale)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   SHAPE_FUNC   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%hatshape(i),i=1,PAW%irc+ishift)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   VLOCFUN      ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%vloc(i),i=1,PAW%irc+ishift)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,&
&    '("   VLOCION      ",3x,i8," #ionic vloc for abinit in Ryd units")') ivion
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%abinitvloc(i),i=1,ivion)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,&
&    '("   VLOCION_NOHAT",3x,i8," #ionic vlocnohat for abinit in Ryd units")') ivion
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%abinitnohat(i),i=1,ivion)
    WRITE(ifatompaw,'("  END     ")')

    DO ib=1,nbase
       WRITE(ifatompaw,'("   TPROJECTOR",i4," #p(r), for p(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%otp(i,ib),i=1,PAW%irc+ishift)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    DO ib=1,nbase
       WRITE(ifatompaw,'("   PHI",i5," #phi(r), for phi(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%ophi(i,ib),i=1,PAW%irc+ishift)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    DO ib=1,nbase
       WRITE(ifatompaw,'("   TPHI",i5," #tphi(r), for tphi(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%otphi(i,ib),i=1,PAW%irc+ishift)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO


    DO ib=1,nbase
       ttphi=PAW%tphi(:,ib)
       CALL trunk(AEGrid,ttphi(1:n),6.d0,10.d0)
       WRITE(ifatompaw,'("   TPHI_LCAO",i4," #tphi0(r) for tphi0(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (ttphi(j),j=1,lcao_points)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    ! spherical matrix elements
    icount=(nbase*(nbase+1))/2
    ALLOCATE(soij(icount),stij(icount),svij(icount),stat=ok)
    IF (ok /= 0) THEN
       WRITE(6,*) 'Error in allocating dia matrix elements',icount,ok
       STOP
    ENDIF

    icount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             l=PAW%l(ib)
             icount=icount+1
             CALL dqij(AEGrid,PAW,ib,ic,soij(icount))
             !CALL dtij(AEGrid,PAW,ib,ic,stij(icount))
             !CALL altdtij(AEGrid,AEPot,PAW,ib,ic,x)
             !  write(6,*) 'altdtij --', ib,ic,stij(icount),x
             CALL altdtij(AEGrid,PAW,ib,ic,stij(icount))
             CALL dvij(AEGrid,PAW,FC,AEPot%nz,ib,ic,svij(icount))
             PAW%tvij(ib,ic)=stij(icount)+svij(icount)
             PAW%tvij(ic,ib)=stij(icount)+svij(icount)
             PAW%kin(ib,ic)=stij(icount)
             PAW%kin(ic,ib)=stij(icount)
             PAW%oij(ib,ic)=soij(icount)
             PAW%oij(ic,ib)=soij(icount)
          ENDIF
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   OVERLAP_SIZE    ",i10)') icount
    WRITE(ifatompaw,'("   OVERLAP_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (soij(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')
    WRITE(ifatompaw,'("   KINETIC_ENERGY_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (stij(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')
    WRITE(ifatompaw,'("   V_ION_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (svij(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifxml,'("<kinetic_energy_differences>")')
       WRITE(ifxml,'(1p,3e25.17)') ((PAW%kin(ib,ic)/2,ic=1,nbase),ib=1,nbase)
    WRITE(ifxml,'("</kinetic_energy_differences>")')
    WRITE(ifxml,'("</paw_setup>")')

    !
    !  angularly dependent matrix elements
    !
    icount=(nbase*(nbase+1))*(PAW%lmax+1)
    ALLOCATE(shij(icount),snij(icount),stat=ok)
    IF (ok /= 0) THEN
       WRITE(6,*) 'Error in allocating L matrix elements',icount,ok
       STOP
    ENDIF

    icount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             icount=icount+1
             ! write(6,*) 'ib,ic,ll',ib,ic,ll,icount
             DO i=1,irc
                rr=AEGrid%r(i)
                dum(i)=(rr**ll)*(PAW%ophi(i,ib)*PAW%ophi(i,ic) &
&                    -PAW%otphi(i,ib)*PAW%otphi(i,ic))
             ENDDO
             snij(icount)=integrator(AEGrid,dum,1,irc)
             if (ll==0)  then
                PAW%v0ij(ib,ic)=snij(icount)
                PAW%v0ij(ic,ib)=snij(icount)
             endif
             CALL hatpotL(AEGrid,PAW,ll,dum)
             DO i=1,irc
                dum(i)=dum(i)*PAW%otphi(i,ib)*PAW%otphi(i,ic)
             ENDDO
             shij(icount)=integrator(AEGrid,dum,1,irc)
             if (ll==0)  then
                PAW%vhatij(ib,ic)=shij(icount)
                PAW%vhatij(ic,ib)=shij(icount)
             endif
          ENDDO
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   DENVHAT_SIZE   ",i10)') icount
    WRITE(ifatompaw,'("   DENSITY   ")')
    PAW%multipole=0.d0
    icount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             icount=icount+1
             WRITE(ifatompaw,'(3i10,1pe25.17)') ib,ic,ll,snij(icount)
             PAW%multipole(ib,ic,ll+1)=snij(icount)
             PAW%multipole(ic,ib,ll+1)=snij(icount)
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END   ")')
    WRITE(ifatompaw,'("   V_HAT   ")')
    icount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             icount=icount+1
             WRITE(ifatompaw,'(3i10,1pe25.17)') ib,ic,ll,shij(icount)
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END   ")')
    !
    icount=(nbase*(nbase+1))/2
    lcount=(icount**2)*((PAW%lmax+1)*2)
    !ALLOCATE(shartree(lcount),sshartree(icount,icount,2*PAW%lmax+2),&
    !&    scself1(lcount),scself2(icount,icount,2*PAW%lmax+2),&
    !&    wf(n),twf(n),rh(n),rth(n),stat=ok)
    ALLOCATE(shartree(lcount),sshartree(icount,icount,2*PAW%lmax+2),&
&        wf(n),twf(n),rh(n),rth(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'Error in hartree allocation', icount,lcount,irc,ok
       STOP
    ENDIF
    !
    !  Hartree matrix elements
    !
    !! Due to precision error in apoisson solver, some equivalent Hartree
    !!  matrix elements are not equal with an error of e-5; take average
    !!  these elements in order to remove assymmetry errors
    !!  Thanks to Francois Jollet for pointing out this problem
    sshartree=0
    shartree=0
    !scself1=0; scself2=0
    icount=0
    lcount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          icount=icount+1
          DO i=1,irc
             wf(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)
             twf(i)=PAW%otphi(i,ib)*PAW%otphi(i,ic)
          ENDDO
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             CALL apoisson(AEGrid,ll,irc,wf,rh)
             CALL apoisson(AEGrid,ll,irc,twf,rth)
             jcount=0
             DO id=1,nbase
                DO ie=id,nbase
                   jcount=jcount+1
                   lp=llmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      !WRITE(6,*) 'ib,ic,ll',ib,ic,id,ie,ll,lcount,icount,jcount
                      dum=0;dum1=0
                      DO i=2,irc
                         rr=AEGrid%r(i)
                         dum(i)=PAW%ophi(i,id)*PAW%ophi(i,ie)*rh(i)/rr     &
&                             -PAW%otphi(i,id)*PAW%otphi(i,ie)*rth(i)/rr
                         dum1(i)=PAW%ophi(i,id)*PAW%ophi(i,ie)*rh(i)/rr
                      ENDDO
                      shartree(lcount)=integrator(AEGrid,dum(1:irc),1,irc)
                      sshartree(icount,jcount,ll+1)=shartree(lcount)
                      !scself1(lcount)=integrator(AEGrid,dum1(1:irc),1,irc)
                      !scself2(icount,jcount,ll+1)=scself1(lcount)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    icount=0
    lcount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          icount=icount+1
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             jcount=0
             DO id=1,nbase
                DO ie=id,nbase
                   jcount=jcount+1
                   lp=llmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      shartree(lcount)=                      &
&                          0.5d0*(sshartree(icount,jcount,ll+1) + &
&                          sshartree(jcount,icount,ll+1))
                      !scself1(lcount)=                      &
                      !&    0.5d0*(scself2(icount,jcount,ll+1) + &
                      !&    scself2(jcount,icount,ll+1))
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    WRITE(6,*) lcount,'   ; # Hartree matrix elements'
    WRITE(ifatompaw,'("   HARTREE_SIZE   ",i10)') lcount
    WRITE(ifatompaw,'("   V_HARTREE   ")')
    !WRITE(ifself,'("   V_SCSELF   ",a2,i10)') AEPot%sym,lcount
    lcount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             DO id=1,nbase
                DO ie=id,nbase
                   lp=llmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      WRITE(ifatompaw,'(5i5,1pe25.17)')ib,ic,id,ie,ll,shartree(lcount)
                      !WRITE(ifself,'(5i5,1pe25.17)')ib,ic,id,ie,ll,scself1(lcount)

                      IF (ll==0) THEN
                        PAW%vhijkl(ib,ic,id,ie)=shartree(lcount)
                        PAW%vhijkl(ic,ib,id,ie)=shartree(lcount)
                        PAW%vhijkl(ib,ic,ie,id)=shartree(lcount)
                        PAW%vhijkl(ic,ib,ie,id)=shartree(lcount)
                        PAW%vhijkl(id,ie,ib,ic)=shartree(lcount)
                        PAW%vhijkl(id,ie,ic,ib)=shartree(lcount)
                        PAW%vhijkl(ie,id,ib,ic)=shartree(lcount)
                        PAW%vhijkl(ie,id,ic,ib)=shartree(lcount)
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END         ")')
    !WRITE(ifself,'("   END         ")')

    WRITE(ifatompaw,'("   HAT_SELF-ENERGY    ",i10)') 2*PAW%lmax
    DO ll=0,2*PAW%lmax
       CALL selfhatpot(AEGrid,PAW,ll,eself)
       WRITE(ifatompaw,'(i10,1pe25.17)') ll,eself
       IF (ll==0) storeself=eself
    ENDDO
    WRITE(ifatompaw,'("   END   ")')

    lcount=(nbase*(nbase+1))/2
    ALLOCATE(bm(nbase),cm(nbase),dm(lcount,lcount),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'Error in allocation before evale', nbase,ok
       STOP
    ENDIF

    !
    !  using matrix elements to re-evaluate evale --
    !
    !
    tildekin=0.d0
    PAW%den=0; PAW%tden=0
    DO ib=1,nbase
       l=PAW%l(ib)
       CALL kinetic(AEGrid,PAW%l(ib),PAW%tphi(:,ib),ekin)
       tildekin=tildekin+PAW%occ(ib)*ekin
       PAW%den=PAW%den+PAW%occ(ib)*(PAW%phi(:,ib))**2
       PAW%tden=PAW%tden+PAW%occ(ib)*(PAW%tphi(:,ib))**2
    ENDDO

    DO i=1,irc
       dum(i)=PAW%den(i)+FC%coreden(i)-PAW%tden(i)-PAW%tcore(i)
    ENDDO
    q00=-AEPot%nz+integrator(AEGrid,dum,1,irc)
    WRITE(6,*) 'q00 = ',q00
    DO i=1,n
       dum(i)=(PAW%tcore(i)+q00*PAW%hatden(i))
    ENDDO
    tq=integrator(AEGrid,dum)
    WRITE(6,*) 'tqcore  = ',tq
    CALL poisson(AEGrid,tq,dum,ttphi,x)
    do i=2,n
       ttphi(i)=ttphi(i)/r(i)+PAW%vloc(i)
    enddo
    call extrapolate(AEGrid,ttphi)
    tildepot=overlap(AEGrid,ttphi,PAW%tden)
    tq=integrator(AEGrid,PAW%tden)
    WRITE(6,*) 'tq  = ',tq
    CALL poisson(AEGrid,tq,PAW%tden,ttphi,x)
    tildepot=tildepot+x

    oneenergy=0.d0
    DO io=1,nbase
       l=PAW%l(io)
       bm=0.d0
       DO ib=1,nbase
          IF (PAW%l(ib).EQ.l)  THEN
             bm(ib)=overlap(AEGrid,PAW%otp(:,ib),PAW%tphi(:,io),1,irc)
             !IF (io==ib) THEN
             !   IF (ABS(bm(ib)-1).GT.1.d-11)&
             !&       WRITE(6,*) 'Warning in pdot',io,ib,bm(ib)
             !ELSE
             !   IF (ABS(bm(ib)).GT.1.d-11)&
             !&       WRITE(6,*) 'Warning in pdot',io,ib,bm(ib)
             !ENDIF
          ENDIF
       ENDDO
       lcount=0
       DO ib=1,nbase
          DO ic=ib,nbase
             IF (PAW%l(ib).EQ.PAW%l(ic)) THEN
                lcount=lcount+1
                fac=PAW%occ(io)
                IF(ic.GT.ib) fac=2*fac
                oneenergy=oneenergy+                                         &
&                    fac*bm(ib)*bm(ic)*(stij(lcount)+svij(lcount))
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    onehat=0.d0
    DO io=1,nbase
       l=PAW%l(io)
       bm=0.d0
       DO ib=1,nbase
          IF (PAW%l(ib).EQ.l)                                       &
&              bm(ib)=overlap(AEGrid,PAW%otp(:,ib),PAW%tphi(:,io),1,irc)
       ENDDO
       lcount=0
       DO ib=1,nbase
          DO ic=ib,nbase
             llmin=ABS(PAW%l(ib)-PAW%l(ic))
             llmax=PAW%l(ib)+PAW%l(ic)
             fac=PAW%occ(io)
             IF (ic.GT.ib) fac=2*fac
             DO ll=llmin,llmax,2
                lcount=lcount+1
                IF (ll.EQ.0) THEN
                   onehat=onehat+q00*fac*bm(ib)*bm(ic)*shij(lcount)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !onehat=onehat + storeself*(q00**2) -- self term excluded here

    onehartree=0.d0
    ! load array dm
    dm=0
    lcount=0
    icount=0
    DO ib=1,nbase
       DO ic=ib,nbase
          icount=icount+1
          llmin=ABS(PAW%l(ib)-PAW%l(ic))
          llmax=PAW%l(ib)+PAW%l(ic)
          DO ll=llmin,llmax,2
             jcount=0
             DO id=1,nbase
                DO ie=id,nbase
                   lp=llmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   jcount=jcount+1
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      IF (ll.EQ.0) dm(icount,jcount)=shartree(lcount)
                   ENDIF
                ENDDO !ie
             ENDDO   !id
          ENDDO    !ll
       ENDDO ! ic
    ENDDO  ! ib

    DO io=1,nbase
       l=PAW%l(io)
       bm=0.d0
       DO ib=1,nbase
          IF (PAW%l(ib).EQ.l)                                       &
&              bm(ib)=overlap(AEGrid,PAW%otp(:,ib),PAW%tphi(:,io),1,irc)
       ENDDO
       icount=0
       DO ib=1,nbase
          DO ic=ib,nbase
             icount=icount+1
             stuff=PAW%occ(io)*bm(ic)*bm(ib)
             IF (ic.GT.ib) stuff=2*stuff
             term=0.d0
             DO iop=1,nbase
                lp=PAW%l(iop)
                cm=0.d0
                DO id=1,nbase
                   IF (PAW%l(id).EQ.lp)                                      &
&                     cm(id)=overlap(AEGrid,PAW%otp(:,id),PAW%tphi(:,iop),1,irc)
                ENDDO
                jcount=0
                DO id=1,nbase
                   DO ie=id,nbase
                      jcount=jcount+1
                      fac=PAW%occ(iop)*cm(id)*cm(ie)
                      IF (ie.GT.id) fac=2*fac
                      term=term+fac*dm(icount,jcount)
                   ENDDO !ie
                ENDDO   !id
             ENDDO  !iop
             onehartree=onehartree+stuff*term
          ENDDO ! ic
       ENDDO  ! ib
    ENDDO  !io

    !  exchange-correlation energy
    !
    CALL exch(AEGrid,PAW%tden,dum,stuff,vtexc)
    CALL exch(AEGrid,PAW%tden,dum,stuff,vtexc1,irc)
    CALL exch(AEGrid,FC%coreden,dum,stuff,cc)
    CALL exch(AEGrid,FC%coreden,dum,stuff,cexc1,irc)
    cc=cc-cexc1
    DO i=1,irc
       ttphi(i)=FC%coreden(i)+PAW%den(i)
    ENDDO
    CALL exch(AEGrid,ttphi,wf,stuff,vexc1,irc)

    texc=vtexc+vexc1-vtexc1-cexc1
    WRITE(6,*) ' exchange-correlation energy', texc
    WRITE(6,*) ' core tail exchange error ', cc
    ttphi=PAW%tcore+PAW%tden
    CALL exch(AEGrid,ttphi,dum,stuff,vtexc)
    CALL exch(AEGrid,ttphi,twf,stuff,vtexc1,irc)
    ttphi=FC%coreden+PAW%den
    call exch(AEGrid,ttphi,wf,stuff,vexc1,irc)
    texc=vtexc+vexc1-vtexc1
    write(6,*) ' New form of exchange-correlation energy', texc
    !
    !  local potential contribution -- already in matrix element -- not added
    !
    vlocal=overlap(AEGrid,PAW%vloc(1:irc),PAW%tden(1:irc),1,irc)
    WRITE(6,*) 'vlocal contribution ', vlocal

    evale=tildekin+tildepot+oneenergy-onehat                        &
&        +0.5d0*onehartree+texc
    WRITE (6,*) ' evale from matrix elements ', evale
    WRITE(6,*) ' evale from AE calculation ',FC%evale
    WRITE(6,*) ' difference in evale results ', evale-FC%evale
    WRITE(6,*) 'tildekin ', tildekin
    WRITE(6,*) 'tildepot ' , tildepot
    WRITE(6,*) 'oneenergy ', oneenergy
    WRITE(6,*) 'onehat ', onehat
    WRITE(6,*) 'onehartree ', onehartree
    WRITE(6,*) 'texc ',texc
    WRITE(6,*) 'vlocal contribution ', vlocal
    !
    call coretailselfenergy(AEGrid,PAW,ctctse,cthatse)
    WRITE(ifatompaw,'("   CORETAILSELFENERGY   ",1pe25.17,"  END")') ctctse
    WRITE(ifatompaw,'("   CORETAILHATENERGY   ",1pe25.17,"  END")') cthatse

    WRITE(ifen,'("evale from matrix elements",1pe25.17)') evale
    WRITE(6,'("evale from matrix elements",1pe25.17)') evale


    WRITE(ifatompaw,'("   ENERGY   ",1pe25.17,"  END")') evale

    CLOSE(ifatompaw)
    CLOSE(ifxml)

    CALL ftprod(AEGrid,PAW)
    CALL SCFPAW(AEGrid,AEPOT%nz,PAW,FC,AEOrbit,PSOrbit,Etotal)
    CALL logderiv(AEGrid,AEPot,PAW)
    CALL fthatpot(AEGrid,PAW)
    CALL ftkin(AEGrid,PAW)
    CALL ftvloc(AEGrid,PAW)
    !
    Do
       WRITE(6,*) 'Enter 0 to end program'
       WRITE(6,*) 'Enter 1 to run SCFPAW'
       WRITE(6,*) 'Enter 2 to run atompaw2abinit'
       WRITE(6,*) 'Enter 3 to run atompaw2pwscf'
       READ(5,*) i
       if (i==0) then
          if(scalarrelativistic) call Deallocate_Scalar_Relativistic
          exit
       elseif (i==1) then
          CALL SCFPAW(AEGrid,AEPOT%nz,PAW,FC,AEOrbit,PSOrbit,Etotal,.true.)
           WRITE(ifen,*) ' '
           WRITE(ifen,*) ' PAW results for new configuration'
           WRITE(6,*) ' '
           WRITE(6,*) ' PAW results for new configuration'
           DO io=1,PSOrbit%norbit
              WRITE(ifen,'(3i5,1p,2e15.7)') io,PSOrbit%np(io),PSOrbit%l(io),&
&               PSOrbit%occ(io),PSOrbit%eig(io)
              WRITE(6,'(3i5,1p,2e15.7)') io,PSOrbit%np(io),PSOrbit%l(io),&
&               PSOrbit%occ(io),PSOrbit%eig(io)
           ENDDO
           WRITE(ifen,'("evale from matrix elements",1pe25.17)') Etotal
           WRITE(6,'("evale from matrix elements",1pe25.17)') Etotal
           !CALL flush(ifen)
       elseif (i==2) then
           CALL atompaw2abinit(AEOrbit,AEPot,PAW,FC,AEGrid)
       elseif (i==3) then
           CALL atompaw2pwscf(AEGrid,AEPot,FC,PAW,ifinput)
       else
          write(6,*) 'Option not recognized -- pgm terminating'
          stop
       endif
    Enddo

!   Free memory
    call DestroyFC(FC)
    call DestroyGrid(AEGrid)
    call DestroyOrbit(AEOrbit)
    call DestroyOrbit(PSOrbit)
    call DestroyPotential(AEPot)
    call DestroyPseudopot(PAW)
    if (have_libxc) call libxc_end()

  END PROGRAM atompaw
