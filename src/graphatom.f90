PROGRAM graphatom  
  !*************************************************************
  !  program to calculate the self-consistent density functional
  !    atom ground state for atom with atomic number nz
  !************************************************************
  USE GlobalMath
  USE atomdata
  USE aeatom
  USE gridmod
  IMPLICIT NONE
  INTEGER, PARAMETER :: ifen=9, ifden=7,ifwfn=8  
  CHARACTER (len=4) :: flnm  
  CHARACTER (len=20) :: nm  
  CHARACTER (len=2) :: sym  
  CHARACTER (len=1) :: syml  
  REAL(8), POINTER :: r(:),den(:),rv(:),wfn(:,:)
  INTEGER, POINTER :: n,norbit,nps,npp,npd,npf,npg

  INTEGER :: i,j,io,many,l,istart
  TYPE (GridInfo), TARGET :: AEGrid
  TYPE (PotentialInfo), TARGET :: AEPot
  TYPE (OrbitInfo), TARGET :: AEOrbit
  TYPE (OrbitInfo), TARGET :: FCOrbit
  TYPE (SCFInfo), TARGET :: AESCF
  TYPE (FCInfo), TARGET :: FC

  CALL Init_GlobalConstants
  CALL iSCFatom(AEGrid,AEPot,AEOrbit,AESCF)
  
  OPEN (unit = ifen, file=TRIM(AEPot%sym)//'.GA', form='formatted')

   i=1
  do
     WRITE(ifen,*) 'Completed calculations for ',TRIM(AEPOT%sym)
     call reportgrid(AEGrid,ifen)
     if (scalarrelativistic) then
         WRITE(ifen,*) 'Scalar relativistic calculation'
     else
         WRITE(ifen,*) 'Non-relativistic calculation'
     endif
     if (i==1)WRITE(ifen,*) '  aeatom converged in',AESCF%iter,' iterations'
     if (i==2)WRITE(ifen,*) '  FCatom converged in',AESCF%iter,' iterations'
     WRITE(ifen,*) '     for nz = ',AEPot%nz
     WRITE(ifen,*) '    delta(density)  = ', AESCF%delta
     WRITE(ifen,*) '  Orbital energies:         '
     WRITE(ifen,*) ' n  l     occupancy       energy'
     If (i==1) then
       DO io=1,AEorbit%norbit
             WRITE(ifen,'(i2,1x,i2,4x,1p,2e15.7)') &
&            AEOrbit%np(io),AEOrbit%l(io),&
&            AEOrbit%occ(io),AEOrbit%eig(io)
       ENDDO
     Else if (i==2) then
       DO io=1,FC%norbit
             IF (.NOT.FC%iscore(io))WRITE(ifen,'(i2,1x,i2,4x,1p,2e15.7)') &
&            AEOrbit%np(io),AEOrbit%l(io),&
&            AEOrbit%occ(io),AEOrbit%eig(io)
       ENDDO
     Endif
     WRITE(ifen,*)
     WRITE(ifen,*) ' Total energy'
     WRITE(ifen,*) '    Total                    :  ',AESCF%etot
     If (i==2) then
       WRITE(ifen,*) '    Valence                  :  ',FC%evale
     EndIf

     write(6,*) ' Input 0  for plotting results and completing program'
     write(6,*) ' Input 1  for changing configuration (all-electron)'
     write(6,*) ' Input 2  for changing configuration (frozen-core)'

     read(5,*) i

     if (i == 0) exit

     if (i == 1) CALL cSCFatom(AEGrid,AEPot,AEOrbit,AESCF)
     if (i == 2) CALL FCSCFatom(AEGrid,AEPot,AEOrbit,FCOrbit,AESCF,FC)
  enddo

  n=>AEGrid%n
  r=>AEGrid%r
  den=>AEPot%den
  rv=>AEPot%rv
  norbit=>AEOrbit%norbit
  nps=>AEOrbit%nps
  npp=>AEOrbit%npp
  npd=>AEOrbit%npd
  npf=>AEOrbit%npf
  npg=>AEOrbit%npg
  wfn=>AEOrbit%wfn

  !
  ! write density and wavefunctions
  !


  CLOSE(ifen)
  OPEN (unit = ifden, file = 'density.GA', form = 'formatted')  
  DO i = 1, n  
     IF (r (i) .LE.6.d0) THEN  
        WRITE (ifden, '(1p,6e12.4)') r (i), den (i)  
     ENDIF
  ENDDO
  CLOSE (ifden)  

  OPEN (unit = ifden, file = 'potential.GA', form = 'formatted')  
  DO i = 1, n  
     IF (r (i) .LE.6.d0) THEN  
        WRITE (ifden, '(1p,6e12.4)') r (i), rv (i)  
     ENDIF
  ENDDO
  CLOSE (ifden)  

  OPEN (unit = ifden, file = 'plotdensity.GA', form = 'formatted')  
  nm = 'density.GA'  
  WRITE (ifden, '("gplot -t ""Radial density for ",a2, &
       &""" -tx ""r (bohr)"" -f ",a10, " 1 2 lines")') AEpot%sym, nm
  CLOSE (ifden)  
  OPEN (unit = ifden, file = 'plotpotential.GA', form = 'formatted')  
  nm = 'potential.GA'  
  WRITE (ifden, '("gplot -t ""rxV(r) for ",a2, &
       &""" -tx ""r (bohr)"" -f ",a12, " 1 2 lines")') AEpot%sym, nm
  CLOSE (ifden)  
  DO i = 1, n  
     DO io = 1, norbit  
        IF (dabs (wfn (i, io) ) .LT.1.d-8) wfn (i, io) = 0.d0  
     ENDDO
  ENDDO
  istart = 0  
  DO l = 0, 4  
     IF (l.EQ.0) many = nps  
     IF (l.EQ.1) many = npp - 1  
     IF (l.EQ.2) many = npd-2  
     IF (l.EQ.3) many = npf - 3  
     IF (l.EQ.4) many = npg - 4  
     IF (l.EQ.0) syml = 's'  
     IF (l.EQ.1) syml = 'p'  
     IF (l.EQ.2) syml = 'd'  
     IF (l.EQ.3) syml = 'f'  
     IF (l.EQ.4) syml = 'g'  
     IF (many.GT.0) THEN  
        CALL mkname (l, flnm)  
        OPEN (unit = ifwfn, file = 'GAwfn'//flnm, form = 'formatted')  
        DO i = 1, n  
           IF (r (i) .LE.6.d0) THEN  
              WRITE (ifwfn, '(1p,8e10.2)') r (i) , (wfn (i, j) , j = &
&                  istart + 1, istart + many)
           ENDIF
        ENDDO
        CLOSE (ifwfn)  
        nm = 'GAwfn'//flnm  
        OPEN (unit = ifwfn, file = 'plotGAwfn'//flnm, form = 'formatted')  
        IF (many.EQ.1) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines")') syml, AEpot%sym, nm
        IF (many.EQ.2) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines")') syml &
             , AEpot%sym, nm, nm
        IF (many.EQ.3) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
             &a8," 1 4 lines")') syml, AEpot%sym, nm, nm, nm
        IF (many.EQ.4) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
             &a8," 1 4 lines -f ",a8," 1 5 lines")') syml, AEpot%sym, nm, nm, nm, nm
        IF (many.EQ.5) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
             &a8," 1 4 lines -f ",a8," 1 5 lines -f ",a8," 1 6 lines")') syml,  &
             &AEpot%sym, nm, nm, nm, nm, nm
        IF (many.EQ.6) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
             &a8," 1 4 lines -f ",a8," 1 5 lines -f ",a8," 1 6 lines -f " &
             &,a8," 1 7 lines")') syml, AEpot%sym, nm, nm, nm, nm, nm, nm
        IF (many.EQ.7) WRITE (ifwfn, '("gplot -t ""Radial ",a1, &
             &   " wavefunctions for ",a2, &
             &""" -tx ""r (bohr)"" -f ",a8, " 1 2 lines -f ",a8," 1 3 lines -f ", &
             &a8," 1 4 lines -f ",a8," 1 5 lines -f ",a8," 1 6 lines -f " &
             &,a8," 1 7 lines -f ",a8," 1 8 lines" )') syml, AEpot%sym, nm, nm, nm, nm &
             &, nm, nm, nm
        CLOSE (ifwfn)  
     ENDIF
     istart = istart + many  

  ENDDO
END PROGRAM graphatom
