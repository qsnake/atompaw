!! NAME
!!  libxc_mod
!!
!! FUNCTION
!!  This module contains routines using features from libXC library
!!  (http://www.tddft.org/programs/octopus/wiki/index.php/Libxc)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module libxc_mod

#if defined HAVE_LIBXC
 use xc_f90_types_m
 use libxc_funcs_m
 use xc_f90_lib_m
 use globalmath
#endif

 implicit none

!!=================================================================
!! CONSTANTS
!!=================================================================

#if defined HAVE_LIBXC
 logical,parameter,public :: have_libxc=.true.
#else
 logical,parameter,public :: have_libxc=.false.
#endif

!!=================================================================
!! STRUCTURED DATATYPES
!!=================================================================

#if defined HAVE_LIBXC
 type libxc_functional
  private
  integer :: family ! LDA, GGA, etc.
  integer :: id     ! identifier
  type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
  type(xc_f90_pointer_t) :: info ! information about the functional
 end type libxc_functional
 type(libxc_functional) :: libxc_funcs(2)
#endif

 private

 public :: libxc_init_func,&
&          libxc_print_func,&
&          libxc_getid_fromname,&
&          libxc_getid,&
&          libxc_getvxc,&
&          libxc_isgga,&
&          libxc_end

 CONTAINS

!!=================================================================
!! NAME
!! libxc_getid_fromname
!!
!! FUNCTION
!! From a character string (given in input file), gives the libXC id(s)
!!
!! PARENTS
!! excor
!!
!!=================================================================
 subroutine libxc_getid_fromname(xcname,id,xcname_short)

 implicit none
 integer,intent(inout) :: id(2)
 character*(*),intent(in) :: xcname
 character*(*),intent(out),optional :: xcname_short

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii,i_plus
 character*50 :: xcstrg(2)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 i_plus=index(xcname,'+')
 if (i_plus<=0) then
  xcstrg(1)=trim(xcname)
  xcstrg(2)=""
 else
  xcstrg(1)=trim(xcname(1:i_plus-1))
  xcstrg(2)=trim(xcname(i_plus+1:))
 end if

 do ii=1,2
  id(ii)=0
  call uppercase(xcstrg(ii))

  if (xcstrg(ii)=="") exit

  if (xcstrg(ii)=="XC_LDA_X")              id(ii)=XC_LDA_X
  if (xcstrg(ii)=="XC_LDA_C_WIGNER")       id(ii)=XC_LDA_C_WIGNER
  if (xcstrg(ii)=="XC_LDA_C_RPA")          id(ii)=XC_LDA_C_RPA
  if (xcstrg(ii)=="XC_LDA_C_HL")           id(ii)=XC_LDA_C_HL
  if (xcstrg(ii)=="XC_LDA_C_GL")           id(ii)=XC_LDA_C_GL
  if (xcstrg(ii)=="XC_LDA_C_XALPHA")       id(ii)=XC_LDA_C_XALPHA
  if (xcstrg(ii)=="XC_LDA_C_VWN")          id(ii)=XC_LDA_C_VWN
  if (xcstrg(ii)=="XC_LDA_C_VWN_RPA")      id(ii)=XC_LDA_C_VWN_RPA
  if (xcstrg(ii)=="XC_LDA_C_PZ")           id(ii)=XC_LDA_C_PZ
  if (xcstrg(ii)=="XC_LDA_C_PZ_MOD")       id(ii)=XC_LDA_C_PZ_MOD
  if (xcstrg(ii)=="XC_LDA_C_OB_PZ")        id(ii)=XC_LDA_C_OB_PZ
  if (xcstrg(ii)=="XC_LDA_C_PW")           id(ii)=XC_LDA_C_PW
  if (xcstrg(ii)=="XC_LDA_C_PW_MOD")       id(ii)=XC_LDA_C_PW_MOD
  if (xcstrg(ii)=="XC_LDA_C_OB_PW")        id(ii)=XC_LDA_C_OB_PW
  if (xcstrg(ii)=="XC_LDA_C_2D_AMGB")      id(ii)=XC_LDA_C_2D_AMGB
  if (xcstrg(ii)=="XC_LDA_C_2D_PRM")       id(ii)=XC_LDA_C_2D_PRM
  if (xcstrg(ii)=="XC_LDA_C_VBH")          id(ii)=XC_LDA_C_vBH
  if (xcstrg(ii)=="XC_LDA_C_1D_CSC")       id(ii)=XC_LDA_C_1D_CSC
  if (xcstrg(ii)=="XC_LDA_X_2D")           id(ii)=XC_LDA_X_2D
  if (xcstrg(ii)=="XC_LDA_XC_TETER93")     id(ii)=XC_LDA_XC_TETER93
  if (xcstrg(ii)=="XC_LDA_X_1D")           id(ii)=XC_LDA_X_1D
  if (xcstrg(ii)=="XC_LDA_C_ML1")          id(ii)=XC_LDA_C_ML1
  if (xcstrg(ii)=="XC_LDA_C_ML2")          id(ii)=XC_LDA_C_ML2
  if (xcstrg(ii)=="XC_GGA_X_PBE")          id(ii)=XC_GGA_X_PBE
  if (xcstrg(ii)=="XC_GGA_X_PBE_R")        id(ii)=XC_GGA_X_PBE_R
  if (xcstrg(ii)=="XC_GGA_X_B86")          id(ii)=XC_GGA_X_B86
  if (xcstrg(ii)=="XC_GGA_X_B86_R")        id(ii)=XC_GGA_X_B86_R
  if (xcstrg(ii)=="XC_GGA_X_B86_MGC")      id(ii)=XC_GGA_X_B86_MGC
  if (xcstrg(ii)=="XC_GGA_X_B88")          id(ii)=XC_GGA_X_B88
  if (xcstrg(ii)=="XC_GGA_X_G96")          id(ii)=XC_GGA_X_G96
  if (xcstrg(ii)=="XC_GGA_X_PW86")         id(ii)=XC_GGA_X_PW86
  if (xcstrg(ii)=="XC_GGA_X_PW91")         id(ii)=XC_GGA_X_PW91
  if (xcstrg(ii)=="XC_GGA_X_OPTX")         id(ii)=XC_GGA_X_OPTX
  if (xcstrg(ii)=="XC_GGA_X_DK87_R1")      id(ii)=XC_GGA_X_DK87_R1
  if (xcstrg(ii)=="XC_GGA_X_DK87_R2")      id(ii)=XC_GGA_X_DK87_R2
  if (xcstrg(ii)=="XC_GGA_X_LG93")         id(ii)=XC_GGA_X_LG93
  if (xcstrg(ii)=="XC_GGA_X_FT97_A")       id(ii)=XC_GGA_X_FT97_A
  if (xcstrg(ii)=="XC_GGA_X_FT97_B")       id(ii)=XC_GGA_X_FT97_B
  if (xcstrg(ii)=="XC_GGA_X_PBE_SOL")      id(ii)=XC_GGA_X_PBE_SOL
  if (xcstrg(ii)=="XC_GGA_X_RPBE")         id(ii)=XC_GGA_X_RPBE
  if (xcstrg(ii)=="XC_GGA_X_WC")           id(ii)=XC_GGA_X_WC
  if (xcstrg(ii)=="XC_GGA_X_MPW91")        id(ii)=XC_GGA_X_mPW91
  if (xcstrg(ii)=="XC_GGA_X_AM05")         id(ii)=XC_GGA_X_AM05
  if (xcstrg(ii)=="XC_GGA_X_PBEA")         id(ii)=XC_GGA_X_PBEA
  if (xcstrg(ii)=="XC_GGA_X_MPBE")         id(ii)=XC_GGA_X_MPBE
  if (xcstrg(ii)=="XC_GGA_X_XPBE")         id(ii)=XC_GGA_X_XPBE
  if (xcstrg(ii)=="XC_GGA_X_2D_B86_MGC")   id(ii)=XC_GGA_X_2D_B86_MGC
  if (xcstrg(ii)=="XC_GGA_X_BAYESIAN")     id(ii)=XC_GGA_X_BAYESIAN
  if (xcstrg(ii)=="XC_GGA_X_PBE_JSJR")     id(ii)=XC_GGA_X_PBE_JSJR
  if (xcstrg(ii)=="XC_GGA_X_2D_B88")       id(ii)=XC_GGA_X_2D_B88
  if (xcstrg(ii)=="XC_GGA_X_2D_B86")       id(ii)=XC_GGA_X_2D_B86
  if (xcstrg(ii)=="XC_GGA_X_2D_PBE")       id(ii)=XC_GGA_X_2D_PBE
  if (xcstrg(ii)=="XC_GGA_C_PBE")          id(ii)=XC_GGA_C_PBE
  if (xcstrg(ii)=="XC_GGA_C_LYP")          id(ii)=XC_GGA_C_LYP
  if (xcstrg(ii)=="XC_GGA_C_P86")          id(ii)=XC_GGA_C_P86
  if (xcstrg(ii)=="XC_GGA_C_PBE_SOL")      id(ii)=XC_GGA_C_PBE_SOL
  if (xcstrg(ii)=="XC_GGA_C_PW91")         id(ii)=XC_GGA_C_PW91
  if (xcstrg(ii)=="XC_GGA_C_AM05")         id(ii)=XC_GGA_C_AM05
  if (xcstrg(ii)=="XC_GGA_C_XPBE")         id(ii)=XC_GGA_C_XPBE
  if (xcstrg(ii)=="XC_GGA_C_LM")           id(ii)=XC_GGA_C_LM
  if (xcstrg(ii)=="XC_GGA_C_PBE_JRGX")     id(ii)=XC_GGA_C_PBE_JRGX
  if (xcstrg(ii)=="XC_GGA_X_OPTB88_VDW")   id(ii)=XC_GGA_X_OPTB88_VDW
  if (xcstrg(ii)=="XC_GGA_X_PBEK1_VDW")    id(ii)=XC_GGA_X_PBEK1_VDW
  if (xcstrg(ii)=="XC_GGA_X_OPTPBE_VDW")   id(ii)=XC_GGA_X_OPTPBE_VDW
  if (xcstrg(ii)=="XC_GGA_X_RGE2")         id(ii)=XC_GGA_X_RGE2
  if (xcstrg(ii)=="XC_GGA_C_RGE2")         id(ii)=XC_GGA_C_RGE2
  if (xcstrg(ii)=="XC_GGA_XC_LB")          id(ii)=XC_GGA_XC_LB
  if (xcstrg(ii)=="XC_GGA_XC_HCTH_93")     id(ii)=XC_GGA_XC_HCTH_93
  if (xcstrg(ii)=="XC_GGA_XC_HCTH_120")    id(ii)=XC_GGA_XC_HCTH_120
  if (xcstrg(ii)=="XC_GGA_XC_HCTH_147")    id(ii)=XC_GGA_XC_HCTH_147
  if (xcstrg(ii)=="XC_GGA_XC_HCTH_407")    id(ii)=XC_GGA_XC_HCTH_407
  if (xcstrg(ii)=="XC_GGA_XC_EDF1")        id(ii)=XC_GGA_XC_EDF1
  if (xcstrg(ii)=="XC_GGA_XC_XLYP")        id(ii)=XC_GGA_XC_XLYP
  if (xcstrg(ii)=="XC_GGA_XC_B97")         id(ii)=XC_GGA_XC_B97
  if (xcstrg(ii)=="XC_GGA_XC_B97_1")       id(ii)=XC_GGA_XC_B97_1
  if (xcstrg(ii)=="XC_GGA_XC_B97_2")       id(ii)=XC_GGA_XC_B97_2
  if (xcstrg(ii)=="XC_GGA_XC_B97_D")       id(ii)=XC_GGA_XC_B97_D
  if (xcstrg(ii)=="XC_GGA_XC_B97_K")       id(ii)=XC_GGA_XC_B97_K
  if (xcstrg(ii)=="XC_GGA_XC_B97_3")       id(ii)=XC_GGA_XC_B97_3
  if (xcstrg(ii)=="XC_GGA_XC_PBE1W")       id(ii)=XC_GGA_XC_PBE1W
  if (xcstrg(ii)=="XC_GGA_XC_MPWLYP1W")    id(ii)=XC_GGA_XC_MPWLYP1W
  if (xcstrg(ii)=="XC_GGA_XC_PBELYP1W")    id(ii)=XC_GGA_XC_PBELYP1W
  if (xcstrg(ii)=="XC_GGA_XC_SB98_1A")     id(ii)=XC_GGA_XC_SB98_1a
  if (xcstrg(ii)=="XC_GGA_XC_SB98_1B")     id(ii)=XC_GGA_XC_SB98_1b
  if (xcstrg(ii)=="XC_GGA_XC_SB98_1C")     id(ii)=XC_GGA_XC_SB98_1c
  if (xcstrg(ii)=="XC_GGA_XC_SB98_2A")     id(ii)=XC_GGA_XC_SB98_2a
  if (xcstrg(ii)=="XC_GGA_XC_SB98_2B")     id(ii)=XC_GGA_XC_SB98_2b
  if (xcstrg(ii)=="XC_GGA_XC_SB98_2C")     id(ii)=XC_GGA_XC_SB98_2c
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B3PW91")  id(ii)=XC_HYB_GGA_XC_B3PW91
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B3LYP")   id(ii)=XC_HYB_GGA_XC_B3LYP
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B3P86")   id(ii)=XC_HYB_GGA_XC_B3P86
  if (xcstrg(ii)=="XC_HYB_GGA_XC_O3LYP")   id(ii)=XC_HYB_GGA_XC_O3LYP
  if (xcstrg(ii)=="XC_HYB_GGA_XC_mPW1K")   id(ii)=XC_HYB_GGA_XC_mPW1K
  if (xcstrg(ii)=="XC_HYB_GGA_XC_PBEH")    id(ii)=XC_HYB_GGA_XC_PBEH
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B97")     id(ii)=XC_HYB_GGA_XC_B97
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B97_1")   id(ii)=XC_HYB_GGA_XC_B97_1
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B97_2")   id(ii)=XC_HYB_GGA_XC_B97_2
  if (xcstrg(ii)=="XC_HYB_GGA_XC_X3LYP")   id(ii)=XC_HYB_GGA_XC_X3LYP
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B1WC")    id(ii)=XC_HYB_GGA_XC_B1WC
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B97_K")   id(ii)=XC_HYB_GGA_XC_B97_K
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B97_3")   id(ii)=XC_HYB_GGA_XC_B97_3
  if (xcstrg(ii)=="XC_HYB_GGA_XC_MPW3PW")  id(ii)=XC_HYB_GGA_XC_mPW3PW
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B1LYP")   id(ii)=XC_HYB_GGA_XC_B1LYP
  if (xcstrg(ii)=="XC_HYB_GGA_XC_B1PW91")  id(ii)=XC_HYB_GGA_XC_B1PW91
  if (xcstrg(ii)=="XC_HYB_GGA_XC_MPW1PW")  id(ii)=XC_HYB_GGA_XC_mPW1PW
  if (xcstrg(ii)=="XC_HYB_GGA_XC_MPW3LYP") id(ii)=XC_HYB_GGA_XC_mPW3LYP
  if (xcstrg(ii)=="XC_HYB_GGA_XC_SB98_1A") id(ii)=XC_HYB_GGA_XC_SB98_1a
  if (xcstrg(ii)=="XC_HYB_GGA_XC_SB98_1B") id(ii)=XC_HYB_GGA_XC_SB98_1b
  if (xcstrg(ii)=="XC_HYB_GGA_XC_SB98_1C") id(ii)=XC_HYB_GGA_XC_SB98_1c
  if (xcstrg(ii)=="XC_HYB_GGA_XC_SB98_2A") id(ii)=XC_HYB_GGA_XC_SB98_2a
  if (xcstrg(ii)=="XC_HYB_GGA_XC_SB98_2B") id(ii)=XC_HYB_GGA_XC_SB98_2b
  if (xcstrg(ii)=="XC_HYB_GGA_XC_SB98_2C") id(ii)=XC_HYB_GGA_XC_SB98_2c
  if (xcstrg(ii)=="XC_MGGA_X_LTA")         id(ii)=XC_MGGA_X_LTA
  if (xcstrg(ii)=="XC_MGGA_X_TPSS")        id(ii)=XC_MGGA_X_TPSS
  if (xcstrg(ii)=="XC_MGGA_X_M06L")        id(ii)=XC_MGGA_X_M06L
  if (xcstrg(ii)=="XC_MGGA_X_GVT4")        id(ii)=XC_MGGA_X_GVT4
  if (xcstrg(ii)=="XC_MGGA_X_TAU_HCTH")    id(ii)=XC_MGGA_X_TAU_HCTH
  if (xcstrg(ii)=="XC_MGGA_X_BR89")        id(ii)=XC_MGGA_X_BR89
  if (xcstrg(ii)=="XC_MGGA_X_BJ06")        id(ii)=XC_MGGA_X_BJ06
  if (xcstrg(ii)=="XC_MGGA_X_TB09")        id(ii)=XC_MGGA_X_TB09
  if (xcstrg(ii)=="XC_MGGA_X_RPP09")       id(ii)=XC_MGGA_X_RPP09
  if (xcstrg(ii)=="XC_MGGA_C_TPSS")        id(ii)=XC_MGGA_C_TPSS
  if (xcstrg(ii)=="XC_MGGA_C_VSXC")        id(ii)=XC_MGGA_C_VSXC
  if (xcstrg(ii)=="XC_LCA_OMC")            id(ii)=XC_LCA_OMC
  if (xcstrg(ii)=="XC_LCA_LCH")            id(ii)=XC_LCA_LCH

  if (id(ii)==0.and.xcstrg(ii)(1:6)=="LIBXC_") then
   read(unit=xcstrg(ii)(7:),fmt=*,err=333,end=333) id(ii)
333 continue
  end if

  if (id(ii)==0) then
   id(ii)=-1
   write(6,'(/,2x,a)') "Error in get_id_from_name:"
   write(6,'(2x,3a)')  " Unknown X, C or XC functionnal (", &
&                     trim(xcstrg(ii)),") !"
   stop
  end if
 end do

 if (present(xcname_short)) then
  xcname_short=""
  if (id(1)>0.and.xcstrg(1)(1:3)=="XC_") xcname_short=trim(xcname_short)     //trim(xcstrg(1)(4:))
  if (id(2)>0.and.xcstrg(2)(1:3)=="XC_") xcname_short=trim(xcname_short)//"+"//trim(xcstrg(2)(4:))
 end if

#else
 id(1:2)=-2
#endif

 end subroutine libxc_getid_fromname


!!=================================================================
!! NAME
!! libxc_getid
!!
!! FUNCTION
!! From LibXC datastructure, gives the libXC id(s)
!!
!! PARENTS
!! rdpawps1
!!
!!=================================================================
 subroutine libxc_getid(id)

 implicit none
 integer :: id(2)

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 do ii=1,2
  id(ii)=libxc_funcs(ii)%id
 end do

#else
 id(1:2)=-2
#endif

 end subroutine libxc_getid


!!=================================================================
!! NAME
!! libxc_init_func
!!
!! FUNCTION
!! Initialize libXC functional(s)
!!
!! PARENTS
!! initexch
!!
!!=================================================================
 subroutine libxc_init_func(id,nsp)

 implicit none
 integer,intent(in)  :: id(2)
 integer, intent(in) :: nsp

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 libxc_funcs(1)%id=id(1)
 libxc_funcs(2)%id=id(2)

 do ii=1,2

  if (libxc_funcs(ii)%id==0) then
   libxc_funcs(ii)%family = 0
   cycle
  end if

! Get XC functional family
  libxc_funcs(ii)%family = xc_f90_family_from_id(libxc_funcs(ii)%id)
  select case (libxc_funcs(ii)%family)
   case (XC_FAMILY_LDA, XC_FAMILY_GGA)
    call xc_f90_func_init(libxc_funcs(ii)%conf,libxc_funcs(ii)%info,libxc_funcs(ii)%id,nsp)
   case default
    write(6,'(4a,i3,4a)' ) char(10),&
&    ' Error in libxc_functionals_init:',char(10),&
&    '  The LibXC functional family ',libxc_funcs(ii)%family,char(10),&
&    '  is currently unsupported by AtomPAW !',char(10),&
&    '  (at present only LGA or GGA are supported)'
    stop
  end select

  if (libxc_funcs(ii)%id == XC_LDA_C_XALPHA) then
   call xc_f90_lda_c_xalpha_set_par(libxc_funcs(ii)%conf,0.d0)
  end if

 end do ! loop on functionals

!Print functional(s) information
 call libxc_print_func(6)

#endif
 end subroutine libxc_init_func


!!=================================================================
!! NAME
!! libxc_end
!!
!! FUNCTION
!! Free libXC functional(s)
!!
!! PARENTS
!!
!!=================================================================
 subroutine libxc_end()

  implicit none

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------
 do ii=1,2
  if (libxc_funcs(ii)%id==0) cycle
  call xc_f90_func_end(libxc_funcs(ii)%conf)
 end do

#endif
 end subroutine libxc_end


!!=================================================================
!! NAME
!! libxc_print_func
!!
!! FUNCTION
!! Print libXC functionnal(s) details
!!
!! PARENTS
!! atompaw,libxc_init_func
!!
!!=================================================================
 subroutine libxc_print_func(unt)

 implicit none
 integer :: unt

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii,jj
 character(len=500) :: msg
 type(xc_f90_pointer_t) :: libxc_str

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 do ii=1,2

  select case (xc_f90_info_kind(libxc_funcs(ii)%info))
   case (XC_EXCHANGE)
    write(unt,'(a)') 'Exchange functional (LibXC):'
   case (XC_CORRELATION)
    write(unt,'(a)') 'Correlation functional (LibXC):'
   case (XC_EXCHANGE_CORRELATION)
    write(unt,'(a)') 'Exchange-Correlation functional (LibXC):'
  end select

  call xc_f90_info_name(libxc_funcs(ii)%info,msg)
  write(unt,'(2x,a)') trim(msg)
  jj=0
  call xc_f90_info_refs(libxc_funcs(ii)%info,jj,libxc_str,msg)
  do while (jj>=0)
   write(unt,'(2x,a)') trim(msg)
   call xc_f90_info_refs(libxc_funcs(ii)%info,jj,libxc_str,msg)
  end do
 end do

#endif
 end subroutine libxc_print_func


!!=================================================================
!! NAME
!! libxc_isgga
!!
!! FUNCTION
!! Returns TRUE is LibXC functional is GGA
!!
!! PARENTS
!! exch
!!
!!=================================================================
 function libxc_isgga()

 implicit none
 logical :: libxc_isgga

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------
 libxc_isgga = .false.
#if defined HAVE_LIBXC
 libxc_isgga = (any(libxc_funcs(:)%family == XC_FAMILY_GGA))
#endif

 end function libxc_isgga


!!=================================================================
!! NAME
!! libxc_getvxc
!!
!! FUNCTION
!! Returns TRUE is LibXC functional is GGA
!!
!! NOTES
!!  nsp=1 : rho is total density (not half)
!!          grho is abs(grad(rho))
!!  nsp=2 : rho is [rho^up,rho^dn]
!!          grho is [abs(grad(rho^up)),abs(grad(rho^dn)),abs(grad(rho^tot))]
!! PARENTS
!! exch
!!
!!=================================================================
 subroutine libxc_getvxc(npts,exc,vxc,nsp,rho,grho,vxcgr)

 implicit none
 integer, intent(in)            :: npts,nsp
 real(8),intent(in)             :: rho(npts,nsp)
 real(8),intent(inout)          :: exc(npts),vxc(npts,nsp)
 real(8),intent(in),optional    :: grho(npts,2*nsp-1)
 real(8),intent(inout),optional :: vxcgr(npts,2*nsp-1)

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 real(8),parameter :: tol=1.d-14
 integer :: ii,ipts,izero
 real(8) :: rhotmp(nsp),exctmp
 real(8) :: sigma(3),vsigma(3),vxctmp(nsp)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (libxc_isgga().and.((.not.present(vxcgr).or.(.not.present(grho))))) then
   write(6,'(/,2x,a)') "Bug in libxc_getvxc:"
   write(6,'(2x,3a)')  " GGA called without grho or vxcgr !"
   stop
 end if

!Initializations
 vxc=0.d0;exc=0.d0
 if (libxc_isgga()) vxcgr=0.d0

!Filter density/gradient when density goes to zero
 izero=0
 do ipts=1,npts
  if (rho(ipts,1)>tol) izero=ipts
 end do
 if (nsp==2.and.izero<npts) then
  do ipts=izero+1,npts
   if (rho(ipts,2)>tol) izero=ipts
  end do
 end if

!Loop over points
 do ipts=1,npts

  vxctmp=0.d0;exctmp=0.d0
  if (ipts<=izero) then
   rhotmp(1:nsp)=rho(ipts,1:nsp)
  else
   rhotmp=tol
  end if

  if (libxc_isgga()) then
   if (ipts<=izero) then
    if (nsp==1) then
    !AtomPAW passes |grho| while LibXC needs |grho|^2
     sigma(1)=grho(ipts,1)**2
    else
    !AtomPAW passes |grho_up|, |grho_dn|, and |grho_tot|
    !while Libxc needs |grho_up|^2, grho_up.grho_dn, and |grho_dn|^2
     sigma(1)= grho(ipts,1)**2
     sigma(3)= grho(ipts,2)**2
     sigma(2)=(grho(ipts,3)**2-sigma(1)-sigma(3))*0.5d0
    end if
   else
    sigma=0.d0
   end if
  end if

! Loop over functionals
  do ii=1,2
   if (libxc_funcs(ii)%id==0) cycle

!  Get the potential (and possibly the energy)
   if (iand(xc_f90_info_flags(libxc_funcs(ii)%info),XC_FLAGS_HAVE_EXC)/=0) then
    select case (libxc_funcs(ii)%family)
     case (XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),exctmp,vxctmp(1))
     case (XC_FAMILY_GGA)
      call xc_f90_gga_exc_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),sigma(1),&
&                             exctmp,vxctmp(1),vsigma(1))
    end select
   else
    exctmp=0.d0
    select case (libxc_funcs(ii)%family)
     case (XC_FAMILY_LDA)
      call xc_f90_lda_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),vxctmp(1))
     case (XC_FAMILY_GGA)
      call xc_f90_gga_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))
    end select
   end if

   exc(ipts)      =exc(ipts)               +2.d0*exctmp           ! From Ha to Ry
   vxc(ipts,1:nsp)=vxc(ipts,1:nsp)         +2.d0*vxctmp(1:nsp)    ! From Ha to Ry
   if (libxc_isgga()) then
    if (nsp==1) then
     vxcgr(ipts,2*nsp-1)=vxcgr(ipts,2*nsp-1)+4.d0*vsigma(2*nsp-1) ! From Ha to Ry
                         ! Note: for nsp=1, vsigma(1) contains 1/2 (up part only)
    else
     vxcgr(ipts,2*nsp-1)=vxcgr(ipts,2*nsp-1)+2.d0*vsigma(2*nsp-1) ! From Ha to Ry
    end if
!   if (nsp==1) then
!    vxcgr(ipts,3)=vxcgr(ipts,3)+vsigma(1)*2.d0
!   else
!    vxcgr(ipts,1)=vxcgr(ipts,1)+2.d0*vsigma(1)-vsigma(2)
!    vxcgr(ipts,2)=vxcgr(ipts,2)+2.d0*vsigma(3)-vsigma(2)
!    vxcgr(ipts,3)=vxcgr(ipts,3)+     vsigma(2)
!   end if
   end if

  end do ! loop over functional(s)
 end do  ! loop over points

#endif
 end subroutine libxc_getvxc

end  Module libxc_mod
