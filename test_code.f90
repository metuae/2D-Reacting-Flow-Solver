!
! Programe for computing the chemical time scales.
!
module chemistry_size
  implicit none
  integer(kind=4), PARAMETER  :: nspec_max   = 7
  real(kind=8), dimension(nspec_max) :: X

end module chemistry_size
!
module thermo
   use chemistry_size
   implicit none
   integer(kind=4)                    :: nspec, nspec_i
   real(kind=8), dimension(nspec_max) :: w, DH0, HS_T0, rgaz, Un_sur_MilleW, MilleW,Un_sur_W,   &
                                         CPL0, CPL1, CPL2, CPL3, CPL4, CPL5, CPL6,              &
                                         CPH0, CPH1, CPH2, CPH3, CPH4, CPH5, CPH6,              &
                                         CPL0_prim, CPL1_prim, CPL2_prim, CPL3_prim, CPL4_prim, &
                                         CPH0_prim, CPH1_prim, CPH2_prim, CPH3_prim, CPH4_prim, &
                                         cp_t0
    real(kind=8)                      :: Wo1,CPMA
    character, DIMENSION(nspec_max)   :: spec*4
    real(kind=8), SAVE :: runiv=8314.51D00  ! R in J/(kmol.K)
    real(kind=8)       :: t0                ! reference temperature for the enthalpy computation
!

end module thermo
!
subroutine read_data()
!
use chemistry_size
use thermo
!
implicit none
!
integer(kind=4)     :: e
real(kind =8)       :: dummy
!
open(unit=457, file='data_chem.dat', status='old', form='formatted')
!--------------- Read thermodynamic data for the species ------
!    units
!    molecular weight                : kg/kmol or g/mole
!    enthalpy of formation           : kJ/mol               
read(457,*) 
read(457,*) nspec
read(457,*) nspec_i
read(457,*) t0
 !
 do e = 1, nspec
    read(457,'(A)', advance='no') spec(e)
    read(457,9,advance='no')
    read(457,*) w(e)
    write(6,*) w(e)
    !
    read(457,10) CPH0(e), CPH1(e), CPH2(e), CPH3(e), CPH4(e)
    read(457,10) CPH5(e), CPH6(e), CPL0(e), CPL1(e), CPL2(e)
    read(457,11) CPL3(e), CPL4(e), CPL5(e), CPL6(e)
    read(457,*)  dummy
    !
 enddo
 !
9   format(61x)
10  format(5(E15.8))
11  format(4(E15.8))

 !
close(457) 
end subroutine read_data
!
program begin_prog
!
use thermo
use chemistry_size
implicit none
!
integer(kind=4)     :: e
!
!'H2:0.029,O2:0.209,N2:0.78'
x (1)=0.029
x (2)=0.209
x (3)=0.78
call read_data()
write(6,*) "Enter Reference Temperature"
read(*,*) t0
cpma=0.d0
do e = 1, nspec
   if (t0 .le. 1000.d0) then 
    cp_t0(e)=cpl5(e)*t0
    cp_t0(e)=(cpl4(e)+cp_t0(e))*t0
    cp_t0(e)=(cpl3(e)+cp_t0(e))*t0
    cp_t0(e)=(cpl2(e)+cp_t0(e))*t0
    cp_t0(e)=(cpl1(e)+cp_t0(e))*runiv
   else 
    cp_t0(e)=cph5(e)*t0
    cp_t0(e)=(cph4(e)+cp_t0(e))*t0
    cp_t0(e)=(cph3(e)+cp_t0(e))*t0
    cp_t0(e)=(cph2(e)+cp_t0(e))*t0
    cp_t0(e)=(cph1(e)+cp_t0(e))*runiv
   endif 
   !
   cpma=cpma+(cp_t0(e)*x(e))
enddo
write(6,*) cpma
!26.12.2016 METUAE
!28.12.2016 mk : You were computing enthalp ?
!28.12.2016 METUAE	:	Enthalpy Calculation
h_ma=0.d0
do e = 1, nspec
   if (t0 .le. 1000.d0) then 
    h_t0(e)=(cpl5(e)/5)*t0
    h_t0(e)=((cpl4(e)/4)+cp_t0(e))*t0
    h_t0(e)=((cpl3(e)/3)+cp_t0(e))*t0
    h_t0(e)=((cpl2(e)/2)+cp_t0(e))*t0
    h_t0(e)=((cpl1(e)+cp_t0(e))*runiv*t0
   else 
    h_t0(e)=(cph5(e)/5)*t0
    h_t0(e)=((cph4(e)/4)+cp_t0(e))*t0
    h_t0(e)=((cph3(e)/3)+cp_t0(e))*t0
    h_t0(e)=((cph2(e)/2)+cp_t0(e))*t0
    h_t0(e)=(cph1(e)+cp_t0(e))*runiv*t0
   endif 
   !
   h_ma=h_ma+(h_t0(e)*x(e))
enddo
end program begin_prog

