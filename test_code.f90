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
                                         CPH0_prim, CPH1_prim, CPH2_prim, CPH3_prim, CPH4_prim
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
X (1)=0.029
X (2)=0.209
X (3)=0.78
call read_data()
write(6,*) "Enter Reference Temperature"
read(*,*) t0
CPMA=0
do e = 1, nspec
HS_T0(e)=runiv*(CPL0(e)+CPL1(e)*t0+CPL2(e)*t0**2+CPL3(e)*t0**3+CPL4(e)*t0**4+CPL5(e)*t0**5+CPL6(e)*t0**6)
!
CPMA=CPMA+(HS_T0(e)*X(e))
enddo
write(6,*) CPMA
!26.12.2016 METUAE
end program begin_prog

