module PARAMETER_MODULE
  use SYSTEM_MODULE
  implicit none
  !RYUDBLE CONSTANTS
  real(MYPRE), parameter :: AMU_RYDBERG          =       2.1798723611035d-18
  !Bohr radius
  real(MYPRE), parameter :: BOHR_RADIUS         =       5.29177210903d-11
  !mass of Electorn 
  !real(MYPRE), parameter :: me         =       9.1093837015d-31
  real(MYPRE), parameter :: ELECTRON_MASS         =       9.1093897d-31
  ![C]
  real(MYPRE), parameter :: ELECTRON_CHARGE =   1.602176634d-19
  !angstrom [SI]
  real(MYPRE), parameter :: ANGSTROM =       1d-10
  !Dolton
  real(MYPRE), parameter :: AMU_DALTON =       1.6605390666050d-27
  !ennsyuuritu
  real(MYPRE),parameter  :: PI         =       4.0d0 * atan(1.0d0)
  !complex integer(MYIP) 
  complex(kind(0d0)),parameter    :: CI         =       (0.0d0,1.0d0)
  !theroshold
  real(MYPRE),parameter  :: EPS =       1d-8
  !
  real(MYPRE),parameter  :: DIELECTRIC_CONSTANTS   =       8.8541878128d-12
  !ボルツマン定数
  real(MYPRE),parameter :: BOLTZMAN_CONSTANTS      =       1.380649d-23
  !プランクの定数
  real(MYPRE),parameter :: PLANCK_NUMBER =     6.62607015d-34 
  
end module PARAMETER_MODULE

