module MDALGORITHM_MODULE
  use MATH_MODULE
  use KAHAN_MODULE
  use ENSEMBLE_MODULE
  implicit none
  !
  contains 
  !
  !============================================================!
  !            IMPLEMNT ALLOCATE GRAND ENSEMBLE                !
  !============================================================!
  !
  subroutine Preprocess_GrandEnsemble(sys)
    type(GrandEnsemble) :: sys
    !
    return 
  end subroutine Preprocess_GrandEnsemble
  !
  !
  !
  subroutine allocate_Ensemble(ensem,freenumber)
    type(Ensemble) :: ensem
    integer,intent(in) :: freenumber
    type(Vector_kahan) :: kaha
    integer :: i
    !
    call allocate_Vector_kahan(kaha,freenumber)
    !
    do i = lbound(ensem%Displacement,1),1
      ensem%Displacement(i) = kaha
    end do 
    do i = lbound(ensem%Velocity,1),1
      ensem%Velocity(i) = kaha
    end do
    do i = lbound(ensem%Force,1),1
      ensem%Force(i) = kaha
    end do 
    call allocate_Vector_d(ensem%Mass,freenumber)
    !
    return 
  end subroutine
  !
  !initialize ensemble
  !
  subroutine Set_GrandEnsemble(sys,Lx,Ly,Lz,freenumber)
    type(GrandEnsemble) :: sys
    integer,intent(in) :: Lx,Ly,Lz
    integer,intent(in) :: freenumber
    !
    integer :: ide
    !
    ide = sys%Ensemble_list(Lx,Ly,Lz)
    !
    call allocate_Ensemble(sys%System(ide),freenumber)
    !
    return
  end subroutine Set_GrandEnsemble
  !
  ! Set Mass Vector
  !
  subroutine Set_GrandEnsemble_Mass(sys,mass,Lx,Ly,Lz,idx)
    type(GrandEnsemble) :: sys
    integer,intent(in) :: Lx,Ly,Lz,idx
    double precision,intent(in) :: mass
    integer :: ide
    !
    ide = sys%Ensemble_list(Lx,Ly,Lz)
    !
    call Set_Ensemble_Mass(sys%System(ide),mass,idx)
    !
    return 
  end subroutine Set_GrandEnsemble_Mass
  !
  subroutine Set_Ensemble_Mass(ensem,mass,idx)
    type(Ensemble) :: ensem
    integer,intent(in) :: idx
    double precision,intent(in) ::  mass
    !
    ensem%Mass%vector(idx) = mass
    !
    return 
  end subroutine Set_Ensemble_Mass
  !
  !=========================================================================!
  !              IMPLEMENT RENEW ROUINTE                                    !
  !=========================================================================!
  !
  subroutine Update_Displacement(sys,ide,idx)
    type(GrandEnsemble) :: sys
    integer,intent(in) :: idx,ide
    double precision :: tmp,mass
    !
    !
    mass = sys%System(ide)%Mass%vector(idx)
    !
    tmp = sys%step * sys%System(ide)%Velocity(1)%vector(idx)%sum
    call add_kahan(sys%System(ide)%Displacement(1)%vector(idx), tmp)
    !
    tmp = 0.5d0 * sys%step * sys%step * sys%System(ide)%Force(0)%vector(idx)%sum / mass
    call add_kahan(sys%System(ide)%Displacement(1)%vector(idx), tmp)
    !
    return
  end subroutine  Update_Displacement
  !
  !Update Velocity
  !
  subroutine Update_Velocity(sys,ide,idx)
    type(GrandEnsemble) :: sys
    integer,intent(in) :: idx,ide
    double precision :: mass
    double precision :: tmp
    !
    mass = sys%System(ide)%Mass%vector(idx)
    !
    tmp = sys%step * sys%System(ide)%Force(0)%vector(idx)%sum * 0.5d0 / mass
    call add_kahan(sys%System(ide)%Velocity(1)%vector(idx), tmp)
    !
    tmp = sys%step * sys%System(ide)%Force(1)%vector(idx)%sum * 0.5d0 / mass
    call add_kahan(sys%System(ide)%Velocity(1)%vector(idx), tmp)
    !
    return 
  end subroutine Update_Velocity
  !
  !swap Force
  !
  subroutine Update_Force(sys)
    type(GrandEnsemble) :: sys
    integer :: it
    integer :: ide
    
    
    do ide = 1,sys%Ensemble_Number
      do it = lbound(sys%System(ide)%Force,1),ubound(sys%System(ide)%Force,1)-1
        sys%System(ide)%Force(it)%vector(:) = sys%System(ide)%Force(it+1)%vector(:)
      end do 
      call Vector_kahan_Zero(sys%System(ide)%Force(1))
    end do 
    !
    !
    return 
  end subroutine  Update_Force
  !
  !
  !
! double precision function Calculate_KineticEnergy_GrandEnsemble(sys)
!   type(GrandEnsemble),intent(in) :: sys
!   integer :: il,jl,kl
!   type(KahanSum) :: keg
!   !
!   call kahan_zero(keg)
!   !
!   do kl = lbound(sys%System,3),ubound(sys%System,3)
!   do jl = lbound(sys%System,2),ubound(sys%System,2)
!   do il = lbound(sys%System,1),ubound(sys%System,1)
!     call Calculate_KineticEnergy_Ensemble(sys%System(il,jl,kl),keg)
!   end do
!   end do
!   end do
!   Calculate_KineticEnergy_GrandEnsemble = keg%sum
!   !
!   return 
! end function Calculate_KineticEnergy_GrandEnsemble
! !
! !
! subroutine Calculate_KineticEnergy_Ensemble(ensem,keg)
!   type(Ensemble),intent(in) :: ensem
!   type(KahanSum) :: keg
!   double precision :: tmp
!   integer  :: i
!   !
!   do i = lbound(ensem%Velocity(1)%vector,1),ubound(ensem%Velocity(1)%vector,1)
!     tmp = ensem%Velocity(1)%vector(i)%sum
!     tmp = 0.5d0 * ensem%Mass%vector(i) * tmp * tmp
!     call add_kahan(keg,tmp)
!   end do 
!   !
!   return
! end subroutine Calculate_KineticEnergy_Ensemble
! ! 
! ! 
! ! 
! double precision function Calculate_PotentialEnergy_GrandEnsemble(sys)
!   type(GrandEnsemble),intent(in) :: sys
!   integer :: il,jl,kl
!   type(KahanSum) :: peg
!   !
!   call kahan_zero(peg)
!   !
!   do kl = lbound(sys%System,3),ubound(sys%System,3)
!   do jl = lbound(sys%System,2),ubound(sys%System,2)
!   do il = lbound(sys%System,1),ubound(sys%System,1)
!     call Calculate_PotentialEnergy_Ensemble(sys%System(il,jl,kl),peg)
!   end do
!   end do
!   end do
!   Calculate_PotentialEnergy_GrandEnsemble = peg%sum
!   !
!   return 
! end function Calculate_PotentialEnergy_GrandEnsemble
! ! 
! ! 
! ! 
! subroutine Calculate_PotentialEnergy_Ensemble(ensem,peg)
!   type(Ensemble),intent(in) :: ensem
!   type(KahanSum) :: peg
!   double precision :: tmp
!   integer  :: i
!   !
!   tmp = 0d0
!   !
!   do i = lbound(ensem%Displacement(1)%vector,1),ubound(ensem%Displacement(1)%vector,1)
!     tmp = -ensem%Displacement(1)%vector(i)%sum * ensem%Force(1)%vector(i)%sum * 0.5d0 
!     call add_kahan(peg,tmp)
!   end do 
!   !
!   return
! end subroutine Calculate_PotentialEnergy_Ensemble
! !
  !
  !
end module MDALGORITHM_MODULE


