module MDALGORITHM_MODULE
  use SYSTEM_MODULE
  use VECTOR_MODULE
  use KAHAN_MODULE
  use ENSEMBLE_MODULE
  !
  implicit none
  !
  !
  contains 
  !
  !
  !============================================================!
  !            IMPLEMENT ALLOCATE GRAND ENSEMBLE                !
  !============================================================!
  !
  !
  subroutine Preprocess_GrandEnsemble(sys)
    type(GrandEnsemble) :: sys
    integer(MYIP) :: ide
    integer(MYIP) :: ensem_num
    !
    !
    print *, "THIS MODULE IS IMPLREMENTED BASED ON BEEMAN ALGORITHM "
    ensem_num = sys%Ensemble_Number
    do ide = 1,ensem_num
      sys%System(ide)%state%Force(0)%vector(:) = sys%System(ide)%state%Force(1)%vector(:)
    end do 
    !
    return 
  end subroutine Preprocess_GrandEnsemble
  !
  !============================================================!
  !
  subroutine allocate_Ensemble(ensem,freenumber)
    type(Ensemble) :: ensem
    integer(MYIP),intent(in) :: freenumber
    type(Vector_kahan) :: kaha
    integer(MYIP) :: i
    !
    call allocate_Vector_kahan(kaha,1,freenumber)
    !
    do i = lbound(ensem%state%Displacement,1),1
      ensem%state%Displacement(i) = kaha
    end do 
    do i = lbound(ensem%state%Velocity,1),1
      ensem%state%Velocity(i) = kaha
    end do
    do i = lbound(ensem%state%Force,1),1
      ensem%state%Force(i) = kaha
    end do 
    !
    call allocate_Vector_d(ensem%info%Mass,1,freenumber)
    !
    return 
  end subroutine allocate_Ensemble
  !
  !
  !initialize ensemble
  !
  subroutine Set_GrandEnsemble(sys,Lx,Ly,Lz,freenumber)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: freenumber
    integer(MYIP),intent(in) :: Lx,Ly,Lz
    integer(MYIP) :: ide
    !
    !
    sys%Ensemble_Number = sys%Ensemble_Number + 1
    ide = sys%Ensemble_list(Lx,Ly,Lz)
    !
    call allocate_Ensemble(sys%System(ide),freenumber)
    !
    return
  end subroutine Set_GrandEnsemble
  !
  subroutine Occupy_GrandEnsemble(sys,Nx,Ny,Nz,Mx,My,Mz)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: Nx,Ny,Nz
    integer(MYIP),intent(in) :: Mx,My,Mz
    !
    sys%Ensemble_list(Mx,My,Mz) = sys%Ensemble_list(Nx,Ny,Nz)
    !
    return 
  end subroutine Occupy_GrandEnsemble 
  !
  subroutine Set_GrandEnsemble_Mass(sys,mass,Lx,Ly,Lz,idx)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: Lx,Ly,Lz,idx
    integer(MYIP) :: ide
    real(MYPRE),intent(in) :: mass
    !
    ide = sys%Ensemble_list(Lx,Ly,Lz)
    !
    call Set_Ensemble_Mass(sys%System(ide),mass,idx)
    !
    return 
  end subroutine 
  !
  subroutine Set_Ensemble_Mass(ensem,mass,idx)
    type(Ensemble) :: ensem
    integer(MYIP),intent(in) :: idx
    real(MYPRE),intent(in) ::  mass
    !
    ensem%info%Mass%vector(idx) = mass
    !
    return 
  end subroutine 
  !
  !=========================================================================!
  !              IMPLEMENT RENEW ROUINTE                                    !
  !=========================================================================!
  !
  subroutine Update_Dislacement(sys,ide,idx)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: idx,ide
    real(MYPRE) :: mass
    real(MYPRE) :: tmp
    !
    mass = sys%System(ide)%info%Mass%vector(idx)
    !
    tmp = sys%step * sys%System(ide)%state%Velocity(1)%vector(idx)%sum
    call add_kahan(sys%System(ide)%state%Displacement(1)%vector(idx), tmp)
    !
    tmp = sys%step * sys%step * 2.d0 * sys%System(ide)%state%Force(0)%vector(idx)%sum / 3.0d0  / mass
    call add_kahan(sys%System(ide)%state%Displacement(1)%vector(idx), tmp)
    !
    tmp = -sys%step * sys%step * sys%System(ide)%state%Force(-1)%vector(idx)%sum / 6.0d0  / mass
    call add_kahan(sys%System(ide)%state%Displacement(1)%vector(idx), tmp)
    !
    return

  end subroutine Update_Dislacement
  !
  !Renew state%Velocity
  !
  subroutine Update_Velocity(sys,ide,idx)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: ide,idx
    real(MYPRE) :: mass,tmp
    !
    mass = sys%System(ide)%info%Mass%vector(idx)
    !
    tmp = sys%step * sys%System(ide)%state%Force(1)%vector(idx)%sum / 3.0d0 / mass
    call add_kahan(sys%System(ide)%state%Velocity(1)%vector(idx), tmp)
    !
    tmp = 5.0d0 * sys%step * sys%System(ide)%state%Force(0)%vector(idx)%sum / 6.0d0 / mass
    call add_kahan(sys%System(ide)%state%Velocity(1)%vector(idx), tmp)
    !
    tmp = -sys%step * sys%System(ide)%state%Force(-1)%vector(idx)%sum / 6.0d0 / mass
    call add_kahan(sys%System(ide)%state%Velocity(1)%vector(idx), tmp)
    !
    return
  end subroutine Update_Velocity
  !
  !swap state%Force
  !
  subroutine Swap_Force(sys)
    type(GrandEnsemble) :: sys
    integer(MYIP) :: it
    integer(MYIP) :: ide
    !
    !
    do ide = 1,sys%Ensemble_Number
      do it = lbound(sys%System(ide)%state%Force,1),ubound(sys%System(ide)%state%Force,1)-1
        sys%System(ide)%state%Force(it)%vector(:) = sys%System(ide)%state%Force(it+1)%vector(:)
      end do
      call Vector_kahan_zero(sys%System(ide)%state%Force(it))
    end do 
    !do kl = lbound(sys%Ensemble_list,3),ubound(sys%Ensemble_list,3)
    !do jl = lbound(sys%Ensemble_list,2),ubound(sys%Ensemble_list,2)
    !do il = lbound(sys%Ensemble_list,1),ubound(sys%Ensemble_list,1)
    !  ide = sys%Ensemble_list(il,jl,kl)
    !  do it = lbound(sys%System(ide)%state%Force,1),ubound(sys%System(ide)%state%Force,1)-1
    !    sys%System(ide)%state%Force(it)%vector(:) = sys%System(ide)%state%Force(it+1)%vector(:)
    !  end do 
    !end do 
    !end do
    !end do
    !
    !
    return 
  end subroutine  Swap_Force
  !
  !===================================================================================!
  !     WRITE AND READ MODULE                                                         !
  !===================================================================================!
  !
  subroutine Write_Ensemble_state_Beeman(ensem,fd)
    type(Ensemble),intent(in) :: ensem
    integer(MYIP) , intent(in) :: fd
    !
    integer(MYIP) :: i 
    integer(MYIP) :: freenumber
    !
    freenumber = size(ensem%info%Mass%vector)
    write(fd,*)(ensem%state%Displacement(1)%vector(i)%sum,i=1,freenumber)
    write(fd,*)(ensem%state%Displacement(0)%vector(i)%sum,i=1,freenumber)
    write(fd,*)(ensem%state%Velocity(1)%vector(i)%sum,i=1,freenumber)
    write(fd,*)(ensem%state%Velocity(0)%vector(i)%sum,i=1,freenumber)
    write(fd,*)(ensem%state%Force(1)%vector(i)%sum,i=1,freenumber)
    write(fd,*)(ensem%state%Force(0)%vector(i)%sum,i=1,freenumber)
    !
    return 
  end subroutine  Write_Ensemble_state_Beeman
  !
  subroutine Write_GrandEnsemble_state_Beeman(sys,fd)
    type(GrandEnsemble),intent(in) :: sys
    integer(MYIP),intent(in) :: fd
    !
    integer(MYIP) :: i
    !
    do i = 1,sys%Ensemble_Number
      write(fd,*) i 
      call Write_Ensemble_state_Beeman(sys%System(i),fd)
    end do 
  end subroutine  Write_GrandEnsemble_state_Beeman
  !
  subroutine Read_Ensemble_state_Beeman(ensem,fd)
    type(Ensemble) :: ensem
    integer(MYIP) , intent(in) :: fd
    !
    integer(MYIP) :: i 
    integer(MYIP) :: freenumber
    !
    freenumber = size(ensem%info%Mass%vector)
    read(fd,*)(ensem%state%Displacement(1)%vector(i)%sum,i=1,freenumber)
    read(fd,*)(ensem%state%Displacement(0)%vector(i)%sum,i=1,freenumber)
    read(fd,*)(ensem%state%Velocity(1)%vector(i)%sum,i=1,freenumber)
    read(fd,*)(ensem%state%Velocity(0)%vector(i)%sum,i=1,freenumber)
    read(fd,*)(ensem%state%Force(1)%vector(i)%sum,i=1,freenumber)
    read(fd,*)(ensem%state%Force(0)%vector(i)%sum,i=1,freenumber)
  end subroutine Read_Ensemble_state_Beeman
  !
  subroutine Read_GrandEnsemble_state_Beeman(sys,fd)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: fd
    !
    integer(MYIP) :: i
    !
    do i = 1,sys%Ensemble_Number
      write(fd,*) i 
      call Write_Ensemble_state_Beeman(sys%System(i),fd)
    end do
  end subroutine Read_GrandEnsemble_state_Beeman
  !
end module MDALGORITHM_MODULE
