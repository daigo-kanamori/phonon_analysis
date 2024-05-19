module ENSEMBLE_MODULE
  use SYSTEM_MODULE
  use VECTOR_MODULE
  use KAHAN_MODULE
  use IOTK_MODULE
  !
  !
  implicit none
  !
  !
  type Ensemble_Variables 
    type(Vector_kahan),allocatable :: Displacement(:)
    type(Vector_kahan),allocatable :: Velocity(:)
    type(Vector_kahan),allocatable :: Force(:)
  end type
  !
  type Ensemble_Information
    type(Vector_d) :: Mass 
  end type 
  !
  type Ensemble
    type(Ensemble_Variables) :: state
    type(Ensemble_Information) :: info
  end type
  !
  type GrandEnsemble
    type(Ensemble),allocatable :: System(:)
    integer(MYIP) ,allocatable :: Ensemble_list(:,:,:)
    integer(MYIP) :: Ensemble_Number
    real(MYPRE) :: step
  end type
  !
  !
  contains 
  !
  subroutine allocate_GrandEnsemble(sys,Nx,Ny,Nz,step)
    type(GrandEnsemble) :: sys
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2) 
    real(MYPRE),intent(in) :: step
    integer(MYIP) :: cell_number
    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: ide
    !
    allocate(sys%Ensemble_list(Nx(1):Nx(2),Ny(1):Ny(2),Nz(1):Nz(2)))
    !
    cell_number = (Nx(2)-Nx(1)+1) * (Ny(2)-Ny(1)+1) * (Nz(2)-Nz(1)+1)
    allocate(sys%System(cell_number))
    !
    ide = 1
    !
    do kl = Nz(1),Nz(2)
    do jl = Ny(1),Ny(2)
    do il = Nx(1),Nx(2)
      sys%Ensemble_list(il,jl,kl) = ide
      allocate(sys%System(ide)%state%Displacement(1))
      allocate(sys%System(ide)%state%Velocity(1))
      allocate(sys%System(ide)%state%Force(-1:1))
      ide = ide + 1
    end do
    end do
    end do
    !
    sys%step = step
    sys%Ensemble_Number = 0
    !
    return 
  end subroutine allocate_GrandEnsemble
  !

  !
  !===========================================================!
  !     ROUTINE : WRITE THE STATE OF ENSEMBLE                 !
  !===========================================================!
  !
  !
  subroutine  Write_Ensemble_state(ensem,fd)
    type(Ensemble),intent(in) :: ensem
    integer(MYIP),intent(in) :: fd
    integer(MYIP) :: freenumber
    !
    freenumber = size(ensem%info%Mass%vector)
    !write(fd,*)freenumber
    call Write_Vector_kahan_line(ensem%state%Displacement(1),fd)
    call Write_Vector_kahan_line(ensem%state%Velocity(1),fd)
    call Write_Vector_kahan_line(ensem%state%Force(1),fd)
    !
    return 
  end subroutine  Write_Ensemble_state
  !
  subroutine Write_GrandEnsemble_state(sys,fd)
    type(GrandEnsemble),intent(in) :: sys
    integer(MYIP),intent(in) :: fd 
    integer(MYIP) :: ide
    !
    !
    do ide = 1,ubound(sys%System,1)
      if(allocated(sys%System(ide)%info%Mass%vector)) then 
        write(fd,*) ide
        call Write_Ensemble_state(sys%System(ide),fd)
      end if 
    end do 
    !
    return 
  end subroutine  Write_GrandEnsemble_state
  !
  !
  !
  !=================================================================!
  !          WRITE ENSEMBLE INFO                                    !
  !=================================================================!
  !
  subroutine Write_Ensemble_info(ensem,fd)
    type(Ensemble) ,intent(in) :: ensem
    integer(MYIP) ,intent(in) :: fd
    integer(MYIP) :: freenumber
    !
    freenumber = size(ensem%info%Mass%vector)
    write(fd,*)freenumber
    call Write_Vector_d_line(ensem%info%Mass%vector,fd)
    !
    return 
  end subroutine  Write_Ensemble_info
  !
  subroutine Write_GrandEnsemble_info(sys,fd)
    type(GrandEnsemble), intent(in) :: sys
    integer(MYIP),intent(in) :: fd

    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: kind_cell
    integer(MYIP) :: ide
    !
    kind_cell = 0
    !
     
    kind_cell = sys%Ensemble_Number
    !
    Nx = [lbound(sys%Ensemble_list,1),ubound(sys%Ensemble_list,1)]
    Ny = [lbound(sys%Ensemble_list,2),ubound(sys%Ensemble_list,2)]
    Nz = [lbound(sys%Ensemble_list,3),ubound(sys%Ensemble_list,3)]
    write(fd,*) kind_cell
    !
    !
    do kl = Nz(1),Nz(2)
    do jl = Ny(1),Ny(2)
    do il = Nx(1),Nx(2)
      ide = sys%Ensemble_list(il,jl,kl)
      if(allocated(sys%System(ide)%info%Mass%vector) .eqv. .false.) then 
        exit
      end if
      write(fd,*) il,",",jl,",",kl,","
      call Write_Ensemble_info(sys%System(ide),fd)
    end do 
    end do 
    end do
    !
    return 
  end subroutine Write_GrandEnsemble_info
  !
  !======================================================================!
  !                 READ THE STATE OF ENSEMBLE                           !
  !======================================================================!
  ! DONT COMPLETE
  !
  subroutine Read_Ensemble_state(ensem,fd)
    type(Ensemble) :: ensem
    integer(MYIP),intent(in) :: fd
    integer(MYIP) :: freenumber
    integer(MYIP) :: i
    !
    read(fd,*) freenumber
    read(fd,*) (ensem%state%Displacement(1)%vector(i)%sum,i=1,freenumber)
    read(fd,*) (ensem%state%Velocity(1)%vector(i)%sum,i=1,freenumber)
    read(fd,*) (ensem%state%Force(1)%vector(i)%sum,i=1,freenumber)
    !
    return 
  end subroutine Read_Ensemble_state
  !
  subroutine Read_GrandEnsemble_state(sys,fd)
    type(GrandEnsemble) :: sys
    integer(MYIP) ,intent(in) :: fd
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    integer(MYIP) :: ir,jr,kr
    integer(MYIP) :: ide,idx
    integer(MYIP) ::cell_number
    !
    Nx = [lbound(sys%Ensemble_list,1),ubound(sys%Ensemble_list,1)]
    Ny = [lbound(sys%Ensemble_list,2),ubound(sys%Ensemble_list,2)]
    Nz = [lbound(sys%Ensemble_list,3),ubound(sys%Ensemble_list,3)]
    !
    read(fd,*) cell_number
    !
    do idx = 1,cell_number
      read(fd,*)ir,jr,kr
      ide = sys%Ensemble_list(ir,jr,kr)
      call Read_Ensemble_state(sys%System(ide),fd)
    end do 
    !
    return 
  end subroutine Read_GrandEnsemble_state
  !
  !======================================================================!
  !                   INITIALIZE ENSEMBLE ROUTINE                        ! 
  !======================================================================!
  !
  subroutine Initialize_Ensemble_Uniform_Disp(ensem,amp)
    type(Ensemble) :: ensem
    double precision,intent(in) :: amp
    double precision :: rnd
    integer(MYIP) i 
    !
    !
    do i = 1,size(ensem%state%Displacement(1)%vector)
      call random_number(rnd)
      rnd = rnd - 0.5d0
      rnd = rnd * amp
      call add_kahan(ensem%state%Displacement(1)%vector(i),rnd)
    end do 
    !
    return 
  end subroutine Initialize_Ensemble_Uniform_Disp 
  !
  subroutine Initialize_Ensemble_Uniform_Velo(ensem,amp)
    type(Ensemble) :: ensem
    double precision,intent(in) :: amp
    double precision :: rnd
    integer(MYIP) i 
    !
    !
    do i = 1,size(ensem%state%Velocity(1)%vector)
      call random_number(rnd)
      rnd = rnd - 0.5d0
      rnd = rnd * amp
      call add_kahan(ensem%state%Velocity(1)%vector(i),rnd)
    end do 
    !
    return 
  end subroutine Initialize_Ensemble_Uniform_Velo
  !
  subroutine Initialize_GrandEnsemble_Uniform_Disp(sys,amp)
    type(GrandEnsemble) :: sys
    double precision,intent(in) :: amp
    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: ide
    !
    !
    do kl = lbound(sys%Ensemble_list,3),ubound(sys%Ensemble_list,3)
    do jl = lbound(sys%Ensemble_list,2),ubound(sys%Ensemble_list,2)
    do il = lbound(sys%Ensemble_list,1),ubound(sys%Ensemble_list,1)
      ide = sys%Ensemble_list(il,jl,kl)
      call Initialize_Ensemble_Uniform_Disp(sys%System(ide),amp)
    end do 
    end do 
    end do 
    !
    return 
  end subroutine Initialize_GrandEnsemble_Uniform_Disp
  !
  !
  !==================================================!
  !    ENERGY CALCULATION ROUTINE                    !
  !==================================================!
  !
  double precision function Calculate_KineticEnergy_GrandEnsemble(sys)
    type(GrandEnsemble),intent(in) :: sys
    integer(MYIP) :: ide
    type(KahanSum) :: keg
    !
    call kahan_zero(keg)
    !
    do ide = 1,size(sys%System,1)
      if(allocated(sys%System(ide)%info%Mass%vector)) then
        call Calculate_KineticEnergy_Ensemble(sys%System(ide),keg)
      end if
    end do 
    !
    Calculate_KineticEnergy_GrandEnsemble = keg%sum
    return 
  end function Calculate_KineticEnergy_GrandEnsemble
  !
  !
  subroutine Calculate_KineticEnergy_Ensemble(ensem,keg)
    type(Ensemble),intent(in) :: ensem
    type(KahanSum) :: keg
    double precision :: tmp
    integer(MYIP)  :: i
    !
    do i = lbound(ensem%state%Velocity(1)%vector,1),ubound(ensem%state%Velocity(1)%vector,1)
      tmp = ensem%state%Velocity(1)%vector(i)%sum
      tmp = 0.5d0 * ensem%info%Mass%vector(i) * tmp * tmp
      call add_kahan(keg,tmp)
    end do 
    !
    return
  end subroutine Calculate_KineticEnergy_Ensemble
  ! 
  ! 
  ! 
  double precision function Calculate_PotentialEnergy_GrandEnsemble(sys)
    type(GrandEnsemble),intent(in) :: sys
    integer(MYIP) :: ide
    type(KahanSum) :: peg
    !
    call kahan_zero(peg)
    !
    do ide = 1,size(sys%System,1)
      if(allocated(sys%System(ide)%info%Mass%vector)) then 
        call Calculate_PotentialEnergy_Ensemble(sys%System(ide),peg)
      end if 
    end do 
    !
    Calculate_PotentialEnergy_GrandEnsemble = peg%sum
    !
    return 
  end function Calculate_PotentialEnergy_GrandEnsemble
  ! 
  ! 
  ! 
  subroutine Calculate_PotentialEnergy_Ensemble(ensem,peg)
    type(Ensemble),intent(in) :: ensem
    type(KahanSum) :: peg
    double precision :: tmp
    integer(MYIP)  :: i
    !
    tmp = 0d0
    !
    do i = lbound(ensem%state%Displacement(1)%vector,1),ubound(ensem%state%Displacement(1)%vector,1)
      tmp = -ensem%state%Displacement(1)%vector(i)%sum * ensem%state%Force(1)%vector(i)%sum * 0.5d0 
      call add_kahan(peg,tmp)
    end do 
    !
    return
  end subroutine Calculate_PotentialEnergy_Ensemble
  !
end module ENSEMBLE_MODULE


