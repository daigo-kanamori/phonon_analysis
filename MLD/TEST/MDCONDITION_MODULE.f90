module MDCONDITION_MODULE
  use LATTICE_MODULE
  implicit none
  !
  contains
  !
  subroutine Calculate_Force(Lattice,now,Nx,Ny,Nz)
    type(Lattice_MD)      :: Lattice
    integer,intent(in) :: Nx,Ny,Nz,now
    !
    integer :: Mx,My,Mz
    integer :: Fx,Fy,Fz
    integer :: nide,mide
    integer :: nidx,midx
    integer :: nnat,mnat
    integer :: il,jl,kl
    integer :: nfree,mfree
    integer :: ia,ja
    integer :: i,j
    double precision :: nfcons,mfcons
    logical :: flag_condition
    !
    nide = Lattice%sys%Ensemble_list(Nx,Ny,Nz)
    nidx = Lattice%Index(Nx,Ny,Nz)
    nnat = Lattice%Cell(nidx)%atom_number
    flag_condition = .true.
    !
    do ia = 1,nnat 
    do i = 1,3
    do kl = Lattice%Cz(1),Lattice%Cz(2)
    do jl = Lattice%Cy(1),Lattice%Cy(2)
    do il = Lattice%Cx(1),Lattice%Cx(2)
      !
      Mx = il + Nx 
      My = jl + Ny 
      Mz = kl + Nz 
      !
      call Check_BoundaryCondition(Lattice,now,Nx,Ny,Nz,Mx,My,Mz,flag_condition)
      if(flag_condition .eqv. .false.) exit
      mide = Lattice%sys%Ensemble_list(Mx,My,Mz)
      midx = Lattice%Index(Mx,My,Mz)
      mnat = Lattice%Cell(midx)%atom_number
      !
      do ja = 1,mnat
      do j = 1,3
        call Check_FreeCondition(Lattice,now,Nx,Ny,Nz,i,ia,nfree,nfcons,&
                                             Mx,My,Mz,j,ja,mfree,mfcons,& 
                                             il,jl,kl,&
                                             flag_condition)
        !
        if(flag_condition .eqv. .true. ) then
          call Calculate_Force_harmonic(Lattice,nide,nfree,nfcons,&
                                                mide,mfree,mfcons)
        end if
      end do 
      end do 
      !
    end do 
    end do 
    end do 
    end do 
    end do 
    !
  end subroutine 
  !
  subroutine Check_BoundaryCondition(Lattice,now,Nx,Ny,Nz,Mx,My,Mz,boundary_flag)
    type(Lattice_MD),intent(in) :: Lattice
    integer,intent(in) :: now
    integer,intent(in) :: Nx,Ny,Nz
    integer :: Mx,My,Mz
    logical :: boundary_flag
    !
    boundary_flag = .true.
    !
    if(Mx > Lattice%Nx(2)) then 
      Mx = Mx - (Lattice%Nx(2) - Lattice%Nx(1)+1)
    end if
    if(Mx < Lattice%Nx(1)) then 
      Mx = Mx + (Lattice%Nx(2) - Lattice%Nx(1)+1)
    end if
    !
    if(My > Lattice%Ny(2)) then 
      My = My - (Lattice%Ny(2) - Lattice%Ny(1)+1)
    end if
    if(My < Lattice%Ny(1)) then 
      My = My + (Lattice%Ny(2) - Lattice%Ny(1)+1)
    end if
    !
    if(Mz > Lattice%Nz(2)) then 
      Mz = Mz - (Lattice%Nz(2) - Lattice%Nz(1)+1)
    end if
    if(Mz < Lattice%Nz(1)) then 
      Mz = Mz + (Lattice%Nz(2) - Lattice%Nz(1)+1)
    end if
    !
    !
    return 
  end subroutine Check_BoundaryCondition
  ! 
  subroutine Check_FreeCondition(Lattice,now,Nx,Ny,Nz,i,ia,nfree,nfcons,& 
                                             Mx,My,Mz,j,ja,mfree,mfcons,&
                                             Fx,Fy,Fz,&
                                             time_flag)
    type(Lattice_MD) :: Lattice
    integer,intent(in) :: now
    integer,intent(in) :: Nx,Ny,Nz
    integer,intent(in) :: Mx,My,Mz
    integer,intent(in) :: Fx,Fy,Fz
    integer,intent(in) :: i,j
    integer,intent(in) :: ia,ja
    integer :: nfree,mfree
    double precision :: nfcons,mfcons
    integer  :: midx,nidx
    logical :: time_flag
    !
    time_flag = .true.
    !
    nfree = 3*(ia-1) + i
    nidx = Lattice%Index(Nx,Ny,Nz)
    nfcons = Lattice%Cell(nidx)%fcons_h(i,j,ia,ja,Fx,Fy,Fz)
    !
    mfree = 3*(ja-1) + j
    midx = Lattice%Index(Mx,My,Mz)
    mfcons = Lattice%Cell(midx)%fcons_h(j,i,ja,ia,-Fx,-Fy,-Fz)
    !
    return 
  end subroutine Check_FreeCondition
  !
  subroutine Apply_BindCondition_System(Lattice,now)
    type(Lattice_MD) :: Lattice
    integer,intent(in) :: now
    !
    !
    return 
  end subroutine Apply_BindCondition_System
  !
  subroutine Apply_BindCondition_Displacment(Lattice,now)
    type(Lattice_MD) :: Lattice
    integer ,intent(in) :: now
    !
    return
  end subroutine 
  !
  subroutine Apply_BindCondition_Force(Lattice,now)
    type(Lattice_MD) :: Lattice
    integer,intent(in) :: now
    !
    !
    return 
  end subroutine 
  !
  subroutine Apply_BindCondition_Velocity(Lattice,now)
    type(Lattice_MD) :: Lattice
    integer,intent(in) :: now
    !
    return 
  end subroutine
  !
  subroutine Calculate_Force_harmonic(Lattice,nide,nfree,nfcons,mide,mfree,mfcons)
    type(Lattice_MD) :: Lattice
    integer , intent(in) :: nide,mide
    integer ,intent(in) :: nfree,mfree
    double precision,intent(in) :: nfcons,mfcons
    double precision :: tmp
    !
    tmp = - 0.5d0*(nfcons + mfcons) * Lattice%sys%System(mide)%state%Displacement(1)%vector(mfree)%sum
    call add_kahan(Lattice%sys%System(nide)%state%Force(1)%vector(nfree),tmp)
    !
    return 
  end subroutine  Calculate_Force_harmonic
  !
  !
  !
  subroutine analyze(Lattice,now,fds) 
    type(Lattice_MD),intent(in) :: Lattice
    integer,intent(in) :: now
    integer,intent(in),allocatable :: fds(:)
    integer :: N,i,j
    integer :: ide
    integer :: il,jl,kl
    !
    return 
  end subroutine analyze
end module MDCONDITION_MODULE
