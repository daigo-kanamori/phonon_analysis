!
! LATTICE_MODULE
!
module LATTICE_MODULE
  use SYSTEM_MODULE
  use VECTOR_MODULE
  use UNITCELL_MODULE
  use KAHAN_MODULE
  use ENSEMBLE_MODULE
  use MDALGORITHM_MODULE
  implicit none
  !
  !
  type Lattice_MD
    type(GrandEnsemble) :: sys
    type(UnitCell),allocatable :: Cell(:)
    integer(MYIP),allocatable :: Index(:,:,:)
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    integer(MYIP) :: Cx(2),Cy(2),Cz(2)
    integer(MYIP) :: allocated_number
  end type Lattice_MD
  !
  !
  contains 
  !
  !
  subroutine allocate_Lattice_MD(lattice,Nx,Ny,Nz,Cx,Cy,Cz,cnum,unit_time)
    type(Lattice_MD) :: lattice
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2)
    integer(MYIP),intent(in) :: Cx(2),Cy(2),Cz(2)
    integer(MYIP),intent(in) :: cnum
    real(MYPRE),intent(in) :: unit_time
    !
    !
    allocate(lattice%Cell(cnum))
    allocate(lattice%Index(Nx(1):Nx(2),Ny(1):Ny(2),Nz(1):Nz(2)))
    lattice%Index(:,:,:) = -1_MYIP
    lattice%allocated_number = 0_MYIP
    lattice%Nx(:) = Nx(:)
    lattice%Ny(:) = Ny(:)
    lattice%Nz(:) = Nz(:)
    lattice%Cx(:) = Cx(:)
    lattice%Cy(:) = Cy(:)
    lattice%Cz(:) = Cz(:)
    !
    !
    call allocate_GrandEnsemble(lattice%sys,Nx,Ny,Nz,unit_time)
    !
    return 
    !
  end subroutine allocate_Lattice_MD
  !
  !
  subroutine Set_Lattice_MD(lattice,Cell,Lx,Ly,Lz)
    type(Lattice_MD) :: lattice
    type(UnitCell) :: Cell
    integer(MYIP),intent(in) :: Lx,Ly,Lz
    !
    integer(MYIP) :: idx,nat
    integer(MYIP) :: ide
    integer(MYIP) :: ia,i,aidx
    !
    if(lattice%Index(Lx,Ly,Lz) > 0) then
      print *, "THE CELL IS ALREADY ALLOCATED."
      print *, "IGONRE"
      return 
    end if
    !
    call Check_Cell_MD(lattice,Cell,idx)
    Lattice%Index(Lx,Ly,Lz) = idx
    nat = lattice%Cell(idx)%atom_number
    call Set_GrandEnsemble(lattice%sys,Lx,Ly,Lz,nat*3)
    !
    ide = lattice%sys%Ensemble_list(Lx,Ly,Lz)
    !
    do ia = 1_MYIP,Cell%atom_number
      aidx =  Cell%atom_index(ia)
      do i = 1_MYIP,3_MYIP
        call Set_Ensemble_Mass(lattice%sys%System(ide),Cell%exatom(aidx)%mass,3_MYIP*(ia-1_MYIP)+i)
      end do 
    end do 
    !
    return 
  end subroutine 
  !
  !
  subroutine Check_Cell_MD(lattice,Cell,idx)
    type(Lattice_MD) :: lattice
    type(UnitCell),intent(in) :: Cell
    integer(MYIP) :: idx
    integer(MYIP) :: i
    logical :: isntfind
    !
    isntfind = .true.
    idx = Lattice%allocated_number + 1
    do i = 1,Lattice%allocated_number
      if(Lattice%Cell(i)%name .eq. Cell%name) then 
        idx = i
        isntfind = .false.
      end if
    end do 
    !
    if(isntfind .eqv. .true.) then 
      Lattice%Cell(idx) = Cell
      Lattice%allocated_number = Lattice%allocated_number+1_MYIP
    end if
  end subroutine 
  !
  !
  subroutine MakeBulk_z(Cell,dir)
    type(UnitCell) :: Cell
    integer(MYIP)  :: dir
    integer(MYIP) :: i,j,k
    integer(MYIP) :: xmax,ymax,zmax,xmin,ymin,zmin
    xmax = int(ubound(Cell%fcons_h,5),MYIP)
    xmin = int(lbound(Cell%fcons_h,5),MYIP)
    ymax = int(ubound(Cell%fcons_h,6),MYIP)
    ymin = int(lbound(Cell%fcons_h,6),MYIP)
    zmax = int(ubound(Cell%fcons_h,7),MYIP)
    zmin = int(lbound(Cell%fcons_h,7),MYIP)
    !
    !dir == 1 then cutting +z 
    !
    if(dir == 1_MYIP) then
      !print *, "dir = ",dir
      do k = 1_MYIP,zmax
      do j = ymin,ymax
      do i = xmin,xmax
        Cell%fcons_h(:,:,:,:,i,j,k) = 0.0_MYPRE
      enddo
      enddo
      enddo 
    else if(dir == 2_MYIP) then
      !print *, "dir = ",dir
      do k = zmin,-1_MYIP
      do j = ymin,ymax
      do i = xmin,xmax
        Cell%fcons_h(:,:,:,:,i,j,k) = 0d0
      enddo
      enddo
      enddo
    end if
  end subroutine
  !
  !
  subroutine Debug_Lattice_MD(Lattice)
    type(Lattice_MD),intent(in) :: Lattice

    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: idx
    integer(MYIP) :: ide
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    !
    Nx(:) = Lattice%Nx(:)
    Ny(:) = Lattice%Ny(:)
    Nz(:) = Lattice%Nz(:)
    print *, "%============================================="
    print *, "%               Debug Lattice                ="
    print *, "%============================================="
    do kl = Nz(1),Nz(2)
    do jl = Ny(1),Ny(2)
    do il = Nx(1),Nx(2)
      idx = Lattice%Index(il,jl,kl)
      ide = Lattice%sys%Ensemble_list(il,jl,kl)
      print *, "%%%% Index = ",il,jl,kl 
      print *, "% CeLL = ",Lattice%Cell(idx)%name, "idx = ",idx
      print *, "freenumber = ",size(Lattice%sys%System(ide)%info%Mass%vector)
    end do 
    end do
    end do
    print *, "Lattice%sys%Ensemble_Number = ",Lattice%sys%Ensemble_Number
  end subroutine 
  !
  !
end module  LATTICE_MODULE
