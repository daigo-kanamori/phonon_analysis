module ATOM_MODULE
  use SYSTEM_MODULE
  implicit none
    type atom
      character(100)                           ::           name
      real(MYPRE)                       ::           mass
    end type atom
  !private :: atom
end module ATOM_MODULE
!MODULE 
module UNITCELL_MODULE
  use SYSTEM_MODULE
  use ATOM_MODULE
  use PARAMETER_MODULE
  use IOTK_MODULE
  use VECTOR_MODULE
  use IBRAV_MODULE
  use KAHAN_MODULE
  !
  implicit none
  !
  type UnitCell
    real(MYPRE), allocatable   ::            Position(:,:)
    real(MYPRE) ::            Lattice_Vector(3,3)
    real(MYPRE) ::            Recipro_Vector(3,3)
    real(MYPRE),allocatable    ::            fcons_h(:,:,:,:,:,:,:)
    real(MYPRE), allocatable   ::            Zeu(:,:,:)
    real(MYPRE) :: epsil(3,3)
    !
    integer(MYIP)                         ::            Lattice_Number(3)
    integer(MYIP),allocatable             ::            atom_index(:)
    integer(MYIP)                         ::            atom_type,atom_number
    type(atom),allocatable          ::            exatom(:)
    character(100)                   ::            Name
    logical                         ::            IsPolorize
  end type UnitCell
  !
  !
  contains
  !
  !
  !======================================================================================
  !                                   UNIT CELL MODULE                                  =
  !======================================================================================
  !
  !
  subroutine allocate_UnitCell(Cell,ntype,natom)
    integer(MYIP),intent(in)        ::     ntype,natom
    type(UnitCell)            ::     Cell
    Cell%IsPolorize   = .false.
    Cell%atom_type    = ntype
    Cell%atom_number  = natom
    !
    print *, ""
    !
    allocate(Cell%Position(3,natom))
    allocate(Cell%atom_index(natom))
    allocate(Cell%exatom(ntype))
    allocate(Cell%Zeu(3,3,natom))

  end subroutine allocate_UnitCell
  !
  !
  subroutine allocate_fcons(fcons_h,Lx,Ly,Lz,natom)

    real(MYPRE) , allocatable ::             fcons_h(:,:,:,:,:,:,:)
    integer(MYIP)                        ::             Lx,Ly,Lz,natom
    integer(MYIP)                        ::             x,y,z

    x = Lx/2
    y = Ly/2
    z = Lz/2
    print *, "力の範囲は"
    print *, Lx,-x,Lx-x-1
    print *, Ly,-y,Ly-y-1
    print *, Lz,-z,Lz-z-1
    allocate(fcons_h(1:3,1:3,1:natom,1:natom,-x:(Lx-x-1),-y:(Ly-y-1),-z:(Lz-z-1)))
    !fcons_lattice_size(1) = [-x,(Lx-x-1)]
    !fcons_lattice_size(2) = [-y,(Ly-y-1)]
    !fcons_lattice_size(3) = [-z,(Lz-z-1)]
    print *, "ALLOCATE FCONS DONE"
  end subroutine allocate_fcons
  !
  !
  subroutine input_UnitCell(Cell,filenumber)

    type(UnitCell)               ::         Cell
    integer(MYIP)                      ::         natom,ntype
    character(100)               ::         tmp_name 
    integer(MYIP)                      ::         i,j,k
    integer(MYIP)                      ::         ind
    integer(MYIP)                      ::         Lx(3)
    integer(MYIP)                      ::         Nx(3)
    integer(MYIP)                      ::         dir1,dir2,atom1,atom2
    integer(MYIP)                      ::         filenumber
    real(MYPRE)             ::         f
    real(MYPRE)             ::         Mass
    real(MYPRE)             ::         vector(3)


    read (filenumber,*) ntype,natom  
    print *, "ntype = ",ntype, ", natom = ",natom
    call Allocate_UnitCell(Cell,ntype,natom)
    !
    do i = 1,3
      read (filenumber,*) Cell%Lattice_Vector(1,i),Cell%Lattice_Vector(2,i),Cell%Lattice_Vector(3,i)
    end do 
    !
    !
    !逆格子ベクトルを作る
    call Make_Recipro_Vector(Cell)
    !
    !原子インデックスと原子の番号を一致させる
    do i = 1,ntype
      read (filenumber,*)ind, tmp_name,Mass  
      Cell%exatom(ind)%name = tmp_name
      Cell%exatom(ind)%mass = Mass
    end do 
    !
    !原子位置を取得
    !
    do i = 1,natom 
      read (filenumber,*) tmp_name ,vector(1),vector(2),vector(3)
      Cell%Position(i,1) = (vector(1)) 
      Cell%Position(i,2) = (vector(2)) 
      Cell%Position(i,3) = (vector(3)) 
      ! 
      do j = 1,Cell%atom_type
        if(Cell%exatom(j)%name .eq. tmp_name) then 
          Cell%atom_index(i) =  j
        end if 
      end do 
    end do 

    read (filenumber,*) Nx(1),Nx(2),Nx(3)
    print *, "Nx ; ",Nx(1),Nx(2),Nx(3)
    print *, Nx(1)/2
    Cell%Lattice_Number(:) = Nx(:)
    call Allocate_fcons(Cell%fcons_h,Nx(1),Nx(2),Nx(3),natom) 
    print *, "INPUT IFC"
    print *, "INPUT IFC INDEX"
    do i = 1,3*3*Cell%atom_number*Cell%atom_number
      read (filenumber,*) dir1,dir2,atom1,atom2
      do j = 1,Nx(1)*Nx(2)*Nx(3)
        read (filenumber,*) Lx(1),Lx(2),Lx(3),f
        !print *, Lx(1),Lx(2),Lx(3),f
        do k = 1,3
          Lx(k) = Lx(k) - 1
          if(Lx(k) > (Nx(k)-1)/2) then 
            Lx(k) = Lx(k) - Nx(k)
          end if 
        end do 
        Cell%fcons_h(dir1,dir2,atom1,atom2,Lx(1),Lx(2),Lx(3)) = f
      end do 
      !
    end do 
    !
    print *, "INPUT DONE"
  end subroutine input_UnitCell
  !
  subroutine Set_UnitCell(Cell,ionode)
    type(UnitCell) :: Cell
    integer(MYIP),intent(in) :: ionode 
    character(100) :: atname
    real(MYPRE) :: celldm(6),atmass,f
    integer(MYIP) :: Nx(3),Lx(3)
    integer(MYIP) :: dir1,dir2,atom1,atom2
    integer(MYIP) :: atidx,attype,nat,ibrav,ntype
    logical :: has_zstar
    integer(MYIP) :: i,j,k,ia
    !
    !
    read(ionode,*)ntype,nat,ibrav,(celldm(i),i=1,6)
    print *, "ntype = ",ntype, "nat = ",nat
    call allocate_UnitCell(Cell,ntype,nat)
    if( ibrav .ne. 0 ) then 
      call ibrav_LatticeVector(Cell%Lattice_Vector,celldm,ibrav)
    else 
      !後ほど実装
    end if
    !
    call Make_Recipro_Vector(Cell)
    !
    do ia = 1,ntype
      read(ionode,*) i,atname,atmass
      Cell%exatom(i)%name = TRIM(atname)
      Cell%exatom(i)%mass = atmass
    end do 
    !
    do ia = 1,nat
      read(ionode,*) atidx,attype,(Cell%Position(i,ia),i=1,3)
      Cell%atom_index(ia) = attype
    end do
    !
    read(ionode,*)has_zstar
    if( has_zstar .eqv. .true. ) then 
      Cell%IsPolorize = .true.
      do i = 1,3
        read(ionode,*) ( Cell%epsil(i,j),j=1,3 )
      end do 
      do ia = 1,nat
        read(ionode,*)
        do i = 1,3
          read(ionode,*) (Cell%Zeu(i,j,ia),j=1,3)
        end do
      end do 
    else
      Cell%Zeu(:,:,:) = 0
      Cell%epsil(:,:) = 0
    end if
    !
    !
    read (ionode,*) Nx(1),Nx(2),Nx(3)
    !
    print *, "Nx ; ",Nx(1),Nx(2),Nx(3)
    Cell%Lattice_Number(:) = Nx(:)
    call allocate_fcons(Cell%fcons_h,Nx(1),Nx(2),Nx(3),nat) 
    print *, "INPUT IFC"
    print *, "INPUT IFC INDEX"
    do i = 1,3*3*Cell%atom_number*Cell%atom_number
      read (ionode,*) dir1,dir2,atom1,atom2
      do j = 1,Nx(1)*Nx(2)*Nx(3)
        read (ionode,*) Lx(1),Lx(2),Lx(3),f
        do k = 1,3
          Lx(k) = Lx(k) - 1
          if(Lx(k) > (Nx(k)-1)/2) then 
            Lx(k) = Lx(k) - Nx(k)
          end if 
        end do 
        Cell%fcons_h(dir1,dir2,atom1,atom2,Lx(1),Lx(2),Lx(3)) = f
      end do 
      !
    end do 
    !
    print *, "INPUT DONE"
    !
    return 
  end subroutine Set_UnitCell 

  !
  subroutine Make_Recipro_Vector(Cell)
    type(UnitCell)                     ::         Cell
    integer(MYIP)                            ::         i
    real(MYPRE)                   ::         x
    do i = 0,2
      Cell%Recipro_Vector(i+1,1) = Cell%Lattice_Vector(1+modulo(i+1,3),2) * Cell%Lattice_Vector(1+modulo(i+2,3),3) &
      - Cell%Lattice_Vector(1+modulo(i+1,3),3) * Cell%Lattice_Vector(1+modulo(i+2,3),2)
    end do 
    do i = 0,2
      Cell%Recipro_Vector(i+1,2) = Cell%Lattice_Vector(1+modulo(i+1,3),3) * Cell%Lattice_Vector(1+modulo(i+2,3),1) &
          - Cell%Lattice_Vector(1+modulo(i+1,3),1) * Cell%Lattice_Vector(1+modulo(i+2,3),3)
    end do

    do i = 0,2
      Cell%Recipro_Vector(i+1,3) = Cell%Lattice_Vector(1+modulo(i+1,3),1) * Cell%Lattice_Vector(1+modulo(i+2,3),2) & 
        - Cell%Lattice_Vector(1+modulo(i+1,3),2) * Cell%Lattice_Vector(1+modulo(i+2,3),1)
    end do 
    x = 0d0
    do i = 1,3
      x = x + Cell%Lattice_Vector(i,1) * Cell%Recipro_Vector(i,1) 
    end do 
    write(*,*) "x = ",x
    do i = 1,3
      Cell%Recipro_Vector(:,i) = Cell%Recipro_Vector(:,i) * 2 * PI / x 
    end do 
    !
    !
    return 
  end subroutine Make_Recipro_Vector

  subroutine input_CellName(uCell,Name)
    character(50)      ::       Name 
    type(UnitCell)     ::       uCell

    uCell%Name= Name
    uCell%Name = TRIM(uCell%Name)
  end subroutine input_CellName


  subroutine ChangeUnit_UnitCell(Cell,uenergy,umetre,uforce,umass)
    type(UnitCell)   :: Cell
    real(MYPRE),intent(in) :: uenergy,umetre,umass,uforce

    print *, umetre
    Cell%Lattice_Vector(:,:)   =   Cell%Lattice_Vector(:,:) * umetre
    Cell%Recipro_Vector(:,:)       =   Cell%Recipro_Vector(:,:) / umetre
    Cell%fcons_h          =   Cell%fcons_h* uforce
    Cell%exatom(:)%mass   =   Cell%exatom(:)%mass * umass
    print *, "UNIT CHANGE DONE"
  end subroutine ChangeUnit_UnitCell

  subroutine ChangeUnit_NA(Cell,ucharge,udietensor)
    type(UnitCell) :: Cell
    real(MYPRE) :: ucharge,udietensor
    if(Cell%IsPolorize) then 
      Cell%Zeu(:,:,:) = Cell%Zeu(:,:,:) * ucharge
      Cell%epsil(:,:) = Cell%epsil(:,:) * udietensor
      print *, "誘電率、有効電荷の単位を変更しました"
    end if 
  end subroutine 
  !
  !
  subroutine Debug_UnitCell(Cell)
    type(UnitCell),intent(in) :: Cell
    integer(MYIP) :: i,j
    integer(MYIP) :: ia

    print *, "%================================================%"  
    print *, "%        CELL NAME = " , Cell%name, "            %"
    print *, "%================================================%"  
    print *, "ntyp = ",Cell%atom_type , ", natom = ",Cell%atom_number
    print *, "OUTPUT ATOM INDEX"
    do ia = 1,Cell%atom_number
      print *, ia,"th atom index : ",Cell%atom_index(ia)
    end do
    print *,"OUTPUT EXIST ATOMS"
    do ia = 1,Cell%atom_type
      print *,ia,"th index atom is ",Cell%exatom(ia)%name
      print *, "mass = ",Cell%exatom(ia)%mass
    end do 
    print *,"OUTPUT POSITION"
    do ia = 1,Cell%atom_number
      print *,ia,"th atom :",Cell%Position(1,ia),",",Cell%Position(2,ia),",",Cell%Position(3,ia)
    end do 
    !
    print *, "OUTPUT LATTICE VECTOR"
    do i = 1,3
      print *, i," : ",Cell%Lattice_Vector(1,i),Cell%Lattice_Vector(2,i),Cell%Lattice_Vector(3,i)
    end do 
    print *,"OUTPUT REV VECTOR"
    do i = 1,3
      print *,i," : ",Cell%Recipro_Vector(1,i),Cell%Recipro_Vector(2,i),Cell%Recipro_Vector(3,i)
    end do 
    print *, "OUTPUT FCONS OF SELF TERM"
    do i = 1,Cell%atom_number
      print *, i,"th atom"
      do j = 1,3
        print *, "  ",j," : = ",Cell%fcons_h(j,j,i,i,0,0,0)
      end do 
    end do 
    print *,"OUTPUT PARAMETER OF CELL DONE"
  end subroutine 
  !
  !
  subroutine output_fcons(Cell,Lx,Ly,Lz,a,b,filenumber)
    integer(MYIP) :: Lx,Ly,Lz,a,b
    integer(MYIP) :: i,j,filenumber
    type(UnitCell) :: Cell
    print *, Lx,Ly,Lz,a,b
    do i = 1,3
      do j = 1,3
        write (filenumber,fmt='(E13.5e3, A)',advance='no') Cell%fcons_h(i,j,a,b,Lx,Ly,Lz),","
      end do 
      write (filenumber,*) 
    end do
    write (*,*)
  end subroutine output_fcons
  
  !ユニットセルを開放する
  !メモリを節約したい
  subroutine ImposeASR(Cell,CX,CY,CZ)
    type(UnitCell)     :: Cell
    type(KahanSum) :: asr
    integer(MYIP),intent(in) :: CX(2),CY(2),CZ(2)
    integer(MYIP)            :: i,j,ia,ja,il,jl,kl
    !
    print *, "ASRを適用します"
    print *, "まずshort rangeから"
    !
    do i = 1,3
    do j = 1,3
      do ia = 1,Cell%atom_number
        call kahan_zero(asr)
        do kl = CZ(1) ,CZ(2)
        do jl = CY(1) ,CY(2)
        do il = CX(1) ,CX(2)
        do ja = 1,Cell%atom_number
          call add_kahan(asr,Cell%fcons_h(i,j,ia,ja,il,jl,kl))
        end do 
        end do 
        end do 
        end do 
        Cell%fcons_h(i,j,ia,ia,0,0,0) = Cell%fcons_h(i,j,ia,ia,0,0,0) - asr%sum
      end do
    end do 
    end do 
    print *,"short range 完了"
    !
    !
    return 
  end subroutine ImposeASR
  !
  subroutine Adjust_fcons(Cell)
    type(UnitCell) :: Cell
    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: i,j
    integer(MYIP) :: ia,ja
    integer(MYIP) :: nat
    real(MYPRE) :: tmp
    !
    nat = Cell%atom_number
    do il = lbound(Cell%fcons_h,5),ubound(Cell%fcons_h,5)
    do jl = lbound(Cell%fcons_h,6),ubound(Cell%fcons_h,6)
    do kl = lbound(Cell%fcons_h,7),ubound(Cell%fcons_h,7)
      do ia = 1,nat
      do ja = 1,nat
      do i = 1,3
      do j = 1,3
        tmp = 0.5_MYPRE * ( Cell%fcons_h(i,j,ia,ja,il,jl,kl) + Cell%fcons_h(j,i,ja,ia,-il,-jl,-kl) ) 
        Cell%fcons_h(i,j,ia,ja,il,jl,kl) = tmp
        Cell%fcons_h(j,i,ja,ia,-il,-jl,-kl) = tmp
      end do
      end do
      end do
      end do
    end do 
    end do
    end do
  end subroutine
  !
  !
  subroutine Adjust_Selfterm(Cell,amp)
    type(UnitCell) :: Cell
    real(MYPRE) ,intent(in) :: amp
    !
    type(KahanSum) :: asr
    integer(MYIP) :: ia,ja,i,j
    integer(MYIP) :: nat
    !
    print *, "ADJUST SELF TERM"
    !
    do i = 1,3
    do ia = 1,Cell%atom_number
      Cell%fcons_h(i,i,ia,ia,0,0,0) = Cell%fcons_h(i,i,ia,ia,0,0,0) *amp
    end do 
    end do 
    !
    return 
  end subroutine Adjust_Selfterm
  !
  subroutine Calculate_volume_UnitCell(Cell,volume)
    type(UnitCell),intent(in) :: Cell
    real(MYPRE) :: volume
    !
    real(MYPRE) :: vec(3)
    real(MYPRE) :: tmp
    type(KahanSum) :: dot
    integer(MYIP) :: i
    !
    do i = 0,2
      vec(1+i)  = Cell%Lattice_Vector(1+modulo(i+1,3),1) * Cell%Lattice_Vector(1+modulo(i+2,3),2) &
              - Cell%Lattice_Vector(1+modulo(i+1,3),2) * Cell%Lattice_Vector(1+modulo(i+2,3),1)
    end do
    !
    call kahan_zero(dot)
    do i = 1,3
      tmp = vec(i) * Cell%Lattice_Vector(i,3)
      call add_kahan(dot,tmp)
    end do 
    !
    volume = abs(dot%sum)
    !
    return 
  end subroutine Calculate_volume_UnitCell
  !
end module UNITCELL_MODULE
!
