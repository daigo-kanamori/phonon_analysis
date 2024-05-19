!ここはLattice Dynamicsによってフォノンの分散を求めようという腹です
!必要なものは格子の情報とかを必要とします
module LD_MODULE
  use SYSTEM_MODULE
  use PARAMETER_MODULE
  use IOTK_MODULE
  use UNITCELL_MODULE
  use ZHEEV_MODULE
  implicit none
  !
  !
  contains 
  !
  !
  subroutine LatticeDynamics(PrimitiveCell,Dynmat,Eigvec,Wavevector,qhat,Nx,Ny,Nz)
    type(UnitCell),intent(in)                    ::     PrimitiveCell
    complex(MYPRE), allocatable   ::     Dynmat(:,:)
    real(MYPRE),intent(in)                  ::     Wavevector(3),qhat(3)
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2)
    real(MYPRE) ,  allocatable   ::     Eigvec(:)
    integer(MYIP)                           ::     aatom,batom
    integer(MYIP)                           ::     i,j,k,ix,jx,jz,ia,ja,ka,il,jl,kl,id,jd,jk,ir,jr
    DynMat(:,:) = 0_MYPRE
    !
    call Make_DM(PrimitiveCell,Dynmat,Wavevector,qhat,Nx,Ny,Nz) 
    call MakeMatrixHermite(Dynmat)
    call exec_zheev(Dynmat,EigVec)
    !
    return 
  end subroutine LatticeDynamics
  !  
  !  
  subroutine Make_DM(PrimitiveCell,Dynmat,Wavevector,qhat,Nx,Ny,Nz)  
    type(UnitCell),intent(in)                    ::     PrimitiveCell
    complex(MYPRE),allocatable    ::     Dynmat(:,:)
    real(MYPRE)                  ::     Wavevector(3),qhat(3)
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2)
    !
    integer(MYIP)                           ::     i,j,k,ia,ja,il,jl,kl,id
    real(MYPRE)                  ::     theta
    complex(MYPRE), allocatable   ::     dyn(:,:,:,:)
    integer(MYIP)                           ::     MaxLx(3),MinLx(3)
    integer(MYIP)                           ::     aatom,batom
    !
    allocate(dyn(3,3,PrimitiveCell%atom_number,PrimitiveCell%atom_number))
    dyn(:,:,:,:) = (0d0,0d0)
    !
    call calc_shortpart(PrimitiveCell,dyn,Wavevector,Nx,Ny,Nz)
    if( PrimitiveCell%IsPolorize .eqv. .true.) then 
      call calc_longpart(PrimitiveCell,dyn,Wavevector,qhat,Nx,Ny,Nz)
      !print *, "LOGPART DONE"
    end if
    do ja = 1,PrimitiveCell%atom_number
    do ia = 1,PrimitiveCell%atom_number
    do j = 1,3
    do i = 1,3
       Dynmat(3*(ia-1)+i,3*(ja-1)+j) = dyn(i,j,ia,ja)   
    end do 
    end do 
    end do 
    end do 

  end subroutine Make_DM
  

  subroutine calc_shortpart(PrimitiveCell,dyn,Wavevector,Nx,Ny,Nz)
    type(UnitCell),intent(in)  ::     PrimitiveCell
    complex(MYPRE),allocatable    ::     dyn(:,:,:,:)
    real(MYPRE),intent(in) ::     Wavevector(3)
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2)
    !
    integer(MYIP)  ::     aatom,batom
    real(MYPRE) :: amass,bmass
    real(MYPRE)    ::     theta
    complex(MYPRE) ::     ftmp
    integer(MYIP)        ::     i,j,k,ia,ja,ka,il,jl,kl,dir
    real(MYPRE) :: Rx(3)
    integer(MYIP) :: ida,idb
    ftmp    = 0d0
    theta   = 0d0

    do i = 1,3
    do j = 1,3
    do ia = 1,PrimitiveCell%atom_number
    do ja = 1,PrimitiveCell%atom_number
      ftmp = 0d0
      do kl = Nz(1),Nz(2)
      do jl = Ny(1),Ny(2)
      do il = Nx(1),Nx(2)
        theta = 0d0
        Rx(:) = 0d0
        !
        do dir = 1,3
          Rx(dir) =  il*PrimitiveCell%Lattice_Vector(dir,1) + jl * PrimitiveCell%Lattice_Vector(dir,2)  &
                    + kl * PrimitiveCell%Lattice_Vector(dir,3) 
        end do
        do dir = 1,3
          theta = theta + Rx(dir) * Wavevector(dir)
        end do 
        !
        ftmp = ftmp + PrimitiveCell%fcons_h(i,j,ia,ja,il,jl,kl)  * CMPLX(cos(theta),sin(theta),kind=kind(0d0))
        !
      end do
      end do 
      end do
      dyn(i,j,ia,ja) = dyn(i,j,ia,ja) + ftmp
    end do 
    end do 
    end do 
    end do 
    !
    do i = 1,3
    do j = 1,3
    do ia = 1,PrimitiveCell%atom_number
    do ja = 1,PrimitiveCell%atom_number
      aatom = PrimitiveCell%atom_index(ia)
      batom = PrimitiveCell%atom_index(ja)
      amass = PrimitiveCell%exatom(aatom)%mass
      bmass = PrimitiveCell%exatom(batom)%mass
      dyn(i,j,ia,ja) = dyn(i,j,ia,ja) / sqrt( amass * bmass)
    end do
    end do
    end do
    end do
    !
    return 
  end subroutine calc_shortpart
  !
  !
  subroutine calc_longpart(PrimitiveCell,dyn,Wavevector,qhat,Nx,Ny,Nz)
    type(UnitCell),intent(in) :: PrimitiveCell
    complex(MYPRE),allocatable :: dyn(:,:,:,:)
    real(MYPRE) :: Wavevector(3),qhat(3)
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    !
    integer(MYIP) :: ia,ja,ka
    integer(MYIP) :: i,j
    real(MYPRE) :: zero_vector(3)
    complex(MYPRE) :: A_dd
    !
    zero_vector(:) = 0d0
    !
    !
    do ia = 1,PrimitiveCell%atom_number
    do ja = 1,PrimitiveCell%atom_number
    do i = 1,3
    do j = 1,3
      !
      call calc_sub_longpart(PrimitiveCell,A_dd,Wavevector,i,j,ia,ja,Nx,Ny,Nz)
      dyn(i,j,ia,ja) = dyn(i,j,ia,ja) + A_dd
      !
      if( i .eq. j ) then 
        do ka = 1,PrimitiveCell%atom_number
          call calc_sub_longpart(PrimitiveCell,A_dd,zero_vector,i,j,ia,ka,Nx,Ny,Nz)
          dyn(i,j,ia,ja) = dyn(i,j,ia,ja) - A_dd
        end do 
      end if
      !
    end do 
    end do 
    end do 
    end do 
  end subroutine calc_longpart
  !
  subroutine calc_Gamma_part
  end subroutine calc_Gamma_part
  !
  subroutine calc_sub_longpart(PrimitiveCell,A_dd,Wavevector,ix,jx,ia,ja,Nx,Ny,Nz)
    type(UnitCell),intent(in) :: PrimitiveCell
    complex(MYPRE) :: A_dd
    real(MYPRE),intent(in) :: Wavevector(3)
    integer(MYIP),intent(in) :: ix,jx,ia,ja
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2)
    !
    integer(MYIP) :: iy,jy 
    integer(MYIP) :: i,j 
    complex(MYPRE) :: A_sub
    real(MYPRE) :: zeu_a,zeu_b
    !
    A_dd = CMPLX(0d0,0d0,kind=kind(0d0))
    A_sub = 0
    !
    do iy = 1,3
    do jy = 1,3
      zeu_a = PrimitiveCell%Zeu(ix,iy,ia)
      zeu_b = PrimitiveCell%Zeu(jx,jy,ja)
      !
      A_sub = CMPLX(0.0_MYPRE,0.0_MYPRE,kind=MYPRE)
      call calc_reciprocal_space(PrimitiveCell,A_sub,Wavevector,iy,jy,ia,ja,Nx,Ny,Nz)
      !
      !print *, A_sub
      !
      A_dd = A_dd + zeu_a * A_sub * zeu_b 
    end do 
    end do 
    !
    return 
  end subroutine calc_sub_longpart 
  !
  !
  subroutine Calc_reciprocal_space(PrimitiveCell,A_sub,Wavevector,ix,jx,ia,ja,Nx,Ny,Nz)
    type(UnitCell),intent(in) :: PrimitiveCell
    real(MYPRE) :: Wavevector(3)
    complex(MYPRE) :: A_sub
    integer(MYIP) ,intent(in) :: ix,jx
    integer(MYIP) ,intent(in) :: ia,ja
    integer(MYIP) ,intent(in) :: Nx(2),Ny(2),Nz(2)
    real(MYPRE),allocatable :: Dynmat(:,:,:,:)
    !
    real(MYPRE) :: volume
    real(MYPRE) :: K_Vector(3)
    real(MYPRE) :: Zeta
    real(MYPRE) :: K_e_K
    real(MYPRE) :: Ewald_Gamma
    real(MYPRE) :: phasor
    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: i,j
    !
    real(MYPRE) :: geg
    !
    Ewald_Gamma = 1.0_MYPRE
    !
    call Calculate_volume_UnitCell(PrimitiveCell,volume)
    !
    do il = 2*Nx(1),2*Nx(2)
    do jl = 2*Ny(1),2*Ny(2)
    do kl = 2*Nz(1),2*Nz(2)
      !
      if(il .eq. 0 .and. jl .eq. 0 .and. kl .eq. 0) exit
      do i = 1,3
        K_Vector(i) = Wavevector(i) &
                    + PrimitiveCell%Recipro_Vector(i,1) * il & 
                    + PrimitiveCell%Recipro_Vector(i,2) * jl &
                    + PrimitiveCell%Recipro_Vector(i,3) * kl
      end do 
      !
      K_e_K = 0_MYPRE
      do i = 1,3
      do j = 1,3
        K_e_K = K_e_K + K_Vector(i) * PrimitiveCell%epsil(i,j) * K_Vector(j)
      end do 
      end do 
      !
      Zeta = 4 *  K_Vector(ix) * K_Vector(jx) / (K_e_K * volume)
      !
      phasor = 0.0_MYPRE
      do i = 1,3
        phasor = phasor &
               + K_Vector(i) &
               * (  &
               ( PrimitiveCell%Lattice_Vector(i,ix) *PrimitiveCell%Position(i,ix) ) &
               - ( PrimitiveCell%Lattice_Vector(i,jx) *PrimitiveCell%Position(i,jx) ) &
                  )
                
      end do
      !
      !
      A_sub = A_sub + CMPLX(cos(phasor),sin(phasor),kind=MYPRE) *  &
                     Zeta * &
                     EXP(-K_e_k/(4*Ewald_Gamma*Ewald_Gamma))
      !print *, "A_sub=", A_sub
      !             
    end do
    end do
    end do
    !
    return 
  end subroutine calc_reciprocal_space
  !
  !
  subroutine MakeMatrixHermite(Mat)
    complex(MYPRE),allocatable     :: Mat(:,:)
    integer(MYIP)                            :: i,j,k
    real(MYPRE)                   :: diff,difrel,dif1
    integer(MYIP)                            :: na
    
    na = size(Mat,1)

    diff  = 0.d0
    difrel= 0.d0
    do i = 1,na
      Mat(i,i) = CMPLX( DBLE(Mat(i,i)),0.d0,kind=kind(0d0))
      do j = 1,i - 1
        dif1 = abs(Mat(i,j)-CONJG(Mat(j,i)))
        if ( dif1 > diff .and. &
          max ( abs(Mat(i,j)), abs(Mat(j,i))) > 1.0d-6) then
          diff = dif1
          difrel=diff / min ( abs(Mat(i,j)), abs(Mat(j,i)))
        end if
        Mat(i,j) = 0.5d0* (Mat(i,j)+CONJG(Mat(j,i)))
        Mat(j,i) = CONJG(Mat(i,j))
      end do
    end do

  end subroutine MakeMatrixHermite
  
  

end module 
