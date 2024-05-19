module MD_MODULE
  use SYSTEM_MODULE
  use PARAMETER_MODULE
  use IOTK_MODULE
  use LATTICE_MODULE 
  use ENSEMBLE_MODULE
  use ENSEMBLE_DATA_MODULE
  use KAHAN_MODULE
  use MDCONDITION_MODULE
  use MDALGORITHM_MODULE
  implicit none
  !
  !
  contains 
  !
  !
  subroutine Update_Latiice_MD( Lattice ,CalculateForce,ide,nfree,now,free_number,il,jl,kl)
    interface
      subroutine CalculateForce(Lattice,now,Nx,Ny,Nz)
        use LATTICE_MODULE
        type(Lattice_MD)             :: Lattice
        integer(MYIP),intent(in)   :: Nx,Ny,Nz
        integer(MYIP),intent(in)   :: now
      end subroutine
    end interface 

    type(Lattice_MD) :: Lattice
    integer(MYIP) ::     ide
    integer(MYIP) ::     nfree
    integer(MYIP) ::     now
    integer(MYIP) :: free_number
    integer(MYIP) ::     il,jl,kl

    !
    !
    !$omp parallel 
    !$omp do private(free_number,nfree)
    do ide = 1,Lattice%sys%Ensemble_Number
      free_number = size(Lattice%sys%System(ide)%info%Mass%vector)
      do nfree = 1,free_number
        call Update_Displacement(Lattice%sys,ide,nfree)
      end do 
    end do 
    !$omp end do 
    !$omp end parallel
    !
    call Apply_BindCondition_Displacement(Lattice,now)
    !
    !$omp parallel 
    !$omp do private(jl,il)
    do kl = Lattice%Nz(1),Lattice%Nz(2)
    do jl = Lattice%Ny(1),Lattice%Ny(2)
    do il = Lattice%Nx(1),Lattice%Nx(2)
      call CalculateForce(Lattice,now,il,jl,kl)
    end do
    end do
    end do
    !$omp end do 
    !$omp end parallel
    !
    call Apply_BindCondition_Force(Lattice,now)
    !
    !$omp parallel 
    !$omp do private(jl,il)
    do ide = 1,Lattice%sys%Ensemble_Number
      free_number = size(Lattice%sys%System(ide)%info%Mass%vector)
      do nfree  = 1,free_number
        call Update_Velocity(Lattice%sys,ide,nfree)
      end do 
    end do 
    !$omp end do 
    !$omp end parallel
    !
    call Apply_BindCondition_Velocity(Lattice,now)
    !
    return 
  end subroutine Update_Latiice_MD
  !
  !
  !
  subroutine Check_Force(Lattice,fd)
    type(Lattice_MD),intent(in) :: Lattice
    type(KahanSum) :: ftmp(3)
    integer(MYIP) :: fd
    integer(MYIP) :: i,il,jl,kl,ia
    integer(MYIP) :: ide
    integer(MYIP) :: N
    integer(MYIP) :: idx
    !
    do i = 1,3
      call kahan_zero(ftmp(i))
    end do
    !
    !
    !
    do kl = Lattice%Nz(1),Lattice%Nz(2)
    do jl = Lattice%Ny(1),Lattice%Ny(2)
    do il = Lattice%Nx(1),Lattice%Nx(2)
      idx = Lattice%Index(il,jl,kl)
      N = Lattice%Cell(idx)%atom_number
      ide = Lattice%sys%Ensemble_list(il,jl,kl)
      do ia = 1,N
        do i = 1,3
          call add_kahan(ftmp(i),Lattice%sys%System(ide)%state%Force(1)%vector(3*(ia-1)+i)%sum)
        enddo
      enddo
    enddo
    enddo
    enddo 
    !
    write(fd,*) ftmp(1)%sum, "," , ftmp(2)%sum , "," , ftmp(3)%sum
    !
    return
  end subroutine
  !
end module MD_MODULE


