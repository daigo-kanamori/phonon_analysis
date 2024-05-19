module ZHEEV_MODULE
  use SYSTEM_MODULE
  implicit none

  !!LAPCCK用の計算パラメータと変数
  character,parameter ::                         JOBZ = 'V'
  character,parameter ::                         UPLO = 'U'
  integer ::                                     INFO, LDA, LWORK,Norder
  complex(MYPRE),allocatable    ::               WORK(:) 
  complex(MYPRE), allocatable ::               RWORK(:) 
  contains 
  subroutine allocate_ZHEEV(Matrix,EigenVector,N)
    integer                       :: N
    complex(MYPRE)  , allocatable :: Matrix(:,:)
    real(MYPRE), allocatable :: EigenVector(:)
    Norder = N
    LWORK  = 2*N - 1 
    LDA = N
    allocate(Matrix(N,N))
    allocate(EigenVector(N))
    allocate(WORK(LWORK))
    allocate(RWORK(3*Norder))
  end subroutine 


  subroutine exec_ZHEEV(Matrix,EigenVector)
   complex(MYPRE),allocatable :: Matrix(:,:)
   real(MYPRE),allocatable :: EigenVector(:)

   call zheev(JOBZ,UPLO,Norder,Matrix,LDA,EigenVector,WORK,LWORK,RWORK,INFO)
  end subroutine

  subroutine free_zheev(Matrix,EigenVector)
    complex(MYPRE) , allocatable :: Matrix(:,:)
    complex(MYPRE) , allocatable :: EigenVector(:)
    deallocate(Matrix)
    deallocate(EigenVector)
    deallocate(WORK)
    deallocate(RWORK)
  end subroutine

end module ZHEEV_MODULE
