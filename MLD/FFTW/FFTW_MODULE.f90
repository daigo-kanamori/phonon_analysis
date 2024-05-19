module FFTW_MODULE
  use SYSTEM_MODULE
  implicit none 
  !
  include 'fftw3.f'
  !
  contains 
  !
  subroutine fftw_1d(xin,xout,FFTW_TYPE)
    complex(MYPRE),allocatable :: xin(:),xout(:)
    integer(MYIP), intent(in) :: FFTW_TYPE
    integer(MYIP) :: N,M
    integer(MYIP) :: plan
    !
    M = lbound(xin,1)
    N = ubound(xin,1)
    if(allocated(xout)) then
      deallocate(xout)
    end if
    !
    allocate(xout(M:N))
    !
    call dfftw_plan_dft_1d(plan,(N-M+1),xin,xout,FFTW_TYPE,FFTW_ESTIMATE)
    call dfftw_execute_dft( plan, xin, xout )
    call dfftw_destroy_plan(plan)
    !
    return 
  end subroutine fftw_1d 
  !
  subroutine fftw_2d(xin,xout,FFTW_TYPE)
    complex(MYPRE),allocatable :: xin(:,:),xout(:,:)
    integer(MYIP),intent(in) :: FFTW_TYPE
    integer(MYIP) :: N(2),M(2)
    integer(MYIP) :: plan
    !
    N(1) = lbound(xin,1)
    N(2) = ubound(xin,1)
    !
    M(1) = lbound(xin,2)
    M(2) = ubound(xin,2)
    !
    if(allocated(xout)) then
      deallocate(xout)
    end if
    allocate(xout(N(1):N(2), M(1):M(2)))
    !
    call dfftw_plan_dft_2d(plan,(N(2)-N(1)+1),(M(2)-M(1)+1),xin,xout,FFTW_TYPE,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,xin,xout)
    call dfftw_destroy_plan(plan)
    !
    return 
  end subroutine fftw_2d 
  !
end module FFTW_MODULE
