program fft_2d
  use FFTW_MODULE
  implicit none
  !
  complex(kind(0d0)),allocatable ::  xout(:,:)
  complex(kind(0d0)),allocatable ::  fftw(:,:)
  double precision,allocatable :: power(:,:)
  double precision,allocatable :: phasor(:,:)
  double precision :: xtmp
  integer :: grid(2)
  integer :: gtmp(2)
  integer :: i,j
  !
  !
  !
  print *, "START FFT 2d"
  read *, grid(1),grid(2)
  print *,"FFT GRID : 1dim = ",grid(1), ", 2dim = ",grid(2)
  allocate(xout(0:grid(1),0:grid(2)))
  !
  do i = 1,grid(1)
  do j = 1,grid(2)
    read *, gtmp(1),gtmp(2),xtmp
    xout(gtmp(1),gtmp(2)) = CMPLX(xtmp,0d0,kind=kind(0d0))
  end do 
  end do 
  !
  call fftw_2d(xout,fftw,FFTW_FORWARD)
  !
  allocate(power(0:grid(1),0:grid(2)))
  allocate(phasor(0:grid(1),0:grid(2)))
  !
  do i = 0,grid(2)
  do j = 0,grid(1)
    power(j,i) = abs(fftw(j,i)) / DBLE((1+grid(1))*(grid(2)+1))
    phasor(j,i) = atan2(aimag(fftw(j,i)),real(fftw(j,i)))
  end do 
  end do 
  !
  !
  do i = 0,grid(2)
  do j = 0,grid(1)
    print*, j,", ",i,", ",power(j,i),", ",phasor(j,i)
  end do 
  end do 
  print *, "FFT 2D DONE"
  !
  return 
end program 
