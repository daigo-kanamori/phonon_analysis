program fft_1d
  use FFTW_MODULE
  implicit none
  !
  complex(kind(0d0)),allocatable ::  xout(:)
  complex(kind(0d0)),allocatable ::  fftw(:)
  double precision,allocatable :: power(:)
  double precision,allocatable :: phasor(:)
  double precision :: xtmp
  integer :: row,now
  integer :: i
  !
  !
  
  print *, "START FFT 1d"
  read *, row
  !
  allocate(xout(0:row))
  !
  do i = 0,row
      read *, now ,xtmp 
      xout(i) = CMPLX(xtmp,0,kind=kind(0d0))
  end do 
  !
  call fftw_1d(xout,fftw,FFTW_FORWARD)
  !
  allocate(power(row))
  allocate(phasor(row))
  !
  do i = 0,row
    power(i) = abs(fftw(i)) / DBLE(row)
    phasor(i) = atan2(aimag(fftw(i)),real(fftw(i)))
  end do 
  !
  do i = 0,row
    print *, i,", " ,power(i),"," ,phasor(i)
  end do 
  print *, "FFT DONE"
  !
  return 
end program fft_1d
