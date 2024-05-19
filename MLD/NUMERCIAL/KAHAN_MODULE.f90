module KAHAN_MODULE
  use SYSTEM_MODULE
  !
  implicit none 
  !
  type KahanSum
    real(MYPRE) :: sum
    real(MYPRE) :: cmps
  end type KahanSum
  !
  type Vector_kahan
    type(KahanSum),allocatable :: vector(:)
  end type Vector_kahan
  !
  type KahanSum_coplex
    type(KahanSum) :: real,imag
  end type
  !
  !
  contains 
  !
  !
  subroutine kahan_zero(kahan)
    type(KahanSum) :: kahan
    !
    kahan%sum = 0d0
    kahan%cmps = 0d0
    return 
  end subroutine kahan_zero
  !
  !
  subroutine add_kahan(kahan,x)
    type(KahanSum) :: kahan
    real(MYPRE),intent(in) ::  x
    !
    real(MYPRE) :: tmp1,tmp2
    !
    tmp1 = x - kahan%cmps
    tmp2 = kahan%sum + tmp1
    kahan%cmps = ( tmp2 - kahan%sum) - tmp1
    kahan%sum = tmp2
    !
    return 
  end subroutine add_kahan
  !
  subroutine Vector_kahan_Zero(veck)
    type(Vector_kahan) :: veck
    !
    integer(MYIP) :: i,N
    !
    N = size(veck%vector)
    do i = 1,N
      call kahan_zero(veck%vector(i))
    end do 
    return 
  end subroutine Vector_kahan_Zero
  !
  subroutine allocate_Vector_kahan(kVec,N,M)
    type(Vector_kahan) :: kVec
    integer(MYIP),intent(in) :: N,M
    !
    integer(MYIP) :: i
    !
    allocate(kVec%vector(N:M))
    do i = N,M
      call kahan_zero(kVec%vector(i))
    end do 
    return 
  end subroutine allocate_Vector_kahan
  !
  !
  subroutine Write_Vector_kahan_line(kahan,fd)
    type(Vector_kahan),intent(in) :: kahan
    integer(MYIP),intent(in) :: fd
    !
    integer(MYIP) :: i,N
    !
    N = size(kahan%vector)
    do i = 1, N-1
      write(fd,fmt='(E18.4 , A)',advance='no') kahan%vector(i)%sum, ","
    end do 
    !
    write(fd,fmt='(E18.4)')kahan%vector(N)%sum
    !
    return 
  end subroutine  Write_Vector_kahan_line
  !
  !
  subroutine kahan_comp_zero(kahan)
    type(KahanSum_coplex) :: kahan
    !
    call kahan_zero(kahan%real)
    call kahan_zero(kahan%imag)
    !
    return 
  end subroutine kahan_comp_zero
  !
  subroutine add_kahan_comp(kahan,x)
    type(KahanSum_coplex) :: kahan
    complex(MYPRE),intent(in) :: x
    !
    real(MYPRE) :: tmp
    !
    tmp = real(x)
    call add_kahan(kahan%real,tmp)
    tmp = aimag(x)
    call add_kahan(kahan%imag,tmp)
    !
    return 
  end subroutine add_kahan_comp
end module

