module VECTOR_MODULE
  use SYSTEM_MODULE
  !
  implicit none
  !
  type :: Vector_d
    real(MYPRE),allocatable :: Vector(:)
  end type Vector_d
  !
  type :: Matrix_d
    real(MYPRE),allocatable :: Matrix(:,:)
  end type Matrix_d
  !
  contains 

  subroutine allocate_Vector_d(vector,N,M)
    type(Vector_d)            :: vector
    integer(MYIP),intent(in)        :: N,M
    allocate(vector%vector(N:M))
    !
    return 
  end subroutine
  !
  subroutine deallocate_Vector_d(vector)
    type(Vector_d)            :: vector
    deallocate(vector%vector)
    !
    return 
  end subroutine

  subroutine allocate_Matrix_d(Matrix,N,M)
    type(Matrix_d)            :: Matrix
    integer(MYIP),intent(in)        :: N,M
    allocate(Matrix%Matrix(N,M))
    !
    return 
  end subroutine

  subroutine Vector_Zero(vector)
    type(Vector_d) :: vector
    integer(MYIP) :: i
    do i = 1,int(size(vector%vector),MYIP)
      vector%vector(i) =  0.0d0
    end do 
    !
    return 
  end subroutine 

  real(MYPRE) function DotProduct(vector1,vector2)
    real(MYPRE),allocatable :: vector1(:),vector2(:)
    integer(MYIP)                      :: N,i

    DotProduct = 0d0
    N = int(min(size(vector1,1), size(vector2,1)),MYIP)
    do i = int(lbound(vector1,1),MYIP),int(ubound(vector1,1),MYIP)
      DotProduct = DotProduct + vector1(i) * vector2(i)
    end do 
  end function 

  function CrossProduct(vector1,vector2) result(res)
    real(MYPRE),allocatable,intent(in)  :: vector1(:),vector2(:)
    real(MYPRE),allocatable :: res(:)
    ! 
    integer(MYIP)                       :: i
    res = vector1
    !allocate(res(3))
    res(:) = 0d0
    do i = int(lbound(res,1),MYIP),int(ubound(res,1),MYIP)
      res(i) = vector1(modulo(i+1,3)) * vector2(modulo(i+2,3)) & 
                        - vector1(modulo(i+2,3)) * vector2(modulo(i+1,3))
    end do 
    !
    return 
  end function CrossProduct
  !
  subroutine CalculateAverage(vector,avr)
    real(MYPRE),allocatable,intent(in) :: vector(:)
    real(MYPRE),intent(out) :: avr
    !
    integer(MYIP)                      :: N
    integer(MYIP)                      :: i
    !
    avr = 0d0
    !
    N = int(size(vector),MYIP)
    do i = int(lbound(vector,1),MYIP),int(ubound(vector,1),MYIP)
      avr = avr + vector(i) 
    end do 
    avr = avr / DBLE(N)
    !
    return 
  end subroutine
  !
  subroutine CalculateSD(vector,avr,variance)
    real(MYPRE),allocatable :: vector(:),vector2(:)
    real(MYPRE),intent(out)             :: avr,variance
    !
    integer(MYIP)                      :: i
    !
    avr = 0d0
    variance = 0d0
    call CalculateAverage(vector,avr)
    vector2 = vector
    !
    do i = int(lbound(vector2,1),MYIP),int(ubound(vector2,1),MYIP)
      vector2(i) = (vector2(i) - avr) * (vector2(i) - avr)
    end do 
    call CalculateAverage(vector2,variance)
    !
    return 
  end subroutine  CalculateSD
  !
  real(MYPRE) function CalculateCovariance(dat1,dat2)
    real(MYPRE),allocatable,intent(in) :: dat1(:),dat2(:)
    !
    real(MYPRE) :: cov
    real(MYPRE) :: avr1,avr2
    integer(MYIP) :: i,N
    !
    cov = 0d0
    call CalculateAverage(dat1,avr1)
    call CalculateAverage(dat2,avr2)
    !
    N = int(size(dat1),MYIP)
    do i = int(lbound(dat1,1),MYIP),int(ubound(dat1,1),MYIP)
      cov = cov + (dat1(i) - avr1) * (dat2(i) - avr2)
    end do
    cov = cov / DBLE(N)
    CalculateCovariance = cov
    !
    return 
  end function CalculateCovariance

end module VECTOR_MODULE 
