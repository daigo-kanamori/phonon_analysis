module FITTING_MODULE
  use VECTOR_MODULE
  implicit none
  contains
  !
  !=================================================================!
  !               LINEAR FITTING                                    !
  !=================================================================!
  !
  subroutine LinearFitting(x,y,a,b)
    double precision,allocatable,intent(in) :: x(:),y(:)
    double precision,intent(out) :: a,b
    double precision :: xavr,yavr,cov,xvari
    !    
    call CalculateAverage(y,yavr)
    call CalculateSD(x,xavr,xvari)
    cov = CalculateCovariance(x,y)
    a = cov / xvari
    b = yavr - a * xavr
    !
    return 
  end subroutine LinearFitting
  !
end module
