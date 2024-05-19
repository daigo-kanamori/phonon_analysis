module IOTK_MODULE
  use SYSTEM_MODULE
  implicit none
  !
  integer(MYIP):: FILE_DISCRIPT_NUMBERS = 10
  integer(MYIP),parameter :: STDOUT = 6
  integer(MYIP),parameter :: STDIN  = 5
  integer(MYIP),parameter :: STDERR = 0
  contains
  !
  subroutine open4read(filename,fd)
    integer(MYIP)         ::       fd
    character(100)   ::       filename
    !
    fd = FILE_DISCRIPT_NUMBERS
    FILE_DISCRIPT_NUMBERS = FILE_DISCRIPT_NUMBERS+1_MYIP
    filename = TRIM(filename)
    open(fd,file=filename,status='old')
    !
    return 
  end subroutine open4read
  !
  subroutine open4write(filename,fd)
    integer(MYIP)         ::       fd
    character(100)   ::       filename
    !
    fd = FILE_DISCRIPT_NUMBERS
    FILE_DISCRIPT_NUMBERS = FILE_DISCRIPT_NUMBERS+1_MYIP
    filename = TRIM(filename)
    open(fd,file=filename,status='replace')
    !
    return 
  end subroutine open4write

  subroutine closefile(fld)
    integer(MYIP) :: fld
    close(fld)
  end subroutine closefile
  !
  subroutine close_all
    integer(MYIP) :: i
    do i = 10_MYIP,FILE_DISCRIPT_NUMBERS-1_MYIP
      close(i)
    end do
  end subroutine close_all
  !
  !ファイルの行数を数えるカウントルーチーン
  !
  subroutine count_file(fd,N)
    integer(MYIP),intent(in) :: fd
    integer(MYIP),intent(out) :: N
    N = 0
    do 
      read(fd,*,end = 666)
      N = N + 1_MYIP
      end do 

    666 continue
    rewind(fd)
    !
    return 
  end subroutine 
  !
  !
  !
  subroutine Write_Vector_d_line(x,fd)
    real(MYPRE),allocatable,intent(in) :: x(:)
    integer(MYIP),intent(in) :: fd
    integer(MYIP) :: i
    !
    do i = int(lbound(x,1),MYIP),int(ubound(x,1)-1,MYIP)
      write(fd,fmt='(E20.7, A)',advance='no')x(i) ,", "
    end do 
    write(fd,fmt='(E20.7)')x(i)
    !
    return 
  end subroutine Write_Vector_d_line
  !
  !
  !
end module IOTK_MODULE


