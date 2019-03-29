module integration
  implicit none
  contains
  subroutine trapezoid(func,up,down,int,num)
    implicit none
    real(kind=8), external :: func
    real(kind=8) :: int, up, down
    integer :: num
    real(kind=8) :: step, tmp
    integer :: i
    
    step = (up - down) / num
    tmp = down
    int = 0.0D0
    
    do i = 1, num-1
      tmp = tmp + step
      int = int + func(tmp)
    end do
    int = 0.5D0 * step * (func(down) + func(up) + 2.0D0 * int)
    
  end subroutine
  subroutine Simpson(func,up,down,int,num)
    implicit none
    real(kind=8), external :: func
    real(kind=8) :: int, up, down
    integer :: num
    real(kind=8) :: step, tmp
    integer :: i
    
    step = 0.5D0 * (up - down) / num
    tmp = down
    int = 0.0D0
    
    do i = 1, num-1
      tmp = tmp + step
      int = int + 4.0D0 * func(tmp)
      tmp = tmp + step
      int = int + 2.0D0 * func(tmp)
    end do
    int = step * (func(down) + func(up) + 4.0D0 * func(tmp + step) + int) / 3.0D0
    
  end subroutine
  subroutine Romberg(func,up,down,int,tol)
    implicit none
    real(kind=8), external :: func
    real(kind=8) :: up, down, int, tol
    ! integer :: tag 
    ! tag == 1, return integral and write every items
    ! tag /= 1, only return integral
    integer, parameter :: maxcyc = 10
    real(kind=8) :: R(maxcyc,maxcyc) = 0.0D0
    real(kind=8) :: step
    integer :: i, j, k
    
    R(1,1) = 0.5D0 * (up - down) * (func(up) + func(down))
    outter: do i = 2, maxcyc
      step = (up - down) / 2.0D0**(i-1)
      do j = 1, 2**(i-2)
        R(i,1) = R(i,1) + func(down + (2*j - 1) * step)
      end do
      R(i,1) = R(i,1) * step + 0.5D0 * R(i-1,1)
      do k = 2, i
        R(i,k) = R(i,k-1) + (R(i,k-1) - R(i-1,k-1)) / (4.0D0**(k-1) - 1.0D0)
      end do
      if (abs(R(i,i)-R(i-1,i-1)) < tol ) then
          int = R(i,i)
          exit outter
      end if
    end do outter
    
    do i = 1, maxcyc
      do j = 1, i
        write(*,*) R(i,j)
      end do
    end do
    
  end subroutine
end module
program main
  use integration
  implicit none
  real(kind=8), external :: f
  real(kind=8) :: up = 2, down = 1
  integer :: num = 10
  real(kind=8) :: int
  real(kind=8), parameter :: pres_res = 2.0D0 * log(2.0D0) - 1.0D0
  
  call trapezoid(f,up,down,int,num)
  write(*,"(F15.13,'  ',F15.12)") int, int - pres_res
  
  call Simpson(f,up,down,int,num)
  write(*,"(F15.13,'  ',F15.12)") int, int - pres_res
  
  call Romberg(f,up,down,int,1.0D-7)
  write(*,"(F15.13,'  ',F15.12)") int, int - pres_res
  
end program
function f(x)
  implicit none
  real(kind=8) :: f, x
  f = log(x)
end function