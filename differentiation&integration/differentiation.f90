module differentiation
  implicit none
  contains
  subroutine derivative_order1(func,deriv,x,step)
    implicit none
    real(kind=8), external :: func
    real(kind=8) :: deriv, x, step
    deriv = (func(x+step) - func(x)) / step
  end subroutine
  subroutine derivative_order2(func,deriv,x,step)
    implicit none
    real(kind=8), external :: func
    real(kind=8) :: deriv, x, step
    deriv = 0.5D0 * (func(x+step) - func(x-step)) / step
  end subroutine
end module

program main 
  use differentiation
  implicit none
  real(kind=8), parameter :: a = 1.0
  real(kind=8), parameter :: exp1 = exp(a)
  real(kind=8) :: steps(10) = [1.0D-1, 1.0D-2, 1.0D-3, 1.0D-4, 1.0D-5, 1.0D-6, 1.0D-7, 1.0D-8, 1.0D-9, 1.0D-10]
  real(kind=8) :: deriv
  real(kind=8), external :: f
  integer :: i
  
  do i = 1, 10
   call derivative_order1(f,deriv,a,steps(i))
   write(*,"(D7.1,'  ',F15.13,'  ',F16.13)") steps(i), deriv, deriv-exp1
  end do

  write(*,*)
  
  do i = 1, 10
   call derivative_order2(f,deriv,a,steps(i))
   write(*,"(D7.1,'  ',F15.13,'  ',F16.13)") steps(i), deriv, deriv-exp1
  end do

end program

function f(x)
  implicit none
  real(kind=8) :: f, x
  f = exp(x)
end function