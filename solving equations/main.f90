program main
  implicit none
  real(kind=8) :: low, up, guess1, guess2
  real(kind=8), external :: f, g1, g2, g3, d, dd
  real(kind=8) :: res

  low = 0.0D0
  up = 1.0D0
  call bisect(f, res, low, up)
  write(*,"('Bisection first root: ',F10.7)") res
  
  low = -2.0D0
  up = -1.0D0
  call bisect(f, res, low, up)
  write(*,"('Bisection second root: ',F10.7)") res
  
  guess1 = 5.0D-1
  call fpi(g1, res, guess1)
  write(*,"('FPI first root: ',F10.7)") res
  
  guess1 = 5.0D-1
  call fpi(g2, res, guess1)
  write(*,"('FPI first root: ',F10.7)") res
  
  guess1 = -1.5D0
  call fpi(g3, res, guess1)
  write(*,"('FPI Second root: ',F10.7)") res
  
  guess1 = 1.0D0
  call newton(f, d, res, guess1)
  write(*,"('Newton first root: ',F10.7)") res
  
  guess1 = -1.0D0
  call newton(f, d, res, guess1)
  write(*,"('Newton second root: ',F10.7)") res
  
  guess1 = 1.0D0
  guess2 = 8.0D-1
  call secant(f, res, guess1, guess2)
  write(*,"('Secant first root: ',F10.7)") res

  guess1 = 1.0D0
  call aug_newton(f, d, dd, res, guess1)
  write(*,"('Aug Newton first root: ',F10.7)") res
  
  stop
end program

function f(x)
  implicit none
  real(kind=8) :: f, x
  f = x*x + x - 3.9D-1
end function

function g1(x)
  ! fpi function g(x) = x^2 - 0.39
  implicit none
  real(kind=8) :: g1, x
  g1 = 3.9D-1 - x*x
end function

function g2(x)
  ! fpi function g(x) = 0.39 / (x + 1)
  implicit none
  real(kind=8) :: g2, x
  g2 = 3.9D-1 / (x + 1.0D0)
end function  
    
function g3(x)
  ! fpi function g(x) = (0.39 - x) / x
  implicit none
  real(kind=8) :: g3, x
  g3 = (3.9D-1 - x) / x
end function
    
function d(x)
  implicit none
  real(kind=8) :: d, x
  d = x + x + 1.0D0
end function

function dd(x)
  implicit none
  real(kind=8) :: dd, x
  dd = 2.0D0
end function
    
subroutine bisect(func, res, low, up)
  ! func: function f(x)
  ! res: result
  ! low, up: upper and lower bounds
  ! steps: times of bisection, depending on the precision
  implicit none
  real(kind=8) :: res, low, up, mid
  real(kind=8), external :: func
  integer, parameter :: steps = 20
  integer :: i

  do i = 1, steps
    mid = (low + up) / 2.0D0
    if (func(low)*func(mid) < 0) then
      up = mid
    else
      low = mid
    end if
  end do
  res = (low + up) / 2.0D0

end subroutine

subroutine fpi(func, res, x1)
  ! func: fpi function g(x)
  ! res: result(may not meet convergence)
  ! x1: initial guess
  ! error: absolute termination condition
  ! maxsteps: max iteration times
  implicit none
  real(kind=8) :: res, x1, x2
  real(kind=8),external :: func
  real(kind=8), parameter :: error = 5.0D-7
  integer, parameter :: maxsteps = 129
  integer :: i
  
  i = 1
  do while(.true.)
    x2 = func(x1)
    if (abs(x2-x1) < error .or. i > maxsteps) exit
    i = i + 1
    x1 = x2
  end do
  res = x2

end subroutine

subroutine newton(func, deriv, res, x1)
  ! func: function f(x)
  ! deriv: derivative of f(x)
  ! res: result
  ! x1: initial guess
  implicit none
  real(kind=8) :: res, x1, x2
  real(kind=8),external :: func, deriv
  real(kind=8), parameter :: error = 5.0D-7
  integer, parameter :: maxsteps = 129
  integer :: i

  i = 1
  do while(.true.)
    x2 = x1 - func(x1) / deriv(x1)
    if (abs(x2-x1) < error .or. i > maxsteps) exit
    i = i + 1
    x1 = x2
  end do
  res = x2

end subroutine

subroutine secant(func, res, x1, x2)
  implicit none
  real(kind=8) :: res, x1, x2
  real(kind=8), external :: func
  real(kind=8), parameter :: error = 5.0D-7
  integer, parameter :: maxsteps = 129
  real(kind=8) :: f1, f2, x3
  integer :: i

  i = 1
  f1 = func(x1)
  f2 = func(x2)
  do while(.true.)
    x3 = x2 - (x2 - x1) * f2 / (f2 - f1)
    if (abs(x3-x2) < error .or. i > maxsteps) exit
    i = i + 1
    x1 = x2
    x2 = x3
    f1 = f2
    f2 = func(x2)
  end do
  res = x3

end subroutine

subroutine aug_newton(func, df, ddf, res, x1)
  implicit none
  real(kind=8) :: res, x1, x2
  real(kind=8), external :: func, df, ddf
  real(kind=8), parameter :: error = 5.0D-7
  integer, parameter :: maxsteps = 129
  integer :: i
  real(kind=8) :: temp1, temp2

  i = 1
  do while(.true.)
    temp1 = func(x1)
    temp2 = df(x1)
    x2 = x1 - temp1 * temp2 / (temp2 * temp2 - temp1 * ddf(x1))
    if (abs(x2-x1) < error .or. i > maxsteps) exit
    i = i + 1
    x1 = x2
  end do
  res = x2

end subroutine

subroutine regfal(func, res, low, up)
  implicit none
  real(kind=8) :: res, low, up, mid, fa, fb, fc
  real(kind=8), external :: func
  real(kind=8), parameter :: error = 5.0D-7
  integer, parameter :: maxsteps = 129
  integer :: i

  i = 1
  fa = func(low)
  fb = func(up)
  do while(.true.)
    mid = (up*fa - low*fb)/(fa - fb)
    fc = func(mid)
    if (abs(fc) < error .or. i > maxsteps) exit
    if (fa*fc < 0.0D0 ) then
      up = mid
      fb = fc
    else
      low = mid
      fa = fc
    end if
    i = i + 1
  end do
  res = mid

end subroutine
