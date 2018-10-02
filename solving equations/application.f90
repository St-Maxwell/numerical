module global
  implicit none
  !oxygen
  real(kind=8),parameter :: p = 1.5D1 ! atm
  real(kind=8),parameter :: R = 8.20578D-2
  real(kind=8),parameter :: T = 3.20D2 ! K
  real(kind=8),parameter :: n = 1.0D0 ! mol
  real(kind=8),parameter :: a = 1.36D0
  real(kind=8),parameter :: b = 3.183D-3
  !benzene vapor
  !real(kind=8),parameter :: p = 2.0D1 ! atm
  !real(kind=8),parameter :: R = 8.20578D-2
  !real(kind=8),parameter :: T = 7.0D2 ! K
  !real(kind=8),parameter :: n = 1.0D0 ! mol
  !real(kind=8),parameter :: a = 1.8D1
  !real(kind=8),parameter :: b = 1.154D-1
end module
program main
  implicit none
  real(kind=8) :: guess
  real(kind=8), external :: vdW, d, ideal
  real(kind=8) :: res
  guess = ideal()
  write(*,"(F4.2)") guess
  call newton(vdW, d, res, guess)
  write(*,"(F4.2)") res
  stop
end program

function vdW(V)
  use global
  implicit none
  real(kind=8) :: vdW, V
  vdW = p*V + n*n*a/V * (1.0D0 - n*b/V) - n*(R*T + b*p)
end function
    
function d(V)
  use global
  implicit none
  real(kind=8) :: d, V
  d = p - n*n*a/(V*V) * (1.0D0 - 2.0D0*n*b/V)
end function   

function ideal()
  use global
  implicit none
  real(kind=8) :: ideal
  ideal = n*R*T/p
end function
    
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