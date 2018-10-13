module global
  implicit none
  contains
  subroutine BezierCurve(BzPts,npts,P0,P1,P2,P3)
    implicit none
    real(kind=8),intent(in) :: P0(2),P1(2)
    real(kind=8),intent(in),optional :: P2(2),P3(2)
    integer,intent(in) :: npts
    real(kind=8),intent(out) :: BzPts(:,:)
    real(kind=8) :: step
    real(kind=8) :: t, t2, t3, tp, tp2, tp3
    integer :: i
    logical :: bl1, bl2 ! determine the existence of P2, P3
    
    bl1 = present(P2)
    bl2 = present(P3)
    step = 1.0D0 / npts
    t = 0.0D0
    if (bl1 .and. bl2) then
      ! cubic bezier curve when four points passed
      BzPts(npts+1,:) = P3
      do i = 1, npts
        tp = (1.0D0 - t)
        tp2 = tp * tp
        tp3 = tp2 * tp
        t2 = t * t
        t3 = t2 * t
        BzPts(i,:) = P0*tp3 + 3.0D0*P1*t*tp2 + 3.0D0*P2*tp*t2 + P3*t3
        t = t + step
      end do
    else if (bl1 .or. bl2) then
      ! quadratic bezier curve when three points passed
      BzPts(npts+1,:) = P2
      do i = 1, npts
        tp = (1.0D0 - t)
        tp2 = tp * tp
        t2 = t * t
        BzPts(i,:) = P0*tp2 + 2.0D0*P1*t*tp + P2*t2
        t = t + step
      end do
    else
      ! linear bezier curve when two points passed
      BzPts(npts+1,:) = P1
      do i = 1, npts
        tp = (1.0D0 - t)
        BzPts(i,:) = P0*tp + P1*t
        t = t + step
      end do
    end if
    
  end subroutine
end module
program main
  use global
  implicit none
  real(kind=8) :: P0(2) = (/ 0.0D0, 0.0D0 /)
  real(kind=8) :: P1(2) = (/ 1.0D0, 1.0D0 /)
  real(kind=8) :: P2(2) = (/ 2.0D0,-1.0D0 /)
  real(kind=8) :: P3(2) = (/ 3.0D0, 0.0D0 /)
  integer :: n = 50
  real(kind=8),allocatable :: B(:,:)
  integer :: i
  
  allocate(B(n+1,2))
  call BezierCurve(B,n,P0,P1)
  write(*,*) "Linear Bezier Curve"
  do i = 1, size(B,1)
    write(*,"(F6.3,'  ',F6.3)") B(i,:)
  end do
  call BezierCurve(B,n,P0,P1,P2)
  write(*,*) "Quadratic Bezier Curve"
  do i = 1, size(B,1)
    write(*,"(F6.3,'  ',F6.3)") B(i,:)
  end do
  call BezierCurve(B,n,P0,P1,P2,P3)
  write(*,*) "Cubic Bezier Curve"
  do i = 1, size(B,1)
    write(*,"(F6.3,'  ',F6.3)") B(i,:)
  end do
  
end program
