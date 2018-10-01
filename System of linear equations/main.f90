module global
  implicit none
  contains
  subroutine gauss_elim(mat,res)
    implicit none
    real(kind=8) :: mat(:,:)
    real(kind=8) :: res(:)
    integer :: ncol, nrow
    integer :: i, j, k
    real(kind=8) :: temp

    ncol = size(mat,1)
    nrow = size(mat,2)

    do i = 1, ncol - 1
      do j = i + 1, nrow
        temp = -mat(i,j)/mat(i,i)
        do k = i, ncol
          mat(k,j) = mat(k,j) + temp*mat(k,i)
        end do
      end do
    end do

    do i = nrow, 1, -1
      do j = i + 1, nrow
        mat(ncol,i) = mat(ncol,i) - mat(j,i)*res(j)
      end do
      res(i) = mat(ncol,i)/mat(i,i)
    end do
    
    !mat(nrow:ncol,nrow) = mat(nrow:ncol,nrow)/mat(nrow,nrow)
    !do i = nrow, 2, -1
    !  do j = 1, i - 1
    !    mat(ncol,j) = mat(ncol,j) - mat(ncol,i)*mat(i,j)/mat(i,i)
    !    mat(i,j) = 0
    !  end do
    !end do
    !write(*,"(4F5.1,/,4F5.1,/,4F5.1)") ((mat(i,j),i=1,4),j=1,3)
    !do i = 1, nrow
    !  mat(ncol,i) = mat(ncol,i)/mat(i,i)
    !  mat(i,i) = 1
    !end do
    !write(*,"(4F5.1,/,4F5.1,/,4F5.1)") ((mat(i,j),i=1,4),j=1,3)
  end subroutine
end module
program main
  use lapack95
  use global
  implicit none
  real(kind=8) :: mat(4,3) = reshape( (/ 2,-2,-1,-2,4,1,-2,1,-2,1,-1,-3 /), (/ 4,3 /) )
  real(kind=8) :: res(3)
  real(kind=8) :: A(3,3) = reshape( (/ 2,4,-2,-2,1,1,-1,-2,-1 /), (/ 3,3 /) )
  real(kind=8) :: B(3) = (/ -2, 1, -3 /)
  real(kind=8) :: ipiv(3)
  integer :: i, j, info

  len_res = size(res)
  write(*,"(4F5.1,/,4F5.1,/,4F5.1)") ((mat(i,j),i=1,4),j=1,3)
  call gauss_elim(mat,res)
  do i = 1, size(res)
    write(*,"('x',I2,' = ',F7.3)") i, res(i)
  end do
  call dgesv( 3, 1, A, 3, ipiv, B, 3, info )
  do i = 1, size(B)
    write(*,"('x',I2,' = ',F7.3)") i, B(i)
  end do
end program
