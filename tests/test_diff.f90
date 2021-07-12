! Created by  on 8/7/21.

program test_diff
    use diff
    implicit none
    call test_tridiag
    call test_sparse_tridiag
    contains
        subroutine test_tridiag
            implicit none
            integer :: N = 4, ierr
            real (kind = 8), dimension(:, :), allocatable :: result

            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 0
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 1
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)

            N = 2
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)
        end subroutine test_tridiag

        subroutine test_sparse_tridiag
            implicit none
            integer :: N, ierr
            real(kind=8), allocatable, dimension(:) :: val
            integer, dimension(:), allocatable :: indx, jndx

            N = 5
            call sparse_tridiag(dble(-2), dble(1), dble(-1), N, val, indx, jndx, ierr)

            Print *, val
            Print *, indx
            Print *, jndx

            deallocate(val, indx, jndx)


        end subroutine test_sparse_tridiag

end program test_diff