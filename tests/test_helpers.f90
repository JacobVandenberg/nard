! Created by  on 10/7/21.

program test_helpers
    use helpers
    use diff
    use ogpf
    implicit none
    Print *, 'speye'
    call test_speye
    Print *, 'sparse_kron'
    call test_sparse_kron
    Print *, 'sparse_update'
    call test_sparse_update
    Print *, 'sparse_multiply'
    call test_sparse_multiply
    Print *, "sparse_add"
    call test_sparse_add

    contains
        subroutine test_sparse_kron
            implicit none
            integer :: N, n1, m1, n2, m2, ierr, nout, mout
            integer, allocatable, dimension(:) :: indx1, jndx1, indx2, jndx2, indxout, jndxout
            real(kind=8), allocatable, dimension(:) :: val1, val2, valout

            N = 3
            n1 = 3
            m1 = 3
            n2 = 2
            m2 = 2
            call sparse_tridiag(dble(-2), dble(1), dble(1), N, val1, indx1, jndx1, ierr)
            Print *, ierr
            Print *, val1, indx1, jndx1
            call speye(2, val2, indx2, jndx2, ierr)
            Print *, ierr
            Print *, val2, indx2, jndx2

            call sparse_kron(val1, indx1, jndx1, n1, m1, val2, indx2, jndx2, n2, m2, valout, indxout, jndxout, nout, mout, ierr)
            Print *, ierr
            Print *, valout
            Print *,indxout
            Print *,jndxout
            deallocate(indx1, indx2, jndx1, jndx2, val1, val2, indxout, jndxout, valout)

        end subroutine test_sparse_kron

        subroutine test_speye
            implicit none
            integer :: N, ierr
            real (kind=8), dimension(:), allocatable :: val
            integer, dimension(:), allocatable :: indx, jndx

            N = 4

            call speye(N, val, indx, jndx, ierr)
            Print *, val, indx, jndx
            deallocate(val, indx, jndx)

        end subroutine test_speye

        subroutine test_sparse_update
            implicit none
            integer :: N, ierr
            real (kind=8), dimension(:), allocatable :: val
            integer, dimension(:), allocatable :: indx, jndx

            N = 4
            call sparse_tridiag(dble(-2), dble(1), dble(1), N, val, indx, jndx, ierr)
            call sparse_update(val, indx, jndx, dble(-3), 2, 2, ierr)
            Print *, val
            Print *, indx
            Print *, jndx
            call sparse_update(val, indx, jndx, dble(-3), 1, 4, ierr)
            Print *, val
            Print *, indx
            Print *, jndx

            deallocate (val, indx, jndx)

        end subroutine test_sparse_update

        subroutine test_sparse_multiply
            implicit none
            type(gpf) :: gp
            real (kind=8), dimension(:), allocatable :: x, vals, y_plot
            real (kind=8), dimension(:, :), allocatable :: y, DDy
            integer, dimension(:), allocatable :: indx, jndx
            integer :: N, ierr

            N = 1000

            allocate(x(N), y(N, 2), DDy(N, 2), y_plot(N))

            call my_linspace(dble(-1), dble(1), N, x)

            y(:, 1) = cos(dble(10) * x)
            y(:, 2) = sin(dble(10) * x)

            call double_diff_FD(x, vals, indx, jndx, ierr)
            call sparse_multiply(vals, indx, jndx, y, DDy)

            DDy(1, :) = 0
            DDy(N, :) = 0
            y_plot = DDy(:, 2)
            call gp%plot(x, y_plot)

            deallocate(x, vals, y, DDy, indx, jndx)
        end subroutine test_sparse_multiply

        subroutine test_sparse_add
            implicit none
            real (kind=8), allocatable, dimension(:) :: val1, val2, valout
            integer, allocatable, dimension(:) :: indx1, indx2, indxout, jndx1, jndx2, jndxout
            integer :: N, ierr

            N = 4
            call speye(N, val1, indx1, jndx1, ierr)
            call speye(N, val2, indx2, jndx2, ierr)
            call sparse_update(val1, indx1, jndx1, dble(10), 1, 2, ierr)

            call sparse_add(val1, indx1, jndx1, val2, indx2, jndx2, valout, indxout, jndxout, ierr)
            Print *, valout
            Print *, indxout
            Print *, jndxout
        end subroutine test_sparse_add

end program test_helpers