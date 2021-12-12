! Created by  on 13/7/21.

program test_sparse_matrices
    use sparse_matrices
    use ogpf
    use diff
    implicit none
    Print *, 'speye'
    call test_speye
    Print *, 'set_value'
    call test_set_value
    Print *, 'get_value'
    call test_get_value
    Print *, 'sparse_kron'
    call test_sparse_kron
    Print *, 'coo_multiply'
    call test_coo_multiply
    Print *, 'sparse_add'
    call test_sparse_add
    Print *, 'coo_to_csr'
    call test_coo_to_csr
    Print *, 'csr_multiply'
    call test_csr_multiply
    Print *, 'test_sparse_lu'
    call test_sparse_lu
    Print *, 'test_insert_block'
    call test_insert_block
    Print *, 'test_scale_rows'
    call test_scale_rows
    Print *, 'test_sparse_direct_solve'
    call test_sparse_direct_solve

    contains
        subroutine test_sparse_kron
            implicit none
            type (coo_matrix) :: mat1, mat2, matout
            integer (ip) :: ierr

            mat1 = speye(2_ip, ierr)
            call mat1%set_value(1_ip, 2_ip, dble(40), ierr)
            mat2 = speye(2_ip, ierr)

            matout = sparse_kron(mat1, mat2, ierr)
            Print *, matout%vals - (/1, 1, 1, 1, 40, 40/)
            Print *, matout%indx - (/1, 2, 3, 4, 1, 2/)
            Print *, matout%jndx - (/1, 2, 3, 4, 3, 4/)

        end subroutine test_sparse_kron

        subroutine test_speye
            implicit none
            type(coo_matrix) :: test_mat
            integer (ip) :: ierr

            test_mat = speye(4_ip, ierr)
            Print *, test_mat%vals - dble ( (/1, 1, 1, 1/) )
            Print *, test_mat%indx - (/1, 2, 3, 4/)
            Print *, test_mat%jndx - (/1, 2, 3, 4/)
            Print *, test_mat%n - 4
            Print *, test_mat%m - 4

        end subroutine test_speye

        subroutine test_set_value
            implicit none
            type(coo_matrix) :: test_mat
            integer (ip) :: ierr

            test_mat = speye(4_ip, ierr)
            call test_mat%set_value(1_ip, 2_ip, dble(40), ierr)
            Print *, test_mat%vals - dble ( (/1, 1, 1, 1, 40/) )
            Print *, test_mat%indx -  ( (/1, 2, 3, 4, 1/) )
            Print *, test_mat%jndx -  ( (/1, 2, 3, 4, 2/) )

        end subroutine test_set_value

        subroutine test_get_value
            implicit none
            type(coo_matrix) :: test_mat
            integer (ip) :: ierr

            test_mat = speye(4_ip, ierr)
            call test_mat%set_value(1_ip, 2_ip, dble(40), ierr)

            Print *, dble(1) - test_mat%get_value(2_ip, 2_ip)
            Print *, dble(40) - test_mat%get_value(1_ip, 2_ip)

        end subroutine test_get_value

        subroutine test_coo_multiply
            implicit none
            type(coo_matrix) :: test_mat
            integer (ip) :: ierr
            real (rp), dimension(:, :), allocatable :: x, y

            test_mat = speye(4_ip, ierr)
            test_mat%vals = dble((/1, 2, 3, 4/))
            call test_mat%set_value(1_ip, 2_ip, dble(40), ierr)

            allocate( x(4, 2), y(4, 2) )
            x = reshape((/1, 1, 1, 1, 2, 2, 2, 2/), (/4, 2/))
            call coo_multiply(test_mat, x, y)
            Print *, y(:, 1) - dble((/41, 2, 3, 4/))
            Print *, y(:, 2) - dble((/82, 4, 6, 8/))

        end subroutine test_coo_multiply

        subroutine test_sparse_add
            implicit none
            type (coo_matrix) :: mat1, mat2
            integer (ip) :: ierr

            mat1 = speye(4_ip, ierr)
            mat2 = speye(4_ip, ierr)
            call mat1%set_value(1_ip, 2_ip, dble(40), ierr)
            call sparse_add(mat1, mat2, ierr)
            Print *, mat1%vals - dble((/2, 40, 2, 2, 2/))
            Print *, mat1%indx - (/1, 1, 2, 3, 4/)
            Print *, mat1%jndx - (/1, 2, 2, 3, 4/)


        end subroutine test_sparse_add

        subroutine test_coo_to_csr
            implicit none
            type (coo_matrix) :: mat1
            type (csr_matrix) :: matout
            integer (ip) :: ierr

            mat1 = speye(4_ip, ierr)
            call mat1%set_value(1_ip, 4_ip, dble(40), ierr)
            call mat1%set_value(2_ip, 4_ip, dble(40), ierr)
            call mat1%set_value(3_ip, 4_ip, dble(40), ierr)
            matout = coo_to_csr(mat1, ierr)
            Print *, matout%vals - dble( (/1, 40, 1, 40, 1, 40, 1/) )
            Print *, matout%rows - (/1, 3, 5, 7, 8/)
            Print *, matout%jndx - (/1, 4, 2, 4, 3, 4, 4/)

        end subroutine test_coo_to_csr

        subroutine test_csr_multiply
            implicit none
            type(coo_matrix) :: test_mat
            type (csr_matrix) :: test_mat_csr
            integer (ip) :: ierr
            real (rp), dimension(:, :), allocatable :: x, y

            test_mat = speye(4_ip, ierr)
            test_mat%vals = dble((/1, 2, 3, 4/))
            call test_mat%set_value(1_ip, 2_ip, dble(40), ierr)

            test_mat_csr = coo_to_csr(test_mat, ierr)

            allocate( x(4, 2), y(4, 2) )
            x = reshape((/1, 1, 1, 1, 2, 2, 2, 2/), (/4, 2/))
            call csr_multiply(test_mat_csr, x, y)
            Print *, y(:, 1) - dble((/41, 2, 3, 4/))
            Print *, y(:, 2) - dble((/82, 4, 6, 8/))
        end subroutine test_csr_multiply

        subroutine test_sparse_lu
            implicit none
            type (coo_matrix) :: mat
            type (sparse_lu_decomposition) :: lu_mat
            real (rp), dimension(:), allocatable :: x, y, fx
            integer (ip) :: N, ierr
            type (gpf) :: gp

            N = 50_ip
            x = linspace(dble(0), dble(N-1), N)
            mat = double_diff_FD(x, 0_ip, ierr)
            Print *, mat%vals
            call mat%set_value(N, N-1_ip, dble(0.0001), ierr)
            call mat%set_value(N, N, dble(1), ierr)
            lu_mat = sparse_lu(mat, ierr)

            allocate (fx(N))
            fx(:) = cos(x * 2 * 3.14159265 / dble(N))
            fx(N) = 0

            y = sparse_lu_solve(lu_mat, fx, ierr)
            Print *, fx, y
            call gp%plot(x, reshape(y, (/N/)))

        end subroutine test_sparse_lu

        subroutine test_insert_block
            implicit none
            type (coo_matrix) :: mat1, mat2

            integer (ip) :: ierr

            mat1 = speye(4_ip, ierr)
            mat2 = speye(2_ip, ierr)
            call mat1%insert_block(mat2, 2_ip, 0_ip, ierr)
            Print *, mat1%vals! - real( (/ 1, 1, 1, 1 /), kind=rp )
            Print *, mat1%indx! - int( (/1, 2, 2, 3/), kind=ip )
            Print *, mat1%jndx! - int( (/1, 2, 3, 4/), kind=ip )
            Print *, mat1%m! - 3_ip
            Print *, mat1%n! - 4_ip

        end subroutine test_insert_block

        subroutine test_scale_rows
            implicit none
            type (coo_matrix) :: mat1
            type (csr_matrix) :: mat2, mat3
            integer (ip) :: ierr
            real (rp), dimension(6) :: factors

            mat1 = speye(6_ip, ierr)
            call mat1%set_value( 1_ip, 6_ip, 2.0_rp, ierr)
            mat2 = coo_to_csr(mat1, ierr)

            factors = real((/ 6, 2, 3, 4, 5, 6 /), kind=rp)

            call csr_scale_rows(mat2, factors)
            call coo_scale_rows(mat1, factors)

            Print *, mat1%vals
            mat3 = copy_csr_matrix(mat2, ierr)

        end subroutine test_scale_rows

        subroutine test_sparse_direct_solve
            implicit none
            type (coo_matrix) :: mat1
            type (csr_matrix) :: mat2
            integer (ip) :: ierr
            real (rp), allocatable, dimension(:, :) :: x
            real (rp), allocatable, dimension(:) :: b

            mat1 = speye(6_ip, ierr)
            call mat1%set_value( 1_ip, 6_ip, 2.0_rp, ierr)

            allocate( b(6_ip) )
            b = (/ 3, 1, 1, 1, 1, 1 /)
            allocate (x(6_ip, 1_ip))

            x(:, 1) = sparse_direct_solve(coo_to_csr(mat1, ierr), b, ierr);
            Print *, x
        end subroutine test_sparse_direct_solve



end program test_sparse_matrices