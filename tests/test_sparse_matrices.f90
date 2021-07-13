! Created by  on 13/7/21.

program test_sparse_matrices
    use sparse_matrices
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

    contains
        subroutine test_sparse_kron
            implicit none
            type (coo_matrix) :: mat1, mat2, matout
            integer :: ierr

            mat1 = speye(2, ierr)
            call mat1%set_value(1, 2, dble(40), ierr)
            mat2 = speye(2, ierr)

            matout = sparse_kron(mat1, mat2, ierr)
            Print *, matout%vals - (/1, 1, 1, 1, 40, 40/)
            Print *, matout%indx - (/1, 2, 3, 4, 1, 2/)
            Print *, matout%jndx - (/1, 2, 3, 4, 3, 4/)

        end subroutine test_sparse_kron

        subroutine test_speye
            implicit none
            type(coo_matrix) :: test_mat
            integer :: ierr

            test_mat = speye(4, ierr)
            Print *, test_mat%vals - dble ( (/1, 1, 1, 1/) )
            Print *, test_mat%indx - (/1, 2, 3, 4/)
            Print *, test_mat%jndx - (/1, 2, 3, 4/)
            Print *, test_mat%n - 4
            Print *, test_mat%m - 4

        end subroutine test_speye

        subroutine test_set_value
            implicit none
            type(coo_matrix) :: test_mat
            integer :: ierr

            test_mat = speye(4, ierr)
            call test_mat%set_value(1, 2, dble(40), ierr)
            Print *, test_mat%vals - dble ( (/1, 1, 1, 1, 40/) )
            Print *, test_mat%indx -  ( (/1, 2, 3, 4, 1/) )
            Print *, test_mat%jndx -  ( (/1, 2, 3, 4, 2/) )

        end subroutine test_set_value

        subroutine test_get_value
            implicit none
            type(coo_matrix) :: test_mat
            integer :: ierr

            test_mat = speye(4, ierr)
            call test_mat%set_value(1, 2, dble(40), ierr)

            Print *, dble(1) - test_mat%get_value(2, 2)
            Print *, dble(40) - test_mat%get_value(1, 2)

        end subroutine test_get_value

        subroutine test_coo_multiply
            implicit none
            type(coo_matrix) :: test_mat
            integer :: ierr
            real (kind=8), dimension(:, :), allocatable :: x, y

            test_mat = speye(4, ierr)
            test_mat%vals = dble((/1, 2, 3, 4/))
            call test_mat%set_value(1, 2, dble(40), ierr)

            allocate( x(4, 2), y(4, 2) )
            x = reshape((/1, 1, 1, 1, 2, 2, 2, 2/), (/4, 2/))
            call coo_multiply(test_mat, x, y)
            Print *, y(:, 1) - dble((/41, 2, 3, 4/))
            Print *, y(:, 2) - dble((/82, 4, 6, 8/))

        end subroutine test_coo_multiply

        subroutine test_sparse_add
            implicit none
            type (coo_matrix) :: mat1, mat2
            integer :: ierr

            mat1 = speye(4, ierr)
            mat2 = speye(4, ierr)
            call mat1%set_value(1, 2, dble(40), ierr)
            call sparse_add(mat1, mat2, ierr)
            Print *, mat1%vals - dble((/2, 40, 2, 2, 2/))
            Print *, mat1%indx - (/1, 1, 2, 3, 4/)
            Print *, mat1%jndx - (/1, 2, 2, 3, 4/)


        end subroutine test_sparse_add

        subroutine test_coo_to_csr
            implicit none
            type (coo_matrix) :: mat1
            type (csr_matrix) :: matout
            integer :: ierr

            mat1 = speye(4, ierr)
            call mat1%set_value(1, 4, dble(40), ierr)
            call mat1%set_value(2, 4, dble(40), ierr)
            call mat1%set_value(3, 4, dble(40), ierr)
            matout = coo_to_csr(mat1, ierr)
            Print *, matout%vals - dble( (/1, 40, 1, 40, 1, 40, 1/) )
            Print *, matout%rows - (/1, 3, 5, 7, 8/)
            Print *, matout%jndx - (/1, 4, 2, 4, 3, 4, 4/)

        end subroutine test_coo_to_csr

        subroutine test_csr_multiply
            implicit none
            type(coo_matrix) :: test_mat
            type (csr_matrix) :: test_mat_csr
            integer :: ierr
            real (kind=8), dimension(:, :), allocatable :: x, y

            test_mat = speye(4, ierr)
            test_mat%vals = dble((/1, 2, 3, 4/))
            call test_mat%set_value(1, 2, dble(40), ierr)

            test_mat_csr = coo_to_csr(test_mat, ierr)

            allocate( x(4, 2), y(4, 2) )
            x = reshape((/1, 1, 1, 1, 2, 2, 2, 2/), (/4, 2/))
            call csr_multiply(test_mat_csr, x, y)
            Print *, y(:, 1) - dble((/41, 2, 3, 4/))
            Print *, y(:, 2) - dble((/82, 4, 6, 8/))
        end subroutine test_csr_multiply



end program test_sparse_matrices