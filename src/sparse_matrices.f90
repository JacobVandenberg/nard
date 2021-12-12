! Created by  on 13/7/21.

module sparse_matrices
    use helpers
    use hdf5
    use precision
    implicit none
    public speye

    type coo_matrix
        integer (ip), allocatable, dimension(:) :: indx, jndx
        real (rp), allocatable, dimension(:) :: vals
        integer (ip) :: n, m

        contains

            procedure, public :: set_value
            procedure, public :: get_value
            procedure, public :: insert_block

    end type coo_matrix

    type csr_matrix
        integer (ip), allocatable, dimension(:) :: rows, jndx
        real (rp), allocatable, dimension(:) :: vals
        integer (ip) :: n, m

    end type csr_matrix

    type sparse_lu_decomposition
        integer (ip) :: n, m, mtype
        integer (ip), dimension(64) :: pt, iparm
        integer (ip), allocatable, dimension(:) :: perm
        type (csr_matrix) :: matrix
    end type sparse_lu_decomposition

    contains
        type(coo_matrix) function sparse_kron(matrix1, matrix2, ierr)
            ! computes the kronecker tensor product of two sparse matrices
            !
            ! arguments: matrix1: the first coo_matrix
            !            matrix2: the second coo_matrix
            !
            ! returns kronecker product of both matrices
            implicit none
            ! inputs
            type(coo_matrix) :: matrix1, matrix2
            integer (ip), optional, intent(out) :: ierr

            ! return values
            type(coo_matrix) :: matrixout

            ! subroutine variables
            integer (ip) :: block_i, block_j, element_gi, numel1, numel2, i1, i2, block_gi

            numel1 = size(matrix1%vals, 1_ip)
            numel2 = size(matrix2%vals, 1_ip)

            allocate(matrixout%indx(numel1 * numel2), matrixout%jndx(numel1 * numel2), matrixout%vals(numel1 * numel2), STAT=ierr)
            if (ierr/=0_ip) then
                return
            end if
            do i1 = 1, numel1
                block_i = (matrix1%indx(i1)-1_ip) * matrix2%n
                block_j = (matrix1%jndx(i1)-1_ip) * matrix2%m
                block_gi = numel2 * (i1 - 1_ip)
                do i2 = 1_ip, numel2
                    element_gi = block_gi + i2
                    matrixout%vals(element_gi) = matrix1%vals(i1) * matrix2%vals(i2)
                    matrixout%indx(element_gi) = block_i + matrix2%indx(i2)
                    matrixout%jndx(element_gi) = block_j + matrix2%jndx(i2)
                end do
            end do
            matrixout%n = matrix1%n * matrix2%n
            matrixout%m = matrix1%m * matrix2%m
            sparse_kron = matrixout
            return
        end function sparse_kron

        type(coo_matrix) function speye(N, stat)
            !
            ! returns a sparse matrix in coo format
            !
            ! arguments: N: Size of matrix
            !            stat: internal error status
            implicit none
            ! arguments
            integer (ip), intent(in) :: N
            integer (ip), intent(out) :: stat

            ! return
            type(coo_matrix) :: matrix_result


            ! allocate values
            allocate (matrix_result%vals(N), stat=stat)
            if (stat/=0) then
                return
            end if
            ! assign values
            matrix_result%vals = dble(1)

            ! assign dimensions
            matrix_result%n = N
            matrix_result%m = N
            ! assign and allocate indice
            call arange_int(1_ip, N, matrix_result%indx, stat)
            if (stat/=0_ip) then
                return
            end if
            call arange_int(1_ip, N, matrix_result%jndx, stat)
            if (stat/=0_ip) then
                return
            end if

            ! return sparse matrix
            speye = matrix_result
        end function speye

        subroutine set_value(this, i_index, j_index, val, ierr)
            !
            ! sets the value at index (i, j) to val
            !
            ! Arguments:
            !   this: coo_matrix
            !   i: i index
            !   j: j index
            !   val: the value to set it to

            implicit none
            class(coo_matrix) :: this
            integer (ip), intent(in) :: i_index, j_index
            real (rp), intent(in) :: val
            integer (ip) :: element_i, numel
            integer (ip), intent(out) :: ierr

            integer (ip), dimension(:), allocatable :: indx_new, jndx_new
            real (rp), dimension(:), allocatable :: vals_new

            numel = size(this%vals, 1)
            do element_i = 1, numel
                ! if the element exists in the array
                if (this%indx(element_i) == i_index .and. this%jndx(element_i) == j_index) then
                    this%vals(element_i) = val
                    return
                end if
            end do
            ! otherwise reallocate array
            allocate (indx_new(numel+1_ip), jndx_new(numel+1_ip), vals_new(numel+1_ip), STAT=ierr)

            if (ierr/=0_ip) then
                return
            end if

            indx_new(1_ip:numel) = this%indx
            jndx_new(1_ip:numel) = this%jndx
            vals_new(1_ip:numel) = this%vals

            indx_new(numel+1_ip) = i_index
            jndx_new(numel+1_ip) = j_index
            vals_new(numel+1_ip) = val
            call move_alloc (indx_new, this%indx)
            call move_alloc (jndx_new, this%jndx)
            call move_alloc (vals_new, this%vals)
            return

        end subroutine set_value

        real (rp) function get_value(this, i_index, j_index)
            !
            ! gets the value at index (i_index, j_index) from a coo_matrix
            implicit none
            ! input
            class (coo_matrix) :: this
            integer (ip), intent(in) :: i_index, j_index
            ! internal
            integer (ip) :: current_el

            do current_el = 1_ip,size(this%vals, 1_ip)
                if (this%indx(current_el) == i_index .and. this%jndx(current_el) == j_index) then
                    get_value = this%vals(current_el)
                    return
                end if
            end do
            get_value = dble(0)
            return

        end function get_value

        subroutine coo_multiply(matrix, input_x, output_y)
            !
            ! Computes matrix - matrix multiplication y <- A x, where x is a dense matrix
            !
            ! relatively slow. use CSR for multiplication if performance is needed
            !
            ! arguments:
            !   matrix: coo_matrix
            !   input_x: dense array
            !   output_y: output dense array

            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type(coo_matrix) :: matrix
            real (rp), dimension(:, :), intent(in) :: input_x
            ! outputs
            real (rp), dimension(:, :), intent(out) :: output_y
            ! subroutine variables
            integer (ip) :: numel, i
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            output_y=0 ! initialise y to 0
            numel = size(matrix%vals, 1_ip) ! get number of elements

            ! perform matrix multiplication elemnetwise from the sparse matrix
            do i = 1_ip, numel
                output_y(matrix%indx(i), :) = matrix%vals(i) * input_x(matrix%jndx(i), :) + output_y(matrix%indx(i), :)
            end do
            return
            ! END SUBROUTINE
        end subroutine coo_multiply

        subroutine sparse_add(mat1, mat2, ierr)
            !
            ! adds two sparse coo matrices together
            ! will also rearrange the elements of both matrices
            ! result is stored in mat1
            !
            ! inputs: mat1, mat2: the two matrices to add
            !   ierr: internal error flag
            !
            ! outputs: mat1
            implicit none

            ! BEGIN DECLARATIONS
            ! input
            type(coo_matrix), intent(inout) :: mat1, mat2
            ! output
            integer (ip), intent(out) :: ierr
            ! subroutine variables
            integer (ip) :: duplicates, numel1, numel2, numelout, i, j, current, num_cols
            integer (ip), allocatable, dimension(:) :: key, indxout, jndxout
            real (rp), dimension(:), allocatable :: valout
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            ! append matrices together
            numel1 = size(mat1%vals, 1_ip)
            numel2 = size(mat2%vals, 1_ip)
            numelout = numel1+numel2 ! this is the maximum
            allocate (indxout(numelout), STAT = ierr)
            num_cols = maxval((/maxval( mat1%jndx ), maxval(mat2%jndx)/))
            indxout(1_ip:numel1) = (mat1%indx - 1_ip) * (num_cols) + mat1%jndx - 1_ip
            indxout(numel1+1 : numelout) = (mat2%indx - 1_ip) * (num_cols) + mat2%jndx - 1_ip


            ! find output size
            duplicates = 0_ip
            call quicksort_int(indxout, 1_ip, numelout)
            duplicates = count(indxout(2_ip:numelout) == indxout(1_ip:numelout-1_ip))
            deallocate(indxout)

            ! sort both input matrices.
            call arange_int(1_ip, numel1, key, ierr)
            call order_int((mat1%indx - 1_ip) * (num_cols) + mat1%jndx - 1_ip, 1_ip, numel1, key)

            ! sort into temporary array
            allocate (indxout(numel1), jndxout(numel1), valout(numel1))
            do i= 1_ip,numel1
                indxout(i) = mat1%indx(key(i))
                jndxout(i) = mat1%jndx(key(i))
                valout(i) = mat1%vals(key(i))
            end do
            call move_alloc(indxout, mat1%indx)
            call move_alloc(jndxout, mat1%jndx)
            call move_alloc(valout, mat1%vals)
            deallocate(key)

            call arange_int(1_ip, numel2, key, ierr)
            call order_int((mat2%indx - 1_ip) * (num_cols) + mat2%jndx - 1_ip, 1_ip, numel2, key)
            ! sort into temporary array
            allocate (indxout(numel2), jndxout(numel2), valout(numel2), STAT=ierr)
            do i= 1_ip,numel2
                indxout(i) = mat2%indx(key(i))
                jndxout(i) = mat2%jndx(key(i))
                valout(i) = mat2%vals(key(i))
            end do
            call move_alloc(indxout, mat2%indx)
            call move_alloc(jndxout, mat2%jndx)
            call move_alloc(valout, mat2%vals)
            deallocate(key)

            ! allocate proper number of elements
            numelout = numel1 + numel2 - duplicates
            allocate (valout(numelout), indxout(numelout), jndxout(numelout), STAT=ierr)
            ! iterate over elements in both matrices, dealing with the elements in CSR order.
            i=1_ip; j=1_ip
            current = 1_ip

            do while (current <= numelout)
                if (i > numel1) then
                    valout(current) = mat2%vals(j)
                    indxout(current) = mat2%indx(j)
                    jndxout(current) = mat2%jndx(j)
                    j = j + 1_ip
                elseif (j > numel2) then
                    valout(current) = mat1%vals(i)
                    indxout(current) = mat1%indx(i)
                    jndxout(current) = mat1%jndx(i)
                    i = i + 1_ip
                elseif (mat1%indx(i) == mat2%indx(j) .and. mat1%jndx(i) == mat2%jndx(j)) then
                    valout(current) = mat1%vals(i) + mat2%vals(j)
                    indxout(current) = mat1%indx(i)
                    jndxout(current) = mat1%jndx(i)
                    i = i + 1_ip; j = j + 1_ip
                elseif  ((mat1%indx(i) - 1_ip) * num_cols + mat1%jndx(i) - 1_ip < (mat2%indx(j) - 1_ip) * num_cols +&
                    mat2%jndx(j) - 1_ip) then
                    valout(current) = mat1%vals(i)
                    indxout(current) = mat1%indx(i)
                    jndxout(current) = mat1%jndx(i)
                    i = i + 1_ip
                else
                    valout(current) = mat2%vals(j)
                    indxout(current) = mat2%indx(j)
                    jndxout(current) = mat2%jndx(j)
                    j = j + 1_ip
                end if
                current = current + 1_ip
            end do
            ! move new entries to input matrix

            call move_alloc(indxout, mat1%indx)
            call move_alloc(jndxout, mat1%jndx)
            call move_alloc(valout, mat1%vals)
            return
            ! END SUBROUTINE
        end subroutine sparse_add

        type (csr_matrix) function coo_to_csr(input_matrix, ierr)
            !
            ! converts a sparse matrix in coordinate form into a sparse matrix in
            ! compressed row form
            !
            ! Arguments:
            !   input_matrix: coo_matrix to convert (may be rearranged)
            !   ierr: output error flag
            !
            ! Returns:
            !   csr_matrix

            implicit none
            ! BEGIN DECLARATIONS
            ! arguments
            type (coo_matrix), intent(inout) :: input_matrix
            integer (ip), intent(out) :: ierr

            ! return value
            type (csr_matrix) :: return_matrix

            !internals
            integer (ip), dimension(:), allocatable :: key
            integer (ip) :: num_cols, numel, ii, current_row
            integer (ip), dimension(:), allocatable :: temp_int
            real (rp), dimension(:), allocatable :: temp_real
            ! END DECLARATIONS

            ! BEGIN FUNCTION
            ! first get the order of the elements
            num_cols = maxval(input_matrix%jndx)
            numel = size(input_matrix%vals, 1_ip)
            call arange_int(1_ip, numel, key, ierr)
            call order_int( (input_matrix%indx-1) * num_cols + input_matrix%jndx - 1_ip, 1_ip, numel, key)

            ! reindex the column indices and the values
            allocate ( return_matrix%vals(numel), return_matrix%jndx(numel), STAT=ierr )
            do ii=1,numel
                return_matrix%vals(ii) = input_matrix%vals(key(ii))
                return_matrix%jndx(ii) = input_matrix%jndx(key(ii))
            end do

            ! find row vector
            num_cols = maxval(input_matrix%indx)  ! well, rows now
            allocate (return_matrix%rows(num_cols+1_ip), STAT=ierr)

            ! assign row vector
            current_row = 0_ip
            do ii = 1_ip, numel
                if (input_matrix%indx(key(ii)) > current_row) then
                    return_matrix%rows(current_row+1:input_matrix%indx(key(ii))) = ii
                    current_row = input_matrix%indx(key(ii))
                end if
            end do
            return_matrix%rows(num_cols+1_ip) = numel+1_ip
            return_matrix%n = input_matrix%n
            return_matrix%m = input_matrix%m
            coo_to_csr = return_matrix
            ! END FUNCTION
        end function coo_to_csr

        subroutine csr_multiply(matrix, x, y)
            ! multiples CSR (compressed sparse matrix) A by dense matrix x, returning y.
            !
            ! ARGUMENTS:
            !   matrix: input matrix
            !   x: dense matrix
            !   y: output dense matrix (overwritten)

            implicit none

            ! BEGIN DECLARATIONS
            ! arguments
            type (csr_matrix) :: matrix
            real (rp), dimension(:, :), intent(in) :: x
            real (rp), dimension(:, :), intent(out) :: y

            ! subroutine variables
            integer (ip) :: num_rows, current_row, current_el, col_end, numel
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            num_rows = size(matrix%rows, 1_ip)-1_ip
            numel = size(matrix%vals, 1_ip)

            !$OMP PARALLEL DO
            do current_row = 1_ip, num_rows
                y(current_row, :) = 0_ip
                do current_el = matrix%rows(current_row), matrix%rows(current_row+1_ip)-1_ip
                    y(current_row, :) = y(current_row, :) + matrix%vals(current_el) * x(matrix%jndx(current_el), :)
                end do
            end do
            !$OMP END PARALLEL DO
            return
            ! END SUBROUTINE
        end subroutine csr_multiply

        type(coo_matrix) function copy_coo_matrix(matrix, ierr)
            !
            ! creates a new coo_matrix and copies data from inmatrix
            !
            ! inputs:
            !   matrix: matrix to copy
            !
            ! outputs:
            !   ierr: internal error flag
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (coo_matrix), intent(in) :: matrix

            ! outputs
            type (coo_matrix) :: out_matrix
            integer (ip), intent(out) :: ierr

            ! runtime
            integer (ip) :: numel
            ! END DECLARATIONS

            numel = size(matrix%vals)
            allocate ( out_matrix%vals(numel), out_matrix%indx(numel), out_matrix%jndx(numel), STAT=ierr )
            if (ierr/=0_ip) then
                return
            end if
            out_matrix%vals(:) = matrix%vals(:)
            out_matrix%indx(:) = matrix%indx(:)
            out_matrix%jndx(:) = matrix%jndx(:)
            out_matrix%n = matrix%n
            out_matrix%m = matrix%m
            copy_coo_matrix = out_matrix
        end function copy_coo_matrix

        type (sparse_lu_decomposition) function sparse_lu(matrix, ierr)
            !
            ! generates a sparse lu decomposition of the given matrix
            !
            ! inputs:
            !   matrix: a coo_matrix
            !
            ! outputs:
            !   ierr: internal error flag
            !   return value: sparse lu decomposition
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (coo_matrix), intent(inout) :: matrix

            ! outputs
            integer (ip), intent(out) :: ierr
            type (sparse_lu_decomposition) :: return_value
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            ! unsymmetric real. Might be symmetric for backwards euler, and with only periodic BCS, but not in general
            return_value%mtype = 11_ip
            call pardisoinit(return_value%pt, return_value%mtype, return_value%iparm)

            return_value%matrix = coo_to_csr(matrix, ierr)
            allocate (return_value%perm(return_value%matrix%m), STAT=ierr)

            call pardiso( return_value%pt, 1_ip, 1_ip, return_value%mtype, 12_ip,& ! 12 means that we are analysing and factoring
                    return_value%matrix%m, return_value%matrix%vals, return_value%matrix%rows, return_value%matrix%jndx,&
                    return_value%perm, 1_ip, return_value%iparm, 0_ip, real((/0/), kind=rp), real((/0/), kind=rp), ierr)

            sparse_lu = return_value
            ! END SUBROUTINE
        end function sparse_lu

        function sparse_lu_solve(lu_decompstn, b, ierr)
            !
            ! solves the system A x = b for x, given A is in a sparse lu decomposition
            !
            ! inputs:
            !   lu_decompstn: sparse lu decomposition of A
            !   b: the RHS vector to solve for. Supports solving for multiple right hand sides.
            !   caution: b will be modified
            !
            ! outputs:
            !   x: matrix of column vector solutions
            !   ierr: internal error flag

            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (sparse_lu_decomposition), intent(inout) :: lu_decompstn
            real (rp), dimension(:), intent(inout) :: b

            ! outputs
            real (rp), dimension(:), allocatable :: sparse_lu_solve
            integer (ip), intent(out) :: ierr
            ! runtime
            ! END DECLARATION

            ! api should be changed to reduce reshaping
            ! BEGIN SUBROUTINE
            allocate (sparse_lu_solve(lu_decompstn%matrix%n), STAT=ierr)
            if (ierr /= 0_ip) then
                return
            end if
            call pardiso( lu_decompstn%pt, 1_ip, 1_ip, lu_decompstn%mtype, 33_ip,& ! 33 means that we are solving
                    lu_decompstn%matrix%m, lu_decompstn%matrix%vals, lu_decompstn%matrix%rows, lu_decompstn%matrix%jndx,&
                    lu_decompstn%perm, 1_ip, lu_decompstn%iparm, 0_ip, b,&
                    sparse_lu_solve, ierr)
            return
            ! END SUBROUTINE
        end function sparse_lu_solve

        subroutine insert_block(this, block_input, i_offset, j_offset, ierr)
            !
            ! inserts a matrix into another matrix.
            ! inserts block_input into this, with offset i_offset for the column index, and j_offset for the row index.
            ! an offset of 0, 0 means that the indices of block_input will not be changed.
            ! if there are overlapping entries, they will be added.
            !
            ! INPUTS:
            !   block_input: the matrix to insert
            !   i_offset: the offset of the column indices
            !   j_offset: the offset of the row indices.
            !   this (inout): the matrix to add to.

            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            class( coo_matrix ), intent(inout) :: this
            type (coo_matrix), intent(in) :: block_input
            integer (ip), intent(in) :: i_offset, j_offset

            ! outputs
            integer (ip), intent(out) :: ierr

            ! runtime
            type (coo_matrix) :: temp_matrix

            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            temp_matrix = copy_coo_matrix(block_input, ierr)

            ! offset matrix
            temp_matrix%indx = temp_matrix%indx + i_offset
            temp_matrix%jndx = temp_matrix%jndx + j_offset
            temp_matrix%m = temp_matrix%m + i_offset
            temp_matrix%n = temp_matrix%n + j_offset

            ! add matrices

            call sparse_add(this, temp_matrix, ierr)

            this%n = maxval( (/ temp_matrix%n, this%n /) )
            this%m = maxval( (/ temp_matrix%m, this%m /) )

            deallocate ( temp_matrix%vals, temp_matrix%indx, temp_matrix%jndx )
            return
            ! END SUBROUTINE
        end subroutine insert_block

        subroutine csr_scale_rows(matrix, scaling_factors)
            !
            ! scales the rows of a sparse matrix (in csr format)
            !
            ! inputs:
            !   matrix (modified): the sparse matrix to modify
            !   scaling_factors: the factors to scale each row by (column vector)
            !
            implicit none
            ! BEGIN DECLARATIONS
            type (csr_matrix), intent(inout) :: matrix
            real (rp), intent(in), dimension(:) :: scaling_factors
            ! runtime
            integer (ip) :: current_row
            ! END DECLARAIONS

            !$OMP PARALLEL DO
            do current_row = 1_ip, matrix%m
                matrix%vals(matrix%rows(current_row):matrix%rows(current_row+1_ip)-1_ip) = &
                        matrix%vals(matrix%rows(current_row):matrix%rows(current_row+1_ip)-1_ip) *&
                                scaling_factors(current_row)
            end do
            !$OMP END PARALLEL DO
            return
        end subroutine csr_scale_rows

        subroutine coo_scale_rows(matrix, scaling_factors)
            !
            ! scales the rows of a sparse matrix (in csr format)
            !
            ! inputs:
            !   matrix (modified): the sparse matrix to modify
            !   scaling_factors: the factors to scale each row by (column vector)
            !
            implicit none
            ! BEGIN DECLARATIONS
            type (coo_matrix), intent(inout) :: matrix
            real (rp), intent(in), dimension(:) :: scaling_factors
            ! runtime
            integer (ip) :: current_el
            ! END DECLARAIONS

            !$OMP PARALLEL DO
            do current_el = 1_ip, size(matrix%vals, kind=ip)
                matrix%vals(current_el) = matrix%vals(current_el) * scaling_factors(matrix%indx(current_el))
            end do
            !$OMP END PARALLEL DO
            return
        end subroutine coo_scale_rows

        function copy_csr_matrix(matrix, ierr)
            !
            ! creates a new csr_matrix and copies data from matrix
            !
            ! inputs:
            !   matrix: matrix to copy
            !
            ! outputs:
            !   ierr: internal error flag
            !
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (csr_matrix), intent(in) :: matrix

            ! outputs
            type (csr_matrix) :: copy_csr_matrix
            integer (ip), intent(out) :: ierr

            ! runtime
            integer (ip) :: numel
            ! END DECLARATIONS

            numel = size(matrix%vals)
            allocate ( copy_csr_matrix%vals(numel), copy_csr_matrix%rows(numel), copy_csr_matrix%jndx(numel), STAT=ierr )
            if (ierr/=0_ip) then
                return
            end if
            copy_csr_matrix%vals(:) = matrix%vals(:)
            copy_csr_matrix%rows(:) = matrix%rows(:)
            copy_csr_matrix%jndx(:) = matrix%jndx(:)
            copy_csr_matrix%n = matrix%n
            copy_csr_matrix%m = matrix%m
            return
        end function copy_csr_matrix

        function sparse_direct_solve(matrix, b, ierr)
            !
            ! solves a system of sparse linear equations described by a csr matrix
            !
            ! inputs: matrix: the csr_matrix
            !   b: the RHS to solve for
            !
            ! outputs: ierr: internal error flag
            implicit none
            ! BEGIN DECLARATIONS
            ! inputs
            type (csr_matrix), intent(in) :: matrix
            real (rp), dimension(:), intent(inout) :: b

            ! outputs
            integer (ip) :: ierr
            real (rp), dimension(:), allocatable :: sparse_direct_solve

            ! runtime
            integer (ip) :: mtype
            integer (ip), dimension(64) :: iparm, pt
            integer (ip), allocatable, dimension(:) :: perm

            ! END DECLARATIONS

            allocate( sparse_direct_solve(matrix%n), perm(matrix%n) )

            mtype = 11_ip
            call pardisoinit(pt, mtype, iparm)

            call pardiso( pt, 1_ip, 1_ip, mtype, 13_ip,&
                    matrix%m, matrix%vals, matrix%rows, matrix%jndx,&
                    perm, 1_ip, iparm, 0_ip, b, sparse_direct_solve, ierr)
            call pardiso( pt, 1_ip, 1_ip, mtype, -1_ip,&
                    matrix%m, matrix%vals, matrix%rows, matrix%jndx,&
                    perm, 1_ip, iparm, 0_ip, b, sparse_direct_solve, ierr)
        end function sparse_direct_solve

end module sparse_matrices