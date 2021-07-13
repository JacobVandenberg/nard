! Created by  on 13/7/21.

module sparse_matrices
    use helpers
    implicit none
    public speye

    type coo_matrix
        integer, allocatable, dimension(:) :: indx, jndx
        real (kind=8), allocatable, dimension(:) :: vals
        integer :: n, m

        contains

            procedure, public :: set_value
            procedure, public :: get_value

    end type coo_matrix

    type csr_matrix
        integer, allocatable, dimension(:) :: rows, jndx
        real (kind=8), allocatable, dimension(:) :: vals
        integer :: n, m

    end type csr_matrix

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
            integer, optional, intent(out) :: ierr

            ! return values
            type(coo_matrix) :: matrixout

            ! subroutine variables
            integer :: block_i, block_j, element_gi, numel1, numel2, i1, i2, block_gi

            numel1 = size(matrix1%vals, 1)
            numel2 = size(matrix2%vals, 1)

            allocate(matrixout%indx(numel1 * numel2), matrixout%jndx(numel1 * numel2), matrixout%vals(numel1 * numel2), STAT=ierr)
            if (ierr/=0) then
                return
            end if
            do i1 = 1, numel1
                block_i = (matrix1%indx(i1)-1) * matrix2%n
                block_j = (matrix1%jndx(i1)-1) * matrix2%m
                block_gi = numel2 * (i1 - 1)
                do i2 = 1, numel2
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
            integer, intent(in) :: N
            integer, intent(out) :: stat

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
            call arange_int(1, N, matrix_result%indx, stat)
            if (stat/=0) then
                return
            end if
            call arange_int(1, N, matrix_result%jndx, stat)
            if (stat/=0) then
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
            integer, intent(in) :: i_index, j_index
            real (kind=8), intent(in) :: val
            integer :: element_i, numel
            integer, intent(out) :: ierr

            integer, dimension(:), allocatable :: indx_new, jndx_new
            real (kind=8), dimension(:), allocatable :: vals_new

            numel = size(this%vals, 1)
            do element_i = 1, numel
                ! if the element exists in the array
                if (this%indx(element_i) == i_index .and. this%jndx(element_i) == j_index) then
                    this%vals(element_i) = val
                    return
                end if
            end do
            ! otherwise reallocate array
            allocate (indx_new(numel+1), jndx_new(numel+1), vals_new(numel+1), STAT=ierr)

            if (ierr/=0) then
                return
            end if

            indx_new(1:numel) = this%indx
            jndx_new(1:numel) = this%jndx
            vals_new(1:numel) = this%vals

            indx_new(numel+1) = i_index
            jndx_new(numel+1) = j_index
            vals_new(numel+1) = val
            call move_alloc (indx_new, this%indx)
            call move_alloc (jndx_new, this%jndx)
            call move_alloc (vals_new, this%vals)
            return

        end subroutine set_value

        real (kind=8) function get_value(this, i_index, j_index)
            !
            ! gets the value at index (i_index, j_index) from a coo_matrix
            implicit none
            ! input
            class (coo_matrix) :: this
            integer, intent(in) :: i_index, j_index
            ! internal
            integer :: current_el

            do current_el = 1,size(this%vals, 1)
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
            real (kind=8), dimension(:, :), intent(in) :: input_x
            ! outputs
            real (kind=8), dimension(:, :), intent(out) :: output_y
            ! subroutine variables
            integer :: numel, i
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            output_y=0 ! initialise y to 0
            numel = size(matrix%vals, 1) ! get number of elements

            ! perform matrix multiplication elemnetwise from the sparse matrix
            do i = 1, numel
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
            integer, intent(out) :: ierr
            ! subroutine variables
            integer :: duplicates, numel1, numel2, numelout, i, j, current, num_cols
            integer, allocatable, dimension(:) :: key, indxout, jndxout
            real (kind=8), dimension(:), allocatable :: valout
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            ! append matrices together
            numel1 = size(mat1%vals, 1)
            numel2 = size(mat2%vals, 1)
            numelout = numel1+numel2 ! this is the maximum
            allocate (indxout(numelout), STAT = ierr)
            num_cols = maxval((/maxval( mat1%jndx ), maxval(mat2%jndx)/))
            indxout(1:numel1) = (mat1%indx - 1) * (num_cols) + mat1%jndx - 1
            indxout(numel1+1 : numelout) = (mat2%indx - 1) * (num_cols) + mat2%jndx - 1


            ! find output size
            duplicates = 0
            call quicksort_int(indxout, 1, numelout)
            duplicates = count(indxout(2:numelout) == indxout(1:numelout-1))
            deallocate(indxout)

            ! sort both input matrices.
            call arange_int(1, numel1, key, ierr)
            call order_int((mat1%indx - 1) * (num_cols) + mat1%jndx - 1, 1, numel1, key)

            ! sort into temporary array
            allocate (indxout(numel1), jndxout(numel1), valout(numel1))
            do i= 1,numel1
                indxout(i) = mat1%indx(key(i))
                jndxout(i) = mat1%jndx(key(i))
                valout(i) = mat1%vals(key(i))
            end do
            call move_alloc(indxout, mat1%indx)
            call move_alloc(jndxout, mat1%jndx)
            call move_alloc(valout, mat1%vals)
            deallocate(key)

            call arange_int(1, numel2, key, ierr)
            call order_int((mat2%indx - 1) * (num_cols) + mat2%jndx - 1, 1, numel2, key)
            ! sort into temporary array
            allocate (indxout(numel2), jndxout(numel2), valout(numel2), STAT=ierr)
            do i= 1,numel2
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
            i=1; j=1
            current = 1

            do while (current <= numelout)
                if (i > numel1) then
                    valout(current) = mat2%vals(j)
                    indxout(current) = mat2%indx(j)
                    jndxout(current) = mat2%jndx(j)
                    j = j + 1
                elseif (j > numel2) then
                    valout(current) = mat2%vals(j)
                    indxout(current) = mat2%indx(j)
                    jndxout(current) = mat2%jndx(j)
                    j = j + 1
                elseif (mat1%indx(i) == mat2%indx(j) .and. mat1%jndx(i) == mat2%jndx(j)) then
                    valout(current) = mat1%vals(i) + mat2%vals(j)
                    indxout(current) = mat1%indx(i)
                    jndxout(current) = mat1%jndx(i)
                    i = i + 1; j = j + 1
                elseif  ((mat1%indx(i) - 1) * num_cols + mat1%jndx(i) - 1 < (mat2%indx(j) - 1) * num_cols + mat2%jndx(j) - 1) then
                    valout(current) = mat1%vals(i)
                    indxout(current) = mat1%indx(i)
                    jndxout(current) = mat1%jndx(i)
                    i = i + 1
                else
                    valout(current) = mat2%vals(j)
                    indxout(current) = mat2%indx(j)
                    jndxout(current) = mat2%jndx(j)
                    j = j + 1
                end if
                current = current + 1
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
            integer, intent(out) :: ierr

            ! return value
            type (csr_matrix) :: return_matrix

            !internals
            integer, dimension(:), allocatable :: key
            integer :: num_cols, numel, ii, current_row
            integer, dimension(:), allocatable :: temp_int
            real (kind=8), dimension(:), allocatable :: temp_real
            ! END DECLARATIONS

            ! BEGIN FUNCTION
            ! first get the order of the elements
            num_cols = maxval(input_matrix%jndx)
            numel = size(input_matrix%vals, 1)
            call arange_int(1, numel, key, ierr)
            call order_int( (input_matrix%indx-1) * num_cols + input_matrix%jndx - 1, 1, numel, key)

            ! reindex the column indices and the values
            allocate ( return_matrix%vals(numel), return_matrix%jndx(numel), STAT=ierr )
            do ii=1,numel
                return_matrix%vals(ii) = input_matrix%vals(key(ii))
                return_matrix%jndx(ii) = input_matrix%jndx(key(ii))
            end do

            ! find row vector
            num_cols = maxval(input_matrix%indx)  ! well, rows now
            allocate (return_matrix%rows(num_cols+1), STAT=ierr)

            ! assign row vector
            current_row = 0
            do ii = 1, numel
                if (input_matrix%indx(key(ii)) > current_row) then
                    return_matrix%rows(current_row+1:input_matrix%indx(key(ii))) = ii
                    current_row = input_matrix%indx(key(ii))
                end if
            end do
            return_matrix%rows(num_cols+1) = numel+1
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
            real (kind=8), dimension(:, :), intent(in) :: x
            real (kind=8), dimension(:, :), intent(out) :: y

            ! subroutine variables
            integer :: num_rows, current_row, current_el, col_end, numel
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            num_rows = size(matrix%rows, 1)-1
            numel = size(matrix%vals, 1)

            !$OMP PARALLEL DO
            do current_row = 1, num_rows
                y(current_row, :) = 0
                do current_el = matrix%rows(current_row), matrix%rows(current_row+1)-1
                    y(current_row, :) = y(current_row, :) + matrix%vals(current_el) * x(matrix%jndx(current_el), :)
                end do
            end do
            !$OMP END PARALLEL DO
            return
            ! END SUBROUTINE
        end subroutine csr_multiply



end module sparse_matrices