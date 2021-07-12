! Created by  on 6/7/21.

module helpers
    implicit none
    contains
    subroutine my_linspace(start, finish, N, x)
        implicit none
        real (kind=8), intent(in) :: start, finish
        integer, intent(in) :: N
        real (kind=8), intent(out) :: x(N)
        integer :: i
        real (kind=8) :: range

        range = finish - start

        if (N == 0) return

        if (N ==1) then
            x(1) = start
            return
        end if

        do i=1, N
            x(i) = start + range * (i-1) / (N-1)
        end do
    end subroutine my_linspace

    subroutine meshgrid3(x, y, z, nx, ny, nz, xxx, yyy, zzz)
        implicit none
        real (kind=8), intent(in) :: x(nx), y(ny), z(nz)
        real (kind=8), intent(out), dimension(:, :, :) :: xxx, yyy, zzz
        integer (kind=8) :: i, j, k, nx, ny, nz

        do i=1, nx
            do j=1, ny
                do k=1,nz
                    xxx(j, i, k) = x(i)
                    yyy(j, i, k) = y(j)
                    zzz(j, i, k) = z(k)
                end do
            end do
        end do

    end subroutine meshgrid3

    subroutine sparse_kron(val1, indx1, jndx1, n1, m1, val2, indx2, jndx2, n2, m2, valout, indxout, jndxout, nout, mout, ierr)
        ! computes the kronecker tensor product of two sparse matrices
        !
        ! arguments: val1, val2, valout (output): value vector for sparse matrices 1, 2 and output
        !         indx1, indx2, indxout (output): i indices for first two input matrices and output matrix
        !         jndx1, jndx2, jndxout (output): j indices for forst two input matrices and output matrix
        !         n1, n2, nout (output): number of rows in input matrices and output matrix
        !         m1, m2, mout (output): number of columns in input matrices and output matrix
        !
        !
        implicit none
        real(kind=8), dimension(:), intent(in) :: val1, val2
        integer, dimension(:), intent(in) :: indx1, indx2, jndx1, jndx2
        integer, intent(in) :: n1, m1, n2, m2
        real(kind=8), allocatable, dimension(:), intent(out) :: valout
        integer, dimension(:), allocatable, intent(out) :: indxout, jndxout
        integer, intent(out) :: nout, mout, ierr
        integer :: numel1, numel2, i1, i2, block_i, block_j, block_gi, element_gi

        numel1 = size(val1, 1)
        numel2 = size(val2, 1)

        allocate(indxout(numel1 * numel2), jndxout(numel1 * numel2), valout(numel1 * numel2), STAT=ierr)

        if (ierr/=0) then
            return
        end if

        do i1 = 1, numel1
            block_i = (indx1(i1)-1) * n2
            block_j = (jndx1(i1)-1) * m2
            block_gi = numel2 * (i1 - 1)
            do i2 = 1, numel2
                element_gi = block_gi + i2
                valout( element_gi) = val1(i1) * val2(i2)
                indxout(element_gi) = block_i + indx2(i2)
                jndxout(element_gi) = block_j + jndx2(i2)
            end do
        end do
        return
    end subroutine sparse_kron

    subroutine speye(N, val, indx, jndx, ierr)
        !
        ! sparse identity matrix (double)
        ! input: N: size of matrix
        !        val: matrix values
        !        indx: i indices
        !        jndx: j indices
        !        ierr: internal error flag
        !
        implicit none
        integer, intent(in) :: N
        integer, dimension(:), allocatable, intent(out) :: indx, jndx
        integer, intent(out) :: ierr
        real (kind=8), intent(out), dimension(:), allocatable :: val
        integer :: i

        allocate(val(N), indx(N), jndx(N), STAT=ierr)
        if (ierr/=0) then
            return
        end if

        val = dble(1)
        do i=1,N
            indx(i) = i
            jndx(i) = i
        end do
        return
    end subroutine speye

    subroutine sparse_update(vals, indx, jndx, val, i, j, ierr)
        ! updates sparse matrix at entry i, j with value val.
        ! if value is zero, the array will be reallocated with the update value
        !
        ! inputs: vals: value vector of sparse matix
        !         indx: i indices
        !         jndx: j indices
        !         val: the value for the entry to update
        !         i: the i index of the entry to update
        !         j: the j index of the entry to update

        implicit none
        integer, intent(in) :: i, j
        real (kind=8), intent(in) :: val
        real (kind=8), intent(inout), dimension(:), allocatable :: vals
        integer, intent(inout), dimension(:), allocatable :: indx, jndx
        integer :: element_i, numel
        integer, intent(out) :: ierr
        integer, dimension(:), allocatable :: indx_new, jndx_new
        real (kind=8), dimension(:), allocatable :: vals_new

        numel = size(vals, 1)
        do element_i = 1, numel
            ! if the element exists in the array
            if (indx(element_i) == i .and. jndx(element_i) == j) then
                vals(element_i) = val
                return
            end if
        end do
        ! otherwise reallocate array
        allocate (indx_new(numel+1), jndx_new(numel+1), vals_new(numel+1), STAT=ierr)

        if (ierr/=0) then
            return
        end if

        indx_new(1:numel) = indx
        jndx_new(1:numel) = jndx
        vals_new(1:numel) = vals

        indx_new(numel+1) = i
        jndx_new(numel+1) = j
        vals_new(numel+1) = val
        call move_alloc (indx_new, indx)
        call move_alloc (jndx_new, jndx)
        call move_alloc (vals_new, vals)
        return
    end subroutine sparse_update

    subroutine sparse_multiply(A_vals, A_indx, A_jndx, x, y)
        !
        ! Computes matrix - matrix multiplication y <- A x, where x is a dense matrix, and alpha is a vector
        ! of column by column scaling factors
        !
        ! inputs: A_vals: values of A matrix
        !         A_indx: i indices of A matrix
        !         A_jndx: j indices of A matrix
        !         alpha: vector of columnwise scaling factors

        ! inputs
        real (kind=8), dimension(:), intent(in) :: A_vals
        real (kind=8), dimension(:, :), intent(in) :: x
        integer, dimension(:), intent(in) :: A_indx, A_jndx
        ! outputs
        real (kind=8), dimension(:, :), intent(out) :: y

        ! subroutine variables
        integer :: numel, i

        y=0
        numel = size(A_vals, 1)

        do i = 1, numel
            y(A_indx(i), :) = A_vals(i) * x(A_jndx(i), :) + y(A_indx(i), :)
        end do
        return
    end subroutine sparse_multiply

    subroutine scale_columns(x, alpha)
        ! scales the columns of x according to alpha.
        ! inputs: x: the matrix to rescale (overwritten)
        !         alpha: the factors to scale the columns of x by

        implicit none
        ! inputs
        real (kind=8), intent(inout), dimension(:,:) :: x
        real (kind=8), intent(in), dimension(:) :: alpha

        ! subroutine variables
        integer :: numalpha, i

        numalpha = size(alpha, 1)
        do i = 1, numalpha
            x(:, i) = alpha(i) * x(:, i)
        end do
    end subroutine scale_columns

    subroutine sparse_add(val1, indx1, jndx1, val2, indx2, jndx2, valout, indxout, jndxout, ierr)
        !
        ! adds two sparse matrices together
        !
        ! inputs: val1, indx1, jndx1: values, and i, j indices of first sparse matrix
        !         val2, indx2, jndx2: values, and i, j indices of second sparse matrix
        !
        ! outputs: valout, indxout, jndxout: values, i and j indices of output matrix
        implicit none
        ! input
        real (kind=8), dimension(:), intent(in) :: val1, val2
        integer, dimension(:), intent(in) :: indx1, indx2, jndx1, jndx2
        ! output
        real (kind=8), dimension(:), intent(out), allocatable :: valout
        integer, dimension(:), intent(out), allocatable :: indxout, jndxout
        integer, intent(out) :: ierr
        ! subroutine variables
        integer :: duplicates, numel1, numel2, numelout, i, j, current

        ! find output size
        duplicates = 0
        numel1 = size(val1, 1)
        numel2 = size(val2, 1)

        do i=1,numel1
            do j=1,numel2
                if (indx1(i) == indx2(j) .and. jndx1(i) == jndx2(j)) then
                    duplicates = duplicates + 1
                end if
            end do
        end do
        numelout = numel1 + numel2 - duplicates
        allocate (valout(numelout), indxout(numelout), jndxout(numelout), STAT=ierr)
        valout(1:numel1) = val1
        indxout(1:numel1) = indx1
        jndxout(1:numel1) = jndx1

        current = numel1 ! number of current elements
        do j=1,numel2
            do i=1,current
                if (indxout(i) == indx2(j) .and. indxout(i) == jndx2(j)) then
                    valout(i) = valout(i) + val2(j)
                    exit
                elseif (i == current) then
                    current = current + 1
                    valout(current) = val2(j)
                    indxout(current) = indx2(j)
                    jndxout(current) = jndx2(j)
                end if
            end do
        end do

    end subroutine sparse_add

    function sparse_get(val, indx, jndx, i, j, ierr)
        !
        ! retrieves a value from a sparse matrix
        ! inputs:
        ! val: value array
        ! indx: i index array
        ! jndx: j index array
        ! i: i index of retrieved value
        ! j: j index of retrieved value
        ! ierr: internal error flag
        implicit none
        ! input
        real (kind=8), dimension(:), intent(in) :: val
        integer, dimension(:), intent(in) :: indx, jndx
        integer, intent(in) :: i, j
        ! output
        integer, intent(out) :: ierr
        ! internal
        integer :: current_el

        do current_el = 1,size(val, 1)
            if (indx(current_el) == i .and. jndx(current_el) == j) then
                ierr = 0
                return val(current_el)
            end if
        end do
        ierr = -1
        return 0
    end function sparse_get

end module helpers