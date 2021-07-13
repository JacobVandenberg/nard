real (kind=8) function sparse_get(val, indx, jndx, i, j, ierr)
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
                sparse_get = val(current_el)
                return
            end if
        end do
        ierr = -1
        sparse_get = dble(0)
        return
    end function sparse_get


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

subroutine speye_(N, val, indx, jndx, ierr)
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
    end subroutine speye_

    subroutine sparse_kron_(val1, indx1, jndx1, n1, m1, val2, indx2, jndx2, n2, m2, valout, indxout, jndxout, nout, mout, ierr)
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
    end subroutine sparse_kron_

subroutine test_sparse_kron_
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
            call speye_(2, val2, indx2, jndx2, ierr)
            Print *, ierr
            Print *, val2, indx2, jndx2

            call sparse_kron_(val1, indx1, jndx1, n1, m1, val2, indx2, jndx2, n2, m2, valout, indxout, jndxout, nout, mout, ierr)
            Print *, ierr
            Print *, valout
            Print *,indxout
            Print *,jndxout
            deallocate(indx1, indx2, jndx1, jndx2, val1, val2, indxout, jndxout, valout)

        end subroutine test_sparse_kron_

subroutine test_speye_
            implicit none
            integer :: N, ierr
            real (kind=8), dimension(:), allocatable :: val
            integer, dimension(:), allocatable :: indx, jndx

            N = 4

            call speye_(N, val, indx, jndx, ierr)
            Print *, val, indx, jndx
            deallocate(val, indx, jndx)

        end subroutine test_speye_

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

            allocate(y(N, 2), DDy(N, 2), y_plot(N))

            call my_linspace(dble(-1), dble(1), N, x, ierr)

            y(:, 1) = cos(dble(10) * x)
            y(:, 2) = sin(dble(10) * x)

            call double_diff_FD(x, vals, indx, jndx, ierr)
            call sparse_multiply(vals, indx, jndx, y, DDy)

            DDy(1, :) = 0
            DDy(N, :) = 0
            y_plot = DDy(:, 2)
            call gp%plot(x, y_plot)

            deallocate(x, vals, y, DDy, indx, jndx)

            call speye_(4, vals, indx, jndx, ierr)
            call sparse_update(vals, indx, jndx, dble(1), 4, 1, ierr)
            allocate(y(4, 1), DDy(4, 1))
            y = 1
            DDy = 100
            call sparse_multiply(vals, indx, jndx, y, DDy)
            Print *, DDy(:, 1)


        end subroutine test_sparse_multiply

subroutine euler (conf)
            !
            ! evolves configuration using eulers method
            ! inputs: conf: configuration structure
            !
            ! outputs: write to file specified by conf
            !

            implicit none
            type(gpf) :: gp
            type (config) :: conf
            real (kind = 8), dimension(:, :), allocatable :: u, u_update, xx, yy, plot_uu, test
            integer, dimension(:), allocatable :: DDx_indx, DDx_jndx, DDy_indx, DDy_jndx, L_indx, L_jndx, eye_indx, eye_jndx
            real (kind = 8), dimension(:), allocatable :: DDy_vals, DDx_vals, L_vals, x, y, eye_vals
            integer :: ierr, N, mout, nout
            real (kind=8) :: t, t_max, dt, plot_time, plot_interval, current_time

            integer, dimension(2) :: u_shape

            N = 100
            plot_interval = 5 * 16

            call cpu_time(plot_time)

            ! allocate the result to the size of the initial conditions
            u_shape = shape(conf%IC)
            allocate (u_update(N*N, 2), u(N*N, 2), plot_uu(N, N), STAT=ierr)

            if (ierr /= 0) then
                Print *, 'allocation error for initial solution array'
                return
            end if
            ! generate grid lines
            call my_linspace(dble(-1) + dble(2)/dble(N), dble(1), N, x,ierr)
            call my_linspace(dble(-1)+ dble(2)/dble(N), dble(1), N, y, ierr)


            call meshgrid(xx,yy,x,y, ierr)

            ! generate laplacian matrix
            ! diff wrt x
            Print *, 'starting'
            call double_diff_FD(x, DDx_vals, DDx_indx, DDx_jndx, ierr)
            ! periodic BCS
            call sparse_update(DDx_vals, DDx_indx, DDx_jndx, sparse_get(DDx_vals, DDx_indx, DDx_jndx, 1, 2, ierr), 1, N, ierr)
            call sparse_update(DDx_vals, DDx_indx, DDx_jndx, sparse_get(DDx_vals, DDx_indx, DDx_jndx, N, N-1, ierr), N, 1, ierr)

            call double_diff_FD(y, DDy_vals, DDy_indx, DDy_jndx, ierr)
            ! periodic BCS
            call sparse_update(DDy_vals, DDy_indx, DDy_jndx, sparse_get(DDy_vals, DDy_indx, DDy_jndx, 1, 2, ierr), 1, N, ierr)
            call sparse_update(DDy_vals, DDy_indx, DDy_jndx, sparse_get(DDy_vals, DDy_indx, DDy_jndx, N, N-1, ierr), N, 1, ierr)

            ! identity matrix
            call speye_(N, eye_vals, eye_indx, eye_jndx, ierr)

            Print *, 'generated diff matrices'
            ! kron tensor product
            !
            !DDx = kron(DDx, eye(length(mesh.y)));
            !DDy = kron(eye(length(mesh.x)), DDy);
            call sparse_kron_(DDx_vals, DDx_indx, DDx_jndx, N, N, eye_vals, eye_indx, eye_jndx, N, N, &
                    L_vals, L_indx, L_jndx, nout, mout, ierr)
            call move_alloc (L_vals, DDx_vals)
            call move_alloc (L_indx, DDx_indx)
            call move_alloc (L_jndx, DDx_jndx)
            Print *, 'kron DDx done'

            call sparse_kron_(eye_vals, eye_indx, eye_jndx, N, N, DDy_vals, DDy_indx, DDy_jndx, N, N, &
                    L_vals, L_indx, L_jndx, nout, mout, ierr)
            call move_alloc (L_vals, DDy_vals)
            call move_alloc (L_indx, DDy_indx)
            call move_alloc (L_jndx, DDy_jndx)
            Print *, 'kron DDy done'



            ! L = DDx + DDy
            call sparse_add(DDx_vals, DDx_indx, DDx_jndx, DDy_vals, DDy_indx, DDy_jndx, L_vals, L_indx, L_jndx, ierr)
            Print *, 'sparse add done'
            call coo_to_csr(L_vals, L_indx, L_jndx, ierr)

            t = 0
            dt = dble(0.5) * dble(10)**dble(-6)
            t_max = 5
            u = conf%IC

            do while (t < t_max)

                call csr_multiply(L_vals, L_indx, L_jndx, u, u_update)
                call scale_columns(u_update, (/dble(1), dble(30)/))

                u_update(:, 1) = u_update(:, 1) + dble(500) *(dble(0.2) - u(:, 1) + u(:, 1)*u(:, 1)*u(:, 2))
                u_update(:, 2) = u_update(:, 2) + dble(500) * (dble(1.3) - u(:, 1)*u(:, 1)*u(:, 2))

                call cpu_time(current_time)
                if (current_time - plot_time > plot_interval) then
                    !call gp%reset()
                    Print *, t
                    plot_uu = reshape(u(:, 1), (/N, N /))
                    call gp%title('u')
                    call gp%options("set terminal png size 1920,1080; set output 'tests/test_evolvePDE_euler.png'")
                    call gp%contour(xx,yy,plot_uu, palette='jet')
                    call cpu_time(plot_time)
                end if
                u = u + dt * u_update
                t = t + dt
            end do

            deallocate(u)


        end subroutine euler

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
        real (kind=8), dimension(:), intent(inout), allocatable :: val1, val2
        integer, dimension(:), intent(inout), allocatable :: indx1, indx2, jndx1, jndx2
        ! output
        real (kind=8), dimension(:), intent(out), allocatable :: valout
        integer, dimension(:), intent(out), allocatable :: indxout, jndxout
        integer, intent(out) :: ierr
        ! subroutine variables
        integer :: duplicates, numel1, numel2, numelout, i, j, current, num_cols
        integer, allocatable, dimension(:) :: key

        ! find output size
        duplicates = 0
        numel1 = size(val1, 1)
        numel2 = size(val2, 1)
        numelout = numel1+numel2

        num_cols = maxval((/maxval(jndx1), maxval(jndx2)/))
        allocate (indxout(numelout), STAT = ierr)
        indxout(1:numel1) = (indx1 - 1) * (num_cols) + jndx1 - 1
        indxout(numel1+1 : numelout) = (indx2 - 1) * (num_cols) + jndx2 - 1

        call quicksort_int(indxout, 1, numelout)
        duplicates = count(indxout(2:numelout) == indxout(1:numelout-1))
        deallocate(indxout)

        ! sort both input matrices.
        call arange_int(1, numel1, key, ierr)
        call order_int((indx1 - 1) * (num_cols) + jndx1 - 1, 1, numel1, key)
        ! sort into temporary array
        allocate (indxout(numel1), jndxout(numel1), valout(numel1))
        do i= 1,numel1
            indxout(i) = indx1(key(i))
            jndxout(i) = jndx1(key(i))
            valout(i) = val1(key(i))
        end do
        call move_alloc(indxout, indx1)
        call move_alloc(jndxout, jndx1)
        call move_alloc(valout, val1)
        deallocate(key)

        call arange_int(1, numel2, key, ierr)
        call order_int((indx2 - 1) * (num_cols) + jndx2 - 1, 1, numel2, key)
        ! sort into temporary array
        allocate (indxout(numel2), jndxout(numel2), valout(numel2), STAT=ierr)
        do i= 1,numel2
            indxout(i) = indx2(key(i))
            jndxout(i) = jndx2(key(i))
            valout(i) = val2(key(i))
        end do
        call move_alloc(indxout, indx2)
        call move_alloc(jndxout, jndx2)
        call move_alloc(valout, val2)
        deallocate(key)

        numelout = numel1 + numel2 - duplicates
        allocate (valout(numelout), indxout(numelout), jndxout(numelout), STAT=ierr)

        i=1; j=1
        current = 1
        do while (current <= numelout)
            if (i > numel1) then
                valout(current) = val2(j)
                indxout(current) = indx2(j)
                jndxout(current) = jndx2(j)
                j = j + 1
            elseif (j > numel2) then
                valout(current) = val2(j)
                indxout(current) = indx2(j)
                jndxout(current) = jndx2(j)
                j = j + 1
            elseif (indx1(i) == indx2(j) .and. jndx1(i) == jndx2(j)) then
                valout(current) = val1(i) + val2(j)
                indxout(current) = indx1(i)
                jndxout(current) = jndx1(i)
                i = i + 1; j = j + 1
            elseif  ((indx1(i) - 1) * num_cols + jndx1(i) - 1 < (indx2(j) - 1) * num_cols + jndx2(j) - 1) then
                valout(current) = val1(i)
                indxout(current) = indx1(i)
                jndxout(current) = jndx1(i)
                i = i + 1
            else
                valout(current) = val2(j)
                indxout(current) = indx2(j)
                jndxout(current) = jndx2(j)
                j = j + 1
            end if
            current = current + 1
        end do

    end subroutine sparse_add


        subroutine test_sparse_add
            implicit none
            real (kind=8), allocatable, dimension(:) :: val1, val2, valout
            integer, allocatable, dimension(:) :: indx1, indx2, indxout, jndx1, jndx2, jndxout
            integer :: N, ierr

            N = 4
            !call speye_(N, val1, indx1, jndx1, ierr)
            !call speye_(N, val2, indx2, jndx2, ierr)
            !call sparse_update(val1, indx1, jndx1, dble(10), 1, 2, ierr)
            !call sparse_update(val1, indx1, jndx1, dble(2), 1, 1, ierr)

            Print *, val1
            Print *, indx1
            Print *, jndx1
            Print *, val2
            Print *, indx2
            Print *, jndx2


            call sparse_add(val1, indx1, jndx1, val2, indx2, jndx2, valout, indxout, jndxout, ierr)
            Print *, valout
            Print *, indxout
            Print *, jndxout
        end subroutine test_sparse_add

    subroutine sparse_multiply(A_vals, A_indx, A_jndx, x, y)
        !
        ! Computes matrix - matrix multiplication y <- A x, where x is a dense matrix, and alpha is a vector
        ! of column by column scaling factors
        !
        ! inputs: A_vals: values of A matrix
        !         A_indx: i indices of A matrix
        !         A_jndx: j indices of A matrix

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

subroutine coo_to_csr(valin, indxin, jndxin, ierr)
        ! converts sparse matrix from coo (coordinate) to
        ! CSR (compressed sparse row)
        !
        ! ARGUMENTS:
        !   valin (inout): the values in the matrix
        !   indxin (inout): the i indices (this will change to the row pointer array)
        !   jndxin (inout): the j indices (or which column the entry is in)
        !   ierr (out): error flag

        implicit none
        ! arguments
        real (kind=8), intent(inout), dimension(:), allocatable :: valin
        integer, intent(inout), dimension(:), allocatable :: indxin, jndxin
        integer, intent(out) :: ierr
        ! internals
        integer, dimension(:), allocatable :: key
        integer :: num_cols, numel, ii, current_row
        integer, dimension(:), allocatable :: temp_int
        real (kind=8), dimension(:), allocatable :: temp_real



        ! first get the order of the elements
        num_cols = maxval(jndxin)
        numel = size(valin, 1)
        call arange_int(1, numel, key, ierr)
        call order_int( (indxin-1) * num_cols + jndxin - 1, 1, numel, key)

        ! reindex the column indices and the values
        allocate ( temp_real(numel),temp_int(numel), STAT=ierr )
        do ii=1,numel
            temp_real(ii) = valin(key(ii))
            temp_int(ii) = jndxin(key(ii))
        end do
        call move_alloc(temp_real, valin)
        call move_alloc(temp_int, jndxin)

        ! find row vector
        num_cols = maxval(indxin)  ! well, rows now
        allocate (temp_int(num_cols+1), STAT=ierr)

        current_row = 0
        do ii = 1, numel
            if (indxin(key(ii)) > current_row) then
                temp_int(current_row+1:indxin(key(ii))) = ii
                current_row = indxin(key(ii))
            end if
        end do
        temp_int(num_cols+1) = numel+1
        call move_alloc(temp_int, indxin)
    end subroutine coo_to_csr

subroutine test_coo_to_csr
            implicit none
            integer, allocatable, dimension(:) :: indx, jndx
            real (kind=8), allocatable, dimension(:) :: vals
            integer :: ierr

            allocate (indx(3), jndx(3), vals(3), STAT=ierr)
            indx = (/ 3, 1, 1, 1, 1/)
            jndx = (/ 2, 1, 3, 4, 5 /)
            vals = (/dble(3), dble(1), dble(2), dble(4), dble(5)/)

            call coo_to_csr(vals, indx, jndx, ierr)
            Print *, vals
            Print *, indx
            Print *, jndx

        end subroutine test_coo_to_csr

subroutine test_csr_multiply
            implicit none
            integer, allocatable, dimension(:) :: indx, jndx
            real (kind=8), allocatable, dimension(:) :: vals
            integer :: ierr
            real (kind=8), allocatable, dimension(:, :) :: x, y

            allocate (indx(3), jndx(3), vals(3), STAT=ierr)
            indx = (/ 3, 1, 1, 1, 1/)
            jndx = (/ 2, 1, 3, 4, 5 /)
            vals = (/dble(3), dble(1), dble(2), dble(4), dble(5)/)

            call coo_to_csr(vals, indx, jndx, ierr)

            allocate(x(5, 2), y(3, 2), STAT = ierr)
            x = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /), (/5, 2/))
            Print *, x(:, 1)
            call csr_multiply(vals, indx, jndx, x, y)
            Print *, y(:, 1)
            Print *, y(:, 2)

        end subroutine test_csr_multiply