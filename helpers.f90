! Created by  on 6/7/21.

module helpers
    implicit none
    public my_linspace, meshgrid3
    contains
    subroutine my_linspace(start, finish, N, x, ierr)
        implicit none
        real (kind=8), intent(in) :: start, finish
        integer, intent(in) :: N
        real (kind=8), intent(out), allocatable, dimension(:) :: x
        integer :: i
        real (kind=8) :: range
        integer, intent(out) :: ierr

        range = finish - start

        allocate (x(N), STAT=ierr)
        if (ierr/=0) then
            return
        end if

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

    recursive subroutine quicksort(a, first, last)
        implicit none
        real (kind=8)  a(*), x, t
        integer first, last
        integer i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
               i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
                t = a(i);  a(i) = a(j);  a(j) = t
                i=i+1
                j=j-1
        end do
        if (first < i-1) call quicksort(a, first, i-1)
        if (j+1 < last)  call quicksort(a, j+1, last)
    end subroutine quicksort

    recursive subroutine quicksort_int(a, first, last)
        implicit none
        integer  a(*), x, t
        integer first, last
        integer i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
               i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
                t = a(i);  a(i) = a(j);  a(j) = t
                i=i+1
                j=j-1
        end do
        if (first < i-1) call quicksort_int(a, first, i-1)
        if (j+1 < last)  call quicksort_int(a, j+1, last)
    end subroutine quicksort_int

    recursive subroutine order_int(a, first, last, key)
        ! gives the order of a list of integers
        ! a: list of integers to get order of
        ! first: the index to start sorting from
        ! last: the last index to sort from
        ! key (in/out): the order of the list a. input should be (/1,2,3,4,5/)
        implicit none
        integer  a(*), x, t, key(*)
        integer first, last
        integer i, j

        x = a( key((first+last) / 2) )
        i = first
        j = last
        do
            do while (a(key(i)) < x)
               i=i+1
            end do
            do while (x < a(key(j)))
                j=j-1
            end do
            if (i >= j) exit
            t = key(i);  key(i) = key(j);  key(j) = t
            i=i+1
            j=j-1
        end do
        if (first < i-1) call order_int(a, first, i-1, key)
        if (j+1 < last)  call order_int(a, j+1, last, key)
    end subroutine order_int

    subroutine arange_int(low, high, array, ierr)
        ! returns an array (/low, low+1, low+2, ..., high/)
        ! ierr: internal error flag
        implicit none
        integer, intent(in) :: low, high
        integer, intent(out) :: ierr
        integer, allocatable, intent(out), dimension(:) :: array
        integer :: ii

        if (high - low < 0) then
            ierr = -1
            return
        end if
        allocate(array(high-low+1), STAT=ierr)
        do ii = low,high
            array(ii+1-low) = ii
        end do
    end subroutine arange_int


end module helpers