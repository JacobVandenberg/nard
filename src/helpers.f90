! Created by  on 6/7/21.

module helpers
    use precision
    implicit none
    public my_linspace, meshgrid3
    contains
    subroutine my_linspace(start, finish, N, x, ierr)
        implicit none
        real (rp), intent(in) :: start, finish
        integer (ip), intent(in) :: N
        real (rp), intent(out), allocatable, dimension(:) :: x
        integer (ip) :: i
        real (rp) :: range
        integer (ip), intent(out) :: ierr

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

    function linspace(a,b,n_elements)
        !..............................................................................
        !   returns a linearly spaced vector with n points in [a, b]
        !   if n is omitted, 100 points will be considered
        !..............................................................................

        real(rp), intent(in)           :: a
        real(rp), intent(in)           :: b
        integer (ip),  intent(in), optional :: n_elements
        real(rp), allocatable          :: linspace(:)

        !   Local vars
        real(rp) :: dx
        integer (ip)  :: i
        integer (ip)  :: n
        integer (ip)  :: ierr

        if (present(n_elements)) then
            if (n_elements <=1 ) then
                print*, "linspace procedure: Error: wrong value of n_elements, use an n_elements > 1"
                stop
            end if
            n=n_elements
        else
            n=100
        end if

        allocate(linspace(n), stat=ierr)
        if (ierr /= 0) then
            print*, "linspace procedure: Fatal Error, Allocation failed in linspace function"
            stop
        end if

        dx=(b-a)/real((n-1),rp)
        linspace=[(i*dx+a, i=0,n-1)]

    end function linspace

    subroutine meshgrid3(x, y, z, xxx, yyy, zzz, ierr)
        implicit none
        real (rp), intent(in), dimension(:) :: x, y, z
        real (rp), intent(out), dimension(:, :, :), allocatable :: xxx, yyy, zzz
        integer (ip), intent(out) :: ierr
        integer (ip) :: i, j, k, nx, ny, nz

        nx = size(x)
        ny = size(y)
        nz = size(z)

        allocate(xxx(ny, nx, nz), yyy(ny, nx, nz), zzz(ny, nx, nz), STAT=ierr)

        if (ierr/=0) then
            return
        end if

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
        real (rp), intent(inout), dimension(:,:) :: x
        real (rp), intent(in), dimension(:) :: alpha

        ! subroutine variables
        integer (ip) :: numalpha, i

        numalpha = size(alpha, 1)
        do i = 1, numalpha
            x(:, i) = alpha(i) * x(:, i)
        end do
    end subroutine scale_columns

    recursive subroutine quicksort(a, first, last)
        implicit none
        real (rp)  a(*), x, t
        integer (ip) first, last
        integer (ip) i, j

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
        integer (ip)  a(*), x, t
        integer (ip) first, last
        integer (ip) i, j

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
        ! gives the order of a list of integer (ips
        ! a: list of integer (ips to get order of
        ! first: the index to start sorting from
        ! last: the last index to sort from
        ! key (in/out): the order of the list a. input should be (/1,2,3,4,5/)
        implicit none
        integer (ip)  a(*), x, t, key(*)
        integer (ip) first, last
        integer (ip) i, j

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
        integer (ip), intent(in) :: low, high
        integer (ip), intent(out) :: ierr
        integer (ip), allocatable, intent(out), dimension(:) :: array
        integer (ip) :: ii

        if (high - low < 0) then
            ierr = -1
            return
        end if
        allocate(array(high-low+1), STAT=ierr)
        do ii = low,high
            array(ii+1-low) = ii
        end do
    end subroutine arange_int

    subroutine meshgrid(x,y,xgv,ygv, ierr)
        !..............................................................................
        !meshgrid generate mesh grid over a rectangular domain of [xmin xmax, ymin, ymax]
        ! Inputs:
        !     xgv, ygv are grid vectors in form of full grid data
        ! Outputs:
        !     X and Y are matrix each of size [ny by nx] contains the grid data.
        !     The coordinates of point (i,j) is [X(i,j), Y(i,j)]
        !     ierr: The error flag
        !     """
        !     # Example
        !     # call meshgrid(X, Y, [0.,1.,2.,3.],[5.,6.,7.,8.])
        !     # X
        !     # [0.0, 1.0, 2.0, 3.0,
        !     #  0.0, 1.0, 2.0, 3.0,
        !     #  0.0, 1.0, 2.0, 3.0,
        !     #  0.0, 1.0, 2.0, 3.0]
        !     #
        !     #Y
        !     #[ 5.0, 5.0, 5.0, 5.0,
        !     #  6.0, 6.0, 6.0, 6.0,
        !     #  7.0, 7.0, 7.0, 7.0,
        !     #  8.0, 8.0, 8.0, 8.0]
        !..............................................................................
        ! Rev 0.2, Feb 2018
        ! New feature added: xgv and ygv as full grid vector are accepted now

        ! Arguments
        real(rp), intent(out), allocatable  :: x(:,:)
        real(rp), intent(out), allocatable  :: y(:,:)
        real(rp), intent(in)                :: xgv(:) ! x grid vector [start, stop, step] or [start, stop]
        real(rp), intent(in),  optional     :: ygv(:) ! y grid vector [start, stop, step] or [start, stop]
        integer(ip),  intent(out), optional     :: ierr   ! the error value

        ! Local variables
        integer(ip):: sv
        integer(ip):: nx
        integer(ip):: ny
        logical:: only_xgv_available

        ! Initial setting
        only_xgv_available  = .false.
        sv=0 !Assume no error

        nx=size(xgv, dim=1)

        if (present(ygv)) then
            ny = size(ygv, dim=1)
        else
            only_xgv_available=.true.
            ny=nx
        end if

        allocate(x(ny,nx),y(ny,nx),stat=sv)
        if (sv /=0) then
            print*, "allocataion erro in meshgrid"
            stop
        end if

        x(1,:)    = xgv
        x(2:ny,:) = spread(xgv, dim=1, ncopies=ny-1)

        if (only_xgv_available) then
            y=transpose(x)
        else
            y(:,1)    = ygv
            y(:,2:nx) = spread(ygv,dim=2,ncopies=nx-1)
        end if

        if (present(ierr)) then
            ierr=sv
        end if

    end subroutine meshgrid


end module helpers