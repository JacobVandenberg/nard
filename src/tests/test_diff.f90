! Created by  on 8/7/21.

program test_diff
    use diff
    use ogpf
    use sparse_matrices
    use precision
    implicit none
    call test_tridiag
    call test_sparse_tridiag
    Print *, 'test_double_diff_FD'
    call test_double_diff_FD
    Print *, 'test_single_diff_FD'
    call test_single_diff_FD
    Print *, 'test_laplacian2D'
    call test_laplacian2D
    Print *, 'test_advection2D'
    call test_advection2D
    Print *, 'kron_order'
    !call kron_order
    !call test_laplacian3D
    contains
        subroutine test_tridiag
            implicit none
            integer (ip) :: N = 4_ip, ierr
            real (rp), dimension(:, :), allocatable :: result

            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 0_ip
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate(result)

            N = 1_ip
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)

            N = 2_ip
            allocate (result(N, N), STAT=ierr)
            if (ierr /= 0_ip) then
                Print *, 'Allocation error for array of size (', N, 'x',N,'). Error code:', ierr
            end if

            call tridiag(dble(-2.), dble(1), dble(1), result)
            Print *, result
            deallocate (result)
        end subroutine test_tridiag

        subroutine test_sparse_tridiag
            implicit none
            integer (ip) :: N, ierr
            real(rp), allocatable, dimension(:) :: val
            integer (ip), dimension(:), allocatable :: indx, jndx

            N = 5_ip
            call sparse_tridiag(dble(-2), dble(1), dble(-1), N, val, indx, jndx, ierr)

            Print *, val
            Print *, indx
            Print *, jndx

            deallocate(val, indx, jndx)


        end subroutine test_sparse_tridiag

        subroutine test_double_diff_FD
            implicit none
            type (coo_matrix) :: mat
            real (rp), dimension(:), allocatable :: x
            integer (ip) :: ierr

            x = linspace(dble(0), dble(4), 5_ip)
            mat = double_diff_FD(x, 1_ip, ierr)
            Print *, mat%vals
            deallocate(x)

            x = linspace(dble(0), dble(4), 5_ip)
            mat = double_diff_FD(x, 0_ip, ierr)
            Print *, mat%vals

        end subroutine test_double_diff_FD

        subroutine test_single_diff_FD
            implicit none
            type (coo_matrix) :: mat
            real (rp), dimension(:), allocatable :: x
            integer (ip) :: ierr

            x = linspace(dble(0), dble(4), 5_ip)
            mat = single_diff_FD(x, 1_ip, ierr)
            Print *, mat%vals
            Print *, mat%indx
            Print *, mat%jndx

            deallocate(x)

            x = linspace(dble(0), dble(4), 5_ip)
            mat = single_diff_FD(x, 0_ip, ierr)
            Print *, mat%vals

        end subroutine test_single_diff_FD

        subroutine test_laplacian2D
            implicit none
            type (coo_matrix) :: mat
            type (gpf) :: gp
            real (rp), dimension(:), allocatable :: x, y
            real (rp), dimension(:, :), allocatable :: xx, yy, Luu, uu_plot
            integer (ip) :: ierr, N
            N = 50_ip
            x = linspace(dble(-1), dble(1), N)
            y = linspace(dble(-1), dble(1), N)

            call meshgrid(xx, yy, x, y, ierr)

            mat = laplacian2D(x, y, 0_ip, 0_ip, ierr)
            allocate(Luu(N * N, 1), uu_plot(N, N))

            Print *, mat%n, mat%m

            call coo_multiply(mat, reshape(cos(3.14159 * xx) + cos(3.14159 * yy), (/N*N, 1_ip/)), Luu)
            uu_plot = reshape(Luu, (/N, N/))
            call gp%surf(xx, yy, uu_plot)

        end subroutine test_laplacian2D

        subroutine test_advection2D
            implicit none
            type (coo_matrix) :: mat
            type (gpf) :: gp
            real (rp), dimension(:), allocatable :: x, y
            real (rp), dimension(:, :), allocatable :: xx, yy, Luu, uu_plot
            integer (ip) :: ierr, N
            N = 50_ip
            x = linspace(dble(-1), dble(1), N)
            y = linspace(dble(-1), dble(1), N)

            call meshgrid(xx, yy, x, y, ierr)

            mat = advection2D(x, y, 1_ip, 1_ip, ierr)
            allocate(Luu(N * N, 1), uu_plot(N, N))

            Print *, mat%n, mat%m

            call coo_multiply(mat, reshape(cos(3.14159 * xx) + cos(3.14159 * yy), (/N*N, 1_ip/)), Luu)
            uu_plot = reshape(Luu, (/N, N/))
            call gp%surf(xx, yy, uu_plot)

        end subroutine test_advection2D

        subroutine kron_order
            ! which order should i kron the double derivative matrix in order to make a 3d laplacian?
            !

            implicit none

            ! declarations
            integer (ip) :: ierr
            real (rp), dimension(:), allocatable :: x, y, z, u_plot
            real (rp), dimension(:, :), allocatable :: u, Lu
            real (rp), dimension(:, :, :), allocatable :: xxx, yyy, zzz, uuu, Luuu
            type (coo_matrix) :: DDz, DDx, DDy, ex, ey, ez, temp
            type (gpf) :: gp


            x = linspace(-1.0_rp, 1.0_rp, 10_ip)
            y = linspace(-1.0_rp, 1.0_rp, 20_ip)
            z = linspace(-1.0_rp, 1.0_rp, 30_ip)

            call meshgrid3(x, y, z, xxx, yyy, zzz, ierr)
            allocate(uuu(10, 20, 30), u(10*20*30, 1), Luuu(10, 20, 30), Lu(10*20*30, 1))


            DDx = double_diff_FD(x, 0_ip, ierr)
            DDy = double_diff_FD(y, 0_ip, ierr)
            DDz = double_diff_FD(z, 0_ip, ierr)

            ex = speye(size(x, 1, kind=ip), ierr)
            ey = speye(size(y, 1, kind=ip), ierr)
            ez = speye(size(z, 1, kind=ip), ierr)


            ! Z
            uuu = xxx**2_rp + 17_rp*yyy**2_rp + 19_rp*zzz**2_rp

            temp = sparse_kron(ex, ey, ierr)
            DDz = sparse_kron(DDz, temp, ierr)
            Lu = 0
            call coo_multiply(DDz, reshape(uuu, (/10*20*30, 1/)), Lu)
            Luuu = reshape(Lu, (/10, 20, 30/))

            allocate(u_plot(30))
            u_plot(:) = Luuu(5, 10, :)
            call gp%title('z')
            call gp%plot(z(2:29), u_plot(2:29))
            deallocate(u_plot)

            ! Y
            uuu = xxx**2_rp + 17_rp*yyy**2_rp + 19_rp*zzz**2_rp

            temp = sparse_kron(ez, ex, ierr)
            DDy = sparse_kron(temp, DDy, ierr)
            !Print *, DDy%vals
            !Print *, DDy%indx
            !Print *, DDy%jndx

            Lu = 0
            call coo_multiply(DDy, reshape(uuu, (/10*20*30, 1/)), Lu)
            Luuu = reshape(Lu, (/10, 20, 30/))

            allocate(u_plot(20))
            u_plot(:) = Luuu(5, :, 15)
            call gp%title('y')
            call gp%plot(y(2:19), u_plot(2:19))
            deallocate(u_plot)

            ! X
            uuu = xxx**2_rp + 17_rp*yyy**2_rp + 19_rp*zzz**2_rp

            temp = sparse_kron(DDx, ey, ierr)

            DDx = sparse_kron(ez, temp, ierr)

            call coo_multiply(DDx, reshape(uuu, (/10*20*30, 1/)), Lu)
            Luuu = reshape(Lu, (/10, 20, 30/))

            allocate(u_plot(10))
            u_plot(:) = Luuu(:, 10, 15)
            call gp%title('x')
            call gp%plot(x(2:9), u_plot(2:9))
            deallocate(u_plot)
            return
        end subroutine kron_order

        subroutine test_laplacian3D
            implicit none
            ! declarations
            integer (ip) :: ierr
            real (rp), dimension(:), allocatable :: x, y, z, u_plot
            real (rp), dimension(:, :), allocatable :: u, Lu
            real (rp), dimension(:, :, :), allocatable :: xxx, yyy, zzz, uuu, Luuu
            type (coo_matrix) :: L
            type (gpf) :: gp
            real (rp), parameter :: pi=3.1415926535


            x = linspace(-1.0_rp + (2.0_rp/10.0_rp), 1.0_rp, 10_ip)
            y = linspace(-1.0_rp+ (2.0_rp/20.0_rp), 1.0_rp, 20_ip)
            z = linspace(-1.0_rp+ (2.0_rp/30.0_rp), 1.0_rp, 30_ip)

            call meshgrid3(x, y, z, xxx, yyy, zzz, ierr)
            allocate(uuu(20, 10, 30), u(10*20*30, 1), Luuu(20, 10, 30), Lu(10*20*30, 1))

            uuu = cos(xxx * pi) * cos(yyy * 2 * pi) * cos(zzz * 3 * pi)

            u = reshape(uuu, (/10*20*30, 1/))

            L = laplacian3D(x, y, z, 1_ip,1_ip, 1_ip, ierr)
            call coo_multiply(L, u, Lu)
            Luuu = reshape(Lu, (/20, 10, 30/))


            call gp%title('x')
            call gp%plot(x, Luuu(20, :, 30))
            call gp%title('y')
            call gp%plot(y, Luuu(:, 10, 30))
            call gp%title('x')
            call gp%plot(z, Luuu(20, 10, :))


        end subroutine test_laplacian3D
end program test_diff