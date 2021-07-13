! Created by  on 11/7/21.

program test_evolvePDE
    use config_m
    use evolvePDE
    implicit none
    call test_euler
    contains
        subroutine test_euler
            implicit none
            type (config2d) :: conf
            integer :: N, ierr
            real (kind=8), allocatable, dimension(:) :: x, y
            real (kind=8), allocatable, dimension(:, :) :: xx, yy
            real (kind=8) :: pi
            pi = dble(3.14159265358979323846264338327950288419716939937510)
            N = 60

            conf%BCx = 1
            conf%BCy = 1
            conf%t_max = 5
            conf%dt = 0.000001
            conf%plot_interval = 5*16

            allocate(conf%IC(N*N, 2))


            ! generate grid lines
            call my_linspace(dble(-1) + dble(2)/dble(N), dble(1), N, x,ierr)
            call my_linspace(dble(-1)+ dble(2)/dble(N), dble(1), N, y, ierr)

            call meshgrid(xx, yy, x, y, ierr)
            call random_number(conf%IC)

            call move_alloc( x, conf%x )
            call move_alloc( y, conf%y )
            conf%IC = conf%IC * 0.00001
            conf%IC(:, 1) = reshape(conf%IC(:, 1) + 0.2+1.3, (/N*N/))
            conf%IC(:, 2) = reshape(conf%IC(:, 2) + 1.3 / (1.5*1.5), (/N*N/))

            conf%explicit_rhs => rhs


            call EX_evolve_euler_2D(conf, ierr)

        end subroutine test_euler

        subroutine rhs(u_in, x_in, t, u_out)
            implicit none
            real (kind=8), dimension(:, :), intent(in) :: u_in, x_in
            real (kind=8), intent(in) :: t
            real (kind=8), intent(out), dimension(:, :) :: u_out

            real (kind=8) :: a, b, gamma

            a = 0.2
            b = 1.3
            gamma = 500

            u_out(:, 1) = gamma * (a - u_in(:, 1) + u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            u_out(:, 2) = gamma * (b -  u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            return
        end subroutine rhs
end program test_evolvePDE