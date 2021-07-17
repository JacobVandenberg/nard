! Created by  on 17/7/21.

program test_bigdata
    use config_m
    use evolvePDE
    implicit none
    call test_euler
    contains
        subroutine test_euler
            implicit none
            type (config2d) :: conf
            integer (ip) :: N, ierr
            real (rp), allocatable, dimension(:) :: x, y
            real (rp), allocatable, dimension(:, :) :: xx, yy
            real (rp) :: pi
            pi = dble(3.14159265358979323846264338327950288419716939937510)
            N = 200_ip

            conf%BCx = 1_ip
            conf%BCy = 1_ip
            conf%t_max = 1_ip
            conf%dt = 0.000001_rp
            conf%plot_interval = 5_ip*16_ip

            conf%savenum = 50000_ip
            conf%max_save_size = 33000000000_ip ! 33 gb
            !allocate(conf%savefilename(7))
            conf%savefilename = "test.h5"

            allocate(conf%IC(N*N, 2_ip))


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

            allocate (conf%diffusion_consts(2))
            conf%diffusion_consts = dble((/ 1, 30 /))


            call IMEX_evolve_impliciteuler_euler_2D(conf, ierr)

        end subroutine test_euler

        subroutine rhs(u_in, x_in, t, u_out)
            implicit none
            real (rp), dimension(:, :), intent(in) :: u_in, x_in
            real (rp), intent(in) :: t
            real (rp), intent(out), dimension(:, :) :: u_out

            real (rp) :: a, b, gamma

            a = 0.2
            b = 1.3
            gamma = 500_ip

            u_out(:, 1) = gamma * (a - u_in(:, 1) + u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            u_out(:, 2) = gamma * (b -  u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            return
        end subroutine rhs

end program test_bigdata