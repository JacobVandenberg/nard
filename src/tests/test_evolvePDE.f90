! Created by  on 11/7/21.

program test_evolvePDE
    use config_m
    use evolvePDE
    implicit none
    !call test_euler
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
            N = 64_ip

            conf%BCx = 1_ip
            conf%BCy = 1_ip
            conf%t_max = 1_ip
            conf%dt = 0.0001_rp
            conf%plot_interval = 5_ip*16_ip

            conf%savenum = 1000_ip
            conf%max_save_size = 10000000000_ip ! 10 gb
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

        subroutine test_euler3D
            implicit none
            type (config3d) :: conf
            integer (ip) :: N, ierr
            real (rp), allocatable, dimension(:) :: x, y, z
            real (rp), allocatable, dimension(:, :, :) :: xxx, yyy, zzz
            real (rp) :: pi
            pi = dble(3.14159265358979323846264338327950288419716939937510)
            N = 30_ip

            conf%BCx = 1_ip
            conf%BCy = 1_ip
            conf%BCz = 1_ip
            conf%t_max = 5_ip
            conf%dt = 0.00001_rp
            conf%plot_interval = 5_ip*16_ip

            allocate(conf%IC(N*N*N, 2_ip))


            ! generate grid lines
            call my_linspace(dble(-1) + dble(2)/dble(N), dble(1), N, x,ierr)
            call my_linspace(dble(-1)+ dble(2)/dble(N), dble(1), N, y, ierr)
            call my_linspace(dble(-1)+ dble(2)/dble(N), dble(1), N, z, ierr)

            call meshgrid3(x, y, z, xxx, yyy, zzz, ierr)
            call random_number(conf%IC)

            call move_alloc( x, conf%x )
            call move_alloc( y, conf%y )
            call move_alloc( z, conf%z )

            conf%IC = conf%IC * 0.00001
            conf%IC(:, 1) = reshape(conf%IC(:, 1) + 0.2+1.3, (/N*N*N/))
            conf%IC(:, 2) = reshape(conf%IC(:, 2) + 1.3 / (1.5*1.5), (/N*N*N/))

            conf%explicit_rhs => rhs

            allocate (conf%diffusion_consts(2))
            conf%diffusion_consts = dble((/ 1, 30 /))


            call IMEX_evolve_impliciteuler_euler_3D(conf, ierr)

        end subroutine test_euler3D

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
end program test_evolvePDE