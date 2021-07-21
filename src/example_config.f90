! Created by  on 20/7/21.

module user_config
    use config_m
    use helpers
    implicit none

    contains
        function get_config()
            implicit none
            type(config2d) :: get_config
            type (config2d) :: conf
            integer (ip) :: N, ierr
            real (rp), allocatable, dimension(:) :: x, y
            real (rp), allocatable, dimension(:, :) :: xx, yy
            real (rp) :: pi
            pi = dble(3.14159265358979323846264338327950288419716939937510)
            N = 1000_ip

            conf%BCx = 0_ip
            conf%BCy = 0_ip
            allocate (conf%DBCy_plus(2), conf%DBCy_minus(2), conf%DBCx_plus(2), conf%DBCx_minus(2))
            conf%DBCy_plus = real((/ 1, 0 /), kind=rp)
            conf%DBCy_minus = real((/ 0, 0 /), kind=rp)
            conf%DBCx_plus = real((/ 0, 0 /), kind=rp)
            conf%DBCx_minus = real((/ 0, 0 /), kind=rp)

            conf%DBCy_plus_mask = (/ .FALSE., .FALSE./)
            conf%DBCy_minus_mask = (/ .FALSE., .FALSE. /)
            conf%DBCx_plus_mask = (/ .FALSE., .FALSE. /)
            conf%DBCx_minus_mask = (/ .FALSE., .FALSE. /)


            conf%t_max = 1_ip
            conf%dt = 0.00001_rp
            conf%plot_interval = 40_ip

            conf%savenum = 100_ip
            conf%max_save_size = 10000000000_ip ! 10 gb
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
            get_config = conf
        end function get_config

        subroutine rhs(u_in, x_in, t, u_out)
            implicit none
            real (rp), dimension(:, :), intent(in) :: u_in, x_in
            real (rp), intent(in) :: t
            real (rp), intent(out), dimension(:, :) :: u_out
            real (rp) :: a, b, gamma

            a = 0.2
            b = 1.3
            gamma = 1000_ip

            u_out(:, 1) = gamma * (a - u_in(:, 1) + u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            u_out(:, 2) = gamma * (b -  u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            return
        end subroutine rhs

end module user_config