
module evolvePDE
    use config_m
    use ogpf
    use diff
    use helpers
    use sparse_matrices
    implicit none

    contains
        subroutine IMEX_evolve_impliciteuler_euler_2D (conf, ierr)
            ! UNDER DEVELOPMENT
            implicit none
            ! arguments
            type(gpf) :: gp
            type (config2d) :: conf
            integer, intent(out) :: ierr

            ! subroutine variables
            real (kind=8) :: t, dt, t_max, current_time
            real (kind=8), dimension(:, :), allocatable :: u
            integer :: timestep

            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), STAT=ierr)
            if (ierr /= 0)then
                return
            end if
            u = conf%IC

            timestep = 0
            t = 0
            dt = conf%dt
            t_max = conf%t_max

            do while (t < t_max)
                timestep = timestep + 1
                t = dble(timestep) * dt
            end do


        end subroutine IMEX_evolve_impliciteuler_euler_2D

        subroutine EX_evolve_euler_2D(conf, ierr)
            !
            ! evolves the specified configuration using the explicit euler method
            !
            ! inputs:
            !   conf: config2d type
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config2d), intent(inout) :: conf
            integer :: ierr

            ! runtime variables
            real (kind=8), allocatable, dimension(:, :) :: u, x_pass_through, u_update1, u_update2, plot_uu
            type (gpf) :: gp
            real (kind=8) :: t_max, t, dt, plot_time, current_time
            integer :: timestep
            type (csr_matrix) :: Laplacian_CSR
            type (coo_matrix) :: Laplacian_COO
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise grid
            call conf%make_mesh(ierr)
            allocate (x_pass_through(size(conf%xx), 2))
            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)

            ! initialise runge kutta variables
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_update1( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (u_update2( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (plot_uu( size(conf%xx, 1), size(conf%xx, 2) ))

            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC

            ! plotting parameters

            call cpu_time(plot_time)
            do while (t < t_max)

                call csr_multiply(Laplacian_CSR, u, u_update1)
                call scale_columns(u_update1, dble((/1, 30/)))
                call conf%explicit_rhs(u, x_pass_through, t, u_update2)

                ! BEGIN PLOTTING
                ! TODO: move this to a subroutine
                call cpu_time(current_time)
                if (current_time - plot_time > conf%plot_interval) then
                    Print *, t
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%title('u')
                    call gp%options("set terminal png size 1920,1080; set output 'tests/test_evolvePDE_euler.png'")
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    call cpu_time(plot_time)
                end if
                ! END PLOTTING

                ! euler step
                u = u + dt * (u_update1 + u_update2)

                timestep = timestep + 1
                t = dt * timestep
            end do

            ! END SUBROUTINE
            return
        end subroutine EX_evolve_euler_2D



end module evolvePDE


