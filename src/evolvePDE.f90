
module evolvePDE
    use config_m
    use ogpf
    use diff
    use data
    use helpers
    use sparse_matrices
    use precision
    use hdf5
    implicit none

    contains
        subroutine IMEX_evolve_impliciteuler_euler_2D (conf, ierr)
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config2d), intent(inout) :: conf
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_update1, u_update2, plot_uu,&
                    u_update_temp1, u_update_temp2, adv_const_field, u_adv1, dudxx
            real (rp), allocatable, dimension(:) :: t_vec_out
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time, update_time
            integer (ip) :: timestep, dummy_i, max_timesteps, savenum, prev_timestep, plot_number
            type (csr_matrix) :: Laplacian_CSR, Advection_CSR
            type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp, Advection_COO
            type (sparse_lu_decomposition), dimension(:), allocatable :: implicit_euler_LU
            character (len=1024) :: timestep_string
            logical :: advection

            ! HDF5 variables
            type (h5file) :: save_file
            type (chunked_3D_space) :: result_space
            integer (ip), dimension(3) :: file_dimensions
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise grid
            call conf%make_mesh(ierr)
            if (ierr /= 0) then
                Print *, "Error making mesh. Error code: ", ierr
                call exit()
            end if

            allocate (x_pass_through(size(conf%xx), 2), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating x pass through. Error code: ", ierr
                call exit()
            end if

            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)
            Advection_COO = advection2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Advection_CSR = coo_to_csr(Advection_COO, ierr)

            ! initialise runge kutta variables
            Print *, "Allocating u and u update variables", size(conf%IC, 1), size(conf%IC, 2), size(conf%xx, 1), size(conf%xx, 2)
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_update1( size(conf%IC, 1), size(conf%IC, 2) ),&
                adv_const_field( size(conf%IC, 1), size(conf%IC, 2) ), u_adv1( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating u and u update variables:" , ierr
                call exit()
            end if
            allocate (u_update2( size(conf%IC, 1), size(conf%IC, 2) ),&
                    u_update_temp1( size(conf%IC, 1), 1 ), u_update_temp2( size(conf%IC, 1), 1 ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating u and u update variables:" , ierr
                call exit()
            end if
            allocate (plot_uu( size(conf%xx, 1), size(conf%xx, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating u and u update variables:" , ierr
                call exit()
            end if
            allocate (dudxx( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating dudxx pass through:" , ierr
                call exit()
            end if


            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC
            savenum = conf%savenum
            advection = conf%advection

            ! initialise saving variables
            max_timesteps = int( t_max / dt)

            ! estimate output size
            Print *, "Initialising save file"
            if ( 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1) > conf%max_save_size) then
                Print *, "ERROR: Desired output file exceeds maximum file size."
                Print *, "Desired bytes: ", 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1)
                Print *, "Maximum bytes: ", conf%max_save_size
                Print *, "increase conf%max_save_size or reduce conf%savenum to resolve"
                ierr = 2
                return
            end if
            call save_file%new_file(conf%savefilename, ierr)

            file_dimensions = (/ size(conf%IC, 1, kind=HSIZE_T), size(conf%IC, 2, kind=HSIZE_T),&
                    INT(savenum+1, KIND=HSIZE_T)/)
            allocate (t_vec_out(INT(savenum+1, KIND=HSIZE_T)) )

            t_vec_out = -1.0_rp
            t_vec_out(1) = 0.0_rp

            ! allocate space for result
            result_space = save_file%allocate_chunked_3D_space("uu", file_dimensions, ierr)
            ! save grid information
            call save_file%save_real_vector(conf%x, 'x', ierr)
            call save_file%save_real_vector(conf%y, 'y', ierr)


            ! implicit euler matrices
            Print *, "Allocating LU decompositions"
            allocate (implicit_euler_LU( size(conf%IC, 2) ), STAT=ierr)
            if (ierr /= 0) then
                Print *, "error allocating LU deocmpositions: ", ierr
            end if

            ! IMPLICIT EULER / EXPLICIT EULER SPECIFIC and autonomous diffusion specific.
            do dummy_i = 1, size(conf%IC, 2)
                Print *, Laplacian_COO%indx(size(Laplacian_COO%indx)), Laplacian_COO%jndx(size(Laplacian_COO%indx))
                Laplacian_COO_temp = copy_coo_matrix(Laplacian_COO, ierr)  ! potential memory leak here

                Print *, Laplacian_COO_temp%indx(size(Laplacian_COO_temp%indx)),&
                        Laplacian_COO_temp%jndx(size(Laplacian_COO_temp%indx))
                Laplacian_COO_temp%vals = Laplacian_COO_temp%vals * (dble(-1)*conf%dt*conf%diffusion_consts(dummy_i))
                
                eye_temp = speye(Laplacian_COO_temp%n, ierr)
                call sparse_add(Laplacian_COO_temp, eye_temp, ierr)
                implicit_euler_LU(dummy_i) = sparse_lu(Laplacian_COO_temp, ierr)
            end do

            ! plotting parameters

            prev_timestep = 0_ip
            plot_number = 0_ip
            plot_time = time()
            update_time = time()

            Print *, "Evolving in time"
            do while (t <= t_max)

                ! BEGIN SAVING
                if ((timestep * savenum) / max_timesteps > result_space%voffset(3)) then ! sketchy
                    call result_space%insert_page(u, ierr)
                    t_vec_out(result_space%voffset(3) + 1) = t
                end if
                ! END SAVING

                ! calculate k values
                call csr_multiply(Laplacian_CSR, u, dudxx)
                u_update1 = dudxx
                call scale_columns(u_update1, conf%diffusion_consts)

                do dummy_i = 1, size(conf%IC, 2)
                    u_update1(:, dummy_i) = sparse_lu_solve(implicit_euler_LU(dummy_i), u_update1(:, dummy_i), ierr)
                end do

                call conf%explicit_rhs(u, x_pass_through, t, u_update2, dudxx, conf%user_params)
                
                if (advection) then
                    call csr_multiply( Advection_CSR, u, u_adv1 )
                    call conf%advection_coefficient(u, x_pass_through, t, adv_const_field, dudxx, conf%user_params)
                    u_update2 = u_update2 + u_adv1 * adv_const_field
                end if
                
                ! euler step
                u = u + dt * (u_update1 + u_update2)

                ! impose boundary conditions
                ! could be improved by masking the flattened u, and reducing branching
                do dummy_i = 1, size(conf%IC, 2)
                    plot_uu = reshape(u(:, dummy_i), shape(conf%xx))

                    if (conf%DBCx_plus_mask(dummy_i)) then
                        plot_uu(:, size(conf%x, kind=ip)) = conf%DBCx_plus(dummy_i)
                    end if
                    if (conf%DBCx_minus_mask(dummy_i)) then
                        plot_uu(:, 1) = conf%DBCx_minus(dummy_i)
                    end if
                    if (conf%DBCy_plus_mask(dummy_i)) then
                        plot_uu(size(conf%y, kind=ip), :) = conf%DBCy_plus(dummy_i)
                    end if
                    if (conf%DBCy_minus_mask(dummy_i)) then
                        plot_uu(1, :) = conf%DBCy_minus(dummy_i)
                    end if
                    u(:, dummy_i) = reshape(plot_uu, (/ size(conf%xx, kind=ip) /))
                end do

                timestep = timestep + 1
                t = dt * timestep

                ! PLOTTING
                ! TODO: move this to a subroutine
                current_time = time()
                if (current_time - plot_time > conf%plot_interval) then
                    plot_number = plot_number + 1
                    write (timestep_string, "(I0)")  plot_number
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%reset()
                    call gp%title('u')
                    call gp%options('set terminal png size 4000,4000 font "Helvetica,45"')
                    call gp%options('set size square')
                    call gp%options('set output "' // trim(conf%plotfilename) // '"')
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    plot_time = time()
                end if

                if (current_time - update_time > 10.0_rp) then
                    Print *, t
                    Print *, "timesteps/second: ", real(timestep - prev_timestep, kind=rp) / (current_time - update_time)
                    prev_timestep = timestep
                    update_time = time()
                end if 

                ! END PLOTTING
            end do

            Print *, "saving"
            call save_file%save_real_vector(t_vec_out, 't', ierr)
            if (ierr /= 0) then
                Print *, "Error saving time values: ", ierr
            end if
            Print *, "closing"
            call save_file%close(ierr)
            if (ierr /= 0) then
                Print *, "Error closing file: ", ierr
            end if
            Print *, "closed"

            ! END SUBROUTINE
            return


        end subroutine IMEX_evolve_impliciteuler_euler_2D


        subroutine IMEX_evolve_CN_heun_2D (conf, ierr)
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config2d), intent(inout) :: conf
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_heun1, u_heun2, u_cn1, u_cn2,&
                plot_uu, u_adv1, u_adv2, adv_const_field, dudxx1, dudxx2
            real (rp), allocatable, dimension(:) :: t_vec_out
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time, update_time
            integer (ip) :: timestep, dummy_i, max_timesteps, savenum, prev_timestep, plot_number
            type (csr_matrix) :: Laplacian_CSR, Advection_CSR
            type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp, Advection_COO
            type (sparse_lu_decomposition), dimension(:), allocatable :: implicit_euler_LU
            character (len=1024) :: timestep_string
            logical :: advection

            ! HDF5 variables
            type (h5file) :: save_file
            type (chunked_3D_space) :: result_space
            integer (ip), dimension(3) :: file_dimensions
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise grid
            call conf%make_mesh(ierr)
            if (ierr /= 0) then
                Print *, "Error making mesh. Error code: ", ierr
                call exit()
            end if

            allocate (x_pass_through(size(conf%xx), 2), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating x pass through. Error code: ", ierr
                call exit()
            end if
            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            if (ierr /= 0) then
                Print *, "Error generating Laplacian COO.", ierr
            end if
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)
            Advection_COO = advection2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Advection_CSR = coo_to_csr(Advection_COO, ierr)

            ! initialise runge kutta variables
            Print *, "Allocating u and u update variables", size(conf%IC, 1), size(conf%IC, 2), size(conf%xx, 1), size(conf%xx, 2)
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_heun1( size(conf%IC, 1), size(conf%IC, 2) ), &
                    u_heun2( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (u_cn1(size(conf%IC, 1), size(conf%IC, 2) ), u_cn2(size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (u_adv1(size(conf%IC, 1), size(conf%IC, 2) ), u_adv2(size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (plot_uu( size(conf%xx, 1), size(conf%xx, 2) ), adv_const_field( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (dudxx1( size(conf%IC, 1), size(conf%IC, 2) ), dudxx2( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating dudxx pass through:" , ierr
                call exit()
            end if

            ! gather info from conf
            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC
            savenum = conf%savenum
            advection = (conf%advection)

            ! initialise saving variables
            max_timesteps = int( t_max / dt)

            ! estimate output size
            Print *, "Initialising save file"
            if ( 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1) > conf%max_save_size) then
                Print *, "ERROR: Desired output file exceeds maximum file size."
                Print *, "Desired bytes: ", 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1)
                Print *, "Maximum bytes: ", conf%max_save_size
                Print *, "increase conf%max_save_size or reduce conf%savenum to resolve"
                ierr = 2
                return
            end if
            Print *, "opening new file: ", conf%savefilename
            call save_file%new_file(conf%savefilename, ierr)
            if (ierr /=0_ip) then
                Print *, "Error in opening file, error code: ", ierr
                return
            end if

            file_dimensions = (/ size(conf%IC, 1, kind=HSIZE_T), size(conf%IC, 2, kind=HSIZE_T),&
                    INT(savenum+1, KIND=HSIZE_T)/)
            allocate (t_vec_out(INT(savenum+1, KIND=HSIZE_T)) )

            t_vec_out = -1.0_rp
            t_vec_out(1) = 0.0_rp

            ! allocate space for result
            result_space = save_file%allocate_chunked_3D_space("uu", file_dimensions, ierr)
            ! save grid information
            call save_file%save_real_vector(conf%x, 'x', ierr)
            call save_file%save_real_vector(conf%y, 'y', ierr)

            ! implicit CN matrices
            Print *, "Allocating LU decompositions"
            allocate (implicit_euler_LU( size(conf%IC, 2) ), STAT=ierr)

            do dummy_i = 1, size(conf%IC, 2)
                Laplacian_COO_temp = copy_coo_matrix(Laplacian_COO, ierr)  ! potential memory leak here
                Laplacian_COO_temp%vals = Laplacian_COO_temp%vals * (dble(-1.0_rp/2.0_rp)*conf%dt*conf%diffusion_consts(dummy_i))

                eye_temp = speye( Laplacian_COO_temp%n, ierr )
                call sparse_add(eye_temp, Laplacian_COO_temp, ierr)
                implicit_euler_LU(dummy_i) = sparse_lu(eye_temp, ierr)
            end do

            ! plotting parameters

            prev_timestep = 0_ip
            plot_number = 0_ip
            plot_time = time()
            update_time = time()
            Print *, "Evolving in time"
            do while (t <= t_max)
                call csr_multiply(Laplacian_CSR, u, dudxx1)
                u_cn1 = dudxx1
                call scale_columns(u_cn1, conf%diffusion_consts)

                call csr_multiply(Laplacian_CSR, u + u_cn1 * 0.5_rp*dt, dudxx2)
                u_cn2 = dudxx2
                call scale_columns(u_cn2, conf%diffusion_consts)

                do dummy_i = 1, size(conf%IC, 2)
                    u_cn2(:, dummy_i) = sparse_lu_solve(implicit_euler_LU(dummy_i), u_cn2(:, dummy_i), ierr)
                end do

                call conf%explicit_rhs(u, x_pass_through, t, u_heun1, dudxx1, conf%user_params)

                ! estimate derivative
                if (advection) then
                    call csr_multiply( Advection_CSR, u, u_adv1 )
                    call conf%advection_coefficient(u, x_pass_through, t, adv_const_field, dudxx1, conf%user_params)
                    u_heun1 = u_heun1 + u_adv1 * adv_const_field

                    u_adv1 = u + dt * u_heun1 ! used here to store the argument for the second heun iteration.
                    call conf%explicit_rhs(u_adv1, x_pass_through, t+conf%dt, u_heun2, dudxx2, conf%user_params)

                    call csr_multiply( Advection_CSR, u_adv1, u_adv2 )
                    call conf%advection_coefficient(u_adv1, x_pass_through, t + dt, adv_const_field, dudxx2, conf%user_params)
                    u_heun2 = u_heun2 + u_adv2 * adv_const_field
                else
                    
                    call conf%explicit_rhs(u + conf%dt * u_heun1, x_pass_through, t+conf%dt, u_heun2, dudxx2, conf%user_params)
                end if

                ! BEGIN SAVING
                if ((timestep * savenum) / max_timesteps > result_space%voffset(3)) then ! sketchy
                    call result_space%insert_page(u, ierr)
                    t_vec_out(result_space%voffset(3) + 1) = t
                end if
                ! END SAVING

                ! euler step
                u = u + dt/2.0_rp * (u_heun1 + u_heun2 + u_cn1 + u_cn2)

                ! impose boundary conditions
                ! could be improved by masking the flattened u, and reducing branching
                do dummy_i = 1, size(conf%IC, 2)
                    plot_uu = reshape(u(:, dummy_i), shape(conf%xx))

                    if (conf%DBCx_plus_mask(dummy_i)) then
                        plot_uu(:, size(conf%x, kind=ip)) = conf%DBCx_plus(dummy_i)
                    end if
                    if (conf%DBCx_minus_mask(dummy_i)) then
                        plot_uu(:, 1) = conf%DBCx_minus(dummy_i)
                    end if
                    if (conf%DBCy_plus_mask(dummy_i)) then
                        plot_uu(size(conf%y, kind=ip), :) = conf%DBCy_plus(dummy_i)
                    end if
                    if (conf%DBCy_minus_mask(dummy_i)) then
                        plot_uu(1, :) = conf%DBCy_minus(dummy_i)
                    end if
                    u(:, dummy_i) = reshape(plot_uu, (/ size(conf%xx, kind=ip) /))
                end do

                timestep = timestep + 1
                t = dt * timestep

                ! PLOTTING
                ! TODO: move this to a subroutine
                current_time = time()
                if (current_time - plot_time > conf%plot_interval) then
                    plot_number = plot_number + 1
                    write (timestep_string, "(I0)")  plot_number
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%reset()
                    call gp%title('u')
                    call gp%options('set terminal png size 4000,4000 font "Helvetica,45"')
                    call gp%options('set size square')
                    call gp%options('set output "' // trim(conf%plotfilename) // '"')
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    plot_time = time()
                end if

                if (current_time - update_time > 10.0_rp) then
                    Print *, t
                    Print *, "timesteps/second: ", real(timestep - prev_timestep, kind=rp) / (current_time - update_time)
                    prev_timestep = timestep
                    update_time = time()
                end if 
            end do

            call save_file%save_real_vector(t_vec_out, 't', ierr)
            call save_file%close(ierr)

            ! END SUBROUTINE
            return


        end subroutine IMEX_evolve_CN_heun_2D

        subroutine IMEX_evolve_impliciteuler_euler_2D_NAD (conf, ierr)
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config2d), intent(inout) :: conf
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_update1, u_update2, plot_uu,&
                    u_update_temp1, u_update_temp2, diff_const_field, adv_const_field, u_adv1, dudxx
            real (rp), allocatable, dimension(:) :: t_vec_out
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time, update_time
            integer (ip) :: timestep, dummy_i, max_timesteps, savenum, prev_timestep, plot_number
            type (csr_matrix) :: Laplacian_CSR, Laplacian_CSR_temp, Advection_CSR
            type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp, Advection_COO
            type (csr_matrix), dimension(:), allocatable :: implicit_system
            character (len=1024) :: timestep_string
            logical :: advection

            ! HDF5 variables
            type (h5file) :: save_file
            type (chunked_3D_space) :: result_space
            integer (ip), dimension(3) :: file_dimensions
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise grid
            Print *, "Initialising grid"
            call conf%make_mesh(ierr)
            allocate (x_pass_through(size(conf%xx), 2))
            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)
            Advection_COO = advection2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)

            ! initialise runge kutta variables
            Print *, "Allocating u and u update variables", size(conf%IC, 1), size(conf%IC, 2), size(conf%xx, 1), size(conf%xx, 2)
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_update1( size(conf%IC, 1), size(conf%IC, 2) ),&
                u_adv1( size(conf%IC, 1), size(conf%IC, 2) ), adv_const_field( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating u and u update variables:" , ierr
                call exit()
            end if
            allocate (u_update2( size(conf%IC, 1), size(conf%IC, 2) ),&
                    u_update_temp1( size(conf%IC, 1), 1 ), u_update_temp2( size(conf%IC, 1), 1 ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating u and u update variables:" , ierr
                call exit()
            end if
            allocate (plot_uu( size(conf%xx, 1), size(conf%xx, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating u and u update variables:" , ierr
                call exit()
            end if
            allocate (diff_const_field( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating diff const field:" , ierr
                call exit()
            end if
            allocate (dudxx( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating dudxx pass through:" , ierr
                call exit()
            end if

            ! gather info from conf
            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC
            savenum = conf%savenum
            advection = conf%advection

            ! initialise saving variables
            max_timesteps = int( t_max / dt)

            ! estimate output size
            Print *, "Initialising save file"
            if ( 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1) > conf%max_save_size) then
                Print *, "ERROR: Desired output file exceeds maximum file size."
                Print *, "Desired bytes: ", 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1)
                Print *, "Maximum bytes: ", conf%max_save_size
                Print *, "increase conf%max_save_size or reduce conf%savenum to resolve"
                ierr = 2
                return
            end if
            call save_file%new_file(conf%savefilename, ierr)

            file_dimensions = (/ size(conf%IC, 1, kind=HSIZE_T), size(conf%IC, 2, kind=HSIZE_T),&
                    INT(savenum+1, KIND=HSIZE_T)/)
            allocate (t_vec_out(INT(savenum+1, KIND=HSIZE_T)) )

            t_vec_out = -1.0_rp
            t_vec_out(1) = 0.0_rp

            ! allocate space for result
            result_space = save_file%allocate_chunked_3D_space("uu", file_dimensions, ierr)
            ! save grid information
            call save_file%save_real_vector(conf%x, 'x', ierr)
            call save_file%save_real_vector(conf%y, 'y', ierr)


            ! implicit euler matrices
            Print *, "Allocating linear systems"
            !allocate (implicit_euler_LU( size(conf%IC, 2) ), STAT=ierr)


            ! plotting parameters

            prev_timestep = 0_ip
            plot_number = 0_ip
            plot_time = time()
            update_time = time()
            Print *, "Evolving in time"
            do while (t <= t_max)

                ! BEGIN SAVING
                if ((timestep * savenum) / max_timesteps > result_space%voffset(3)) then ! sketchy
                    call result_space%insert_page(u, ierr)
                    t_vec_out(result_space%voffset(3) + 1) = t
                end if
                ! END SAVING

                call csr_multiply(Laplacian_CSR, u, dudxx)
                u_update1 = dudxx

                ! calculate k values
                call conf%diffusivity(u, x_pass_through, t, diff_const_field, dudxx, conf%user_params)

                u_update1 = u_update1 * diff_const_field

                do dummy_i = 1, size(conf%IC, 2)
                    Laplacian_COO_temp = copy_coo_matrix(Laplacian_COO, ierr)  ! potential memory leak here

                    call coo_scale_rows(Laplacian_COO_temp, diff_const_field(:, dummy_i) * (-1.0_rp) * conf%dt)

                    eye_temp = speye(Laplacian_COO_temp%n, ierr)
                    call sparse_add(Laplacian_COO_temp, eye_temp, ierr)
                    Laplacian_CSR_temp = coo_to_csr(Laplacian_COO_temp, ierr)
                    u_update1(:, dummy_i) = sparse_direct_solve(Laplacian_CSR_temp, u_update1(:, dummy_i), ierr)

                    deallocate (Laplacian_COO_temp%vals, Laplacian_COO_temp%indx, Laplacian_COO_temp%jndx)
                    deallocate (Laplacian_CSR_temp%vals, Laplacian_CSR_temp%rows, Laplacian_CSR_temp%jndx)
                    deallocate (eye_temp%vals, eye_temp%indx, eye_temp%jndx)
                end do

                call conf%explicit_rhs(u, x_pass_through, t, u_update2, dudxx, conf%user_params)

                if (advection) then
                    call csr_multiply(Advection_CSR, u, u_adv1)
                    call conf%advection_coefficient(u, x_pass_through, t, adv_const_field, dudxx, conf%user_params)
                    u_update2 = u_update2 + u_adv1 * adv_const_field
                end if

                ! euler step
                u = u + dt * (u_update1 + u_update2)

                                ! impose boundary conditions
                ! could be improved by masking the flattened u, and reducing branching
                do dummy_i = 1, size(conf%IC, 2)
                    plot_uu = reshape(u(:, dummy_i), shape(conf%xx))

                    if (conf%DBCx_plus_mask(dummy_i)) then
                        plot_uu(:, size(conf%x, kind=ip)) = conf%DBCx_plus(dummy_i)
                    end if
                    if (conf%DBCx_minus_mask(dummy_i)) then
                        plot_uu(:, 1) = conf%DBCx_minus(dummy_i)
                    end if
                    if (conf%DBCy_plus_mask(dummy_i)) then
                        plot_uu(size(conf%y, kind=ip), :) = conf%DBCy_plus(dummy_i)
                    end if
                    if (conf%DBCy_minus_mask(dummy_i)) then
                        plot_uu(1, :) = conf%DBCy_minus(dummy_i)
                    end if
                    u(:, dummy_i) = reshape(plot_uu, (/ size(conf%xx, kind=ip) /))
                end do

                timestep = timestep + 1
                t = dt * timestep

                ! PLOTTING
                ! TODO: move this to a subroutine
                current_time = time()
                if (current_time - plot_time > conf%plot_interval) then
                    plot_number = plot_number + 1
                    write (timestep_string, "(I0)")  plot_number
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%reset()
                    call gp%title('u')
                    call gp%options('set terminal png size 4000,4000 font "Helvetica,45"')
                    call gp%options('set size square')
                    call gp%options('set output "' // trim(conf%plotfilename) // '"')
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    plot_time = time()
                end if

                if (current_time - update_time > 10.0_rp) then
                    Print *, t
                    Print *, "timesteps/second: ", real(timestep - prev_timestep, kind=rp) / (current_time - update_time)
                    prev_timestep = timestep
                    update_time = time()
                end if 
                ! END PLOTTING
            end do

            call save_file%save_real_vector(t_vec_out, 't', ierr)
            call save_file%close(ierr)

            ! END SUBROUTINE
            return


        end subroutine IMEX_evolve_impliciteuler_euler_2D_NAD

        subroutine IMEX_evolve_CN_heun_2D_NAD (conf, ierr)
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config2d), intent(inout) :: conf
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_heun1, u_heun2, u_cn1, u_cn2, u_adv1, u_adv2, plot_uu,&
                    diff_const_field, adv_const_field, dudxx1, dudxx2
            real (rp), allocatable, dimension(:) :: t_vec_out
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time, update_time
            integer (ip) :: timestep, dummy_i, max_timesteps, savenum, prev_timestep, plot_number
            type (csr_matrix) :: Laplacian_CSR, Laplacian_CSR_temp, Advection_CSR
            type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp, Advection_COO
            type (csr_matrix), dimension(:), allocatable :: implicit_system
            character (len=1024) :: timestep_string
            logical :: advection

            ! HDF5 variables
            type (h5file) :: save_file
            type (chunked_3D_space) :: result_space
            integer (ip), dimension(3) :: file_dimensions
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise grid
            Print *, "Initialising grid"
            call conf%make_mesh(ierr)
            allocate (x_pass_through(size(conf%xx), 2))
            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)
            Advection_COO = advection2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Advection_CSR = coo_to_csr(Advection_COO, ierr)

            ! initialise runge kutta variables
            Print *, "Allocating u and u update variables", size(conf%IC, 1), size(conf%IC, 2), size(conf%xx, 1), size(conf%xx, 2)
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_heun1( size(conf%IC, 1), size(conf%IC, 2) ), &
                    u_heun2( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate ( u_cn1(size(conf%IC, 1), size(conf%IC, 2)), u_cn2(size(conf%IC, 1), size(conf%IC, 2)) )
            allocate ( u_adv1(size(conf%IC, 1), size(conf%IC, 2)), u_adv2(size(conf%IC, 1), size(conf%IC, 2)) )
            allocate (plot_uu( size(conf%xx, 1), size(conf%xx, 2) ))

            allocate (diff_const_field( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (adv_const_field( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (dudxx1( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating dudxx pass through:" , ierr
                call exit()
            end if
            allocate (dudxx2( size(conf%IC, 1), size(conf%IC, 2) ), stat=ierr)
            if (ierr /= 0) then
                Print *, "Error allocating dudxx pass through:" , ierr
                call exit()
            end if


            ! gather info from conf
            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC
            savenum = conf%savenum
            advection = (conf%advection)

            ! initialise saving variables
            max_timesteps = int( t_max / dt)

            ! estimate output size
            Print *, "Initialising save file"
            if ( 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1) > conf%max_save_size) then
                Print *, "ERROR: Desired output file exceeds maximum file size."
                Print *, "Desired bytes: ", 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1)
                Print *, "Maximum bytes: ", conf%max_save_size
                Print *, "increase conf%max_save_size or reduce conf%savenum to resolve"
                ierr = 2
                return
            end if
            call save_file%new_file(conf%savefilename, ierr)

            file_dimensions = (/ size(conf%IC, 1, kind=HSIZE_T), size(conf%IC, 2, kind=HSIZE_T),&
                    INT(savenum+1, KIND=HSIZE_T)/)
            allocate (t_vec_out(INT(savenum+1, KIND=HSIZE_T)) )

            t_vec_out = -1.0_rp
            t_vec_out(1) = 0.0_rp

            ! allocate space for result
            result_space = save_file%allocate_chunked_3D_space("uu", file_dimensions, ierr)
            ! save grid information
            call save_file%save_real_vector(conf%x, 'x', ierr)
            call save_file%save_real_vector(conf%y, 'y', ierr)


            ! implicit euler matrices


            ! plotting parameters

            prev_timestep = 0_ip
            plot_number = 0_ip
            plot_time = time()
            update_time = time()
            Print *, "Evolving in time IMEX_evolve_CN_heun_2D_NAD"
            ! calculate initial diffusivity values
            call csr_multiply(Laplacian_CSR, u, dudxx1)
            call conf%diffusivity(u, x_pass_through, t, diff_const_field, dudxx1, conf%user_params)
            do while (t <= t_max)

                ! BEGIN SAVING
                if ((timestep * savenum) / max_timesteps > result_space%voffset(3)) then ! sketchy
                    call result_space%insert_page(u, ierr)
                    t_vec_out(result_space%voffset(3) + 1) = t
                end if
                ! END SAVING
                

                call csr_multiply(Laplacian_CSR, u, dudxx1)
                u_cn1 = dudxx1
                u_cn1 = u_cn1 * diff_const_field
                call csr_multiply(Laplacian_CSR, u + u_cn1 * 0.5_rp*dt, dudxx2)
                u_cn2 = dudxx2

                do dummy_i = 1, size(conf%IC, 2)
                    ! update diffusivity
                    call conf%diffusivity(u, x_pass_through, (timestep+1_ip)*conf%dt, diff_const_field, dudxx2, &
                            conf%user_params)
                    Laplacian_COO_temp = copy_coo_matrix(Laplacian_COO, ierr)  ! potential memory leak here
                    call coo_scale_rows(Laplacian_COO_temp, diff_const_field(:, dummy_i) * (-1.0_rp/2.0_rp) * conf%dt)

                    eye_temp = speye( Laplacian_COO_temp%n, ierr )
                    call sparse_add(Laplacian_COO_temp, eye_temp, ierr)
                    Laplacian_CSR_temp = coo_to_csr(Laplacian_COO_temp, ierr)

                    u_cn2(:, dummy_i) = diff_const_field(:, dummy_i) * u_cn2(:, dummy_i)
                    u_cn2(:, dummy_i) = sparse_direct_solve(Laplacian_CSR_temp, u_cn2(:, dummy_i), ierr)

                    deallocate (Laplacian_COO_temp%vals, Laplacian_COO_temp%indx, Laplacian_COO_temp%jndx)
                    deallocate (Laplacian_CSR_temp%vals, Laplacian_CSR_temp%rows, Laplacian_CSR_temp%jndx)
                    deallocate(eye_temp%vals, eye_temp%indx, eye_temp%jndx)
                end do

                call conf%explicit_rhs(u, x_pass_through, t, u_heun1, dudxx1, conf%user_params)

                ! estimate derivative
                if (advection) then
                    call csr_multiply( Advection_CSR, u, u_adv1 )
                    call conf%advection_coefficient(u, x_pass_through, t, adv_const_field, dudxx1, conf%user_params)
                    u_heun1 = u_heun1 + u_adv1 * adv_const_field

                    u_adv1 = u + dt * u_heun1 ! used here to store the argument for the second heun iteration.
                    call conf%explicit_rhs(u_adv1, x_pass_through, t+conf%dt, u_heun2, dudxx2, conf%user_params)

                    call csr_multiply( Advection_CSR, u_adv1, u_adv2 )
                    call conf%advection_coefficient(u_adv1, x_pass_through, t + dt, adv_const_field, dudxx2, conf%user_params)
                    u_heun2 = u_heun2 + u_adv2 * adv_const_field
                else
                    
                    call conf%explicit_rhs(u + conf%dt * u_heun1, x_pass_through, t+conf%dt, u_heun2, dudxx2, conf%user_params)
                end if

                ! runge kutta step
                u = u + dt/2.0_rp * (u_heun1 + u_heun2 + u_cn1 + u_cn2)
                

                ! impose boundary conditions
                ! could be improved by masking the flattened u, and reducing branching
                do dummy_i = 1, size(conf%IC, 2)
                    plot_uu = reshape(u(:, dummy_i), shape(conf%xx))

                    if (conf%DBCx_plus_mask(dummy_i)) then
                        plot_uu(:, size(conf%x, kind=ip)) = conf%DBCx_plus(dummy_i)
                    end if
                    if (conf%DBCx_minus_mask(dummy_i)) then
                        plot_uu(:, 1) = conf%DBCx_minus(dummy_i)
                    end if
                    if (conf%DBCy_plus_mask(dummy_i)) then
                        plot_uu(size(conf%y, kind=ip), :) = conf%DBCy_plus(dummy_i)
                    end if
                    if (conf%DBCy_minus_mask(dummy_i)) then
                        plot_uu(1, :) = conf%DBCy_minus(dummy_i)
                    end if
                    u(:, dummy_i) = reshape(plot_uu, (/ size(conf%xx, kind=ip) /))
                end do

                timestep = timestep + 1
                t = dt * timestep

                ! PLOTTING
                ! TODO: move this to a subroutine
                current_time = time()
                if (current_time - plot_time > conf%plot_interval) then
                    plot_number = plot_number + 1
                    write (timestep_string, "(I0)")  plot_number
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%reset()
                    call gp%title('u')
                    call gp%options('set terminal png size 4000,4000 font "Helvetica,45"')
                    call gp%options('set size square')
                    call gp%options('set output "' // trim(conf%plotfilename) // '"')
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    plot_time = time()
                end if

                if (current_time - update_time > 10.0_rp) then
                    Print *, t
                    Print *, "timesteps/second: ", real(timestep - prev_timestep, kind=rp) / (current_time - update_time)
                    prev_timestep = timestep
                    update_time = time()
                end if 
                ! END PLOTTING
            end do

            call save_file%save_real_vector(t_vec_out, 't', ierr)
            call save_file%close(ierr)

            ! END SUBROUTINE
            return


        end subroutine IMEX_evolve_CN_heun_2D_NAD
end module evolvePDE


