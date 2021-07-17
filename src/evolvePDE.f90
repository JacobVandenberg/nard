
module evolvePDE
    use config_m
    use ogpf
    use diff
    use helpers
    use sparse_matrices
    use precision
    use hdf5
    implicit none

    contains
        subroutine IMEX_evolve_impliciteuler_euler_2D (conf, ierr)
            ! UNDER DEVELOPMENT
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config2d), intent(inout) :: conf
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_update1, u_update2, plot_uu,&
                    u_update_temp1, u_update_temp2
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time
            integer (ip) :: timestep, dummy_i, max_timesteps, savenum
            type (csr_matrix) :: Laplacian_CSR
            type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp
            type (sparse_lu_decomposition), dimension(:), allocatable :: implicit_euler_LU

            ! HDF5 variables
            integer :: hdf5err
            integer (HID_T) :: outfile_id
            integer (HID_T) :: dsetv_id
            integer (HID_T) :: dataspacev_id
            integer (HID_T) :: memspace
            integer (HSIZE_T), dimension(3) :: file_dimensions, chunk_dimensions, voffset, vchunkcount
            integer :: rank = 3
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise runge kutta variables
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_update1( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (u_update2( size(conf%IC, 1), size(conf%IC, 2) ),&
                    u_update_temp1( size(conf%IC, 1), 1 ), u_update_temp2( size(conf%IC, 1), 1 ))
            allocate (plot_uu( size(conf%xx, 1), size(conf%xx, 2) ))

            ! initialise grid
            call conf%make_mesh(ierr)
            allocate (x_pass_through(size(conf%xx), 2))
            x_pass_through(:, 1) = reshape(conf%xx, (/size(conf%xx)/))
            x_pass_through(:, 2) = reshape(conf%yy, (/size(conf%yy)/))
            Laplacian_COO = laplacian2d(conf%x, conf%y, conf%BCx, conf%BCy, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)

            ! gather info from conf
            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC
            savenum = conf%savenum

            ! initialise saving variables
            max_timesteps = int( t_max / dt)

            ! estimate output size
            if ( (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1) > conf%max_save_size) then
                Print *, "ERROR: Desired output file exceeds maximum file size."
                Print *, "Desired bytes: ", (size(conf%IC, kind=ip) + 1_ip) * savenum
                Print *, "Maximum bytes: ", conf%max_save_size
                Print *, "increase conf%max_save_size or reduce conf%savenum to resolve"
                ierr = 2
                return
            end if

            file_dimensions = (/ size(conf%IC, 1, kind=HSIZE_T), size(conf%IC, 2, kind=HSIZE_T),&
                    INT(savenum+1, KIND=HSIZE_T)/)
            chunk_dimensions = file_dimensions
            chunk_dimensions(3) = 1

            call h5open_f(hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error code ", hdf5err
            end if
            ! Create a new file using the default properties.
            call h5fcreate_f(conf%savefilename, H5F_ACC_TRUNC_F, outfile_id, hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error code ", hdf5err
            end if

            ! Create the data space for the binned dataset.
            call h5screate_simple_f(rank, file_dimensions, dataspacev_id,  hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error code ", hdf5err
            end if

            ! Create the chunked dataset.
            call h5dcreate_f(outfile_id, "data", H5T_NATIVE_DOUBLE, dataspacev_id, dsetv_id, hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error code ", hdf5err
            end if

            ! Create the memory space for the selection
            call h5screate_simple_f(rank, chunk_dimensions, memspace,  hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error code ", hdf5err
            end if

            voffset(1:2) = 0
            vchunkcount = chunk_dimensions

            voffset(3) = -1

            ! implicit euler matrices
            allocate (implicit_euler_LU( size(conf%IC, 2) ), STAT=ierr)

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

            call cpu_time(plot_time)
            do while (t < t_max)

                call csr_multiply(Laplacian_CSR, u, u_update1)
                call scale_columns(u_update1, conf%diffusion_consts)
                do dummy_i = 1, size(conf%IC, 2)
                    u_update_temp1(:, 1) = u_update1(:, dummy_i)
                    call sparse_lu_solve(implicit_euler_LU(dummy_i), u_update_temp1, u_update_temp2, ierr)
                    u_update1(:, dummy_i) = u_update_temp2(:, 1)
                end do

                call conf%explicit_rhs(u, x_pass_through, t, u_update2)

                ! BEGIN SAVING
                if ((timestep * savenum) / max_timesteps > voffset(3)) then
                    voffset(3) = ((timestep * savenum) / max_timesteps)

                    call h5sselect_hyperslab_f(dataspacev_id, H5S_SELECT_SET_F, voffset, vchunkcount, hdf5err)
                    if (hdf5err /=0) then
                        Print *, "HDF5 ERROR: error in h5sselect_hyperslab_f. error code ", hdf5err
                        return
                    end if

                    ! Write the data to the dataset.
                    call h5dwrite_f(dsetv_id, H5T_NATIVE_DOUBLE, reshape(u, chunk_dimensions), chunk_dimensions,&
                            hdf5err, memspace, dataspacev_id)
                    if (hdf5err /=0) then
                        Print *, "HDF5 ERROR: error in h5dwrite_f. error code ", hdf5err
                        return
                    end if
                    call h5sclose_f(dataspacev_id, hdf5err)
                    if (hdf5err /=0) then
                        Print *, "HDF5 ERROR: error in h5sclose_f. error code ", hdf5err
                        return
                    end if
                    call h5dget_space_f(dsetv_id, dataspacev_id, hdf5err)
                    if (hdf5err /=0) then
                        Print *, "HDF5 ERROR: error in h5dget_space_f. error code ", hdf5err
                        return
                    end if
                end if
                ! END SAVING

                ! BEGIN PLOTTING
                ! TODO: move this to a subroutine
                call cpu_time(current_time)
                if (current_time - plot_time > conf%plot_interval) then
                    Print *, t
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%title('u')
                    call gp%options('set terminal png size 1920,1080; set output "src/tests/test_evolvePDE_euler.png"')
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    call cpu_time(plot_time)
                end if
                ! END PLOTTING

                ! euler step
                u = u + dt * (u_update1 + u_update2)

                timestep = timestep + 1
                t = dt * timestep
            end do

            call h5dclose_f(dsetv_id, hdf5err)
            call h5sclose_f(memspace, hdf5err)
            call h5fclose_f(outfile_id, hdf5err)

            call h5close_f(hdf5err)

            ! END SUBROUTINE
            return


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
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_update1, u_update2, plot_uu
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time
            integer (ip) :: timestep
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


        subroutine IMEX_evolve_impliciteuler_euler_3D (conf, ierr)
            ! UNDER DEVELOPMENT
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            type (config3d), intent(inout) :: conf
            integer (ip) :: ierr

            ! runtime variables
            real (rp), allocatable, dimension(:, :, :) :: uuu
            real (rp), allocatable, dimension(:, :) :: u, x_pass_through, u_update1, u_update2, plot_uu,&
                    u_update_temp1, u_update_temp2, xx, yy
            type (gpf) :: gp
            real (rp) :: t_max, t, dt, plot_time, current_time
            integer (ip) :: timestep, dummy_i
            type (csr_matrix) :: Laplacian_CSR
            type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp
            type (sparse_lu_decomposition), dimension(:), allocatable :: implicit_euler_LU
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            ! initialise grid
            call conf%make_mesh(ierr)
            allocate (x_pass_through(size(conf%xxx), 3))
            x_pass_through(:, 1) = reshape(conf%xxx, (/size(conf%xxx)/))
            x_pass_through(:, 2) = reshape(conf%yyy, (/size(conf%yyy)/))
            x_pass_through(:, 3) = reshape(conf%zzz, (/size(conf%zzz)/))
            Laplacian_COO = laplacian3d(conf%x, conf%y, conf%z, conf%BCx, conf%BCy, conf%BCz, ierr)
            Laplacian_CSR = coo_to_csr(Laplacian_COO, ierr)

            allocate (xx(size(conf%y), size(conf%x)),&
                    yy(size(conf%y), size(conf%x)), plot_uu(size(conf%y), size(conf%x)),&
                    uuu(size(conf%y), size(conf%x), size(conf%z)), STAT=ierr)
            call meshgrid(xx, yy, conf%x, conf%y, ierr)

            ! initialise runge kutta variables
            allocate (u( size(conf%IC, 1), size(conf%IC, 2) ), u_update1( size(conf%IC, 1), size(conf%IC, 2) ))
            allocate (u_update2( size(conf%IC, 1), size(conf%IC, 2) ),&
                    u_update_temp1( size(conf%IC, 1), 1 ), u_update_temp2( size(conf%IC, 1), 1 ))

            ! implicit euler matrices
            allocate (implicit_euler_LU( size(conf%IC, 2) ), STAT=ierr)

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



            timestep = 0
            t = dble(0)
            dt = conf%dt
            t_max = conf%t_max
            u = conf%IC

            ! plotting parameters

            call cpu_time(plot_time)
            do while (t < t_max)

                call csr_multiply(Laplacian_CSR, u, u_update1)
                call scale_columns(u_update1, conf%diffusion_consts)
                do dummy_i = 1, size(conf%IC, 2)
                    u_update_temp1(:, 1) = u_update1(:, dummy_i)
                    call sparse_lu_solve(implicit_euler_LU(dummy_i), u_update_temp1, u_update_temp2, ierr)
                    u_update1(:, dummy_i) = u_update_temp2(:, 1)
                end do

                call conf%explicit_rhs(u, x_pass_through, t, u_update2)

                ! BEGIN PLOTTING
                ! TODO: move this to a subroutine
                call cpu_time(current_time)
                if (current_time - plot_time > conf%plot_interval) then
                    Print *, t
                    uuu = reshape(u, shape(conf%xxx))
                    plot_uu(:, :) = uuu(:, :, size(conf%z)/2)
                    call gp%title('u')
                    call gp%options("set terminal png size 1920,1080; set output 'tests/test_evolvePDE_euler3d.png'")
                    call gp%contour(xx,yy,plot_uu, palette='jet')
                    call cpu_time(plot_time)
                end if
                ! END PLOTTING

                ! RK step
                u = u + dt * (u_update1 + u_update2)

                timestep = timestep + 1
                t = dt * timestep
            end do

            ! END SUBROUTINE
            return


        end subroutine IMEX_evolve_impliciteuler_euler_3D



end module evolvePDE


