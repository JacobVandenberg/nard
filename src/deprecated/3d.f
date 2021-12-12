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
    real (rp), allocatable, dimension(:) :: t_vec_out
    type (gpf) :: gp
    real (rp) :: t_max, t, dt, plot_time, current_time
    integer (ip) :: timestep, dummy_i, savenum, max_timesteps
    type (csr_matrix) :: Laplacian_CSR
    type (coo_matrix) :: Laplacian_COO, Laplacian_COO_temp, eye_temp
    type (sparse_lu_decomposition), dimension(:), allocatable :: implicit_euler_LU

    ! HDF5 variables
    type (h5file) :: save_file
    type (chunked_3D_space) :: result_space
    integer (ip), dimension(3) :: file_dimensions
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

    ! gather info from conf
    timestep = 0
    t = dble(0)
    dt = conf%dt
    t_max = conf%t_max
    u = conf%IC
    savenum = conf%savenum

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

    ! initialise saving variables
    max_timesteps = int( t_max / dt)

    ! estimate output size
    if ( 8_ip * (size(conf%IC, kind=ip) + 1_ip) * (savenum + 1) > conf%max_save_size) then
        Print *, "ERROR: Desired output file exceeds maximum file size."
        Print *, "Desired bytes: ", (size(conf%IC, kind=ip) + 1_ip) * savenum
        Print *, "Maximum bytes: ", conf%max_save_size
        Print *, "increase conf%max_save_size or reduce conf%savenum to resolve"
        ierr = 2
        return
    end if
    call save_file%new_file(conf%savefilename, ierr)
    if (ierr /= 0) then
        Print *, "Error in creating file", ierr
        return
    end if

    file_dimensions = (/ size(conf%IC, 1, kind=HSIZE_T), size(conf%IC, 2, kind=HSIZE_T),&
            INT(savenum+1, KIND=HSIZE_T)/)
    allocate (t_vec_out(INT(savenum+1, KIND=HSIZE_T)), STAT=ierr)
    if (ierr /= 0) then
        Print *, "Error in allocating t_vec_out", ierr
        return
    end if

    t_vec_out = -1.0_rp
    t_vec_out(1) = 0.0_rp

    ! allocate space for result
    result_space = save_file%allocate_chunked_3D_space("uu", file_dimensions, ierr)
    if (ierr /= 0) then
        Print *, "Error in allocating chunked 3d space", ierr
        return
    end if
    ! save grid information
    call save_file%save_real_vector(conf%x, 'x', ierr)
    if (ierr /= 0) then
        Print *, "Error in saving x", ierr
        return
    end if
    call save_file%save_real_vector(conf%y, 'y', ierr)
    if (ierr /= 0) then
        Print *, "Error in saving y", ierr
        return
    end if
    call save_file%save_real_vector(conf%z, 'z', ierr)
    if (ierr /= 0) then
        Print *, "Error in saving z", ierr
        return
    end if

    call cpu_time(plot_time)
    do while (t <= t_max)
        call csr_multiply(Laplacian_CSR, u, u_update1)
        call scale_columns(u_update1, conf%diffusion_consts)
        do dummy_i = 1, size(conf%IC, 2)
            u_update1(:, dummy_i) = sparse_lu_solve(implicit_euler_LU(dummy_i), u_update1(:, dummy_i), ierr)
        end do

        call conf%explicit_rhs(u, x_pass_through, t, u_update2, conf%user_params)

        ! BEGIN SAVING
        if ((timestep * savenum) / max_timesteps > result_space%voffset(3)) then ! sketchy
            call result_space%insert_page(u, ierr)
            if (ierr /= 0) then
                Print *, "Error in inserting page", ierr
                return
            end if
            t_vec_out(result_space%voffset(3) + 1) = t
        end if
        ! END SAVING

        ! BEGIN PLOTTING
        ! TODO: move this to a subroutine
        call cpu_time(current_time)
        if (current_time - plot_time > conf%plot_interval) then
            Print *, t
            uuu = reshape(u, shape(conf%xxx))
            plot_uu(:, :) = uuu(:, :, size(conf%z)/2)
            call gp%title('u')
            call gp%options("set terminal png size 1920,1080; set output 'src/tests/test_evolvePDE_euler3d.png'")
            call gp%contour(xx,yy,plot_uu, palette='jet')
            call cpu_time(plot_time)
        end if
        ! END PLOTTING

        ! RK step
        u = u + dt * (u_update1 + u_update2)

        timestep = timestep + 1
        t = dt * timestep
    end do

    call save_file%save_real_vector(t_vec_out, 't', ierr)
    call save_file%close(ierr)

    ! END SUBROUTINE
    return


end subroutine IMEX_evolve_impliciteuler_euler_3D