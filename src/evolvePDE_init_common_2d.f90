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