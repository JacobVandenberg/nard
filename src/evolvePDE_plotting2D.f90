! BEGIN PLOTTING
                ! TODO: move this to a subroutine
                current_time = time()
                if (current_time - plot_time > conf%plot_interval) then
                    Print *, t
                    Print *, "timesteps/second: ", real(timestep - prev_timestep, kind=rp) / (current_time - plot_time)
                    plot_number = plot_number + 1
                    write (timestep_string, "(I0)")  plot_number
                    plot_uu = reshape(u(:, 1), shape(conf%xx))
                    call gp%reset()
                    call gp%title('u')
                    call gp%options('set terminal png size 4000,4000 font "Helvetica,45"')
                    call gp%options('set size square')
                    call gp%options('set output "' // trim(conf%plotfilename) // '"')
                    call gp%contour(conf%xx,conf%yy,plot_uu, palette='jet')
                    prev_timestep = timestep
                    plot_time = time()
                end if