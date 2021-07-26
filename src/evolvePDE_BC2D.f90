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