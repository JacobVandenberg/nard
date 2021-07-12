
module evolvePDE
    use config_m
    use ogpf
    use diff
    implicit none

    contains
        subroutine evolve (conf)
            implicit none
            type (config) :: conf
            real (kind = 8), dimension(:, :, :, :) :: uuu
            real (kind = 8), dimension(:) :: x, y, z


        end subroutine evolve

        subroutine euler (conf)
            !
            ! evolves configuration using eulers method
            ! inputs: conf: configuration structure
            !
            ! outputs: write to file specified by conf
            !

            implicit none
            type (config) :: conf
            real (kind = 8), dimension(:, :), allocatable :: u
            integer, dimension(:), allocatable :: DDx_indx, DDx_jndx, DDy_indx, DDy_jndx, L_indx, L_jndx, eye_indx, eye_jndx
            real, dimension(:), allocatable :: DDy_vals, DDx_vals, L_vals, x, y, eye_vals
            integer :: ierr, N

            N = 30

            ! allocate the result to the size of the initial conditions
            allocate ( u(shape(conf%IC)), STAT=ierr)
            if (ierr /= 0)
                Print *, 'allocation error for initial solution array'
                exit
            end if

            ! generate grid lines
            call my_linspace(dble(-1), dble(1), N, x)
            call my_linspace(dble(-1), dble(1), N, y)

            ! generate laplacian matrix
            ! diff wrt x
            call double_diff_FD(x, DDx_vals, DDx_indx, DDx_jndx, ierr)
            ! periodic BCS
            call sparse_update(DDx_vals, DDx_indx, DDx_jndx, sparse_get(DDx_vals, DDx_indx, DDx_jndx, 1, 2), 1, N, ierr)
            call sparse_update(DDx_vals, DDx_indx, DDx_jndx, sparse_get(DDx_vals, DDx_indx, DDx_jndx, 1, 2), N, 1, ierr)

            call double_diff_FD(y, DDy_vals, DDy_indx, DDy_jndx, ierr)
            ! periodic BCS
            call sparse_update(DDy_vals, DDy_indx, DDy_jndx, sparse_get(DDy_vals, DDy_indx, DDy_jndx, 1, 2), 1, N, ierr)
            call sparse_update(DDy_vals, DDy_indx, DDy_jndx, sparse_get(DDy_vals, DDy_indx, DDy_jndx, 1, 2), N, 1, ierr)

            ! identity matrix
            call speye(N, eye_vals, eye_indx, eye_jndx, ierr)

            ! kron tensor product
            !
            !DDx = kron(DDx, eye(length(mesh.y)));
            !DDy = kron(eye(length(mesh.x)), DDy);
            call sparse_kron(DDx_vals, DDx_indx, DDx_jndx, N, N, eye_vals, eye_indx, eye_jndx, N, N, L_vals, L_indx, L_jndx, nout, mout, ierr)
            call move_alloc (L_vals, DDx_vals)
            call move_alloc (L_indx, DDx_indx)
            call move_alloc (L_jndx, DDx_jndx)

            call sparse_kron(eye_vals, eye_indx, eye_jndx, N, N, DDy_vals, DDy_indx, DDy_jndx, N, N, L_vals, L_indx, L_jndx, nout, mout, ierr)
            call move_alloc (L_vals, DDy_vals)
            call move_alloc (L_indx, DDy_indx)
            call move_alloc (L_jndx, DDy_jndx)

            ! L = DDx + DDy
            call sparse_add(DDx_vals, DDx_indx, DDx_jndx, DDy_vals, DDy_indx, DDy_jndx, L_vals, L_indx, L_jndx, ierr)


            deallocate(u)


        end subroutine euler

end module evolvePDE


