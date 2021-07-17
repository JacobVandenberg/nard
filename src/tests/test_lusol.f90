! Created by  on 14/7/21.

program test_lusol
    use lusol
    use ogpf
    use diff
    use sparse_matrices
    Print *, 'test_1'
    call test_1
    contains

        subroutine test_1
            implicit none

            type (coo_matrix) :: mat

            real (rp) , dimension(:), allocatable :: x, vals_input, parmlu, w, fx, u
            integer (ip), dimension(:), allocatable :: indx_input, jndx_input, luparm,&
                    p, q, lenc, lenr, locc, locr, iploc, iqloc, ipinv, iqinv
            integer (ip) :: N, nelem, inform, ierr
            type (gpf) :: gp

            N = 50_ip
            x = linspace(dble(0), dble(N-1), N)
            mat = double_diff_FD(x, 0_ip, ierr)
            call mat%set_value(N, N-1_ip, dble(0), ierr)
            call mat%set_value(N, N, dble(1), ierr)

            nelem = size(mat%vals)
            allocate (luparm(30_ip), parmlu(30_ip))
            luparm = 0_ip
            parmlu = dble(0)

            luparm(1) = 6_ip
            luparm(3) = 5_ip
            luparm(8) = 1_ip
            parmlu(1) = dble(100)
            parmlu(2) = dble(10)
            parmlu(3) = dble(3.0 * 10**(-13))
            parmlu(4) = dble(0.67 * 10**(-11))
            parmlu(5) = dble(0.67 * 10**(-11))
            parmlu(6) = dble(3.0)
            parmlu(7) = dble(0.3)
            parmlu(8) = dble(0.5)



            allocate (vals_input (nelem * 3), indx_input(nelem*3), jndx_input(nelem * 3))
            vals_input(1:nelem) = mat%vals
            indx_input(1:nelem) = mat%indx
            jndx_input(1:nelem) = mat%jndx

            allocate(p(N), q(N), lenc(N), lenr(N), locc(N), locr(N), iploc(N), iqloc(N), ipinv(N), iqinv(N), w(N))

            call lu1fac (N, N, nelem, nelem*3, luparm, parmlu, vals_input, indx_input, jndx_input, p, q,&
                    lenc, lenr, locc, locr, iploc, iqloc, ipinv, iqinv, w, inform)

            allocate (fx(N), u(N))
            fx = cos(x * 2 * 3.14159265 / dble(N))
            fx(N) = 0

            call lu6sol(5_ip, N, N, fx, u, 3* nelem, luparm, parmlu, vals_input, indx_input, jndx_input, p, q,&
            lenc, lenr, locc, locr, inform)

            call gp%plot(x, u)

        end subroutine test_1

end program test_lusol