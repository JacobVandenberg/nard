program test_mkl
    use precision
    use sparse_matrices
    use ogpf
    use diff
    implicit none

    integer (ip), dimension(64) :: pt
    integer (ip) :: mtype
    integer (ip), dimension(64) :: iparm
    type (coo_matrix) :: mat
    type (csr_matrix) :: mat_csr
    type (sparse_lu_decomposition) :: lu_mat
    real (rp), dimension(:), allocatable :: x
    real (rp), dimension(:), allocatable :: fx, y
    integer (ip), dimension(:), allocatable :: perm
    integer (ip) :: N, ierr
    type (gpf) :: gp

    mtype = 11
    call pardisoinit(pt, mtype, iparm)

    Print *, pt
    Print *, mtype
    Print *, iparm

    N = 1000000_ip
    x = linspace(dble(0), dble(N-1), N)
    mat = double_diff_FD(x, 0_ip, ierr)

    !mat = laplacian2D(x, x, 0_ip, 0_ip, ierr)

    !Print *, mat%vals
    call mat%set_value(N, N-1_ip, dble(0.0001), ierr)
    call mat%set_value(N, N, dble(1), ierr)
    mat_csr = coo_to_csr(mat, ierr)
    allocate (perm(N))



    call pardiso( pt, 1_ip, 1_ip, mtype, 12_ip, mat_csr%m, mat_csr%vals, mat_csr%rows, mat_csr%jndx,&
            perm, 1_ip, iparm, 1_ip, fx, y, ierr)
    Print *, ierr
    allocate (fx(N), y(N))
    fx(:) = cos(x * 2 * 3.14159265 / dble(N))
    fx(N) = 0
    call pardiso( pt, 1_ip, 1_ip, mtype, 33_ip, mat_csr%m, mat_csr%vals, mat_csr%rows, mat_csr%jndx,&
            perm, 1_ip, iparm, 1_ip, fx, y, ierr)


    call gp%plot(x, reshape(y, (/N/)))


end program test_mkl