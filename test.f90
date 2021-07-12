! Created by  on 6/7/21.

program test
    use config_m
    use ogpf
    implicit none
    integer (kind=8) :: N = 30, spacing = 0, BC=0, i, j, k
    !real (kind=8), dimension(:), allocatable :: x, y, z
    real (kind=8), dimension(2) :: domain = (/ 0.0, 1.0 /)
    type (config) :: conf
    real (kind=8), dimension(:, :), allocatable :: xx, yy, zz
    type(gpf) :: gp

    allocate (conf%x(N))
    allocate (conf%y(N))
    allocate (conf%z(N))

    allocate (conf%xxx(N,N,N))
    allocate (conf%yyy(N,N,N))
    allocate (conf%zzz(N,N,N))

    conf%nx = N
    conf%ny = N
    conf%nz = N

    call gridline(N, BC, (domain *2)-1, spacing, conf%x)
    call gridline(N, BC, (domain *4)-2, spacing, conf%y)
    call gridline(N, BC, (domain *6)-3, spacing, conf%z)

    call make_mesh(conf)

    allocate(xx(N, N), yy(N, N), zz(N, N) )

    xx = conf%xxx(:, :, 1)
    yy = conf%yyy(:, :, 1)

    zz = cos(xx) * sin(yy)
    call gp%surf(xx, yy, zz)

end program test