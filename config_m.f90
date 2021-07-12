! Created by  on 6/7/21.

module config_m
    use helpers, only: linspace, meshgrid3
    implicit none
    !private
    public gridline
    type config
        real (kind=8), dimension(:), allocatable :: x
        real (kind=8), dimension(:), allocatable :: y
        real (kind=8), dimension(:), allocatable :: z
        real (kind=8), dimension(:, :, :), allocatable :: xxx
        real (kind=8), dimension(:, :, :), allocatable :: yyy
        real (kind=8), dimension(:, :, :), allocatable :: zzz
        real (kind=8), dimension(:, :), allocatable :: IC
        real (kind=8), dimension(:, :, :), allocatable :: plot_interval
        integer (kind=8) :: nx, ny, nz

    contains
        procedure :: make_mesh

    end type config

    !type Gridline

    !end type Gridline

    contains
        subroutine gridline(N, BC, domain, spacing, g)
            implicit none
            integer (kind=8), intent(in) :: N
            integer (kind=8), intent(in) :: BC
            real (kind=8), intent(in) :: domain(2)
            integer (kind=8), intent(in) :: spacing
            real (kind=8), intent(out) :: g(N)

            call linspace (domain(1), domain(2), N, g)
            return
        end subroutine gridline

        subroutine make_mesh(this)
            implicit none
            class(config), intent(inout) :: this
            call meshgrid3(this%x, this%y,this%z, this%nx, this%ny, this%nz, this%xxx, this%yyy, this%zzz)

        end subroutine make_mesh


end module config_m