! Created by  on 6/7/21.

module config_m
    use helpers, only: my_linspace, meshgrid3
    use ogpf
    implicit none
    !private
    type config2d
        real (kind=8), dimension(:), allocatable :: x
        real (kind=8), dimension(:), allocatable :: y
        real (kind=8), dimension(:, :), allocatable :: IC, xx, yy
        real (kind=8) :: dt, t_max, plot_interval
        integer :: BCx, BCy
        procedure (example_explicit_rhs), pointer, nopass :: explicit_rhs


    contains
        procedure, pass, public :: make_mesh => make_mesh2D

    end type config2d

    interface
        subroutine example_explicit_rhs(u_in, x_in, t_in, u_out)
            ! example right hand side
            implicit none
            real (kind=8), dimension(:, :), intent(in) :: u_in, x_in
            real (kind=8), intent(in) :: t_in
            real (kind=8), dimension(:, :), intent(out) :: u_out
        end subroutine example_explicit_rhs
    end interface

    contains

        subroutine make_mesh2D(this, ierr)
            ! generates the 2d mesh from the 1d gridlines

            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            class (config2d), intent(inout) :: this

            ! outputs
            integer, intent(out) :: ierr

            ! END DECLARATIONS

            ! BEGIN FUNCTION
            call meshgrid(this%xx, this%yy, this%x, this%y, ierr)
            return
            ! END FUNCTION
        end subroutine make_mesh2D

end module config_m