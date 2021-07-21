! Created by  on 6/7/21.

module config_m
    use helpers
    use ogpf
    use precision
    implicit none
    !private
    type config2d
        real (rp), dimension(:), allocatable :: x, y, diffusion_consts, DBCx_plus, DBCx_minus, DBCy_plus, DBCy_minus
        logical, dimension(:), allocatable :: DBCx_plus_mask, DBCx_minus_mask, DBCy_plus_mask, DBCy_minus_mask
        real (rp), dimension(:, :), allocatable :: IC, xx, yy
        real (rp) :: dt, t_max, plot_interval
        integer (ip) :: BCx, BCy, savenum, max_save_size
        procedure (example_explicit_rhs), pointer, nopass :: explicit_rhs
        character (:), allocatable :: savefilename, plotfilename


    contains
        procedure, pass, public :: make_mesh => make_mesh2D

    end type config2d

    type config3d
        real (rp), dimension(:), allocatable :: x, y, z, diffusion_consts
        real (rp), dimension(:, :), allocatable :: IC
        real (rp), dimension(:, :, :), allocatable :: xxx, yyy, zzz
        real (rp) :: dt, t_max, plot_interval
        integer (ip) :: BCx, BCy, BCz, savenum, max_save_size
        procedure (example_explicit_rhs), pointer, nopass :: explicit_rhs
        character (:), allocatable :: savefilename, plotfilename


    contains
        procedure, pass, public :: make_mesh => make_mesh3D

    end type config3d

    interface
        subroutine example_explicit_rhs(u_in, x_in, t_in, u_out)
            import :: rp
            ! example right hand side
            implicit none
            real (rp), dimension(:, :), intent(in) :: u_in, x_in
            real (rp), intent(in) :: t_in
            real (rp), dimension(:, :), intent(out) :: u_out
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
            integer (ip), intent(out) :: ierr

            ! END DECLARATIONS

            ! BEGIN FUNCTION
            call meshgrid(this%xx, this%yy, this%x, this%y, ierr)
            return
            ! END FUNCTION
        end subroutine make_mesh2D

        subroutine make_mesh3D(this, ierr)
            ! generates the 2d mesh from the 1d gridlines

            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            class (config3d), intent(inout) :: this

            ! outputs
            integer (ip), intent(out) :: ierr

            ! END DECLARATIONS

            ! BEGIN FUNCTION
            call meshgrid3(this%x, this%y, this%z, this%xxx, this%yyy, this%zzz, ierr)
            return
            ! END FUNCTION
        end subroutine make_mesh3D

end module config_m