! Created by  on 6/7/21.

module config_m
    use helpers
    use ogpf
    use precision
    use data
    implicit none
    !private
    type config2d
        real (rp), dimension(:), allocatable :: x, y, diffusion_consts, DBCx_plus, DBCx_minus, DBCy_plus, DBCy_minus,&
                user_params
        logical, dimension(:), allocatable :: DBCx_plus_mask, DBCx_minus_mask, DBCy_plus_mask, DBCy_minus_mask
        real (rp), dimension(:, :), allocatable :: IC, xx, yy
        real (rp) :: dt, t_max, plot_interval
        integer (ip) :: BCx, BCy, savenum, max_save_size, timestepping_method
        procedure (example_explicit_rhs), pointer, nopass :: explicit_rhs, diffusivity,&
         advection_coefficient, diffusivity_derivative
        character (:), allocatable :: savefilename, plotfilename
        integer (ip), dimension(64) :: iparams
        logical :: advection


    contains
        procedure, pass, public :: make_mesh => make_mesh2D
        procedure, pass, public :: read_config => read_config2d

    end type config2d

    type config3d
        real (rp), dimension(:), allocatable :: x, y, z, diffusion_consts, user_params
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
        subroutine example_explicit_rhs(u_in, x_in, t_in, u_out, dudxx, user_params)
            import :: rp
            ! example right hand side
            implicit none
            real (rp), dimension(:, :), intent(in) :: u_in, x_in, dudxx
            real (rp), intent(in) :: t_in
            real (rp), dimension(:, :), intent(out) :: u_out
            real (rp), dimension(:), intent(in) :: user_params
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

        subroutine read_config2d(this, filename, ierr)
            ! reads a config from an h5 file
            !
            ! inputs: filename: the file to read the config from
            ! outputs: ierr: the internal error flag
            !
            implicit none
            ! BEGIN DECLARATIONS
            ! arguments
            class (config2d), intent(inout) :: this
            character (len=*), intent(in) :: filename
            integer (ip), intent(out) :: ierr

            ! runtime
            type (h5file) :: file_o
            real (rp), dimension(64) :: rparams
            integer (ip), dimension(64) :: iparams
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            real (rp), dimension(:, :), allocatable :: IC
            character (:), allocatable :: savefilename, plotfilename

            call file_o%open_file_readonly(trim(filename), ierr)
            if (ierr /= 0) then
                Print *, "Error in opening file ", trim(filename), "  ", ierr
            end if
            this%x = file_o%read_real_vector("x", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading x"
            end if
            this%y = file_o%read_real_vector("y", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading y"
            end if
            this%diffusion_consts = file_o%read_real_vector("diffusion_consts", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading diffusion constants"
            end if
            this%DBCx_plus = file_o%read_real_vector("DBCx_plus", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCx_plus"
            end if
            this%DBCx_minus = file_o%read_real_vector("DBCx_minus", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCx_minus"
            end if
            this%DBCy_plus = file_o%read_real_vector("DBCy_plus", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCy_plus"
            end if
            this%DBCy_minus = file_o%read_real_vector("DBCy_minus", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCy_minus"
            end if

            this%IC = file_o%read_real_matrix("IC", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading IC"
            end if

            rparams = file_o%read_real_vector("rparams", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading rparams"
            end if
            this%dt = rparams(1)
            this%t_max = rparams(2)
            this%plot_interval = rparams(3)

            this%user_params = file_o%read_real_vector("user_params", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading user_params"
            end if

            iparams = file_o%read_integer_vector("iparams", ierr)
            this%iparams = iparams
            if (ierr /= 0) then
                Print *, "Error in reading iparams"
            end if
            this%BCx = iparams(3)
            this%BCy = iparams(4)
            this%savenum = iparams(1)
            this%max_save_size = iparams(2)
            this%timestepping_method = iparams(6)
            this%advection = iparams(8)

            this%DBCx_plus_mask = file_o%read_integer_vector("DBCx_plus_mask", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCx_plus_mask"
            end if
            this%DBCx_minus_mask = file_o%read_integer_vector("DBCx_minus_mask", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCx_minus_mask"
            end if
            this%DBCy_plus_mask = file_o%read_integer_vector("DBCy_plus_mask", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCy_plus_mask"
            end if
            this%DBCy_minus_mask = file_o%read_integer_vector("DBCy_minus_mask", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading DBCy_minus_mask"
            end if

            this%savefilename = file_o%read_string("savefilename", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading savefilename"
            end if
            this%plotfilename = file_o%read_string("plotfilename", ierr)
            if (ierr /= 0) then
                Print *, "Error in reading plotfilename"
            end if

            return
            ! END SUBROUTINE
        end subroutine read_config2d

end module config_m