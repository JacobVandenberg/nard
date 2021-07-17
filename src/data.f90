! Created by  on 17/7/21.

module data
    use hdf5
    use precision
    implicit none

    type h5file
        integer (HID_T) :: file_id

        contains
        procedure, pass, public :: new_file
        procedure, pass, public :: save_real_vector
        procedure, pass, public :: close
        procedure, pass, public :: allocate_chunked_3D_space

    end type h5file

    type chunked_3D_space
        integer (HID_T) :: file_id
        integer (HID_T) :: dsetv_id
        integer (HID_T) :: dataspacev_id
        integer (HID_T) :: memspace
        integer (HSIZE_T), dimension(3) :: space_dim, chunk_dim, voffset, vchunkcount

        contains
        procedure, pass, public :: insert_page


    end type chunked_3D_space

    contains

        subroutine new_file(this, filename, ierr)
            !
            ! creates a new file
            !
            ! output: ierr
            integer (ip), intent(out) :: ierr
            class (h5file), intent(inout) :: this
            character (len=*) :: filename
            integer :: hdf5err

            call h5open_f(hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, this%file_id, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if
            ierr = int(hdf5err, kind=ip)

        end subroutine new_file

        subroutine save_real_vector(this, vector, name, ierr)
            !
            ! saves a vector of reals to the file
            !
            ! input: vector: 1D vector of reals
            !   name: name for the data
            !
            ! output: ierr
            implicit none

            ! BEGIN DECLARATIONS
            ! inputs
            class (h5file), intent(in) :: this
            real (rp), dimension(:), intent(in), target :: vector
            character (len=*), intent(in) :: name
            ! outputs
            integer (ip), intent(out) :: ierr
            ! runtime
            integer :: hdf5err
            integer (HID_T)  :: space, dset
            type (C_PTR) :: f_ptr
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            ! Create dataspace.  Setting maximum size to be the current size.
            call h5screate_simple_f(1, (/ size(vector, kind=ip) /), space, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            call h5dcreate_f(this%file_id, name, H5T_NATIVE_DOUBLE, space, dset, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            f_ptr = C_LOC(vector(1))
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            CALL h5dclose_f(dset , hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            CALL h5sclose_f(space, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            return
            ! END SUBROUTINE

        end subroutine save_real_vector

        subroutine close(this, ierr)
            implicit none

            integer (ip), intent(out) :: ierr
            class (h5file), intent(inout) :: this
            integer :: hdf5err

            CALL h5fclose_f(this%file_id , hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if
        end subroutine close

        type (chunked_3D_space) function allocate_chunked_3D_space(this, name, space_dim, ierr)
            !
            ! allocates 3d space of size space dim, with chunk size of chunk dim.
            ! chunks must be pages such that space_dim(1:2) == chunk_dim(1:2)
            !
            ! arguments:
            !   inputs: this: the h5file
            !   space_dim: 3-vector specifying size of array
            !
            ! outputs:
            !   ierr: internal error flag
            !   return value: chunked_3D_space object
            implicit none
            ! BEGIN DECLARATIONS
            ! inputs
            class (h5file), intent(in) :: this
            integer (ip), dimension(3), intent(in) :: space_dim

            ! outputs
            integer (ip), intent(out) :: ierr
            type (chunked_3D_space) :: result
            character (len=*), intent(in) :: name

            ! runtime
            integer :: hdf5err
            integer (ip), dimension(3) :: chunk_dim
            integer (HSIZE_T), dimension(3) :: chunk_dimensions, voffset, vchunkcount

            ! END DECLARATIONS

            ! BEGIN SUBROUTINE

            chunk_dim = space_dim
            chunk_dim(3) = 1_ip

            result%file_id = this%file_id
            result%space_dim = int( space_dim, kind=HSIZE_T )
            result%chunk_dim = int( chunk_dim, kind=HSIZE_T )
            result%voffset = int( (/0, 0, -1/), kind=HSIZE_T )
            result%vchunkcount = result%chunk_dim


            ! Create the data space for the binned dataset.
            call h5screate_simple_f(3, result%space_dim, result%dataspacev_id,  hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            ! Create the chunked dataset.
            call h5dcreate_f(result%file_id, name, H5T_NATIVE_DOUBLE, result%dataspacev_id, result%dsetv_id, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if

            ! Create the memory space for the selection
            call h5screate_simple_f(3, result%chunk_dim, result%memspace, hdf5err)
            if (hdf5err /= 0) then
                ierr = int(hdf5err, kind=ip)
                return
            end if
            allocate_chunked_3D_space = result
            ! END SUBROUTINE
        end function allocate_chunked_3D_space

        subroutine insert_page(this, new_page, ierr)
            !
            ! inserts a new page into the allocated space
            !
            ! inputs:
            !   this: chunked_3d_space
            !   new_page: the new page to put into the 3d space
            !
            ! outputs: ierr: internal error flag
            !
            implicit none

            ! BEGIN DELARATIONS
            class (chunked_3D_space), intent(inout) :: this
            real (rp), dimension(:, :), intent(in) :: new_page
            integer (ip), intent(out) :: ierr

            integer :: hdf5err

            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            this%voffset(3) = this%voffset(3)+1

            call h5sselect_hyperslab_f(this%dataspacev_id, H5S_SELECT_SET_F, this%voffset, this%vchunkcount, hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error in h5sselect_hyperslab_f. error code ", hdf5err
                return
            end if

            ! Write the data to the dataset.
            call h5dwrite_f(this%dsetv_id, H5T_NATIVE_DOUBLE, reshape(new_page, this%chunk_dim), this%chunk_dim,&
                    hdf5err, this%memspace, this%dataspacev_id)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error in h5dwrite_f. error code ", hdf5err
                return
            end if
            call h5sclose_f(this%dataspacev_id, hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error in h5sclose_f. error code ", hdf5err
                return
            end if
            call h5dget_space_f(this%dsetv_id, this%dataspacev_id, hdf5err)
            if (hdf5err /=0) then
                Print *, "HDF5 ERROR: error in h5dget_space_f. error code ", hdf5err
                return
            end if
            ! END SUBROUTINE
        return
        end subroutine insert_page

end module data