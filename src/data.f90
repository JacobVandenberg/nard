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
        procedure, pass, public :: open_file_readonly
        procedure, pass, public :: read_real_vector
        procedure, pass, public :: read_integer_vector
        procedure, pass, public :: read_string
        procedure, pass, public :: read_real_matrix

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
            integer (kind=4) :: hdf5err

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
            integer (kind=4) :: hdf5err
            integer (HID_T)  :: space, dset
            type (C_PTR) :: f_ptr
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            ! Create dataspace.  Setting maximum size to be the current size.
            call h5screate_simple_f(int(1, kind=4), (/ size(vector, kind=ip) /), space, hdf5err)
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
            integer (kind=4) :: hdf5err

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
            integer (kind=4) :: hdf5err
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
            call h5screate_simple_f(int(3, kind=4), result%space_dim, result%dataspacev_id,  hdf5err)
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
            call h5screate_simple_f(int(3, kind=4), result%chunk_dim, result%memspace, hdf5err)
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

            integer (kind=4) :: hdf5err

            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            this%voffset(3) = this%voffset(3)+1_ip
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

        function read_real_vector(this, identifier, ierr)
            ! reads a real vector from h5 file
            !
            ! arguments:
            !   this: h5file object
            !   identifier: string telling us what to read from file
            !   ierr: internal error flag
            !
            implicit none
            ! BEGIN DECLARATIONS
            ! arguments
            class (h5file) :: this
            character (*) :: identifier
            integer (ip) :: ierr
            ! output
            real (rp), dimension(:), allocatable, target :: read_real_vector
            ! runtime
            integer (HID_T)  :: space, dset
            integer (kind=4) :: hdferr
            integer (HSIZE_T), dimension(1) :: dims, maxdims
            TYPE(C_PTR) :: f_ptr
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            call h5open_f(hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening h5 api"
                ierr = int( hdferr,kind=ip )
                return
            end if

            CALL h5dopen_f (this%file_id, identifier, dset, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if
            !
            ! Get dataspace and allocate memory for read buffer.
            !
            CALL h5dget_space_f(dset, space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in getting dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL h5sget_simple_extent_dims_f (space, dims, maxdims, hdferr)
            if (hdferr == int(-1, kind=4)) then
                Print *, "Error in getting dimensions"
                Print *, dims, maxdims, hdferr
                ierr = int( hdferr,kind=ip )
                return
            end if

            ALLOCATE(read_real_vector(dims(1)), stat=ierr)
            if (ierr /= 0) then
                Print *, "Allocation error. Code:", ierr
            end if
            !
            ! Read the data.
            !
            f_ptr = C_LOC(read_real_vector(1))
            CALL h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in reading data. Code:", hdferr
                ierr = int( hdferr,kind=ip)
                return
            end if
            !
            ! Close and release resources.
            !
            CALL h5dclose_f(dset , hdferr)
            if (hdferr /= 0) then
                Print *, "Error in closing dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL h5sclose_f(space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in closing dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            return
            ! END SUBROUTINE
        end function read_real_vector

        function read_integer_vector(this, identifier, ierr)
            ! reads a real vector from h5 file
            !
            ! arguments:
            !   this: h5file object
            !   identifier: string telling us what to read from file
            !   ierr: internal error flag
            !
            implicit none
            ! BEGIN DECLARATIONS
            ! arguments
            class (h5file) :: this
            character (*) :: identifier
            integer (ip) :: ierr
            ! output
            integer (ip), dimension(:), allocatable, target :: read_integer_vector
            ! runtime
            integer (HID_T)  :: space, dset
            integer (kind=4) :: hdferr
            integer (HSIZE_T), dimension(1) :: dims, maxdims
            TYPE(C_PTR) :: f_ptr
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            CALL h5dopen_f (this%file_id, identifier, dset, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if
            !
            ! Get dataspace and allocate memory for read buffer.
            !
            CALL h5dget_space_f(dset, space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in getting dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL h5sget_simple_extent_dims_f (space, dims, maxdims, hdferr)
            if (hdferr == int(-1, kind=4)) then
                Print *, "Error in getting dimensions"
                Print *, dims, maxdims, hdferr
                ierr = int( hdferr,kind=ip )
                return
            end if

            ALLOCATE(read_integer_vector(dims(1)))
            !
            ! Read the data.
            !
            f_ptr = C_LOC(read_integer_vector(1))
            CALL h5dread_f( dset, H5T_STD_I64LE, f_ptr, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in reading data"
                ierr = int( hdferr,kind=ip )
                return
            end if
            !
            ! Close and release resources.
            !
            CALL h5dclose_f(dset , hdferr)
            if (hdferr /= 0) then
                Print *, "Error in closing dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL h5sclose_f(space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in closing dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            return
            ! END SUBROUTINE
        end function read_integer_vector

        function read_string(this, identifier, ierr)
            ! reads a string from an h5 file
            !
            ! arguments:
            !   this: h5file object
            !   identifier: what to get from the h5file
            !   ierr: internal error flag
            implicit none
            ! BEGIN DECLARATIONS
            INTEGER(SIZE_T), parameter :: sdim = 7
            ! inputs
            class (h5file) :: this
            character (*) :: identifier
            ! outputs
            integer (ip) :: ierr
            CHARACTER (LEN=100, kind=c_char) :: read_string

            ! runtime
            INTEGER(HSIZE_T), DIMENSION(1:1) :: dims
            INTEGER(HSIZE_T), DIMENSION(1:1) :: maxdims
            INTEGER (kind=4):: hdferr
            INTEGER(HID_T)  :: filetype, memtype, space, dset
            INTEGER(ip) :: str_len
            TYPE(C_PTR) :: f_ptr
            type(C_PTR), dimension(:), allocatable, target :: rdata
            CHARACTER (LEN=100, kind=c_char), pointer :: read_string_pointer
            ! END DECLARATIONS

            ! BEGIN FUNCTION

            CALL h5dopen_f(this%file_id, identifier, dset, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if

            !
            ! Get the datatype.
            !
            CALL H5Dget_type_f(dset, filetype, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening filetype"
                ierr = int( hdferr,kind=ip )
                return
            end if
            !CALL H5Tget_size_f(filetype, str_size, hdferr)
            !if (hdferr /= 0) then
            !    Print *, "Error in getting size"
            !    ierr = int( hdferr,kind=ip )
            !    return
            !end if
            !Print *, str_size
            !
            ! Get dataspace and allocate memory for read buffer.
            !
            CALL H5Dget_space_f(dset, space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in getting dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

            if (hdferr == int(-1, kind=ip)) then
                Print *, "Error in getting dimensions", dims, maxdims, hdferr
                ierr = int( hdferr,kind=ip )
                return
            end if
            ALLOCATE(rdata(dims(1)))

            !CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
            !if (hdferr /= 0) then
            !    Print *, "Error in copying data"
            !    ierr = int( hdferr,kind=ip )
            !    return
            !end if
            !CALL H5Tset_size_f(memtype, str_size, hdferr)
            !if (hdferr /= 0) then
            !    Print *, "Error in setting size"
            !    ierr = int( hdferr,kind=ip )
            !    return
            !end if
            !
            ! Read the data.
            !
            f_ptr = C_LOC(rdata(1))
            CALL H5Dread_f(dset, filetype, f_ptr, hdferr, space)
            if (hdferr /= 0) then
                Print *, "Error in reading data"
                ierr = int( hdferr,kind=ip )
                return
            end if
            call C_F_POINTER( rdata(1), read_string_pointer )
            str_len = 0_ip
            do
                if (read_string_pointer(str_len+1_ip:str_len+1_ip) == C_NULL_CHAR .or. str_len >= 99_ip) exit
                str_len = str_len + 1_ip
            end do
            write (read_string, '(A)') read_string_pointer(1_ip:str_len)

            CALL h5dclose_f(dset , hdferr)
            CALL h5sclose_f(space, hdferr)
            CALL H5Tclose_f(filetype, hdferr)


            return
            ! END FUNCTION

        end function read_string

        function read_real_matrix(this, identifier, ierr)
            ! reads a real vector from h5 file
            !
            ! arguments:
            !   this: h5file object
            !   identifier: string telling us what to read from file
            !   ierr: internal error flag
            !
            implicit none
            ! BEGIN DECLARATIONS
            ! arguments
            class (h5file) :: this
            character (*) :: identifier
            integer (ip) :: ierr
            ! output
            real (rp), dimension(:, :), allocatable, target :: read_real_matrix
            ! runtime
            integer (HID_T)  :: space, dset
            integer (kind=4) :: hdferr
            integer (HSIZE_T), dimension(2) :: dims, maxdims
            TYPE(C_PTR) :: f_ptr
            ! END DECLARATIONS

            ! BEGIN SUBROUTINE
            call h5open_f(hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening h5 api"
                ierr = int( hdferr,kind=ip )
                return
            end if

            CALL h5dopen_f (this%file_id, identifier, dset, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in opening dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if
            !
            ! Get dataspace and allocate memory for read buffer.
            !
            CALL h5dget_space_f(dset, space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in getting dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL h5sget_simple_extent_dims_f (space, dims, maxdims, hdferr)
            if (hdferr == int(-1, kind=4)) then
                Print *, "Error in getting dimensions"
                Print *, dims, maxdims, hdferr
                ierr = int( hdferr,kind=ip )
                return
            end if

            ALLOCATE(read_real_matrix(dims(1), dims(2)), stat=ierr)
            if (ierr /= 0) then
                Print *, "Allocation error. Code:", ierr
            end if
            !
            ! Read the data.
            !
            f_ptr = C_LOC(read_real_matrix(1, 1))
            CALL h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in reading data. Code:", hdferr
                ierr = int( hdferr,kind=ip)
                return
            end if
            !
            ! Close and release resources.
            !
            CALL h5dclose_f(dset , hdferr)
            if (hdferr /= 0) then
                Print *, "Error in closing dataset"
                ierr = int( hdferr,kind=ip )
                return
            end if
            CALL h5sclose_f(space, hdferr)
            if (hdferr /= 0) then
                Print *, "Error in closing dataspace"
                ierr = int( hdferr,kind=ip )
                return
            end if
            return
            ! END SUBROUTINE
        end function read_real_matrix

        subroutine open_file_readonly(this, filename, ierr)
            ! opens a file for read only
            !
            ! inputs:
            !   filename (string): the name of the file to open
            !   this: h5file
            !
            ! outputs:
            !   ierr: internal error flag
            !
            implicit none
            ! BEGIN DECLARATIONS
            ! arguments
            class (h5file) :: this
            character (*) :: filename
            integer (ip), intent(out) :: ierr
            ! runtime
            integer (kind=4) :: hdferr
            ! END DECLARATIONS
            ierr = 0_ip;
            call h5open_f(hdferr)
            if (hdferr /= 0) then
                ierr = int( hdferr,kind=ip )
                return
            end if
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, this%file_id, hdferr)
            if (hdferr /= 0) then
                ierr = int( hdferr,kind=ip )
                return
            end if
            return
        end subroutine open_file_readonly

end module data