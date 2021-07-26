! Created by  on 17/7/21.

program test_data
    use data
    Print *, 'test_new_file'
    call test_new_file
    Print *, 'test_save_real_vector'
    call test_save_real_vector
    Print *, 'test_allocate_chunked_3D_space'
    call test_allocate_chunked_3D_space
    Print *, 'test_chunked_3D_space'
    call test_chunked_3D_space
    Print *, 'test_read_real_vector'
    call test_read_real_vector
    Print *, 'test_read_integer_vector'
    call test_read_integer_vector
    Print *, 'test_read_string'
    call test_read_string

    contains
        subroutine test_new_file
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            call file%new_file( "src/tests/test1.h5", ierr )
            Print *, ierr
            call file%close(ierr)

        end subroutine test_new_file

        subroutine test_save_real_vector
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            real (rp), dimension(:), allocatable :: test_vec

            allocate( test_vec(10) )

            test_vec = real( (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/) ,kind=rp )

            call file%new_file( "src/tests/test2.h5", ierr )

            call file%save_real_vector( test_vec, "testname", ierr)
            Print *, ierr

            call file%close(ierr)

        end subroutine test_save_real_vector

        subroutine test_allocate_chunked_3D_space
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            type (chunked_3D_space) :: space_3d


            call file%new_file( "src/tests/test3.h5", ierr )

            space_3d =  file%allocate_chunked_3D_space( "3Dspace", int((/ 100, 2, 40 /), kind=ip), ierr )

            call file%close(ierr)


        end subroutine test_allocate_chunked_3D_space

        subroutine test_chunked_3D_space
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            type (chunked_3D_space) :: space_3d

            real (rp), dimension(100, 2) :: next_page

            integer (ip) :: i


            call file%new_file( "src/tests/test4.h5", ierr )

            space_3d =  file%allocate_chunked_3D_space( "3Dspace", int((/ 100, 2, 40 /), kind=ip), ierr )

            do i=1, 40
                call random_number (next_page)
                call space_3d%insert_page(next_page, ierr)
            end do

            call file%close(ierr)


        end subroutine test_chunked_3D_space

        subroutine test_read_real_vector
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            real (rp), dimension(:), allocatable :: rdata

            call file%open_file_readonly("src/tests/test_read.h5", ierr)
            rdata = file%read_real_vector("rparams", ierr)
            Print *, rdata
        end subroutine test_read_real_vector
        subroutine test_read_integer_vector
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            integer (ip), dimension(:), allocatable :: idata

            call file%open_file_readonly("src/tests/test_read.h5", ierr)
            idata = file%read_integer_vector("iparams", ierr)
            Print *, idata
        end subroutine test_read_integer_vector

        subroutine test_read_string
            implicit none
            type (h5file) :: file
            integer (ip) :: ierr
            CHARACTER (LEN=100):: sdata
            !character (len=100,kind=c_char), dimension(:), allocatable :: sdata

            call file%open_file_readonly("src/tests/test_read.h5", ierr)
            sdata = file%read_string("sparams5", ierr)
            Print *, sdata
        end subroutine test_read_string

end program test_data