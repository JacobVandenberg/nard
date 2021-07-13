! Created by  on 10/7/21.

program test_helpers
    use helpers
    use ogpf
    implicit none
    Print *, "order_int"
    call test_order_int
    Print * , 'arange_int'
    call test_arange_int


    contains
        subroutine test_order_int
            implicit none
            integer, dimension(:), allocatable :: numbers, key

            allocate(numbers(10))
            numbers = (/ 1, 2, 3, 5, 4, 6, 7, 8, 9, 0/)
            key = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
            call order_int(numbers, 1, 10, key)
            Print *, key

        end subroutine test_order_int

        subroutine test_arange_int
            implicit none
            integer, allocatable, dimension(:) :: key
            integer :: ierr

            call arange_int(4, 10, key, ierr)
            Print *, key
        end subroutine test_arange_int

end program test_helpers