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
            integer(ip), dimension(:), allocatable :: numbers, key

            allocate(numbers(10_ip))
            numbers = int((/ 1, 2, 3, 5, 4, 6, 7, 8, 9, 0/), kind=ip)
            key = int((/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/), kind=ip)
            call order_int(numbers, 1_ip, 10_ip, key)
            Print *, key

        end subroutine test_order_int

        subroutine test_arange_int
            implicit none
            integer(ip), allocatable, dimension(:) :: key
            integer(ip) :: ierr

            call arange_int(4_ip, 10_ip, key, ierr)
            Print *, key
        end subroutine test_arange_int

end program test_helpers