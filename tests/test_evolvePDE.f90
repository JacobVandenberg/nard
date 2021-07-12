! Created by  on 11/7/21.

program test_evolvePDE
    use config_m
    implicit none
    call test_euler
    contains
        subroutine test_euler
            implicit none
            type (config) :: conf
            integer :: N
            real (kind=8), allocatable, dimension(:) :: x, y
            N = 30

            allocate(conf%IC(N*N, 2))

            ! generate grid lines
            call my_linspace(dble(-1), dble(1), N, x)
            call my_linspace(dble(-1), dble(1), N, y)


        end subroutine test_euler
end program test_evolvePDE