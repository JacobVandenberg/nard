! Created by  on 25/7/21.

program test_config_m
    use config_m
    implicit none

    Print *, 'test_read_config2d'
    call test_read_config2d

    contains

        subroutine test_read_config2d
            implicit none
            type(config2d) :: conf
            integer (ip) :: ierr

            call conf%read_config("/home/jacob/coding/RD_project/nard/src/tests/test_config.h5", ierr)

            Print *, conf%DBCx_plus_mask

        end subroutine test_read_config2d

end program test_config_m