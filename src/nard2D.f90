! Created by  on 20/7/21.

program nard2D
    use user_functions
    use evolvePDE
    use config_m
    use data
    implicit none

    ! BEGIN DECLARATIONS
    integer (ip) :: ierr
    type (config2d) :: conf
    character (len=128) :: config_filename
    ! END DECLARATIONS

    ! BEGIN PROGRAM
    Print *, "getting user RHS function"
    conf%explicit_rhs => reaction_term

    Print *, "getting config"
    IF (COMMAND_ARGUMENT_COUNT() /= 1_ip)THEN
        Print *, 'Please provide config filename'
        STOP
    ENDIF
    call get_command_argument(1_ip, config_filename)
    call conf%read_config(config_filename, ierr)

    Print *, "evolving PDE"
    if (conf%timestepping_method == 2) then
        Print *, "Using IMEX_evolve_CN_heun_2D"
        call IMEX_evolve_CN_heun_2D(conf, ierr)
    elseif (conf%timestepping_method == 1) then
        Print *, "Using IMEX_evolve_impliciteuler_euler_2D"
        call IMEX_evolve_impliciteuler_euler_2D(conf, ierr)
    else
        Print *, "Invalid timestepping method ", conf%timestepping_method,&
                "Please choose 1 for implicit Euler / Euler or 2 for Crank-Nicholson / Heun"
    end if
    ! END PROGRAM
end program nard2D