! Created by  on 20/7/21.

program nard2D
    use user_config, only: get_config
    use evolvePDE
    use config_m
    implicit none

    integer (ip) :: ierr
    type (config2d) :: conf
    Print *, "getting user config"
    conf = get_config()
    Print *, "evolving PDE"
    call IMEX_evolve_CN_heun_2D(conf, ierr)


end program nard2D