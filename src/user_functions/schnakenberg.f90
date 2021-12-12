! Created by  on 24/7/21.
module user_functions
    use precision
    implicit none
    contains

    subroutine reaction_term(u_in, x_in, t, u_out, user_params)
        implicit none
        real (rp), dimension(:, :), intent(in) :: u_in, x_in
        real (rp), intent(in) :: t
        real (rp), intent(out), dimension(:, :) :: u_out
        real (rp), dimension(:), intent(in) :: user_params
        ! user_params = (/a, b, gamma/)
            u_out(:, 1) = user_params(3) * (user_params(1) - u_in(:, 1) + u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
            u_out(:, 2) = user_params(3) * (user_params(2) -  u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
        return
    end subroutine reaction_term

    subroutine diffusivity(u_in, x_in, t, diff_out, user_params)
        implicit none
        real (rp), dimension(:, :), intent(in) :: u_in, x_in
        real (rp), intent(in) :: t
        real (rp), intent(out), dimension(:, :) :: diff_out
        real (rp), dimension(:), intent(in) :: user_params

    end subroutine diffusivity

end module user_functions