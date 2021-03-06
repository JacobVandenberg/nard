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

        u_out(:, 1) = user_params(5) * (u_in(:, 1) - u_in(:, 1)**3_rp - u_in(:, 2))
        u_out(:, 2) = user_params(5) * (user_params(3) * (u_in(:, 1) -&
            user_params(2) * u_in(:, 2) - user_params(1)))
        
        return
    end subroutine reaction_term

    subroutine diffusivity(u_in, x_in, t, diff_out, user_params)
        implicit none
        real (rp), dimension(:, :), intent(in) :: u_in, x_in
        real (rp), intent(in) :: t
        real (rp), intent(out), dimension(:, :) :: diff_out
        real (rp), dimension(:), intent(in) :: user_params

        diff_out(:, 1_ip) = user_params(5) * 1_rp
        diff_out(:, 2_ip) = user_params(5) * 4_rp

    end subroutine diffusivity

    subroutine advection_coefficient(u_in, x_in, t, adv_out, user_params)
        implicit none
        real (rp), dimension(:, :), intent(in) :: u_in, x_in
        real (rp), intent(in) :: t
        real (rp), intent(out), dimension(:, :) :: adv_out
        real (rp), dimension(:), intent(in) :: user_params

        adv_out(:, 1) = user_params(4) * 1/200 * (100 + x_in(:, 2))
        adv_out(:, 2) = user_params(4) * 1/200 * (100 + x_in(:, 2))

    end subroutine advection_coefficient

end module user_functions