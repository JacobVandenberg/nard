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
        u_out = 0_rp
        return
    end subroutine reaction_term

    subroutine diffusivity(u_in, x_in, t, diff_out, user_params)
        implicit none
        real (rp), dimension(:, :), intent(in) :: u_in, x_in
        real (rp), intent(in) :: t
        real (rp), intent(out), dimension(:, :) :: diff_out
        real (rp), dimension(:), intent(in) :: user_params

        diff_out(:, 1) = t**2.0_rp

    end subroutine diffusivity

end module user_functions