! Created by  on 24/7/21.
module user_functions
    use precision
    implicit none
    contains

    subroutine reaction_term(u_in, x_in, t, u_out)
        implicit none
        real (rp), dimension(:, :), intent(in) :: u_in, x_in
        real (rp), intent(in) :: t
        real (rp), intent(out), dimension(:, :) :: u_out
        real (rp) :: a, b, gamma

        a = 0.2
        b = 1.3
        gamma = 1000_ip

        u_out(:, 1) = gamma * (a - u_in(:, 1) + u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
        u_out(:, 2) = gamma * (b -  u_in(:, 1) * u_in(:, 2) * u_in(:, 1))
        return
    end subroutine reaction_term

end module user_functions

