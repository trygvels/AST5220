module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function

contains

  subroutine initialize_rec_mod
    implicit none

    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx
    real(dp)     :: f, n_e0, X_e0, xstart, xstop, yp1, ypn, eps, hmin, step
    real(dp)     :: X_econst, C_r

    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopo

    ! Spline variables
    yp1 = 1.d30
    ypn = 1.d30

    ! Integration variables
    eps = 1.d-10
    hmin = 0.d0


    ! Allocating arrays
    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    ! Uniform x-grid with 1000 points
    dx = (xstop-xstart)/(n-1)
    x_rec(1) = xstart
    do i=2,n
      x_rec(i) = x_rec(i-1) + dx
    end do

    step = abs(1.d-3*(x_rec(1)-x_rec(2))) ! Step length for ODE

    ! Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1, n
       n_b = Omega_b*rho_c/(m_H*exp(x_rec(i))**3)

       if (use_saha) then
         ! Use the Saha equation
          T_b = T_0/exp(x_rec(i))
          X_econst = ((m_e*k_b*T_b)/(2.d0*pi*hbar**2))**1.5d0*exp(-epsilon_0/(k_b*T_b))/n_b
          X_e(i) = (-X_econst + sqrt(X_econst**2 +4.d0*X_econst))/2.d0
      if (X_e(i) < saha_limit) use_saha = .false.
      else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, step, hmin, dX_edx, bsstep, output)
       end if
       n_e(i) = X_e(i)*n_b ! Electron density
    end do

    ! Write to file - x_rec, X_e
    open(1, file="Xe.dat", action="write",status="replace")
    do i=1, n
       write(1,*) x_rec(i), X_e(i)
    end do
    close(1)

    ! Task: Compute splined (log of) electron density function
    call spline(log(n_e), eta,yp1,ypn,n_e2) !n_e2 is now log

    ! Task: Compute optical depth at all grid points
    tau(n) = 0.d0
    do i=n-1,1,-1
      tau(i) = tau(i+1)
      call odeint(tau(i:i),x_rec(i+1),x_rec(i),eps,step,hmin,dtaudx,bsstep,output)
    end do

    ! Task: Compute splined (log of) optical depth
    call spline(x_rec, tau, yp1, ypn,tau2)
    ! Task: Compute splined second derivative of (log of) optical depth
    call spline(x_rec,tau2,yp1,ypn,tau22)

    ! Write to file - x_rec, X_e
    open(2, file="tau.dat", action="write",status="replace")
    do i=1, n
       write(2,*) tau(i),tau2(i),tau22(i)
    end do
    close(2)

    ! Task: Compute splined visibility function
    ! Task: Compute splined second derivative of visibility function


  end subroutine initialize_rec_mod

  ! Saha equation for integration
  subroutine dX_edx(x, X_e, dydx)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: X_e
        real(dp), dimension(:), intent(out) :: dydx
        real(dp) :: T_b,n_b,phi2,alpha2,beta,beta2,n1s,lambda_alpha,C_r, Xe, a, H
        Xe = X_e(1)
        a  = exp(x)
        H  = get_H(x)
        T_b          = T_0/a
        n_b          = Omega_b*rho_c/(m_H*a**3)
        phi2         = 0.448d0*log(epsilon_0/(k_b*T_b))
        alpha2       = 64.d0*pi/sqrt(27.d0*pi)*(alpha/m_e)**2*sqrt(epsilon_0/(k_b*T_b))*phi2*hbar**2/c
        beta         = alpha2 *((m_e*k_b*T_b)/(2.d0*pi*hbar**2))**1.5*exp(-epsilon_0/(k_b*T_b))

        !This part is needed since the exponent
        !in beta2 becomes so large that the computer
        !sets it to infinity. However beta goes to zero before that
        !so it should be 0 even if the exponent is enormous.
        if(T_b <= 169.d0) then
            beta2    = 0.d0
        else
            beta2    = beta*exp((3.d0*epsilon_0)/(4.d0*k_b*T_b))
        end if

        n1s          = (1.d0-Xe)*n_b
        lambda_alpha = H*(3.d0*epsilon_0)**3/((8.d0*pi)**2*n1s) /(c*hbar)**3
        C_r          = (lambda_2s1s +lambda_alpha)/(lambda_2s1s+lambda_alpha+beta2)
        dydx         = C_r/H*(beta*(1.d0-Xe)- n_b*alpha2*Xe**2)

    end subroutine dX_edx

    subroutine dtaudx(x,tau, dydx)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: tau
        real(dp), dimension(:), intent(out) :: dydx
        real(dp)                            :: n_e
        real(dp)                            :: H
        n_e  = get_n_e(x)
        H    = get_H(x)
        dydx = -n_e*sigma_T/H*c
    end subroutine dtaudx


  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x,
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

  end function get_ddg


end module rec_mod
