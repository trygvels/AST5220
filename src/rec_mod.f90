module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b)                        :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:) :: x_rec             ! Grid
  real(dp), allocatable, dimension(:) :: X_e               ! Fractional electron density, n_e / n_H
  real(dp), allocatable, dimension(:) :: tau, tau2, tau22  ! tau, tau'', (tau'')''
  real(dp), allocatable, dimension(:) :: dtau              ! First derivative of tau: tau'
  real(dp), allocatable, dimension(:) :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:) :: g, g2, g22        ! Visibility: g, g'',(g'')''
  real(dp), allocatable, dimension(:) :: dg                ! First derivative of g: g'
  real(dp)                            :: yp1, ypn, eps, hmin
contains

  subroutine initialize_rec_mod
    implicit none

    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx
    real(dp)     :: f, n_e0, X_e0, xstart, xstop, step
    real(dp)     :: C_r

    logical(lgt) :: use_saha

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
    allocate(dtau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(dg(n))
    allocate(g2(n))
    allocate(g22(n))

    !---------------------- Time-grid ----------------------

    ! Uniform x-grid with 1000 points from a = 10^-10 to 1
    dx = (xstop-xstart)/(n-1)
    x_rec(1) = xstart
    do i=2,n
      x_rec(i) = x_rec(i-1) + dx
    end do

    step = abs(1.d-3*(x_rec(1)-x_rec(2))) ! Step length for ODE

    !---------------------- X_e calculation ----------------------

    ! Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1, n
       n_b = Omega_b*rho_c/(m_H*exp(x_rec(i))**3)
       if (use_saha) then
         ! Use the Saha equation
          T_b = T_0/exp(x_rec(i))
          X_e0= ((m_e*k_b*T_b)/(2.d0*pi*hbar**2))**1.5d0*exp(-epsilon_0/(k_b*T_b))/n_b
          X_e(i) = (-X_e0 + sqrt(X_e0**2 +4.d0*X_e0))/2.d0
      if (X_e(i) < saha_limit) use_saha = .false.
      else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, step, hmin, dX_edx, bsstep, output)
       end if
       n_e(i) = X_e(i)*n_b ! Electron density
    end do


    !---------------------- Electron density ----------------------
    !  Compute splined (log of) electron density function
    n_e = log(n_e)
    call spline(x_rec,n_e,yp1,ypn,n_e2)
    ! ---------------------- Optical depth ----------------------
    !  Compute optical depth at all grid points (Reverse integration)
    tau(n) = 0.d0
    do i=n-1,1,-1
      tau(i) = tau(i+1)
      call odeint(tau(i:i),x_rec(i+1),x_rec(i),eps,step,hmin,dtaudx,bsstep,output)
    end do

    ! Compute splined (log of) optical depth
    call spline(x_rec, tau, yp1, ypn,tau2)

    ! Compute splined second derivative of (log of) optical depth
    call spline(x_rec,tau2,yp1,ypn,tau22)

    ! Saving values of dtau
    do i = 1,n
      dtau(i) = get_dtau(x_rec(i))
    end do

    !---------------------- Visibility function ----------------------
    ! Computing visibility function gg
    do i=1,n
      g(i) = -get_dtau(x_rec(i))*exp(-tau(i))
    end do

    !  Compute splined visibility function
    call spline(x_rec,g,yp1,ypn,g2)

    !  Compute splined second derivative of visibility function
    call spline(x_rec,g2,yp1,ypn,g22)

    ! Saving values of dg
    do i = 1,n
      dg(i) = get_dg(x_rec(i))
    end do

    ! ---------------------- Generalization test  ----------------------
    do i = 1,n
      if ((exp(n_e(i)) - get_n_e(x_rec(i))) /= 0.d0) write(*,*) "n_e not zero at i = ", i
      if ((tau(i) - get_tau(x_rec(i))) /= 0.d0) write(*,*) "tau not zero at i = ", i
      if ((tau2(i) - get_ddtau(x_rec(i))) /= 0.d0) write(*,*) "tau2 not zero at i = ", i
      if ((g(i) - get_g(x_rec(i))) /= 0.d0) write(*,*) "g not zero at i = ", i
      if ((g2(i) - get_ddg(x_rec(i))) /= 0.d0) write(*,*) "g2 not zero at i = ", i
    end do

  end subroutine initialize_rec_mod

  ! ---------------------- Saha equation for integration ----------------------
  subroutine dX_edx(x, X_e, dydx)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: X_e
        real(dp), dimension(:), intent(out) :: dydx
        real(dp) :: T_b,n_b,phi2,alpha2,beta,beta2,n1s,lambda_alpha,C_r,H
        H      = get_H(x)
        T_b    = T_0/exp(x)
        n_b    = Omega_b*rho_c/(m_H*exp(x)**3)
        phi2   = 0.448d0*log(epsilon_0/(k_b*T_b))
        alpha2 = 64.d0*pi/sqrt(27.d0*pi)*(alpha/m_e)**2*sqrt(epsilon_0/(k_b*T_b))*phi2*hbar**2/c
        beta   = alpha2 *((m_e*k_b*T_b)/(2.d0*pi*hbar**2))**1.5*exp(-epsilon_0/(k_b*T_b))

        ! To avoid beta2 going to infinity, set it to 0
        if(T_b <= 169.d0) then
            beta2    = 0.d0
        else
            beta2    = beta*exp((3.d0*epsilon_0)/(4.d0*k_b*T_b))
        end if

        n1s          = (1.d0-X_e(1))*n_b
        lambda_alpha = H*(3.d0*epsilon_0)**3/((8.d0*pi)**2*n1s) /(c*hbar)**3
        C_r          = (lambda_2s1s +lambda_alpha)/(lambda_2s1s+lambda_alpha+beta2)
        dydx         = C_r/H*(beta*(1.d0-X_e(1)) - n_b*alpha2*X_e(1)**2)

    end subroutine dX_edx

    !---------------------- Optical thickness for ODE ----------------------
    subroutine dtaudx(x,tau, dydx)
        use healpix_types
        implicit none
        real(dp),               intent(in)  :: x
        real(dp), dimension(:), intent(in)  :: tau
        real(dp), dimension(:), intent(out) :: dydx
        real(dp)                            :: n_e
        real(dp)                            :: H

        dydx = -get_n_e(x)*sigma_T*exp(x)*c/get_H_p(x)
    end subroutine dtaudx

  !---------------------- Functions for generalization ----------------------
  ! Complete routine for computing n_e at arbitrary x, using precomputed information
  function get_n_e(x_in)
      implicit none
      real(dp), intent(in) :: x_in
      real(dp)             :: get_n_e
      !Spline integration with precalculated logarithmic values
      ! TODO: Should n_e be log'ed?
      get_n_e = exp(splint(x_rec, n_e, n_e2, x_in))
  end function get_n_e

  !  Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

    get_tau = splint(x_rec,tau,tau2,x)
  end function get_tau

  function get_dtau(x)
       implicit none
       real(dp), intent(in) :: x
       real(dp)             :: get_dtau
       real(dp)             :: n_e,H_p
       !H_p = get_H_p(x)
       !n_e = get_n_e(x)
       !get_dtau = -n_e*sigma_T*exp(x)*c/H_p
       get_dtau =  splint_deriv(x_rec, tau, tau2, x)
   end function get_dtau

  ! : Complete routine for computing the second derivative of tau at arbitrary x,
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    get_ddtau = splint(x_rec,tau2, tau22, x)
  end function get_ddtau

  !  Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g
    get_g = splint(x_rec, g, g2, x)
  end function get_g

  !  Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg
    get_dg = splint_deriv(x_rec,g,g2,x)
  end function get_dg

  !  Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
    get_ddg = splint(x_rec,g2,g22,x)
  end function get_ddg


end module rec_mod
