module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

contains

subroutine initialize_time_mod
  implicit none

  integer(i4b) :: i, n, n1, n2
  real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, eta_init, a_init, step, eps, stepmin, yp1,ypn, rho_cc, rho_m, rho_b, rho_r, rho_lambda, z

  ! Define two epochs, 1) during and 2) after recombination.
  n1          = 200                       ! Number of grid points during recombination
  n2          = 300                       ! Number of grid points after recombination
  n_t         = n1 + n2                   ! Total number of grid points
  z_start_rec = 1630.4d0                  ! Redshift of start of recombination
  z_end_rec   = 614.2d0                   ! Redshift of end of recombination
  z_0         = 0.d0                      ! Redshift today
  x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
  x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
  x_0         = 0.d0                      ! x today
  
  n_eta       = 1000                      ! Number of eta grid points (for spline)
  a_init      = 1.d-10                    ! Start value of a for eta evaluation
  x_eta1      = log(a_init)               ! Start value of x for eta evaluation
  x_eta2      = 0.d0                      ! End value of x for eta evaluation

  yp1 = 1.d30
  ypn = 1.d30
  stepmin = 0
  eps =  1.d-10
  eta_init = a_init/(H_0*sqrt(Omega_r))

  ! Task: Fill in x and a grids
  allocate(x_t(n_t))
  allocate(a_t(n_t))

 	x_t(1) = x_start_rec 					! initial x
	a_t(1) = exp(x_t(1))					! initial a

	do i = 2,n_t							! filling arrays
		if (i < n1 + 1) then
			dx = (x_end_rec - x_start_rec)/(n1-1)
		else 
			dx = (x_0 - x_end_rec)/(n2-1)
		end if
		x_t(i) = x_t(i-1) + dx 
		a_t(i) = exp(x_t(i))
	end do
  

  ! Task: 1) Compute the conformal time at each eta time step
  !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
  allocate(x_eta(n_eta))
  allocate(eta(n_eta))
  allocate(eta2(n_eta))

  dx = (x_eta2-x_eta1)/(n_eta-1)
	x_eta(1) = x_eta1	  ! Uniformly spaced x-grid
	do i = 2, n_eta
		x_eta(i) = x_eta(i-1) + dx
	end do

! Integrating for eta        
	step =  abs(1.d-2*(x_eta(1)-x_eta(2)))      ! step length
  eta(1) = eta_init                          ! initial value of et
	do i = 2, n_eta
		eta(i) = eta(i-1)
		call odeint(eta(i:i), x_eta(i-1),x_eta(i), eps, step, stepmin, derivs, bsstep, output) 	
  end do
  
  ! Write to file - Eta, x_eta
  open(1, file="eta.dat", action="write")
  do i=1,n_eta
     write(1,"(3f8.3)",advance="no") eta(i), x_eta(i)
  end do
  close(1)
 

  ! Splining eta
  call spline(x_eta, eta,yp1,ypn,eta2)
  
  ! Spline + Interplolation write to file
  open (2,file="etasplint.dat",action="write")
  do i=1,n_t
     write (2,"(3f8.3)",advance="no") get_eta(x_t(i)), x_t(i)
  end do
  close(2)

  ! Calculating Omegas
  open(3, file="omegas.dat", action="write")
  do i=1, n_t
    rho_cc = 3*get_H(x_t(i))/(8*pi*G_grav)

    rho_m = Omega_m*rho_c*exp(x_t(i))**-3
    rho_b = Omega_b*rho_c*exp(x_t(i))**-3
    rho_r = Omega_r*rho_c*exp(x_t(i))**-4
    rho_lambda = Omega_lambda*rho_c
    write(3,"(3f8.3)",advance="no") rho_m/rho_cc, rho_m/rho_cc,rho_b/rho_cc,rho_lambda/rho_cc
  end do
  close(3)

  ! H values write - H(x), H(z)
  open(4, file="HxHz.dat", action="write")
  do i=1,n_t
    z = 1/x_t(i) - 1
    write(4,"(3f8.3)",advance="no") get_H(x_t(i)), z, get_H(z)
  end do
  close(4)  
  end subroutine initialize_time_mod

! dnu/dx=c/H_p
  subroutine derivs(x,eta, derivative)
    use healpix_types
    implicit none
    real(dp), intent(in)                :: x
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: derivative
    derivative = c/(get_H_p(x))
    end subroutine derivs  

 ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
	get_H = H_0*sqrt((Omega_b+Omega_m)*exp(-3*x)+(Omega_r+Omega_nu)*exp(-4*x)+Omega_lambda)
  end function get_H


  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
	get_H_p = H_0*sqrt((Omega_b+Omega_m)*exp(-x)+(Omega_r+Omega_nu)*exp(-2*x)+Omega_lambda*exp(2*x))
  end function get_H_p



  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
	get_dH_p = -0.5*H_0*((Omega_b+Omega_m)*exp(-x)+2*(Omega_r+Omega_nu)*exp(-2*x)-2*Omega_lambda*exp(2*x))/sqrt((Omega_b+Omega_m)*exp(-x)+(Omega_r+Omega_nu)*exp(-2*x)+Omega_lambda*exp(2*x))
  end function get_dH_p



  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    get_eta =  splint(x_eta, eta, eta2, x_in)
  end function get_eta


  subroutine output(x, y)
     implicit none
     real(dp),               intent(in)  :: x
     real(dp), dimension(:), intent(in)  :: y
  end subroutine output

end module time_mod
