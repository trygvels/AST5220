module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

  real(dp),    allocatable, dimension(:) :: uniform_a          ! uniform conformal time

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init

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

    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))

 	x_t(1) = x_start_rec 					! initial x
	a_t(1) = exp(x_t(1))					! initial a

	do i = 1,n_t-1							! filling arrays
		if (i < n1 + 1) then
			dx = (x_end_rec - x_start_rec)/n1
		else 
			dx = (x_0 - x_end_rec)/n2	
		end if
		x_t(i+1) = x_t(i) + dx 
		a_t(i+1) = exp(x_t(i+1))
	end do
	! LAST NUMBER WRONG
	write(*,*) x_t(1), x_t(500)
	write(*,*) a_t(1), a_t(500)

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))
	allocate(uniform_a(n_eta))    

	uniform_a(1) = a_init ! initial values
	x_eta(1) = x_eta1	  ! initial values
	do i = 1, n_eta-1
		uniform_a(i+1) = uniform_a(i) + (1-a_init)/(n_eta)
		x_eta(i+1) = log(uniform_a(i+1))
	end do
	
	! write(*,*) uniform_a(1), uniform_a(1000)
	! write(*,*) x_eta(1), x_eta(1000)

  end subroutine initialize_time_mod


  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
	get_H = H_0*sqrt((omega_b+omega_m)*exp**(-3*x)+(omega_r+omega_nu)*exp**(-4*x)+omega_lambda)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
	get_H_p = H_*sqrt((omega_b+omega_m)*exp**(-x)+(omega_r+omega_nu)*exp**(-2*x)+omega_lambda*exp(2*x))
  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
	get_dH_p = -0.5*H_0((omega_b+omega_m)*e**(-x)+2*(omega_r+omega_nu)*e**(-2*x)-2*omega_lambda*e**(2*x))/sqrt((omega_b+omega_m)*exp**(-x)+(omega_r+omega_nu)*exp**(-2*x)+omega_lambda*exp(2*x))
  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

  end function get_eta

end module time_mod
