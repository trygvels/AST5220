module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-10
  real(dp),     parameter, private :: x_init   = log(a_init)
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current, ck
  integer(i4b), private :: npar = 6+lmax_int

contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:) :: S_lores

    ! TASK: Output a pre-computed 2D array (over k and x) for the
    !       source function, S(k,x). Remember to set up (and allocate) output
    !       k and x arrays too.
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays

  end subroutine get_hires_source_function


  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, k, i_tc
    ! DONE: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1,n_k
        ks(k) = k_min +(k_max - k_min)*((k-1.d0)/(n_k-1.d0))**2
    end do

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))

    ! Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1d0
    delta(0,:)   = 1.5d0*Phi(0,:)
    delta_b(0,:) = delta(0,:)
    Theta(0,0,:) = 0.5d0*Phi(0,:) !Sped up when not in loop
    do k = 1, n_k
       v(0,k)       = c*k/(2.d0*get_H_p(x_init))*Phi(0,k)
       v_b(0,k)     = v(0,k)
       Theta(0,1,k) = -c*k/(6.d0*get_H_p(x_init))*Phi(0,k)
       Theta(0,2,k) = -20.d0*c*k/(45.d0*get_H_p(x_init)*get_dtau(x_init))*Theta(0,1,k)
       do l = 3, lmax_int
          Theta(0,l,k) = -l*c*k*Theta(0,l-1,k)/((2*l+1)*get_H_p(x_init)*get_dtau(x_init))
       end do
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, k, l, i_tc
    real(dp)     :: x1, x2
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx

    h1     = 1.d-5
    eps    = 1.d-8
    hmin   = 0.d0
    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    ! Propagate each k-mode independently
    do k = 1, n_k

       k_current = ks(k)  ! Store k_current as a global module variable
       ck = k_current*c   ! One calculation for each iteration

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)

       ! Find the time to which tight coupling is assumed,
       ! and integrate equations to that time
       ! Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       x_tc = get_tight_coupling_time(k_current) !x value at start of tc

       !################# 1 Grid to rule them all
       do i=1,n_t !700 points
         if (x_t(i)<x_tc) then
           !Integrate with tight coupling
           call odeint(y_tight_coupling, x_t(i),x_t(i+1), eps,h1,hmin,dytc, bsstep, output)
           !Trenger vi egentlig bare ett punkt? Siden verdien er lik under TC??
          !Save variables one value at a time
           delta(i,k)   = y_tight_coupling(1)
           delta_b(i,k) = y_tight_coupling(2)
           v(i,k)       = y_tight_coupling(3)
           v_b(i,k)     = y_tight_coupling(4)
           Phi(i,k)     = y_tight_coupling(5)
           Theta(i,0,k) = y_tight_coupling(6)
           Theta(i,1,k) = y_tight_coupling(7)
           Theta(i,2,k) = -20.d0*ck/45.d0/get_H_p(x_t(i))/get_dtau(x_t(i))*Theta(i,1,k)
           do l = 3, lmax_int
              Theta(i,l,k) = -l/(2.d0*l+1.d0)*ck/get_H_p(x_t(i))/get_dtau(x_t(i))*Theta(i,l-1,k)
           end do

           ! STORE DERIVATIVES
           call dytc(x_t(i),y_tight_coupling,dydx) !Call subroutine which calculates dydx array
           dPhi(i,k)     = dydx(4)
           dv_b(i,k)     = dydx(5)
           dTheta(i,:,k) = dydx(6)
           dPsi(i,k)     = dydx(7)
           dTheta(i,2,k) = 2.d0/5.d0*ck/get_H_p(x_t(i))*Theta(i,1,k) - 3.d0/5.d0*ck/get_H_p(x_t(i))*Theta(i,3,k)+get_dtau(x_t(i))*0.9d0*Theta(i,2,k)

           do l=3,lmax_int-1
           dTheta(i,l,k) = l/(2.d0*l+1.d0)*ck/get_H_p(x_t(i))*dTheta(i,l-1,k) - &
                       (l+1.d0)/(2.d0*l+1.d0)*ck/get_H_p(x_t(i))*dTheta(i,l+1,k) +get_dtau(x_t(i))*Theta(i,l,k)
           end do
         else
           i_tc = i
           exit
         end if
       end do
       !Save initital conditions for next integration
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(i_tc-1,2,k) !Save last index of tight coupling
       do l = 3, lmax_int
          y(6+l) = Theta(i_tc-1,l,k)
       end do

       !#### DO INTEGRATION FOR AFTER TC ####
      do i = i_tc, n_t
           ! Integrate equations from tight coupling to today
           Call odeint(y,x_t(i-1),x_t(i), eps,h1,hmin,dy, bsstep, output)

           ! Store variables at time step i in gloabl variables
           delta(i,k)   = y(1)
           delta_b(i,k) = y(2)
           v(i,k)       = y(3)
           v_b(i,k)     = y(4)
           Phi(i,k)     = y(5)
           do l = 0, lmax_int
              Theta(i,l,k) = y(6+l)
           end do
           Psi(i,k)     = - Phi(i,k) - (12.d0*H_0**2.d0)/(ck*a_t(i))**2.d0*Omega_r*Theta(i,2,k)

            ! Store derivatives that are required for C_l estimation
            call dy(x_t(i),y,dydx)
            dv_b(i,k)     = dydx(4)
            dPhi(i,k)     = dydx(5)
            do l=0,lmax_int
                dTheta(i,l,k) = dydx(6+l)
            end do
            dPsi(i,k)     = -dPhi(i,k) - 12.d0*H_0**2.d0/(ck*a_t(i))**2.d0*&
                             Omega_r*(-2.d0*Theta(i,2,k)+dTheta(i,2,k))
       end do ! Ends i
    end do ! Ends k

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns

  ! Derivatives for integration of y_tight_coupling
  subroutine dytc(x,y_tc, dydx)
    use healpix_types
          implicit none
          real(dp),               intent(in)  :: x
          real(dp), dimension(:), intent(in)  :: y_tc
          real(dp), dimension(:), intent(out) :: dydx

          real(dp) :: d_delta
          real(dp) :: d_delta_b
          real(dp) :: d_v
          real(dp) :: q,R

          real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2
          real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1
          real(dp) :: dtau,ddtau,a,H_p,dH_p,ckH_p


          delta   = y_tc(1)
          delta_b = y_tc(2)
          v       = y_tc(3)
          v_b     = y_tc(4)
          Phi     = y_tc(5)
          Theta0  = y_tc(6)
          Theta1  = y_tc(7)

          dtau  = get_dtau(x)
          ddtau = get_ddtau(x)
          a     = exp(x)
          H_p   = get_H_p(x)
          dH_p  = get_dH_p(x)
          ckH_p = ck/H_p

          Theta2    = -20.d0*ckH_p/(45.d0*dtau)*Theta1
          R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
          Psi       = -Phi - 12.d0*(H_0/(ck*a))**2.d0*Omega_r*Theta2
          dPhi      = Psi -(ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p**2.d0)*(Omega_m/a*delta+Omega_b/a*delta_b+4.d0*Omega_r*Theta0/a**2.d0)
          dTheta0   = -ckH_p*Theta1 - dPhi
          d_delta   = ckH_p*v   - 3.d0*dPhi
          d_delta_b = ckH_p*v_b - 3.d0*dPhi
          d_v       = -v -ckH_p*Psi
          q         = ( -((1.d0-2.d0*R)*dtau + (1.d0+R)*ddtau) *&
                      (3.d0*Theta1+v_b) - ckH_p*Psi +(1.d0-dH_p/H_p)*&
                      ckH_p*(-Theta0 + 2.d0*Theta2) - ckH_p*dTheta0) / &
                      ((1.d0+R)*dtau+dH_p/H_p -1.d0)
          dv_b      = (1.d0/(1.d0+R)) *(-v_b - ckH_p*Psi + &
                      R*(q+ckH_p*(-Theta0 + 2.d0*Theta2)-ckH_p*Psi))
          dTheta1   = (1.d0/3.d0)*(q-dv_b)

          dydx(1) = d_delta
          dydx(2) = d_delta_b
          dydx(3) = d_v
          dydx(4) = dv_b
          dydx(5) = dPhi
          dydx(6) = dTheta0
          dydx(7) = dTheta1
  end subroutine dytc
  ! Derivatives for integration of y
  subroutine dy(x,y, dydx)
    use healpix_types
     implicit none
     real(dp),               intent(in)  :: x
     real(dp), dimension(:), intent(in)  :: y
     real(dp), dimension(:), intent(out) :: dydx

     real(dp) :: d_delta
     real(dp) :: d_delta_b
     real(dp) :: d_v
     real(dp) :: R
     integer(i4b) :: l
     real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2,Theta3,Theta4,Theta5,Theta6
     real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1,dTheta2
     real(dp) :: a,H_p,ckH_p,dtau

     delta   = y(1)
     delta_b = y(2)
     v       = y(3)
     v_b     = y(4)
     Phi     = y(5)
     Theta0  = y(6)
     Theta1  = y(7)
     Theta2  = y(8)
     Theta3  = y(9)
     Theta4  = y(10)
     Theta5  = y(11)
     Theta6  = y(12)

     a = exp(x)
     H_p = get_H_p(x)
     ckH_p = ck/H_p
     dtau = get_dtau(x)


     R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
     Psi       = -Phi - 12.d0*(H_0/ck/a)**2.d0*Omega_r*Theta2

     dPhi      = Psi -(ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p**2.d0)*(Omega_m/a*delta+Omega_b/a*delta_b+4.d0*Omega_r*Theta0/a**2.d0)

     dTheta0   = -ckH_p*Theta1 - dPhi
     d_delta   = ckH_p*v   - 3.d0*dPhi
     d_delta_b = ckH_p*v_b - 3.d0*dPhi
     d_v       = -v -ckH_p*Psi

     dv_b      = -v_b -ckH_p*Psi +dtau*R*(3.d0*Theta1+v_b)
     dTheta1   = ckH_p/3.d0*Theta0 -2.d0/3.d0*ckH_p*Theta2 + &
                 ckH_p/3.d0*Psi +dtau*(Theta1+v_b/3.d0)

     dTheta2   = 2.d0/5.d0*ckH_p*Theta1 - 3.d0/5.d0*ckH_p*Theta3+dtau*0.9d0*Theta2
     do l=3,lmax_int-1
         dydx(6+l) = l/(2.d0*l+1.d0)*ckH_p*y(5+l) - &
                     (l+1.d0)/(2.d0*l+1.d0)*ckH_p*y(7+l) +dtau*y(6+l)
     end do

     dydx(6+lmax_int) = ckH_p*y(6+lmax_int-1) -c*(lmax_int+1.d0)/H_p/get_eta(x)*y(6+lmax_int) +dtau*y(6+lmax_int)

     dydx(1) = d_delta
     dydx(2) = d_delta_b
     dydx(3) = d_v
     dydx(4) = dv_b
     dydx(5) = dPhi
     dydx(6) = dTheta0
     dydx(7) = dTheta1
     dydx(8) = dTheta2

  end subroutine dy


  !       Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none
    integer(i4b)          :: i,n
    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time, x, x_start_rec
    x_start_rec = -log(1.d0 + 1630.4d0 )  ! x of start of recombination
    n = 1d4
    do i = 0,n
      x = x_init - i*x_init/n
      if (x < x_start_rec .and. abs(c*k/(get_H_p(x)*get_dtau(x))) <= 0.1d0 .and. abs(get_dtau(x)) > 10.d0) then
        get_tight_coupling_time = x
      end if
    end do
  end function get_tight_coupling_time

end module evolution_mod
