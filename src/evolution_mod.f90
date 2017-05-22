module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
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
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int

  ! Effectivisation variables
  real(dp),     private :: ck, ckH_p, dt, a !dt = dtau

  !Hires source function variables
  integer(i4b), parameter             :: k_num = 5000
  integer(i4b), parameter             :: x_num = 5000

contains


  subroutine get_hires_source_function(x_hires, k_hires, S)
    implicit none

    real(dp), allocatable, dimension(:),   intent(out) :: x_hires, k_hires
    real(dp), allocatable, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, k
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi
    real(dp), allocatable, dimension(:,:)     :: S_lores
    real(dp), allocatable, dimension(:,:,:,:) :: S_coeff

    ! Output a pre-computed 2D array (over k and x) for the source function, S(k,x).
    allocate(x_hires(x_num))
    allocate(k_hires(k_num))

    ! Generate hires grid for x and k
    do i = 1, x_num !x_num and k_num the same length
      x_hires(i) = x_init - x_init*(i-1.d0)/(x_num-1.d0)
      k_hires(i)= k_min  + (k_max - k_min)*((i-1.d0)/(k_num-1.d0))**2
    end do

    allocate(S_lores(n_t,n_k))  ! Lores source function
    allocate(S(x_num,k_num))    ! Hires source function
    allocate(S_coeff(4,4,n_t,n_k))

    do i=1,n_t
      g     = get_g(x_t(i))
      dg    = get_dg(x_t(i))
      ddg   = get_ddg(x_t(i))
      tau   = get_tau(x_t(i))
      dt    = get_dtau(x_t(i))
      ddt   = get_ddtau(x_t(i))
      H_p   = get_H_p(x_t(i))
      dH_p  = get_dH_p(x_t(i))

      do k=1,n_k
        ck = c*ks(k)
        
        Pi    = Theta(i,2,k)
        dPi   = dTheta(i,2,k)

        ! Source function with low resolution (Note preduct rule on derivatives)
        ddPi  = 2.d0*ck/(5.d0*H_p)*(-dH_p/H_p*Theta(i,1,k) + dTheta(i,1,k)) &
                +0.3d0*(ddt*Pi+dt*dPi) &
                -3.d0*ck/(5.d0*H_p)*(-dH_p/H_p*Theta(i,3,k) + dTheta(i,3,k))

        S_lores(i,k) = g*(Theta(i,0,k) +Psi(i,k) + .25d0*Pi) &
                       +exp(-tau)*(dPsi(i,k)-dPhi(i,k)) &
                       -1.d0/ck*(H_p*(g*dv_b(i,k) + v_b(i,k)*dg) + g*v_b(i,k)*dH_p) &
                       +.75d0/ck**2*((H_0**2/2.d0*((Omega_m+Omega_b)/exp(x_t(i)) &
                       +4.d0*Omega_r/exp(2.d0*x_t(i)) +4.d0*Omega_lambda*exp(2.d0*x_t(i))))*&
                       g*Pi +3.d0*H_p*dH_p*(dg*Pi+g*dPi)+H_p**2* &
                       (ddg*Pi +2.d0*dg*dPi+g*ddPi))
      end do
    end do

    !Spline results
    call splie2_full_precomp(x_t,ks,S_lores,S_coeff)

    ! Make hires Source function 5kx5k grid
    do k=1, k_num
        do i=1, x_num
            S(i,k) = splin2_full_precomp(x_t, ks, S_coeff, x_hires(i), k_hires(k))
        end do
    end do
  end subroutine get_hires_source_function

  !######### MILESTONE 3 #########

  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, k, i_tc

    ! Allocate arrays for perturbation quantities
    allocate(Theta(n_t, 0:lmax_int, n_k))
    allocate(delta(n_t, n_k))
    allocate(delta_b(n_t, n_k))
    allocate(v(n_t, n_k))
    allocate(v_b(n_t, n_k))
    allocate(Phi(n_t, n_k))
    allocate(Psi(n_t, n_k))
    allocate(dPhi(n_t, n_k))
    allocate(dPsi(n_t, n_k))
    allocate(dv_b(n_t, n_k))
    allocate(dTheta(n_t, 0:lmax_int, n_k))

    ! Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1,n_k
        ks(k) = k_min +(k_max - k_min)*((k-1.d0)/(n_k-1.d0))**2
    end do

    ! Initial contitions for B-E-equations
    Phi(1,:)     = 1d0
    delta(1,:)   = 1.5d0*Phi(1,:)
    delta_b(1,:) = delta(1,:)
    Theta(1,0,:) = 0.5d0*Phi(1,:)
    do k = 1, n_k

      ! Effectivisation
       ckH_p        = c*ks(k)/get_H_p(x_init)
       dt           = get_dtau(x_init)

       v(1,k)       = ckH_p/(2.d0)*Phi(1,k)
       v_b(1,k)     = v(1,k)
       Theta(1,1,k) = -ckH_p/(6.d0)*Phi(1,k)
       Theta(1,2,k) = -20.d0*ckH_p/(45.d0*dt)*Theta(1,1,k)
       do l = 3, lmax_int
          Theta(1,l,k) = -l*ckH_P*Theta(1,l-1,k)/((2.d0*l+1)*dt)
       end do
       Psi(1,k)     = -Phi(1,k) - 12.d0*H_0**2/(ks(k)*c*a_t(1))**2*Omega_r*Theta(1,2,k)
    end do

  end subroutine initialize_perturbation_eqns


  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, k, l, i_tc
    real(dp)     :: eps, hmin, h1, x_tc, dt

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
       y_tight_coupling(1) = delta(1,k)
       y_tight_coupling(2) = delta_b(1,k)
       y_tight_coupling(3) = v(1,k)
       y_tight_coupling(4) = v_b(1,k)
       y_tight_coupling(5) = Phi(1,k)
       y_tight_coupling(6) = Theta(1,0,k)
       y_tight_coupling(7) = Theta(1,1,k)

       ! Find the time to which tight coupling is assumed
       x_tc = get_tight_coupling_time(k_current)

       ! Integrate from in tight coupling regime
       i_tc = 2
       do while (x_t(i_tc)<x_tc)
           !Integrate with tight coupling
           call odeint(y_tight_coupling, x_t(i_tc-1),x_t(i_tc), eps, h1, hmin, dytc, bsstep, output)

           ! Effectivisation
           ckH_p        = ck*get_H_p(x_t(i_tc))
           dt           = get_dtau(x_t(i_tc))

           !Save variables one value at a time
           delta(i_tc,k)   = y_tight_coupling(1)
           delta_b(i_tc,k) = y_tight_coupling(2)
           v(i_tc,k)       = y_tight_coupling(3)
           v_b(i_tc,k)     = y_tight_coupling(4)
           Phi(i_tc,k)     = y_tight_coupling(5)
           Theta(i_tc,0,k) = y_tight_coupling(6)
           Theta(i_tc,1,k) = y_tight_coupling(7)
           Theta(i_tc,2,k) = -20.d0*ckH_p/45.d0/dt*Theta(i_tc,1,k)
           do l = 3, lmax_int
              Theta(i_tc,l,k) = -l/(2.d0*l+1.d0)*ckH_p/dt*Theta(i_tc,l-1,k)
           end do
           Psi(i_tc,k)      = -Phi(i_tc,k) - 12.d0*H_0**2.d0/(ck*exp(x_t(i_tc)))**2.d0*Omega_r*Theta(i_tc,2,k)

           ! Storing derivatives
           call dytc(x_t(i_tc),y_tight_coupling,dydx) !Call subroutine which calculates dydx array
           dPhi(i_tc,k)     = dydx(4)
           dv_b(i_tc,k)     = dydx(5)
           dTheta(i_tc,:,k) = dydx(6)
           dPsi(i_tc,k)     = dydx(7)
           dTheta(i_tc,2,k) = 2.d0/5.d0*ckH_p*Theta(i_tc,1,k) - 3.d0/5.d0*ckH_p*Theta(i_tc,3,k)+dt*0.9d0*Theta(i_tc,2,k)

           do l=3,lmax_int-1
           dTheta(i_tc,l,k) = l/(2.d0*l+1.d0)*ckH_p*dTheta(i_tc,l-1,k) - &
                       (l+1.d0)/(2.d0*l+1.d0)*ckH_p*dTheta(i_tc,l+1,k) + dt*Theta(i_tc,l,k)
           end do

           dPsi(i_tc,k)     = -dPhi(i_tc,k) - 12.d0*H_0**2.d0/(ck*exp(x_t(i_tc)))**2.d0 *Omega_r*(-2.d0*Theta(i_tc,2,k)+dTheta(i_tc,2,k))

           i_tc = i_tc + 1 !iteration
       end do ! ends while

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
           Psi(i,k)     = - Phi(i,k) - (12.d0*H_0**2)/(ck*exp(x_t(i)))**2.d0*Omega_r*Theta(i,2,k)

            ! Store derivatives that are required for C_l estimation
            call dy(x_t(i),y,dydx) !Call derivatives subroutine
            dv_b(i,k)     = dydx(4)
            dPhi(i,k)     = dydx(5)
            do l=0,lmax_int
                dTheta(i,l,k) = dydx(6+l)
            end do
            dPsi(i,k)     = -dPhi(i,k) - 12.d0*H_0**2.d0/(ck*a_t(i))**2.d0*Omega_r*(-2.d0*Theta(i,2,k)+dTheta(i,2,k))
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
          real(dp) :: dt,a,H_p,dH_p,ckH_p

          ! Save variables for derivation
          delta   = y_tc(1)
          delta_b = y_tc(2)
          v       = y_tc(3)
          v_b     = y_tc(4)
          Phi     = y_tc(5)
          Theta0  = y_tc(6)
          Theta1  = y_tc(7)

          ! Variables for effectivisation
          dt    = get_dtau(x)
          a     = exp(x)
          H_p   = get_H_p(x)
          dH_p  = get_dH_p(x)
          ckH_p = ck/H_p

          ! Calculate derivatives
          Theta2    = -20.d0*ckH_p/(45.d0*dt)*Theta1
          R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
          Psi       = -Phi - 12.d0*(H_0/(ck*a))**2.d0*Omega_r*Theta2
          dPhi      = Psi -(ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p**2.d0)*(Omega_m/a*delta+Omega_b/a*delta_b+4.d0*Omega_r*Theta0/a**2.d0)
          dTheta0   = -ckH_p*Theta1 - dPhi
          d_delta   = ckH_p*v   - 3.d0*dPhi
          d_delta_b = ckH_p*v_b - 3.d0*dPhi
          d_v       = -v -ckH_p*Psi
          q         = ( -((1.d0-2.d0*R)*dt + (1.d0+R)*get_ddtau(x)) *&
                      (3.d0*Theta1+v_b) - ckH_p*Psi +(1.d0-dH_p/H_p)*&
                      ckH_p*(-Theta0 + 2.d0*Theta2) - ckH_p*dTheta0) / &
                      ((1.d0+R)*dt+dH_p/H_p -1.d0)
          dv_b      = (1.d0/(1.d0+R)) *(-v_b - ckH_p*Psi + &
                      R*(q+ckH_p*(-Theta0 + 2.d0*Theta2)-ckH_p*Psi))
          dTheta1   = (1.d0/3.d0)*(q-dv_b)

          ! Store derivatives
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

     integer(i4b) :: l
     real(dp) :: d_delta, d_delta_b, d_v, R
     real(dp) :: delta, delta_b, v, v_b, Phi
     real(dp) :: Theta0, Theta1, Theta2, Theta3, Theta4, Theta5, Theta6
     real(dp) :: Psi, dPhi, dTheta0, dv_b, dTheta1, dTheta2
     real(dp) :: a,H_p,ckH_p,dt

     ! Save variables for derivation
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

     ! Effectivisation variables
     a = exp(x)
     H_p = get_H_p(x)
     ckH_p = ck/H_p
     dt = get_dtau(x)

     ! Derivation
     R         = (4.d0*Omega_r)/(3.d0*Omega_b*a)
     Psi       = - Phi - 12.d0*(H_0/(ck*a))**2.d0*Omega_r*Theta2
     dPhi      = Psi - (ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p**2.d0)*(Omega_m/a*delta+Omega_b/a*delta_b+4.d0*Omega_r*Theta0/a**2.d0)
     dTheta0   = - ckH_p*Theta1 - dPhi
     d_delta   = ckH_p*v   - 3.d0*dPhi
     d_delta_b = ckH_p*v_b - 3.d0*dPhi
     d_v       = - v - ckH_p*Psi
     dv_b      = - v_b - ckH_p*Psi + dt*R*(3.d0*Theta1+v_b)
     dTheta1   = ckH_p/3.d0*Theta0 - 2.d0/3.d0*ckH_p*Theta2 + ckH_p/3.d0*Psi + dt*(Theta1+v_b/3.d0)
     !dTheta2   = 2.d0/5.d0*ckH_p*Theta1 - 3.d0/5.d0*ckH_p*Theta3+dt*0.9d0*Theta2


     do l = 2, lmax_int-1
       dydx(6+l) = l*ckH_p/(2.d0*l+1.d0)*y(6+l-1) - (l+1.d0)*ckH_p/(2.d0*l+1.d0)*y(6+l+1) + dt*(y(6+l) - 1.d0/10.d0*y(6+l)*abs(l==2))
     end do

     dydx(1) = d_delta
     dydx(2) = d_delta_b
     dydx(3) = d_v
     dydx(4) = dv_b
     dydx(5) = dPhi
     dydx(6) = dTheta0
     dydx(7) = dTheta1
     dydx(8) = dTheta2
     dydx(6+l) = ckH_p*y(6+l-1) - c*(l+1.d0)/(H_p*get_eta(x))*y(6+l) + dt*y(6+l)
     !dydx(6+lmax_int) = ckH_p*y(6+lmax_int-1) - c*(lmax_int+1.d0)/H_p/get_eta(x)*y(6+lmax_int) + dt*y(6+lmax_int)


  end subroutine dy

  ! Time at which tight coupling ends.
  ! dt < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)

  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    integer(i4b)          :: i, n_test
    real(dp)              :: get_tight_coupling_time, x, x_start_rec

    n_test = 1000
    x_start_rec = -log(1.d0+1630.4d0)
    x = x_init
    do i = 1, n_test
       x = x + (x_start_rec-x_init)/(n_test-1)
       if ((abs(c*k/(get_H_p(x)*get_dtau(x))) >= 0.1d0 .and. abs(get_dtau(x))<=10.d0) .or. ( x >= x_start_rec)) then
          ! (abs(get_dtau(x)) <= 10.d0 .and.
          get_tight_coupling_time = x
          exit
       end if
    end do
  end function get_tight_coupling_time

end module evolution_mod
