module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none
  real(dp),     allocatable, dimension(:,:)     :: j_l, j_l2
  real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
  real(dp),     allocatable, dimension(:)       :: x_hires, k_hires,  l_hires, cl_hires, x_lores

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none
    CHARACTER(*), PARAMETER :: fileplace = "/uio/hume/student-u68/trygvels/AST5220/src/data/"
    integer(i4b) :: i, k, l, x_num, l_num, k_num,  n_spline, ilo, method, savetrans
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable,     dimension(:)       :: ls
    real(dp),     allocatable,     dimension(:)       :: integrandx, integrandk
    real(dp),     allocatable,     dimension(:)       :: cls, cls2, ls_dp
    real(dp),     allocatable,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable,     dimension(:,:)     :: S, S2
    real(dp),     allocatable,     dimension(:,:)     :: Theta_l, Theta2_l
    real :: start, finish
    real(dp) :: integralx, integralk, h1, h2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Hires array numbers
    x_num = 5000
    k_num = 5000

    ! Allocate Source function and Hires grids
    allocate(S(x_num,k_num))
    allocate(x_hires(x_num))
    allocate(k_hires(k_num))

    ! Calculate Hires source function from evolution_mod
    call get_hires_source_function(x_hires,k_hires,S)
    write(*,*) k_hires(1)*H0/c,k_hires(1000)*H0/c,k_hires(2000)*H0/c,k_hires(3000)*H0/c,k_hires(4000)*H0/c
    n_spline = 5400
    allocate(z_spline(n_spline))    !j_l(z), not redshift
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    ! Calculate bessel functions, with 5400 sampled points between z = 0 to 3500
    do i = 1, n_spline
      z_spline(i) = (i-1.d0)*3500.d0/(n_spline-1.d0)
      do l=1, l_num
        if (z_spline(i)>0.d0) then
          call sphbes(ls(l),z_spline(i),j_l(i,l))
        endif
      end do
    end do

    !Spline bessel functions, get second derivative for later splint
    do l=1,l_num
          call spline(z_spline, j_l(:,l), 1.d30, 1.d30, j_l2(:,l))
    end do


    allocate(Theta_l(l_num,k_num))
    allocate(integrandx(x_num))
    allocate(integrandk(k_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(x_lores(x_num/10))

    ! #### C_l COMPUTATION OVER l's ####
    ! Method 1 = fast, method 2 = slow
    method = 2

    do l = 1, l_num
      if (method == 1) then
        !###############################################
        !################### METHOD 1 ##################
        !###############################################

        integralk = 0 ! Reset integral for each value of cls
        do k = 1, k_num
          integralx = 0 ! Reset integral for each value of Theta_l
          ! Integrate Theta_l
          do i = 1, x_num/10
            ilo = 1 + (i-1)*(x_num-1)/(x_num/10-1) !Speed up integration
            x_lores(i) = x_hires(ilo)
            integrandx(i) = S(ilo,k)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(k)*(get_eta(0.d0)-get_eta(x_hires(ilo))))
            integralx = integralx + integrandx(i)
          end do
          ! Subtract half of first and last integrand for x
          h1 = (x_lores(x_num/10) - x_lores(1))/(x_num/10)
          Theta_l(l,k) = h1*(integralx - 0.5d0*(integrandx(1)+integrandx(x_num/10)))

          ! Integrate C_l
          integrandk(k) = (c*k_hires(k)/H_0)**(n_s-1.d0)*Theta_l(l,k)**2/k_hires(k)
          integralk = integralk + integrandk(k)
        end do

        ! Subtract half of first and last integrand for k
        h2 = (k_hires(k_num) - k_hires(1))/k_num
        integralk = h2*(integralk - 0.5d0*(integrandk(1)+integrandk(k_num)))

        ! Store C_l in an array. Optionally output to file
        cls(l) = integralk*ls(l)*(ls(l)+1.d0)/(2.d0*pi)
      else if (method == 2) then
        !###############################################
        !################### METHOD 2 ##################
        !###############################################

        integralk = 0 ! Reset integral for each value of cls
        do k = 1, k_num
          integralx = 0 ! Reset integral for each value of Theta_l
          ! Integrate Theta_l
          do i = 1, x_num
            integrandx(i) = S(i,k)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(k)*(get_eta(0.d0)-get_eta(x_hires(i))))
            integralx = integralx + integrandx(i)
          end do
          ! Subtract half of first and last integrand for x
          h1 = (x_hires(x_num) - x_hires(1))/(x_num)
          Theta_l(l,k) = h1*(integralx - 0.5d0*(integrandx(1)+integrandx(x_num)))

          ! Integrate C_l
          integrandk(k) = (c*k_hires(k)/H_0)**(n_s-1.d0)*Theta_l(l,k)**2/k_hires(k)
          integralk = integralk + integrandk(k)
        end do

        ! Subtract half of first and last integrand for k
        h2 = (k_hires(k_num) - k_hires(1))/k_num
        integralk = h2*(integralk - 0.5d0*(integrandk(1)+integrandk(k_num)))

        ! Store C_l in an array. Optionally output to file
        cls(l) = integralk*ls(l)*(ls(l)+1.d0)/(2.d0*pi)
      end if


     ! Timer for loop
       call cpu_time(finish)
       write(*,*) "l = ", l
       print '("Time = ",f7.2," seconds.")',finish-start
    end do


    ! Convert ls to double precision
    allocate(ls_dp(l_num))

    do l=1,l_num
        ls_dp(l) = ls(l)
    end do

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    allocate(cl_hires(int(maxval(ls))))
    allocate(l_hires(int(maxval(ls))))

    ! Spline cl, get second derivative for splint
    call spline(ls_dp, cls, yp1, ypn, cls2)

    !Spline to get cl at 1200 l's
    do l = 1, ls(l_num)
      l_hires(l) = ls_dp(1) + (l-1.d0)*(ls_dp(l_num)-ls_dp(1))/(ls_dp(l_num)-2.d0) ! dp hires l
      cl_hires(l) =  splint(ls_dp, cls, cls2, l_hires(l))
    end do

    ! Transfer functions
    allocate(Theta2_l(l_num,5))

    ! Splining transfer function
    savetrans = 0
    if (savetrans = 1) then

      call spline(ls_dp, Theta_l(:,1),yp1,ypn,Theta2_l(:,1))
      call spline(ls_dp, Theta_l(:,1000),yp1,ypn,Theta2_l(:,2))
      call spline(ls_dp, Theta_l(:,2000),yp1,ypn,Theta2_l(:,3))
      call spline(ls_dp, Theta_l(:,3000),yp1,ypn,Theta2_l(:,4))
      call spline(ls_dp, Theta_l(:,4000),yp1,ypn,Theta2_l(:,5))

      open (unit=31, file=fileplace//"transfer1.dat", action="write", status="replace")
      open (unit=32, file=fileplace//"transfer2.dat", action="write", status="replace")
      open (unit=33, file=fileplace//"transfer3.dat", action="write", status="replace")
      open (unit=34, file=fileplace//"transfer4.dat", action="write", status="replace")
      open (unit=35, file=fileplace//"transfer5.dat", action="write", status="replace")

      ! Write splint transfer functions
      write (31,'(*(2X, ES14.6E3))') c*k_hires(1)/H0
      write (32,'(*(2X, ES14.6E3))') c*k_hires(1000)/H0
      write (33,'(*(2X, ES14.6E3))') c*k_hires(2000)/H0
      write (34,'(*(2X, ES14.6E3))') c*k_hires(3000)/H0
      write (35,'(*(2X, ES14.6E3))') c*k_hires(4000)/H0
      do l = 1, 1200
        write (31,'(*(2X, ES14.6E3))') splint(ls_dp, Theta_l(:,1), Theta2_l(:,1), l_hires(l))
        write (32,'(*(2X, ES14.6E3))') splint(ls_dp, Theta_l(:,1000), Theta2_l(:,2), l_hires(l))
        write (33,'(*(2X, ES14.6E3))') splint(ls_dp, Theta_l(:,2000), Theta2_l(:,3), l_hires(l))
        write (34,'(*(2X, ES14.6E3))') splint(ls_dp, Theta_l(:,3000), Theta2_l(:,4), l_hires(l))
        write (35,'(*(2X, ES14.6E3))') splint(ls_dp, Theta_l(:,4000), Theta2_l(:,5), l_hires(l))
      end do
      close(31)
      close(32)
      close(33)
      close(34)
      close(35)
    endif
  end subroutine compute_cls

end module cl_mod
