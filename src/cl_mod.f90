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
    integer(i4b) :: i, k, l, x_num, l_num, k_num,  n_spline, ilo
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable,     dimension(:)       :: ls
    real(dp),     allocatable,     dimension(:)       :: integrandx, integrandk
    real(dp),     allocatable,     dimension(:)       :: cls, cls2, ls_dp
    real(dp),     allocatable,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable,     dimension(:,:)     :: S, S2
    real(dp),     allocatable,     dimension(:,:)     :: Theta
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

    ! Open files to write transfer functions
    !open(unit=123, file=fileplace//"integrand1.dat", action="write", status="replace")
    !open(unit=124, file=fileplace//"integrand2.dat", action="write", status="replace")
    !open(unit=125, file=fileplace//"integrand3.dat", action="write", status="replace")
    !open(unit=126, file=fileplace//"integrand4.dat", action="write", status="replace")
    !open(unit=127, file=fileplace//"integrand5.dat", action="write", status="replace")
    !open(unit=128, file=fileplace//"integrand6.dat", action="write", status="replace")

    !Spline bessel functions, get second derivative for later splint
    do l=1,l_num
          call spline(z_spline, j_l(:,l), 1.d30, 1.d30, j_l2(:,l))
    end do


    allocate(Theta(l_num,k_num))
    allocate(integrandx(x_num))
    allocate(integrandk(k_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(x_lores(x_num/10))

    ! #### C_l COMPUTATION OVER l's ####

    do l = 1, l_num
      !################### METHOD 1 ##################
      !integralk = 0.d0 ! Reset k integral
      !do k = 1, k_num
      !  integralx = 0.d0 ! Reset x integral

      !  ! Compute integrand over x with 1/10th total array
      !  do i = 1, x_num/10.d0
      !     ilo = 1 + (i-1)*(x_num-1)/(x_num/10-1) !Speed up integration
      !     x_lores(i) = x_hires(ilo)
      !     integrandx(i)=S(ilo,k)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(k)*(get_eta(0.d0)-get_eta(x_hires(ilo))))
      !  end do

      !  ! Trapezoidal integration over x
      !  do i=1, x_num/10.d0-1
      !     integralx = integralx + (x_lores(i+1)-x_lores(i))*(integrandx(i+1)+integrandx(i))/2.d0
      !  end do
      !  Theta(l,k) = integralx
      !  ! Compute integrand over k
      !  integrandk(k) = (c*k_hires(k)/H_0)**(n_s-1.d0)*integralx**2/k_hires(k)
      !end do

      !! Trapezoidal integration over k
      !do k=1, k_num-1
      !   integralk = integralk + (k_hires(k+1)-k_hires(k))*(integrandk(k+1)+integrandk(k))/2.d0
      !end do

      !! Task: Store C_l in an array. Optionally output to file
      !cls(l) = integralk*ls(l)*(ls(l)+1.d0)/(2.d0*pi)

      !################### METHOD 2 ##################

      integralk = 0 ! Reset integral for each value of cls
      do k = 1, k_num
        integralx = 0 ! Reset integral for each value of theta

        ! Integrate theta
        do i = 1, x_num/10
          ilo = 1 + (i-1)*(x_num-1)/(x_num/10-1) !Speed up integration
          x_lores(i) = x_hires(ilo)
          integrandx(i) = S(ilo,k)*splint(z_spline,j_l(:,l),j_l2(:,l),k_hires(k)*(get_eta(0.d0)-get_eta(x_hires(ilo))))
          integralx = integralx + integrandx(i)
        end do
        ! Subtract half of first and last integrand for x
        Theta(l,k) = h1*(integralx - 0.5d0*(integrandx(1)+integrandx(x_num)))


        ! Integrate C_l
        integrandk(k) = (c*k_hires(k)/H_0)**(n_s-1.d0)*Theta(l,k)**2/k_hires(k)
        integralk = integralk + integrandk(k)
      end do

      ! Subtract half of first and last integrand for k
      integralk = h2*(integralk - 0.5d0*(integrandk(1)+integrandk(k_num)))

      ! Store C_l in an array. Optionally output to file
      cls(l) = integralk*ls(l)*(ls(l)+1.d0)/(2.d0*pi)
       !write the transfer function to file
      ! if(ls(l)==2) then
      !     do k=1,k_num
      !         write (123,'(*(2X, ES14.6E3))') c*k_hires(k)/H_0 , ls(l)*(ls(l)+1.d0)*&
      !               Theta(l,k)**2/(c*k_hires(k)/H_0)
      !     end do
      ! end if
      ! if(ls(l)==50) then
      !     do k=1,k_num
      !         write (124,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
      !               /(c*k_hires(k)/H_0)
      !     end do
      ! end if
      ! if(ls(l)==200) then
      !     do k=1,k_num
      !         write (125,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
      !               /(c*k_hires(k)/H_0)
      !     end do
      ! end if
      ! if(ls(l)==500) then
      !     do k=1,k_num
      !         write (126,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
      !               /(c*k_hires(k)/H_0)
      !     end do
      ! end if
      ! if(ls(l)==800) then
      !     do k=1,k_num
      !         write (127,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
      !               /(c*k_hires(k)/H_0)
      !     end do
      ! end if
      ! if(ls(l)==1200) then
      !     do k=1,k_num
      !         write (128,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
      !               /(c*k_hires(k)/H_0)
      !     end do
      ! end if

     ! Timer for loop
       call cpu_time(finish)
       write(*,*) "l = ", l
       print '("Time = ",f7.2," seconds.")',finish-start
    end do

    ! Close open data files
    !close(123)
    !close(124)
    !close(125)
    !close(126)
    !close(127)
    !close(128)

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

  end subroutine compute_cls

end module cl_mod
