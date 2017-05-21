module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none
  real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
  real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
  real(dp),     allocatable, dimension(:)       :: x_hires, k_hires,  l_hires, cl_hires

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, k, l, l_num, x_num, n_spline, k_num
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrandx, integrandk
    real(dp),     allocatable,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     allocatable,     dimension(:)       :: x !,k
    real(dp),     allocatable,     dimension(:,:,:,:) :: S_coeff
    real(dp),     allocatable,     dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta
    real :: start, finish
    real(dp)           :: integralx, integralk, h1, h2
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Hires array numbers
    x_num = 5000
    k_num = 5000

    ! Task: Get source function from evolution_mod

    ! Allocate Source function and Hires grids
    allocate(S(x_num,k_num))
    allocate(x_hires(x_num))
    allocate(k_hires(k_num))

    ! Calculate Hires source function from evolution_mod
    write(*,*) "Calculating hires S, x and k"
    call get_hires_source_function(x_hires,k_hires,S) ! S is pointer, LEARN
    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    ! Calculate bessel functions, needed for LOS integration
    do i = 1, n_spline
      z_spline(i) = (i-1)*3400.d0/(n_spline-1.d0)
      ! TODO: WHAT IS HAPPENING HERE
      do l=1, l_num
        if (z_spline(i)>2.d0) then
          call sphbes(ls(l),z_spline(i),j_l(i,l))
        endif
      end do
    end do


    open(unit=123, file="integrand1.dat", action="write", status="replace")
    open(unit=124, file="integrand2.dat", action="write", status="replace")
    open(unit=125, file="integrand3.dat", action="write", status="replace")
    open(unit=126, file="integrand4.dat", action="write", status="replace")
    open(unit=127, file="integrand5.dat", action="write", status="replace")
    open(unit=128, file="integrand6.dat", action="write", status="replace")
    !Spline bessel functions, get second derivative for later splint
    do l=1,l_num
          call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
    end do


    allocate(Theta(l_num,k_num))
    allocate(integrandx(x_num))
    allocate(integrandk(k_num))
    allocate(cls(l_num))
    allocate(cls2(l_num))


    ! #### C_l COMPUTATION OVER l's ####
    ! Constants for trapezoidal integration
    h1 = (x_hires(x_num) - x_hires(1))/x_num
    h2 = (k_hires(k_num) - k_hires(1))/k_num
    do l = 1, l_num

       ! ##### UNIFORM TRAPEZOIDAL INTEGRATION #####
       ! Integrating both functions in same loop!

       integralk = 0 ! Reset integral for each value of cls
       do k = 1, k_num
         integralx = 0 ! Reset integral for each value of theta

         ! Integrate theta
         do i = 1, x_num
           integrandx(i) = S(i,k)*j_lfunc(l,k_hires(k),x_hires(i))
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
       if(ls(l)==2) then
           do k=1,k_num
               write (123,'(*(2X, ES14.6E3))') c*k_hires(k)/H_0 , ls(l)*(ls(l)+1.d0)*&
                     Theta(l,k)**2/(c*k_hires(k)/H_0)
           end do
       end if
       if(ls(l)==50) then
           do k=1,k_num
               write (124,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
                     /(c*k_hires(k)/H_0)
           end do
       end if
       if(ls(l)==200) then
           do k=1,k_num
               write (125,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
                     /(c*k_hires(k)/H_0)
           end do
       end if
       if(ls(l)==500) then
           do k=1,k_num
               write (126,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
                     /(c*k_hires(k)/H_0)
           end do
       end if
       if(ls(l)==800) then
           do k=1,k_num
               write (127,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
                     /(c*k_hires(k)/H_0)
           end do
       end if
       if(ls(l)==1200) then
           do k=1,k_num
               write (128,'(*(2X, ES14.6E3))') ls(l)*(ls(l)+1.d0)*Theta(l,k)**2 &
                     /(c*k_hires(k)/H_0)
           end do
       end if

      ! Timer for loop
      call cpu_time(finish)
      write(*,*) "l = ", l
      print '("Time = ",f7.2," seconds.")',finish-start
    end do

    ! Close open data files
    close(123)
    close(124)
    close(125)
    close(126)
    close(127)
    close(128)


    write(*,*) 'converting ls to double precision'
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
    do l = 1, int(maxval(ls))
      l_hires(l) = l ! dp hires l
      cl_hires(l) =  splint(ls_dp, cls, cls2, l_hires(l))
    end do
  end subroutine compute_cls

  function j_lfunc(l,k,x)
    implicit none
    integer(i4b), intent(in) :: l
    real(dp), intent(in)     :: x, k
    real(dp)                 :: j_lfunc
    ! Spline integration for j_l bessel functions
    j_lfunc = splint(z_spline, j_l(:,l), j_l2(:,l), k*(get_eta(0.d0)-get_eta(x)))
  end function j_lfunc

end module cl_mod
