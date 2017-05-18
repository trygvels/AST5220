module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none
  real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
  real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
  real(dp),     allocatable, dimension(:)       :: x_hires, k_hires

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, k, l, l_num, x_num, n_spline, k_num
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: x !,k
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     pointer,     dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta

    real(dp)           :: t1, t2, integral
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
    call get_hires_source_function(k_hires,x_hires,S) !These are pointers, work this out.

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    ! Calculate bessel functions
    do i = 1, n_spline
      z_spline(i) = (i-1)*3500d0/(n_spline-1)

      if (z_spline(i)>2.d0) then
        do l=1, l_num
          call sphbes(ls(l),z_spline(i),j_l(i,l))
        end do
      endif
    end do

    !Spline bessel functions TODO Should we use 2D spline?
    do l=1,l_num
          call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
    end do


    allocate(Theta(k_num,x_num))
    ! Overall task: Compute the C_l's for each given l
    do l = 1, l_num

    !   ! Task: Compute the transfer function, Theta_l(k)
    !   do k=1, k_num
!
    !     do i=1, x_num
    !       integrand(i) = S(i,k)*j_lfunc(l,k_hires(k),x_hires(i))
    !     end do
!
    !     Theta_l(l,k) = 0.5d0*(integrand(1)+integrand(x_num))

    !     do i = 2, x_num - 1
    !       Theta_l(l,k) = Theta_l(l,k) + integrand(i)
    !     end do
!
    !     Theta_l(l,k) = h1*Theta_l(l,k)
    !   end do
    !   ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
!
    !   do k = 1, k_num
    !     integrand(l,k) = (c*k_hires(k)/H_0)**(n_s-1)*Theta_l(l,k)**2/k_hires(k)
    !   end do

    !   ! Task: Store C_l in an array. Optionally output to file
    !   C_lint = 0.5d0*(integrand2(l,1)+integrand2(l,n_k_highres))
    !   do k=2,n_k_highres-1
               !write(*,*) 'k=',k
    !           C_lint = C_lint + integrand2(l,k)
    !   end do

       !Store C_l in an array.
    !   cls(l) = h2*C_lint *ls(l)*(ls(l)+1.d0)/(2.d0*pi)
    end do


    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    !call spline(ls, cls, yp1, ypn, cls2)

  end subroutine compute_cls

  function j_lfunc(l,k,x)
    implicit none
    integer(i4b), intent(in) :: l
    real(dp), intent(in)     :: x, k
    real(dp)                 :: j_lfunc

    j_lfunc = splint(z_spline, j_l(:,l), j_l2(:,l), k*(get_eta(0.d0)-get_eta(x)))
  end function j_lfunc

end module cl_mod
