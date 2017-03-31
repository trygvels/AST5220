program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  implicit none
  integer(i4b) :: i      ! Number of grid points

  ! Initialize time grids
  call initialize_time_mod

  call initialize_rec_mod

  ! ------ Output to file desired quantities here ------
  ! Write to file - x_rec, X_e
  open(1, file="Xe.dat", action="write",status="replace")
  open(2, file="neLog.dat", action="write",status="replace") ! Eldens
  open(3, file="tau.dat", action="write",status="replace") ! taus
  open(4, file="g.dat", action="write",status="replace") ! visibility

  do i=1, n
     write(1,*) x_rec(i), X_e(i)
     write(2,*) n_e(i), n_e2(i)
     write(3,*) tau(i),dtau(i), tau2(i)
     write(4,*) g(i),dg(i), g2(i)
  end do

  close(1)
  close(2)
  close(3)
  close(4)




end program cmbspec
