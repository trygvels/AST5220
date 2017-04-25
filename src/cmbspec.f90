program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  implicit none
  integer(i4b) :: i      ! Number of grid points

  ! Initialize time grids
  call initialize_time_mod
  ! Initialize recombination module
  call initialize_rec_mod

  ! ------ Output to file desired quantities here ------
  ! Write to file - x_rec, X_e
  !open(1, file="Xe.dat", action="write",status="replace")
  !open(2, file="neLog.dat", action="write",status="replace") ! Eldens
  !open(3, file="tau.dat", action="write",status="replace") ! taus
  !open(4, file="g.dat", action="write",status="replace") ! visibility

  !do i=1, n
  !   write(1,*) x_rec(i), X_e(i)
  !   write(2,*) n_e(i), n_e2(i)
  !   write(3,*) tau(i),dtau(i), tau2(i)
  !   write(4,*) g(i),dg(i), g2(i)
  !end do

  !close(1)
  !close(2)
  !close(3)
  !close(4)

  ! Initialize perturbation module
  call initialize_perturbation_eqns
  call integrate_perturbation_eqns
  !write(*,*) 'Saving perturbations to file'
  open(1,file="delta.dat",action="write",status="replace")
  open(2,file="delta_b.dat",action="write",status="replace")
  open(3,file="v.dat",action="write",status="replace")
  open(4,file="v_b.dat",action="write",status="replace")
  open(5,file="Phi.dat",action="write",status="replace")
  open(6,file="Psi.dat",action="write",status="replace")
  open(7,file="Theta0.dat",action="write",status="replace")
  open(8,file="dPhi.dat",action="write",status="replace")
  open(9,file="dPsi.dat",action="write",status="replace")
  
  do i=0,n_t
      write(1,'(*(2X, ES14.6))') delta(i,1),delta(i,5),delta(i,10),delta(i,40),delta(i,60),delta(i,100)
      write(2,'(*(2X, ES14.6))') delta_b(i,1),delta_b(i,5),delta_b(i,10),delta_b(i,40),delta_b(i,60),delta_b(i,100)
      write(3,'(*(2X, ES14.6))') v(i,1),v(i,5),v(i,10),v(i,40),v(i,60),v(i,100)
      write(4,'(*(2X, ES14.6))') v_b(i,1),v_b(i,5),v_b(i,10),v_b(i,40),v_b(i,60),v_b(i,100)
      write(5,'(*(2X, ES14.6))') Phi(i,1),Phi(i,5),Phi(i,10),Phi(i,40),Phi(i,60),Phi(i,100)
      write(6,'(*(2X, ES14.6))') Psi(i,1),Psi(i,5),Psi(i,10),Psi(i,40),Psi(i,60),Psi(i,100)
      write(7,'(*(2X, ES14.6))') Theta(i,0,1),Theta(i,0,5),Theta(i,0,10),Theta(i,0,40),Theta(i,0,60),Theta(i,0,100)
      write(8,'(*(2X, ES14.6))') dPhi(i,1),dPhi(i,5),dPhi(i,10),dPhi(i,40),dPhi(i,60),dPhi(i,100)
      write(9,'(*(2X, ES14.6))') dPsi(i,1),dPsi(i,5),dPsi(i,10),dPsi(i,40),dPsi(i,60),dPsi(i,100)
  end do

  close(1)
  close(2)
  close(3)
  close(4)
  close(5)
  close(6)
  close(7)
  close(8)
  close(9)

  open(123,file="dTheta2.dat",action="write",status="replace")
  open(124,file="Theta2.dat",action="write",status="replace")
  do i=1,n_t
      write(123,'(*(2X, ES14.6))') dTheta(i,2,1),dTheta(i,2,5),dTheta(i,2,10),dTheta(i,2,40),dTheta(i,2,60),dTheta(i,2,100)
      write(124,'(*(2X, ES14.6))') Theta(i,2,1),Theta(i,2,5),Theta(i,2,10),Theta(i,2,40),Theta(i,2,60),Theta(i,2,100)
  end do
  close(123)
  close(124)



end program cmbspec
