!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  character(LEN=160)::infile
  real(dp)::a1,a2,z0,f,T0,rho0,sigma0
  real(dp)::rx,ry,rz,rr
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2


  ! Conversion factor from user units to cgs units 
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)


  ! Parameters of the model - convert to cgs
  a1=gravity_params(1) ! disc component factor
  a2=gravity_params(2) ! halo component factor
  z0=gravity_params(3) ! scale height of the disc
  T0=gravity_params(4) ! temperature in the mid-plane
  sigma0=gravity_params(5) ! gas surface density
  a1=a1*3.08d21/(1d6*365.*24.*3600.)**2 ! cm/s2
  a2=a2/(1d6*365.*24.*3600.)**2 ! s-2
  z0=z0*3.08d21 ! cm
  T0=T0 ! K
  f=0.6*1.66d-24/1.38d-16/T0 ! s2/cm2
  sigma0=sigma0*2d33/(3.08d21)**2 ! g/cm2
  rho0=sigma0/1.32171d21 ! g/cm3

  do i=1,ncell
     rz=(x(i,3)-zwind)*scale_l
     d(i)=rho0*exp(-(a1*f*(sqrt(rz**2+z0**2)-z0)+a2*f*(rz**2)/2))/scale_d
  end do

end subroutine rho_ana
