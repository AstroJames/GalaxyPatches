!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz

  character(LEN=160)::infile
  real(dp)::a1,a2,z0,fac,T0,rho0,sigma0
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::smoothing,smoothing_l
  real(dp)::grav,potential_true,potential_lim,potential


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Kuijen & Gilmore 1989 disc
  if(gravity_type==2)then 
 
     ! Parameters of the model - convert to cgs
     a1=gravity_params(1) ! disc component factor
     a2=gravity_params(2) ! halo component factor
     z0=gravity_params(3) ! scale height of the disc
     T0=gravity_params(4) ! temperature in the mid-plane
     rho0=gravity_params(5) ! gas surface density
     a1=a1*3.08d21/(1d6*365.*24.*3600.)**2 ! cm/s2
     a2=a2/(1d6*365.*24.*3600.)**2 ! s-2
     z0=z0*3.08d21 ! cm
       
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xwind
#if NDIM>1
        ry=x(i,2)-ywind
#endif
#if NDIM>2
        rz=(x(i,3)-zwind)*scale_l
#endif
        f(i,1)=0.0d0
#if NDIM>1
        f(i,2)=0.0d0
#endif
#if NDIM>2
        smoothing_l=0.25*boxlen/2*scale_l
        smoothing=1.0d0
        f(i,3)=smoothing*(-a1*rz/sqrt(rz**2+z0**2)-a2*rz)*scale_t**2/scale_l
#endif
     end do
  end if

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
