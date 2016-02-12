!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,ind_cell,x_ref,ind_cell_ref,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen,levelmin
  use hydro_parameters, ONLY: nvar,boundary_var
  use hydro_commons, ONLY: uold,gamma
  use poisson_commons, ONLY: f
  use poisson_parameters, ONLY: gravity_params
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp),dimension(1:nvector,1:ndim)::x_ref ! Neighboring cell center position.
  integer,dimension(1:nvector)::ind_cell  ! Cell index
  integer,dimension(1:nvector)::ind_cell_ref  ! Neighboring cell index   
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i
  real(dp)::switch
  real(dp)::mu,XH,a1,a2,z0,T0,rho0,fp
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::dz,rhoprime,fprime,dz_,smoothing,smoothing_l,grav

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! As in courant_file.f90. 
  mu = 0.6
  XH=0.76

  ! Parameters of the stratified medium model - set gravity_type=2 too - convert to cgs
  a1=gravity_params(1) ! disc component factor
  a2=gravity_params(2) ! halo component factor
  z0=gravity_params(3) ! scale height of the disc 
  T0=gravity_params(4) ! temperature in the mid-plane
  rho0=gravity_params(5) ! density in the mid-plane g/cm3
  a1=a1*3.08d21/(1d6*365.*24.*3600.)**2 ! cm/s2
  a2=a2/(1d6*365.*24.*3600.)**2 ! s-2 
  z0=z0*3.08d21 ! cm
  T0=T0 ! K
  fp=0.6*1.66d-24/1.38d-16/T0 ! s2/cm2

  do i=1,ncell
     do ivar=1,nvar
        switch=(x(i,3)-boxlen/2.)/abs(x(i,3)-boxlen/2.)
        if(ivar.eq.1)then
           u(i,ivar)=uold(ind_cell_ref(i),ivar)
        end if
        if(ivar.gt.1.and.ivar.lt.ndim+1)then
           u(i,ivar)=uold(ind_cell_ref(i),ivar)
        end if
        if(ivar.eq.ndim+1)then
           !u(i,ivar)=0.0 ! z velocity at the boundary set to 0
           !u(i,ivar)=uold(ind_cell_ref(i),ivar)
           if(switch*uold(ind_cell_ref(i),ivar)>0)then
              u(i,ivar)=uold(ind_cell_ref(i),ivar)
           else
              u(i,ivar)=0.0
           end if
        end if
        if(ivar.eq.ndim+2)then
           !u(i,ivar)=uold(ind_cell_ref(i),ivar) ! zero pressure gradient
           dz=x(i,3)-x_ref(i,3)
           fprime=(f(ind_cell(i),3)-f(ind_cell_ref(i),3))/dz ! force gradient
           rhoprime=(u(i,1)-uold(ind_cell_ref(i),1))/dz
           u(i,ivar)=uold(ind_cell_ref(i),ivar) !+ &
                !& dz*uold(ind_cell_ref(i),1)*f(ind_cell_ref(i),3)/(gamma-1.0d0)+ & 
                !& dz**2/2*uold(ind_cell_ref(i),1)*fprime/(gamma-1.0d0)+ &
                !& dz**2/2*rhoprime*f(ind_cell_ref(i),3)/(gamma-1.0d0)+ &
                !& dz**3/3*fprime*rhoprime/(gamma-1.0d0)
        end if
        if(ivar.gt.ndim+2)then
           u(i,ivar)=uold(ind_cell_ref(i),ivar)
        end if
     end do
  end do

end subroutine boundana
