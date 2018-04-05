subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg
  
  real(dp)::scale,dx0

  ! Rescaling factor 
  scale=dble(icoarse_max-icoarse_min+1)/boxlen


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  dx0=0.5D0**ilevel
  vol=dx0**ndim

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        
        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do
        
        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,ndim+2)*vol
        end do
        
        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,ndim+2)*vol
        end do
        do ivar=1,ndim
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol
           end do
        end do
        
        ! Compute CFL time-step
        if(nleaf>0)then
           dx=dx0/scale 
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)


           ! DEBUG: much larger than we would expect for dx/vin...
           !write(*,*) 'dt_loc='
	   !write(*,*) dt_loc


        end if
        
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass_loc,mass_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(eint_loc,eint_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(  dt_loc,  dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
#else
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif  
  
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)



  ! Comment out continuous wind injection for SNR case. 
  !if(ilevel==nlevelmax)call wind_fine(ilevel)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine


subroutine wind_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_parameters
  implicit none
  integer::ilevel
  
  integer::igrid,ncache,i,ind,iskip,ngrid
  integer::ix,iy,iz,idim,ivar
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:nvector),save::rr
 

  real(dp),parameter:: pi = 3.14159265
  real(dp),parameter:: kB = 1.3806200e-16 ! Boltzmann constant, cgs 
  character(LEN=160)::infile
  real(dp)::cx, cy, cz, xcell, ycell, zcell, dctr 
  real(dp)::Lcgs, Luser, Luser_dens, mu
  real(dp)::vx_BC, vy_BC, vz_BC, rho_BC, P_BC


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units 
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)


  ! Parameters of energy injection.
  !
  ! Now read in from namelist; expect rc, vin, Mdot_in in user units.


  ! DEBUG: hardcoded mean molecular weight.
  mu = 0.6


  ! Assume grid units = user (code) units, since boxlen=1.
  Lcgs = (1.0D51)/31536000.0 ! 1 SNe per year
  Luser = Lcgs/(scale_d*scale_l**2.0*scale_v**3.0)

  Luser_dens = Luser/(4.0*pi*rc**3.0/3.0)


  ! Mesh size
  dx=0.5d0**ilevel

  ! Rescaling factor
  scale=dble(icoarse_max-icoarse_min+1)/boxlen
  

  ! Location of the point source (grid or code units?)
  !
  ! Need this correction here, but not in condinit. Distinction
  ! between code and user units?
  xwind = xwind*scale+dble(icoarse_min)
  ywind = ywind*scale+dble(jcoarse_min)
  zwind = zwind*scale+dble(kcoarse_min)


  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx !-xwind
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx !-ywind
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx !-zwind
  end do
 

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells (i.e., cell 1, 2, 3, ..., 8 for each oct).
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        
	do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)

           xcell = xg(ind_grid(i),1)+xc(ind,1)
           ycell = xg(ind_grid(i),2)+xc(ind,2)
           zcell = xg(ind_grid(i),3)+xc(ind,3)

           dctr = ((xcell-xwind)**2.0 + (ycell-ywind)**2.0 + (zcell-zwind)**2.0)**0.5

	   ! Spatial condition, then add energy.
           if(dctr<rc)then

             ! Enforce boundary conditions at the sphere boundary.
             
	     ! Wind velocity, pointing outward from the center.
	     
	     vx_BC = (xcell - xwind)*vin/dctr
	     vy_BC = (ycell - ywind)*vin/dctr
	     vz_BC = (zcell - zwind)*vin/dctr

             rho_BC = Mdot_in/(4.0*pi*rc**2.0*vin)
	     rho_BC = rho_BC*(dctr/(alph*rc+(1.0-alph)*dctr))**(-2.0)
             P_BC = (Tin_K/mu)*rho_BC/scale_T2

             uold(ind_cell(i),1) = rho_BC
             uold(ind_cell(i),2) = rho_BC*vx_BC
             uold(ind_cell(i),3) = rho_BC*vy_BC
             uold(ind_cell(i),4) = rho_BC*vz_BC
             uold(ind_cell(i),5) = 0.5*rho_BC*vin**2.0 + P_BC/(gamma-1.0d0)

           end if 

       end do

     end do

  end do


  ! Update boundaries
  do ivar=1,nvar
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do

111 format('   Entering wind_fine for level ',I2)

end subroutine wind_fine
