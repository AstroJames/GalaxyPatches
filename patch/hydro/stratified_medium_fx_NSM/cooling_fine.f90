subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Compute sink accretion rates
  if(sink)call compute_accretion_rate(0)

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if(cooling.and.ilevel==levelmin.and.cosmo)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc,ix,iy,iz
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant,xcell,ycell,zcell,LAGN,mu
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk,r
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost,x_leaf,y_leaf,z_leaf
  real(kind=8),dimension(1:nvector),save::lambda_cool
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(dp)::xwinds,ywinds,zwinds
  real(dp),dimension(1:twotondim,1:ndim)::xc


  ! Variables needed for namelist read in.
  character(LEN=160)::infile


  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc) ! inverse of 'scale' in wind_fine,
                            !       so will divide insead of multiply
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)


  ! Location of the point source (grid or code units?)
  !
  ! Need this correction here, even though not in condinit and
  ! courant_fine. Distinction b/w user and code units? 
  xwinds = xwind/scale+dble(icoarse_min)
  ywinds = ywind/scale+dble(jcoarse_min)
  zwinds = zwind/scale+dble(kcoarse_min)


  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do 


  ! AGN bolometric luminosity in erg/s.
  LAGN = 10d0**logLAGN


  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)

           ! Follow procedure in wind_fine() to extract cell
	   ! coordinates. 
           x_leaf(nleaf) = xg(ind_grid(i),1)+xc(ind,1)
           y_leaf(nleaf) = xg(ind_grid(i),2)+xc(ind,2)
           z_leaf(nleaf) = xg(ind_grid(i),3)+xc(ind,3)

        end if
     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do
     
     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=(uold(ind_leaf(i),imetal)+uold(ind_leaf(i),imetal+1)+uold(ind_leaf(i),imetal+2))/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     
     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute distance from central AGN in cm.
     do i=1,nleaf
        xcell = x_leaf(i)
	ycell = y_leaf(i)
	zcell = z_leaf(i)

        
	! DEBUG
	!write(*,*) '---------------------------------'
	!write(*,*) 'xcell='
	!write(*,*) xcell
	!write(*,*) 'ycell='
	!write(*,*) ycell
	!write(*,*) 'zcell='
	!write(*,*) zcell
	!write(*,*) 'xwind='
	!write(*,*) xwind
	!write(*,*) 'ywind='
	!write(*,*) ywind
	!write(*,*) 'zwind='
	!write(*,*) zwind


        r(i) = ((xcell-xwinds)**2.0 + (ycell-ywinds)**2.0 + (zcell-zwinds)**2.0)**0.5

	!Convert from code units to cgs.
	r(i) = r(i)*scale_l
     end do


     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=exp(-nH(i)/0.01)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

     ! Compute net cooling at constant nH
     if(cooling)then
        ! Compute "thermal" temperature by substracting polytrope
        do i=1,nleaf
           T2(i) = max(T2(i)-T2min(i),T2_min_fix)
        end do
        
	! nH and T2 are cgs here.

	! provide r and LAGN in cgs units here.
	call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf,LAGN,r,lambda_cool)
     endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Compute net energy sink
     if(cooling)then
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch=uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1)
              if(cooling_switch>1d-3)then
                 delta_T2(i)=0
              endif
           end do
        endif
     endif

     ! Compute minimal total energy from polytrope
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0) + ekk(i)
     end do

     ! Update total fluid energy
     do i=1,nleaf
        T2(i) = uold(ind_leaf(i),ndim+2)
     end do
     if(cooling)then
        do i=1,nleaf
           T2(i) = T2(i)+delta_T2(i)
        end do
     endif
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2min(i)
        end do
     else
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = max(T2(i),T2min(i))
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=20d0*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=uold(ind_leaf(i),idelay)*damp_factor
        end do
     endif

     ! Store cooling rates in uold array
     if(cooling)then
        do i=1,nleaf
           ! Store as density-weighted cooling rate
           ! lambda_cool is in erg cm^-3 s^-1 (net cooling rate per unit volume)
           ! nH is in cm^-3 (number density of hydrogen)
           ! Product is stored to maintain conservation properties
           uold(ind_leaf(i),icool) = lambda_cool(i) * nH(i)
        end do
     endif

  end do
  ! End loop over cells

end subroutine coolfine1
