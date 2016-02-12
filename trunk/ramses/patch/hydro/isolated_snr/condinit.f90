! Set NX_LN as preprocessor variable: need to also set in init_flow_fine.f90
#define NX_LN 512

!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn,ln_d)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  
  real*4,dimension(1:NX_LN**3)::ln_d ! hard coding nx_ln
  

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  real(dp)::dctr, cs_h, dx_ln
  real(dp)::xcell, ycell, zcell, dx_, dy_, dz_, dr, sin_theta
  real(dp)::vx_BC, vy_BC, vz_BC, rho_BC, P_BC, disk_fact, lognorm_fact
  real(dp)::vri,vxi,vyi,vzi
  real(dp)::Ei,Vi,Pi_, M_ej, rho_ej, E_SN_th_int, E_SN_kin_int

  real(dp)::mu
  real(dp),parameter:: pi = 3.14159265
  real(dp),parameter:: kB = 1.3806200e-16 ! Boltzmann constant, cgs 

  real(dp)::scale, scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  real*4:: tmpvar
  character(LEN=160)::inj_type
  integer::i,j,k,id,iu,iv,iw,ip, i_ln, j_ln, k_ln, n_ln
  !integer,parameter::nx_ln = NX_LN
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables


  dx_ln = 1.0/real(NX_LN,dp)


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  
  ! Only for AMR runs trick to always inject energy
  if(levelmin.ne.nlevelmax.and.rc<dx)rc=dx  

  ! Convert SN thermal and kinetic energies from cgs to code units.
  E_SN_th_int = E_SN_th/(scale_v**2.0*scale_d*scale_l**3.0) 
  E_SN_kin_int = E_SN_kin/(scale_v**2.0*scale_d*scale_l**3.0) 


  ! Initial pressure in the injection radius, in code units.
  !
  ! Use for injection in thermal form.
  Vi = (4.0*pi*rc**3.0)/3.0 ! code units

  Pi_ = (2.0/3.0)*(E_SN_th_int/Vi)


  ! Density of SN ejecta, assuming uniformly distributed within rc, as
  ! in Thornton et al. (1998).
  M_ej = 5.97*1E33 ! 3 Msun, cgs

  ! Mass from cgs to code units.
  M_ej = M_ej/(scale_d*scale_l**3.0)

  rho_ej = M_ej/Vi


  ! Radial velocity of the gas within rc, so that total kinetic energy
  ! is E_SN_kin_int. To be used for energy injection in kinetic form.
  vri = (3.0*E_SN_kin_int/(2.0*pi*(rhocool + rho_ej)*rc**3.0))**0.5


  ! As in courant_file.f90.
  mu = 0.6


  ! Rescaling factor
  scale=dble(icoarse_max-icoarse_min+1)/boxlen
  

  ! Indexes of physical quantities.

  id = 1
  iu = 2
  ip = 3
#if NDIM>1
  id = 1
  iu = 2
  iv = 3
  ip = 4
#endif
#if NDIM>2
  id = 1
  iu = 2
  iv = 3
  iw = 4
  ip = 5
#endif

  ! Call built-in initial conditions generator: this is for simple
  ! 'square' and 'point' regions.
  call region_condinit(x,q,dx,nn)

  
  ! Dimensionless lognormal field is read in init_flow_fine.f90,
  ! which calls this function.
    
  ! Loop through positions for the current sub-grid.
  do i=1,nn

    ! Distance from wind source.
    xcell = x(i,1) 
#if NDIM>1
    ycell = x(i,2)
#endif
#if NDIM>2
    zcell = x(i,3)
#endif

    dctr = ((xcell-xwind)**2.0)**0.5
#if NDIM>1
    dctr = ((xcell-xwind)**2.0 + (ycell-ywind)**2.0)**0.5
#endif
#if NDIM>2
    dctr = ((xcell-xwind)**2.0 + (ycell-ywind)**2.0 + (zcell-zwind)**2.0)**0.5
#endif

    if (dctr<rc) then

      ! Within r_c, inject energy of the supernova remmant in thermal
      ! form. Must make sure that r_c is well within Sedov radius, to
      ! avoid initial energy losses.


      ! Velocity components for radial kinetic energy injection.
      vxi = (xcell - xwind)*vri/dctr
      vyi = (ycell - ywind)*vri/dctr
      vzi = (zcell - zwind)*vri/dctr

      q(i,id) = rhocool + rho_ej
      q(i,iu) = vxi
#if NDIM>1
      q(i,iu) = vxi
      q(i,iv) = vyi
#endif
#if NDIM>2
      q(i,iu) = vxi
      q(i,iv) = vyi
      q(i,iw) = vzi
#endif
      q(i,ip) = Pi_
      if(E_SN_th_int < 0.0d0) q(i,ip) = Pc*(rhocool + rho_ej)/rhocool ! For runs with no thermal energy

    else 
      ! Outside inner boundary condition.

      ! Smooth or fractal ambient medium?
      if (medium_type=='uni' .or. medium_type=='disk') then
        lognorm_fact = 1.0

      else if (medium_type=='fractal' .or. medium_type=='disk_fractal') then
          	  
        ! Given (x, y, z), find nearest grid point in pre-computed
        ! lognormal field.

        ! Assume coordinates are from 0 to 1 along each direction.
        i_ln = int(x(i,1)*scale/dx_ln)+1
#if NDIM>1
        j_ln = int(x(i,2)*scale/dx_ln)+1
#endif
#if NDIM>2
        k_ln = int(x(i,3)*scale/dx_ln)+1
#endif

        n_ln = (i_ln-1)*int(NX_LN)**2.0 + (j_ln-1)*int(NX_LN) + k_ln

        lognorm_fact = real(ln_d(n_ln),dp)
      end if

      ! Uniform or disk-like ambient medium?
      if (medium_type=='uni' .or. medium_type=='fractal') then
	disk_fact = 1.0
	
      else if (medium_type=='disk' .or. medium_type=='disk_fractal') then
 	! r0, h from namelist 

	! calculate theta here.
	dx_ = xcell - xwind
	dy_ = ycell - ywind
	dz_ = zcell - zwind
        dr = sqrt(dx_**2.0 + dy_**2.0)
	sin_theta = dr/dctr

        
        ! h is really h/r ratio.
        
	! As in Nathan Roth's paper.
	!disk_fact = ((dctr/r0)**(-gam))*exp((h**(-2.0))*(sin_theta-1.0))

        ! Self-gravitating (fg=1), isothermal sheet. Assume nHc is
	! define as rho(100 pc). h is h/r ratio.
	disk_fact = 0.5*((dr/0.1)**(-2.0))
	disk_fact = disk_fact*4.0/((exp(abs(dz_)/(2.0*h*dr)) + exp(-abs(dz_)/(2.0*h*dr)))**2.0)
     
        ! Impose a minimum value on disk_fact, to avoid errors due to
	! very small numbers.
	!
	! Be careful: 10^-6 of number density of ULIRG disk is not
	! that small, so may over estimate halo gas density.
	disk_fact = max(disk_fact, 10.0**(-6.0))
     
     end if 
	  

      ! Density.

      q(i,id) = lognorm_fact*disk_fact*rhocool


      ! Constant pressure.
      q(i,ip) = Pc

      ! Static.
      q(i,iu) = 0.0
#if NDIM>1
      q(i,iv) = 0.0
#endif
#if NDIM >2
      q(i,iw) = 0.0
#endif
    end if ! inner BC or not

    if(metal)q(i,ndim+3)=z_ave*0.02 ! customized metallicity value.

  end do ! loop over cells

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
