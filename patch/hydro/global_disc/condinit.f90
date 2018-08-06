!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use cooling_module, only: mH, kB 
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
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
#if NENER>0 || NVAR>NDIM+2+NENER
  integer::ivar
#endif
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  integer::i
  real(dp)::rcirc,zcirc,rx,ry,rz,Rspher
  real(dp)::rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot
  real(dp)::minusdphidz,phi,vcirc
  real(dp)::phi1,phi0
  
  real(dp)::f_tot_i,f_nfw_i,f_bulge_i,f_disc_i
  real(dp)::phi_tot_i,phi_nfw_i,phi_bulge_i,phi_disc_i,Tref,zef
  real(dp)::rhocell,Tcell,Pcell,vc,vx,vy,vz
  real(dp)::scale_m,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::pi,GG

  pi=ACOS(-1.0D0)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)  
  scale_m=scale_l**3*scale_d
  
  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Davide Martizzi: initial conditions for a gas exponential disc

  ! Some model parameters
  rho_s=gravity_params(1)/scale_d         ! NFW scale density in g/cm3 to code units 
  R_s=gravity_params(2)*3.08d21/scale_l   ! NFW scale radius kpc to code units
  Mb=gravity_params(3)*2.0d33/scale_m     ! Mass of the bulge Msun to code units
  hb=gravity_params(4)*3.08d21/scale_l    ! Scale radius of the bulge kpc to code units
  Mds=gravity_params(5)*2.0d33/scale_m    ! Stellar mass of the disc Msun to code units
  Mdg=gravity_params(6)*2.0d33/scale_m    ! Gas mass of the disc Msun to code units
  hr=gravity_params(7)*3.08d21/scale_l    ! Disc scale radius kpc to code units
  hz=gravity_params(8)*3.08d21/scale_l    ! Disc scale height kpc to code units

  GG=6.67d-8*scale_t**2*scale_d  ! G in code units
  
  if(.not.gas_sg)then
     Mdtot = Mds + Mdg ! simulation with static potential
  else
     Mdtot = Mds ! simulation with static potential plus gas self-gravity
  end if
  
  do i=1,nn
    ! Position of the cell
    rx = (x(i,1)-boxlen/2) ! code units 
    ry = (x(i,2)-boxlen/2) ! code units 
    rz = (x(i,3)-boxlen/2) ! code units 

    rcirc=SQRT(rx*rx+ry*ry)
    zcirc=rz
    Rspher=SQRT(rcirc**2+zcirc**2)
    
    ! Set gas in hydrostatic equilibrium
    Tcell = Tmu0gas/scale_T2
    phi1 = phi(rcirc,zcirc,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot)    
    phi0 = phi(rcirc,0.0d0,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot)
    ! Set density 
    rhocell = rho0gas/scale_d*EXP(-(rcirc/hr)**2)*EXP(-(phi1-phi0)/Tcell)

    if(rhocell < rhoamb/scale_d) then
       ! Ambient density
       rhocell = rhoamb/scale_d
       ! Set circular velocity
       vc = vcirc(rcirc,zcirc,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot,Tcell,1)
    else
       ! Set circular velocity
       vc = vcirc(rcirc,zcirc,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot,Tcell,0)
    end if

    ! Set Pressure
    Pcell = rhocell*Tcell
    
    ! density
    q(i,1)=rhocell
    ! vx
    q(i,2)=-vc*ry/rcirc
#if NDIM>1
    ! vy
    q(i,3)=vc*rx/rcirc
#endif
#if NDIM>2
    ! vz
    q(i,4)=0.0d0
#endif
    ! pressure 
    q(i,ndim+2)=Pcell
#if NENER>0
  ! radiative pressure 
  do ivar=1,nener
     q(i,ndim+2+ivar)=0.0d0
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! background metallicity
  if(metal)q(i,imetal) = z_ave*0.02 
#endif
 end do
  
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
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,ndim+2+ivar)=q(1:nn,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
real*8 FUNCTION phi(rcirc,zcirc,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot)
  implicit none
  real*8::rcirc,zcirc,Rspher,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot
  real*8::phi_tot_i,phi_nfw_i,phi_bulge_i,phi_disc_i
  real*8::pi,GG

  pi=ACOS(-1.0D0)
  GG=1.0

  Rspher=SQRT(rcirc**2+zcirc**2)  
  phi_nfw_i = -4.0*pi*GG*rho_s*(R_s**3)*LOG(1+(Rspher+1.0d-10)/R_s)/(Rspher+1.0d-10)
  phi_bulge_i = -GG*Mb/SQRT(Rspher**2+hb**2)
  phi_disc_i = -GG*(Mdtot)/SQRT(rcirc**2+(hr+SQRT(zcirc**2+hz**2))**2)
  phi_tot_i = phi_nfw_i + phi_bulge_i + phi_disc_i ! code units
  phi = phi_tot_i
  
  return 
end FUNCTION phi
!================================================================
!================================================================
!================================================================
!================================================================
real*8 FUNCTION vcirc(rcirc,zcirc,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot,TT,sw)
  implicit none
  integer::sw
  real*8::rcirc,zcirc,Rspher,rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot,TT
  real*8::fr_tot_i,fr_nfw_i,fr_bulge_i,fr_disc_i,rho_term,T_term
  real*8::phi
  real*8::pi,GG,ww

  pi=ACOS(-1.0D0)
  GG=1.0
 
  ! Term 1: radial gradient of potential
  Rspher = rcirc ! the term for this model is in the midplane
  fr_nfw_i=4.0*pi*GG*rho_s*(R_s**3)*(-rcirc/(Rspher+1.0d-10)**2/(R_s+Rspher)+rcirc*LOG(1+(Rspher+1.0d-10)/R_s)/(Rspher+1.0d-10)**3)
  fr_bulge_i=GG*Mb*rcirc/SQRT(Rspher**2+hb**2)**3
  fr_disc_i=GG*(Mdtot)*rcirc/SQRT(rcirc**2+(hr+hz)**2)**3
  fr_tot_i = fr_nfw_i + fr_bulge_i + fr_disc_i ! code units
    
  ! Term 2: radial gradient of pressure (density term)
  if(sw == 0) rho_term = -2*TT*rcirc/hr**2
  if(sw == 1) rho_term = 0.0d0
  
  ! Term 3: radial gradient of pressure (temperature term)
  T_term = 0.0 ! gradient of temperature is zero
  
  if(fr_tot_i+rho_term+T_term < 0)write(*,*)fr_tot_i,rho_term,T_term
  
  ww = 1.0 
  vcirc = SQRT((fr_tot_i+rho_term+T_term)*rcirc)*ww ! code units

  return
end FUNCTION vcirc
