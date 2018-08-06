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
  integer::nx_loc,ny_loc,nz_loc
  real(dp)::scale_m,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale
  real(dp),dimension(1:3)::skip_loc
  real(dp)::rcirc,zcirc,Rspher,rx,ry,rz,wx,wy,wz
  real(dp)::rho_s,R_s,Mb,hb,Mds,Mdg,hr,hz,Mdtot
  real(dp)::f_tot_i,f_nfw_i,f_bulge_i,f_disc_i
  real(dp)::pi,fourpi,GG
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m = scale_d*scale_l**3

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  scale=boxlen/dble(nx_loc)

  pi=ACOS(-1.0D0)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale

  ! Davide Martizzi: Model of NFW Halo + Bulge + Disc
  ! this section is for a non self-gravitating setup
  if(gravity_type==2)then
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
     
     do i=1,ncell
        rx=(x(i,1)-boxlen/2) ! code units  
        ry=(x(i,2)-boxlen/2) ! code units 
        rz=(x(i,3)-boxlen/2) ! code units 

        wx = 1.0 !TANH(-(ABS(rx)-0.9*boxlen/2))
        wy = 1.0 !TANH(-(ABS(ry)-0.9*boxlen/2))
        wz = 1.0 !TANH(-(ABS(rz)-0.9*boxlen/2))
        !if(ABS(rx)>0.95*boxlen/2) rx = 0.95*boxlen/2*rx/ABS(rx)
        !if(ABS(ry)>0.95*boxlen/2) ry = 0.95*boxlen/2*ry/ABS(ry)
	!if(ABS(rz)>0.95*boxlen/2) rz = 0.95*boxlen/2*rz/ABS(rz)
        
        rcirc=SQRT(rx*rx+ry*ry)
        zcirc=rz
        Rspher=SQRT(rcirc**2+zcirc**2)

        ! Force along x axis
        f_nfw_i=4.0*pi*GG*rho_s*(R_s**3)*(rx/(Rspher+1.0d-10)**2/(R_s+Rspher)-rx*LOG(1+(Rspher+1.0d-10)/R_s)/(Rspher+1.0d-10)**3)
        f_bulge_i=-GG*Mb*rx/SQRT(Rspher**2+hb**2)**3
        f_disc_i=-GG*(Mdtot)*rx/SQRT(rcirc**2+(hr+SQRT(zcirc**2+hz**2))**2)**3
        f_tot_i = f_nfw_i + f_bulge_i + f_disc_i
        f(i,1)=f_tot_i*wx*wy*wz !*(dx) ! code units
        
        ! Force along y axis
        f_nfw_i=4.0*pi*GG*rho_s*(R_s**3)*(ry/(Rspher+1.0d-10)**2/(R_s+Rspher)-ry*LOG(1+(Rspher+1.0d-10)/R_s)/(Rspher+1.0d-10)**3)
        f_bulge_i=-GG*Mb*ry/SQRT(Rspher**2+hb**2)**3
	f_disc_i=-GG*(Mdtot)*ry/SQRT(rcirc**2+(hr+SQRT(zcirc**2+hz**2))**2)**3
        f_tot_i = f_nfw_i + f_bulge_i + f_disc_i
        f(i,2)=f_tot_i*wx*wy*wz !*(dx) ! code units
        
        ! Force along z axis
        f_nfw_i=4.0*pi*GG*rho_s*(R_s**3)*(rz/(Rspher+1.0d-10)**2/(R_s+Rspher)-rz*LOG(1+(Rspher+1.0d-10)/R_s)/(Rspher+1.0d-10)**3)
        f_bulge_i=-GG*Mb*rz/SQRT(Rspher**2+hb**2)**3
        f_disc_i=-GG*(Mdtot)*zcirc*(hr+SQRT(zcirc**2+hz**2))/SQRT(zcirc**2+hz**2)/SQRT(rcirc**2+(hr+SQRT(zcirc**2+hz**2))**2)**3
        f_tot_i = f_nfw_i + f_bulge_i + f_disc_i
        f(i,3)=f_tot_i*wx*wy*wz !*(dx) ! code units
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
