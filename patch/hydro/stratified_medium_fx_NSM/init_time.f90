subroutine init_time
  use amr_commons
  use hydro_commons
  use cooling_module
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Initialize time and cooling parameters
  ! This patch version uses runtime cooling_model
  !--------------------------------------------------
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(kind=8),save::T2_sim
  integer::i_frw
  real(dp)::time_simu,tau_simu

  if(myid==1)write(*,*)'Entering init_time'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(cooling)then
     ! Read cooling model from namelist (cooling_model parameter)
     ! cooling_model = 6 (Courty - default complex chemical network)  
     ! cooling_model = 7 (Schure - bistable ISM with constant UV heating)
     if(myid==1)then
        write(*,*)'Computing cooling model'
        write(*,*)'Using cooling_model = ',cooling_model
        if(cooling_model==7)then
           write(*,*)'Schure bistable ISM cooling'
           write(*,*)'UV heating rate = ',schure_heating_rate,' erg/s'
        else if(cooling_model==6)then
           write(*,*)'Courty chemical network cooling'
        else
           write(*,*)'Standard cooling model ',cooling_model
        endif
     endif
     
     ! Initialize cooling with selected model
     if(cosmo)then
        ! Reionization redshift has to be later than starting redshift
        z_reion=min(1./(1.1*aexp_ini)-1.,z_reion)
        call set_model(cooling_model,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call set_model(cooling_model,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(1.),T2_sim)
        T2_start=T2_sim
     endif
  endif

  ! Initialize time variables
  t=0.0d0
  nstep_coarse=0
  nstep_coarse_old=0
  
  ! Get initial time steps
  if(cosmo)then
     ! Find neighboring times
     i_frw=1
     do while(aexp_frw(i_frw)>aexp)
        i_frw=i_frw+1
     end do
     ! Compute initial time
     time_simu=t_frw(i_frw  )*(aexp-aexp_frw(i_frw-1))/(aexp_frw(i_frw  )-aexp_frw(i_frw-1))+ &
          & t_frw(i_frw-1)*(aexp-aexp_frw(i_frw  ))/(aexp_frw(i_frw-1)-aexp_frw(i_frw  ))
     ! Compute initial conformal time
     tau_simu=tau_frw(i_frw  )*(aexp-aexp_frw(i_frw-1))/(aexp_frw(i_frw  )-aexp_frw(i_frw-1))+ &
          & tau_frw(i_frw-1)*(aexp-aexp_frw(i_frw  ))/(aexp_frw(i_frw-1)-aexp_frw(i_frw  ))
  endif

  if(myid==1)write(*,*)'init_time done'

end subroutine init_time