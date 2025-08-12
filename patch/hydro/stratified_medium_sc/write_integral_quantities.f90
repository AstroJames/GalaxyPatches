!! Function to calculate integral quantities and write to file. 
!! Called by adaptive loop, each (?) coarse timestep
subroutine write_integral_quantities(isFirst, simTime)
!! We need access to the uold array 
    use hydro_commons
    use amr_commons
    use poisson_commons

    implicit none

#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    logical, intent(inout) :: isFirst
    real(dp), intent(in) :: simTime

    ! Grid and level variables for looping over the cells, communicating the with MPI and indexing the grid
    integer:: ilevel, ind, ngrid, iskip, ilev, i, j, igrid, ncache, info, ierr, iz
    integer, dimension(1:nvector), save::ind_grid, ind_cell
    logical, dimension(1:nvector)::ok !to check for leaf cells 

    ! Cell sizes and scale conversions 
    real(dp)::dx, dx_loc, vol_loc, dm, nx_loc, u_sq, dens, z, skip_loc, P, c_s, Int_E
    real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
    real(dp),dimension(1:twotondim)::xc

    ! Counters and arrays for writing to file
    integer:: iq, funit = 99, istat
    integer, parameter :: max_q = 100
    character(len=25), dimension(max_q)::name_out
    real(dp), dimension(max_q)::quant_out
    character(len=25)::tmp_str
    character(len=10000)::header
    character(len=80)::filename

    ! Arrays to store the local and global sums
    real(dp), dimension(1:max_q) :: lsum, gsum 
    real(dp) :: lmin_dens, lmax_dens, lmin_dens_floor, gmin_dens, gmax_dens, gmin_dens_floor


    if(verbose)write(*,*) 'Entering write_integral_quantities'
     ! Conversion factor from user units to cgs units (NB user units =/= code units!!!)
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
 
    nx_loc =(icoarse_max - icoarse_min+1)
    scale = boxlen/nx_loc
    skip_loc=dble(kcoarse_min)

!! Quantities to calculate: total mass, totol energy, total momentum, min, max and average gas density, 
!! average velocity (decomposed?), velocity rms, average sound speed and/or average Mach number,
!! average and rms of passive scalar fields, median of passive scalar fields (?)

     lsum = 0.d0
     lmin_dens = HUGE(1.0)
     lmin_dens_floor = smallr*scale_d
     lmax_dens = TINY(0.0)

     ! Loop over levels 
     do ilev = levelmin, nlevelmax
        dx = 0.5D0**ilev !local cell size 
        dx_loc = dx*scale*scale_l ! in cgs (scale goes from code to user, scale_l goes from user to cgs)
        vol_loc = dx_loc**ndim ! cell volume 

        ! Cell center positions relative to grid center positions (only need z)
        do ind=1,twotondim
            iz=(ind-1)/4
            xc(ind)=(dble(iz)-0.5D0)*dx
         end do

        ! loop over grids 
        ncache = active(ilev)%ngrid
        do igrid =1,ncache,nvector 
            ngrid=MIN(nvector, ncache-igrid+1)
            do i=1,ngrid
                ind_grid(i)=active(ilev)%igrid(igrid+i-1)
            end do

            ! Loop over cells
            !write(*,*)'Preparing to loop over cells, twotondim =', twondim
            do ind = 1, twotondim
                iskip = ncoarse+(ind-1)*ngridmax
                do i =1, ngrid
                    ind_cell(i) = iskip+ind_grid(i)
                end do 

                ! Flag leaf cells 
                do i=1, ngrid
                    ok(i)=son(ind_cell(i))==0
                end do 
                !write(*,*)'Flagged the leaf cells, there are:', count(ok), ' leaf cells at this level'
                do i = 1, ngrid
                    if(ok(i))then ! Only add values from leaf cells 
                        dm = uold(ind_cell(i), 1)*scale_d*vol_loc ! cell mass 
                        dens = uold(ind_cell(i),1)*scale_d

                        lsum(1) = lsum(1) + vol_loc ! volume
                        lsum(2) = lsum(2) + dm ! mass

                        ! Momenta along x, y, z:
                        do j = 1,3
                            lsum(2+j) = lsum(2+j) + uold(ind_cell(i),j+1)*scale_v*scale_d*vol_loc !uold = conserved quantities
                        end do

                            ! Density 
                        lmin_dens = MIN(lmin_dens, dens)
                        lmax_dens = MAX(lmax_dens, dens)
                        
                        ! compute min dens with floor 
                        lmin_dens_floor = MAX(lmin_dens_floor, dens)

                        lsum(6) = lsum(6) + dens*vol_loc
                        lsum(7) = lsum(7) + dens**2*vol_loc

                        ! Energies 
                        ! Velocities
                        u_sq = 0.0
                        do j=1,3 ! Loop on velocities
                            u_sq = u_sq + ((uold(ind_cell(i), j+1)/uold(ind_cell(i),1))*scale_v)**2
                        enddo 

                        ! Cell center position along z
                        z=(xg(ind_grid(i),3)+xc(ind)-skip_loc)*scale*scale_l

                        lsum(8) = lsum(8) + 0.5*dm*u_sq ! Kinetic
                        lsum(9) = lsum(9) + f(ind_cell(i),3)*(scale_l/scale_t**2)*dm*z ! Potential
                        lsum(10) = lsum(10) + uold(ind_cell(i),5)*scale_d*scale_v**2 ! Internal

                        P = uold(ind_cell(i), 5)

                        do j = 1,3 ! subtract kinetic energy
                            lsum(10) = lsum(10) - 0.5*(uold(ind_cell(i), j+1)**2/uold(ind_cell(i),1))*scale_d*scale_v**2
                            P = P - (0.5*uold(ind_cell(i), j+1)**2/uold(ind_cell(i),1))
                        enddo
                        ! Internal energy: total energy minus kinetic energy
                        lsum(10) = lsum(10) + (uold(ind_cell(i),5) - 0.5*(uold(ind_cell(i),2)**2 + uold(ind_cell(i),3)**2 + uold(ind_cell(i),4)**2)/uold(ind_cell(i),1))*scale_d*scale_v**2

                        P = uold(ind_cell(i), 5) - (0.5*(uold(ind_cell(i),2)**2 + uold(ind_cell(i),3)**2 + uold(ind_cell(i),4)**2)/uold(ind_cell(i),1))
                        P = (gamma - 1.0d0)*P*scale_d*scale_v**2 ! P = 2/3*e

                        lsum(11) = lsum(11) + uold(ind_cell(i), 5)*scale_d*scale_v**2 ! Total fluid energy

                        !vx, vy, vz (average + rms)
                        do j = 1,3
                            lsum(11+j) = lsum(11+j) + ((uold(ind_cell(i),j+1)/uold(ind_cell(i),1))*scale_v)*vol_loc
                            lsum(14+j) = lsum(14+j) + ((uold(ind_cell(i),j+1)/uold(ind_cell(i),1))*scale_v)**2*vol_loc
                        enddo

                        ! Velocity squared
                        lsum(18) = lsum(18) + u_sq*vol_loc
                        
                        ! Mach number 
                        c_s = sqrt(gamma*P/dens)

                        lsum(19) = lsum(19) + c_s*vol_loc
                        lsum(20) = lsum(20) + (sqrt(u_sq)/c_s)*vol_loc ! Mach number
                        lsum(21) = lsum(21) + u_sq/c_s**2*vol_loc ! rms Mach number (u/c_s)**2 = u_sq/c_s**2


                        ! Metalicities
                        do j = 1,3 ! Loop on passive scalars 
                            lsum(21+j) = lsum(21+j) + uold(ind_cell(i), 5+j)*vol_loc
                            lsum(24+j) = lsum(24+j) + uold(ind_cell(i), 5+j)**2*vol_loc
                        enddo 

                        lsum(28) = lsum(28) + sqrt(u_sq)*vol_loc
                        lsum(29) = lsum(29) + P*vol_loc

                    endif
                end do ! Loop grids (?)
            end do ! Loop on cells 
        end do ! Loop on grids
    end do ! Loop on levels 

!! Communicate local sums to the root cpu (MPI_REDUCE, not ALLREDUCE) 
#ifndef WITHOUTMPI
    gsum = 0.0d0
    call MPI_REDUCE(lsum, gsum, max_q, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    gmin_dens = 0.0d0
    call MPI_REDUCE(lmin_dens, gmin_dens, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    gmax_dens = 0.0d0
    call MPI_REDUCE(lmax_dens, gmax_dens, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    gmin_dens_floor = 0.0d0
    call MPI_REDUCE(lmin_dens_floor, gmin_dens_floor, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
#endif
!! Calculate the average and 
!! Create and write to file 
!! If this is the very first step we need to create the file, otherwise we want to just write to it. 
if (myid == 0) then !Only root cpu writes to file 
    iq = 0 
    iq = iq+1; name_out(iq) = 'time';           quant_out(iq) = simTime*scale_t
    iq = iq+1; name_out(iq) = 'gas_mass';       quant_out(iq) = gsum(2)
    iq = iq+1; name_out(iq) = 'momentum_x';     quant_out(iq) = gsum(3)
    iq = iq+1; name_out(iq) = 'momentum_y';     quant_out(iq) = gsum(4)
    iq = iq+1; name_out(iq) = 'momentum_z';     quant_out(iq) = gsum(5)

    iq = iq+1; name_out(iq) = 'ave_density';    quant_out(iq) = gsum(6)/gsum(1)
    iq = iq+1; name_out(iq) = 'rms_density';    quant_out(iq) = sqrt(gsum(7)/gsum(1))
    iq = iq+1; name_out(iq) = 'min_density';    quant_out(iq) = gmin_dens
    iq = iq+1; name_out(iq) = 'min_floor_dens'; quant_out(iq) = gmin_dens_floor
    iq = iq+1; name_out(iq) = 'max_density';    quant_out(iq) = gmax_dens

    iq = iq+1; name_out(iq) = 'Kin_energy';     quant_out(iq) = gsum(8)
    iq = iq+1; name_out(iq) = 'Pot_energy';     quant_out(iq) = gsum(9)
    iq = iq+1; name_out(iq) = 'Int_energy';     quant_out(iq) = gsum(10)
    iq = iq+1; name_out(iq) = 'Tot_energy';     quant_out(iq) = gsum(11)

    iq = iq+1; name_out(iq) = 'velocity_x';     quant_out(iq) = gsum(12)/gsum(1)
    iq = iq+1; name_out(iq) = 'velocity_y';     quant_out(iq) = gsum(13)/gsum(1)
    iq = iq+1; name_out(iq) = 'velocity_z';     quant_out(iq) = gsum(14)/gsum(1)
    iq = iq+1; name_out(iq) = 'rms_vx';         quant_out(iq) = sqrt(gsum(15)/gsum(1))
    iq = iq+1; name_out(iq) = 'rms_vy';         quant_out(iq) = sqrt(gsum(16)/gsum(1))
    iq = iq+1; name_out(iq) = 'rms_vz';         quant_out(iq) = sqrt(gsum(17)/gsum(1))
    iq = iq+1; name_out(iq) = 'ave_u_sq';       quant_out(iq) = gsum(18)/gsum(1)
    iq = iq+1; name_out(iq) = 'ave_u';          quant_out(iq) = gsum(28)/gsum(1)
    iq = iq+1; name_out(iq) = 'rms_u';          quant_out(iq) = sqrt(gsum(18)/gsum(1))
    iq = iq+1; name_out(iq) = 'ave_c_s';        quant_out(iq) = gsum(19)/gsum(1)
    iq = iq+1; name_out(iq) = 'Mach_numb';      quant_out(iq) = gsum(20)/gsum(1)
    iq = iq+1; name_out(iq) = 'rms_Mach_numb';  quant_out(iq) = sqrt(gsum(21)/gsum(1))

    iq = iq+1; name_out(iq) = 'ave_P';          quant_out(iq) = gsum(29)/gsum(1)

    iq = iq+1; name_out(iq) = 'Z_ave';          quant_out(iq) = gsum(22)/gsum(1)
    iq = iq+1; name_out(iq) = 'alpha_ave';      quant_out(iq) = gsum(23)/gsum(1)
    iq = iq+1; name_out(iq) = 'rp_ave';         quant_out(iq) = gsum(24)/gsum(1)
    iq = iq+1; name_out(iq) = 'Z_rms';          quant_out(iq) = sqrt(gsum(25)/gsum(1))
    iq = iq+1; name_out(iq) = 'alpha_rms';      quant_out(iq) = sqrt(gsum(26)/gsum(1))
    iq = iq+1; name_out(iq) = 'rp_rms';         quant_out(iq) = sqrt(gsum(27)/gsum(1))    

    filename = 'global_quantities.txt'
    ! Check if the file already exsisits
    open(funit, file=trim(filename), position = 'APPEND', status = 'OLD', iostat = istat)

    ! Check if the file already exists
    open(funit, file=trim(filename), position = 'APPEND', status = 'OLD', iostat = istat)

    if (istat .NE. 0) then ! File didn't exist already 
        open(funit, file = trim(filename), position = 'APPEND')
    endif

    if ((isFirst) .and. (nrestart == 0 .or. istat .NE. 0)) then ! it's not a restart or if the file did not already exist
        ! Create the header 
        header = ''
        do i = 1, iq
            write(tmp_str, '(A)') trim(name_out(i))//',' !quantity name 
            header = trim(header)//trim(tmp_str)
        enddo 

        write(funit, '(A)') trim(header) ! Write the header to file 
    else if (isFirst .and. (nrestart .NE. 0)) then ! it's a restart
        write(funit, '(A)') '# Simulation restarted'
    endif
    
    write(tmp_str, '(I3)')iq-1
    write(funit, '(1x,'//trim(tmp_str)//'(ES20.5, ","), ES20.5)')quant_out(1:iq)
    close(funit) ! close the file 

endif ! Root cpu 

! set isFirst to false to avoid repeating the header 
isFirst = .false. 
end subroutine write_integral_quantities
