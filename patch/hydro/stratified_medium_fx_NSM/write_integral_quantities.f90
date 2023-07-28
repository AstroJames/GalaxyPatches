!! Function to calculate integral quantities and write to file. 
!! Called by adaptive loop, each (?) coarse timestep

subroutine write_integral_quantities(isFirst, simTime)
!! We need access to the uold array 
    use hydro_commons
    use amr_commons

    implicit none

#ifndef WITHOUTMPI
    include 'mpif.h'
!!    integer::info
#endif

    integer, intent(in) :: isFirst
    integer, intent(in) :: simTime

    integer :: n_glob = 25 ! Number of variables to calculate - maybe shouldn't be hardcoded

    ! Grid and level variables for looping over the cells, communicating the with MPI and indexing the grid
    integer:: ilevel, ind, ngrid, iskip, ilev, i, igrid, ncache, info
    integer, dimension(1:nvector), save::ind_grid, ind_cell
    logical, dimension(1:nvector)::ok !to check for leaf cells 

    ! Cell sizes and scale conversions 
    real(dp)::dx, dx_loc, vel_loc, dm
    real(dp)::scale,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
    
    ! Arrays to store the local and global sums
    real(dp), dimension(0:n_glob) :: lsum, gsum 
    real(dp), dimension(1) :: lmin_dens, lmax_dens, gmin_dens, gmax_dens

    ! Counters and arrays for writing to file
    integer:: iq, funit = 99
    integer, parameter :: max_q = 100
    character(len=25), dimension(max_q)::name_out
    real(dp), dimension(max_q)::quant_out
    character(len=25)::tmp_str
    character(len=10000)::header
    character(len=80)::filename


    if(myid==1 .and. verbose)write(*,111) 'Entering write_integral_quantities'
     ! Conversion factor from user units to cgs units
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

!! Quantities to calculate: total mass, totol energy, total momentum, min, max and average gas density, 
!! average velocity (decomposed?), velocity rms, average sound speed and/or average Mach number,
!! average and rms of passive scalar fields, median of passive scalar fields (?)

     lsum = 0.d0

     ! Loop over levels 
     do ilev = levelmin, levelmax
        dx = 0.5D0**ilev !local cell size 
        dx_loc = dx*scale ! cgs units 
        vol_loc = dx_loc**ndim ! cell volume 

        ! loop over grids 
        ncache = active(ilev)%ngrid
        do igrid =1,ncache,nvector 
            ngrid=MIN(nvector, ncache-igrid+1)
            do i=1,ngrid
                ind_grid(i)=active(ilev)%igrid(igrid+i-1)
            end do

            ! Loop over cells 
            do i = 1, twotondim
                iskip = ncoarse+(ind-1)*ngridmax
                do i =1, ngrid
                    ind_cell(i) = iskip+ind_grid(i)
                end do 

                ! Flag leaf cells 
                do i=1, ngrid
                    ok(i)=son(ind_cell(i))==0
                end do 

                do i, ngrid
                    if(ok(i))then ! Only add values from leaf cells 
                        dm = uold(ind_cell(i))*scale_d*vol_loc ! call mass 

                        lsum(0) = lsum(0) + vol_loc ! volume
                        lsum(1) = lsum(1) + dm ! mass
                        ! Momenta along x, y, z:
                        lsum(2) = lsum(2) + uold(ind_cell(i), 2)*scale_v*dm 
                        lsum(3) = lsum(3) + uold(ind_cell(i),3)*scale_v*dm
                        lsum(4) = lsum(4) + uold(ind_cell(i),4)*scale_v*dm 

                        ! Energies 
                        !lsum(5) = lsum(5) + ! Internal 
                        ! Kinetic 
                        ! Potential 
                    endif
                end do ! Loop grids (?)
            end do ! Loop on cells 
        end do ! Loop on grids
    end do ! Loop on levels 

!! Communicate local sums to the root cpu (MPI_REDUCE, not ALLREDUCE) 
#ifndef WITHOUTMPI
    gsum = 0.0d0
    call MPI_Reduce(lsum, gsum, n_glob+1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD)
#endif
!! Calculate the average and 
!! Create and write to file 
!! If this is the very first step we need to create the file, otherwise we want to just write to it. 
if (myid == 1) then !Only root cpu writes to file 
    iq = 0 
    iq = iq+1; name_out(iq) = 'time';           quant_out(iq) = simTime
    iq = iq+1; name_out(iq) = 'gas_mass';       quant_out(iq) = gsum(1)
    iq = iq+1; name_out(iq) = 'momentum_x';     quant_out(iq) = gsum(2)
    iq = iq+1; name_out(iq) = 'momentum_y';     quant_out(iq) = gsum(3)
    iq = iq+1; name_out(iq) = 'momentum_z';     quant_out(iq) = gsum(4)


    filename = 'global_quantities.txt'
    ! Check if the file already exsisits
    open(funit, file=trim(filename), position = 'APPEND', status = 'OLD', iostat = istat)
    if (istat .NE. 0) then ! File didn't exsists already 
        open(funit, file = trim(filename), position = 'APPEND')
    endif
    if (isFirst .EQ. 1) ! .AND. (nrestart .NE. 0 .or. istart .NE. 0)) then 
        ! Create the header 
        header = ''
        do i = 1, iq
            write(tmp_str, '(I2.2)') i ! column number 
            write(tmp_str, '(A24)') '#'//trim(tmp_str)//'_'//trim(name_out(i)) ! quantity name 
            header = trim(header)//trim(tmp_str)
        enddo 
        write(funit, '(A)') trim(header) ! Write the header to file 
    endif
    if (nrestart .NE. 0) then !it's a restart 
        write(funit, '(A)') '# Simulation restarted'
    endif
    write(tmp_str, '(I3)') iq ! how many quantities are we writing 
    write(funit, '('//trim(tmp_str)//'(IX,ES23.16))') quant_out(1:iq)

    close(funit) ! close the file 
endif ! Root cpu 

!! How to deal with restarts from an older checkpoint? 

end subroutine write_integral_quantities