# Generate log-normal density field with prescribed power-law spatial
# power spectrum. Following Lewins & Austin (2002).
#
# http://faculty.uoit.ca/lewis/pubs/lewis_lognorm.pdf

import random
import numpy
import sys
import scipy
import array as array_mod
from math import *
from numpy import *
from scipy.interpolate import *
from scipy.stats import *

# Compute normalization factor for power spectrum; assuming 0 mean.
def norm_Pk(Pk_unnorm, nx, ndim, sigma2_):
  if ndim==1:
    B = nx*(nx-1)*sigma2_
    den = 0.0

    for k in range(1, nx):
      den = den + Pk_unnorm(k)

    B = B/den

  elif ndim==2:
    B = nx**2.0*(nx**2.0-1)*sigma2_
    den = 0.0

    for kx in range(0, nx):
      for ky in range(0, nx):
	if not (kx, ky)==(0, 0):
	  k = sqrt(kx**2.0 + ky**2.0)
	  den = den + Pk_unnorm(k)

    B = B/den

  elif ndim==3:
    B = nx**3.0*(nx**3.0-1)*sigma2_ 
    den = 0.0

    for kx in range(0, nx):
      for ky in range(0, nx):
        for kz in range(0, nx):
	  if not (kx, ky, kz)==(0, 0, 0):
	    k = sqrt(kx**2.0 + ky**2.0 + kz**2.0)
	    den = den + Pk_unnorm(k)

    B = B/den

  return B

def make_gauss(nx, ndim, mean_Gauss, sigma2_Gauss, Pk_unnorm):
  # Fix seed of random number generator.
  numpy.random.seed(13)

  
  ### First, create a Gaussian random field, following Peacock, p. 504 ###
  kspace_arr = numpy.zeros((nx,)*ndim, numpy.complex_)
 
  # Assume mean = 0 for now.


  # Normalize power spectrum.
  B = norm_Pk(Pk_unnorm, nx, ndim, sigma2_Gauss)
  

  Pk_norm = lambda k: B*Pk_unnorm(k)

  # Choose k modes according to Rayleigh distribution, using above
  # power spectrum to define variance in each case.
  if ndim==1:
    for kx in range(0, nx):
      if not kx==0:
        k = kx

        # Assign random phase between 0 and 2*pi.
        Pk_norm_ = float(Pk_norm(k))
	
	if Pk_norm_>0.0:
	  Ak_norm = numpy.random.rayleigh(sqrt(Pk_norm_/2.0))
	  Ak_phase = 2.0*pi*numpy.random.random()
	else:
	  Ak_norm = 0.0
          Ak_phase = 0.0

	kspace_arr[kx] = complex(Ak_norm*cos(Ak_phase), Ak_norm*sin(Ak_phase))

  elif ndim==2:
    for kx in range(0, nx):
      for ky in range(0, nx):
	if not (kx, ky)==(0, 0):
	  k = sqrt(kx**2.0 + ky**2.0)

          Pk_norm_ = float(Pk_norm(k))
	  
	  if Pk_norm_>0.0:
 	    Ak_norm = numpy.random.rayleigh(sqrt(Pk_norm_/2.0))  
	    Ak_phase = 2.0*pi*numpy.random.random()
	  else: 
  	    Ak_norm = 0.0
            Ak_phase = 0.0
	  
	  kspace_arr[kx][ky] = complex(Ak_norm*cos(Ak_phase), Ak_norm*sin(Ak_phase))
  
  elif ndim==3:
    for kx in range(0, nx):
      for ky in range(0, nx):
        for kz in range(0, nx):
	  if not (kx, ky, kz)==(0, 0, 0):
	    k = sqrt(kx**2.0 + ky**2.0 + kz**2.0)

            Pk_norm_ = float(Pk_norm(k))

	    if Pk_norm_>0.0:
 	      Ak_norm = numpy.random.rayleigh(sqrt(Pk_norm_/2.0))
	      Ak_phase = 2.0*pi*numpy.random.random()
	    else:
	      Ak_norm = 0.0
              Ak_phase = 0.0    
	    
	    kspace_arr[kx][ky][kz] = complex(Ak_norm*cos(Ak_phase), Ak_norm*sin(Ak_phase))

    # Case of (kx, ky, kz)=0 (DC mode). Assume no DC power, i.e. leave
    # 0.
  else:
    print "code missing for specified ndim=",
    print ndim
    sys.exit()


  # Inverse Fourier transform to get Gaussian random field with
  # specified power spectrum.
  xspace_arr = numpy.fft.ifftn(kspace_arr)

  
  # Add desired mean.
  xspace_arr = xspace_arr + mean_Gauss

  return xspace_arr

if __name__=="__main__":
  # Parameters of the numerical grid.
  nx = 128 # number of grid points along each dimension
  ndim = 3 # number of dimensions.


  # New or restart?
  restart = 0
  #restart_file = "/clusterfs/henyey/cgiguere/lognorm_nx_512_not_conv.dat"
  restart_file = ""

  #save_file = "/clusterfs/henyey/cgiguere/lognorm_nx_" + str(nx) + ".dat"
  #save_file = "/clusterfs/henyey/cgiguere/lognorm_nx_" + str(nx) + "_Mach10_beta_2.dat"
  save_file = "/indirect/o/dav.martizzi/lognorm_nx_" + str(nx) + "_Mach30_beta_2_kmax_128.dat"
  #save_file = "/clusterfs/henyey/cgiguere/lognorm_nx_" + str(nx) + "_Mach30_beta_1.dat"
  #save_file = "/clusterfs/henyey/cgiguere/lognorm_nx_" + str(nx) + "_Mach10_beta_1.dat"


  # Desired Mach number and slope of the power spectrum.
  #Mach = 10.0
  Mach = 30.0

  #beta = 5.0/3.0 # Kolmogorov spectrum for incompressible turbulence
  beta = 2.0
  #beta = 1.0


  ### Procedure to get lognormal, from Lewis & Austin ###

  # Lemaster & Stone (2009), eq. 5, fitting formula for mean of
  # ln(rho/rho_mean) distribution for driven supersonic
  # turbulence.
  #
  # Note, this is really the mean of the Gaussian distribution.
  mu_V = -0.36*log(1.0 + 0.5*Mach**2.0) + 0.10 
  mean_Gauss = mu_V

  # Corresponding variance (coming from constraint that expected value
  # of rho/rho_mean must be 1).
  sigma2_Gauss = 2.0*abs(mu_V)


  # Variance of lognorm.
  sigma2_lognorm = exp(2.0*mean_Gauss + sigma2_Gauss)*(exp(sigma2_Gauss) - 1.0)


  # Maximum wavenumber with non-zero power
  kmax = nx/4


  # Target power spectrum of lognormal density field, before
  # normalization.
  Pk_target_unnorm = lambda k: where(k<=kmax, k**-beta, 0.0)

  # First guess for power spectrum in *linear* space.
  Pk_lin_unnorm = lambda k: where(k<=kmax, k**-beta, 0.0)


  # Normalize power spectra.
  B_target = norm_Pk(Pk_target_unnorm, nx, ndim, sigma2_lognorm)
  Pk_target = lambda k: B_target*Pk_target_unnorm(k)
  
  B_lin = norm_Pk(Pk_lin_unnorm, nx, ndim, sigma2_Gauss) 
  Pk_lin = lambda k: B_lin*Pk_lin_unnorm(k)


  for j in range(0, 100):
    if restart==1 and j==0:
      # Restart from previous snapshot.
   
      # Read in lognormal density field from binary.
      f = open(restart_file, 'rb')
      binvalues = array_mod.array('f')
      binvalues.read(f, nx**3)

      lognorm_trial_arr = numpy.array(binvalues, dtype=numpy.float_)
      lognorm_trial_arr = numpy.reshape(lognorm_trial_arr, (nx, nx, nx))

      f.close()

    else:
      # Start from scratch, or continue with loop.

      # Generate complex Gaussian random field, with twice the variance
      # desired for the real Gaussian random field.
      gauss_trial_arr = make_gauss(nx, ndim, mean_Gauss, 2.0*sigma2_Gauss, Pk_lin)

      # Exponentiate real part of Gaussian to get trial lognormal.
      lognorm_trial_arr = numpy.exp(gauss_trial_arr.real)


    # DEBUG
    #print "numpy.mean(lognorm_trial_arr)=",
    #print numpy.mean(lognorm_trial_arr)


    ## Calculate power spectrum of trial lognormal.
    
    # FFT.
    lognorm_FFT_arr = numpy.fft.fftn(lognorm_trial_arr)


    # Assume isotropy, approximate <|Ak|^2> for k=1, ..., nx**(ndim/2.0)
    #
    # nk_arr[i] stores # of k-space samples with i+0.5=<k<i+1.5
    k_arr = numpy.arange(1, int(sqrt(ndim)*nx))
    nk_arr = numpy.zeros(int(sqrt(ndim)*nx)-1, numpy.float_)
    Pk_mean_arr = numpy.zeros(len(nk_arr), numpy.float_)

    if ndim==1:
      for kx in range(0, nx):
      	if not kx==0:
    	  k = kx
          Ak = lognorm_FFT_arr[kx]
          Pk_ = (Ak.real)**2.0 + (Ak.imag)**2.0
   
          k_bin_i = int(round(k))-1

          nk_arr[k_bin_i] = nk_arr[k_bin_i] + 1
          Pk_mean_arr[k_bin_i] = Pk_mean_arr[k_bin_i] + Pk_

    elif ndim==2:
      for kx in range(0, nx):
        for ky in range(0, nx):
    	  if not (kx, ky)==(0, 0):
    	    k = sqrt(kx**2.0 + ky**2.0)
            Ak = lognorm_FFT_arr[kx][ky]
            Pk_ = (Ak.real)**2.0 + (Ak.imag)**2.0
   
            k_bin_i = int(round(k))-1

            nk_arr[k_bin_i] = nk_arr[k_bin_i] + 1
            Pk_mean_arr[k_bin_i] = Pk_mean_arr[k_bin_i] + Pk_

    elif ndim==3:
      for kx in range(0, nx):
        for ky in range(0, nx):
          for kz in range(0, nx):
    	    if not (kx, ky, kz)==(0, 0, 0):
    	      k = sqrt(kx**2.0 + ky**2.0 + kz**2.0)
              Ak = lognorm_FFT_arr[kx][ky][kz]
              Pk_ = (Ak.real)**2.0 + (Ak.imag)**2.0
              
	      # DEBUG
	      #print "Pk_=",
	      #print Pk_

              k_bin_i = int(round(k))-1

              nk_arr[k_bin_i] = nk_arr[k_bin_i] + 1
              Pk_mean_arr[k_bin_i] = Pk_mean_arr[k_bin_i] + Pk_

    Pk_mean_arr = Pk_mean_arr/nk_arr

    #print "nk_arr=",
    #print nk_arr

    #print "Pk_mean_arr="
    #print Pk_mean_arr


    # Fit polynomial to trial Pk vs. k, for lognormal field.
    deg = 5
    #deg = 7
    Pcoeffs = numpy.polyfit(numpy.log(k_arr), numpy.log(Pk_mean_arr), deg)

    
    # Polynomial function.
    Pfunc = lambda k: Pcoeffs[5] + Pcoeffs[4]*k + Pcoeffs[3]*k**2.0 + Pcoeffs[2]*k**3.0 + Pcoeffs[1]*k**4.0 + Pcoeffs[0]*k**5.0
    #Pfunc = lambda k: Pcoeffs[7] + Pcoeffs[6]*k + Pcoeffs[5]*k**2.0 + Pcoeffs[4]*k**3.0 + Pcoeffs[3]*k**4.0 + Pcoeffs[2]*k**5.0 + Pcoeffs[1]*k**6.0 + Pcoeffs[0]*k**7.0


    # Define search direction.
    dk_func = lambda k: Pk_target(k) - numpy.exp(Pfunc(numpy.log(k)))
   

    
    # Compare power spectrum at previous iteration to power spectrum
    # after interation, to check convergence.
    if j>0:
      print "numpy.mean(abs((Pk_mean_arr - Pk_mean_arr_prev)/Pk_mean_arr_prev))=",
      print numpy.mean(abs((Pk_mean_arr - Pk_mean_arr_prev)/Pk_mean_arr_prev))
 

    # Store power spectrum of current iteration.
    Pk_mean_arr_prev = Pk_mean_arr


    # Set energy spectrum for the next iterate Gaussian field.
   
    # Delta_step is some dimensionless number <1 definining the step
    # size.
    if j<1000:
      Delta_step = 0.5
    else:
      Delta_step = 0.1

    # Create an interpolation function for linear-space power spectrum
    # to use in next iteration.
    Pk_lin_prev = Pk_lin
 
    Pk_lin_next_arr = Pk_lin_prev(k_arr) + Delta_step*dk_func(k_arr)
    Pk_lin_floor = 10.0**-10.0*max(Pk_lin_next_arr)

    #print "Pk_lin_floor=",
    #print Pk_lin_floor

    Pk_lin_next_arr = numpy.where(Pk_lin_next_arr>0.0, Pk_lin_next_arr, Pk_lin_floor)
    
    Pk_lin = interp1d(k_arr, Pk_lin_next_arr, copy=True, bounds_error=False, fill_value=Pk_lin_floor)

    print "iteration complete"


    # Write array to disk.
    print "writing array to disk...",

    outvalues = array_mod.array('f')
    outvalues.fromlist(lognorm_trial_arr.flatten().tolist())
  
    f = open(save_file, "wb")
    outvalues.tofile(f)
    f.close()

    print "done"
