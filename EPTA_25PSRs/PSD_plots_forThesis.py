#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import glob
from scipy import interpolate as interp
from scipy import linalg
from scipy import optimize

from astropy import coordinates as coord
from astropy import units as units
from astropy import constants as aconst

colours = ['blue','green','red','cyan','orange','forestgreen','plum','magenta','lime','purple','cyan',           'pink','darkgreen','rosybrown','maroon',          'magenta','olive','peru','yellow','darkblue',           'fuchsia','crimson']

plt.rc('text',usetex=True)
plt.rc('font', size=18)
plt.rc('font', weight='bold')
plt.rc('font', family='serif')          # controls default text sizes
#plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=40)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=35)    # fontsize of the tick labels
plt.rc('ytick', labelsize=35)    # fontsize of the tick labels
plt.rc('legend', fontsize=20)    # legend fontsize
plt.rc('figure', titlesize=15)  # fontsize of the figure title
plt.rc('ytick', direction='in')
plt.rc('xtick', direction='in')

plt.rc('xtick.major',size=8)
plt.rc('ytick.major',size=8)
plt.rc('xtick.minor',size=6)
plt.rc('ytick.minor',size=6)

current_dir = os.getcwd()


# In[ ]:





# In[2]:


#read residuals file
def read_residuals(path='',resfile='residuals.txt'):
    # read times and residuals:
    times = [] #(days), 
    residuals = [] #(seconds), 
    errors = [] #(seconds)
    
    if path[-1] != '/':
        path = path +'/'

    fopen = open(path+resfile)
    flines = fopen.readlines()

    for line in flines:
        split_line = line.split()
        times.append(float(split_line[0]))
        residuals.append(float(split_line[1]))
        errors.append(float(split_line[2])*1e-6)
 
    times = np.array(times)
    residuals = np.array(residuals)
    errors = np.array(errors)
    return times, residuals, errors


# In[3]:


def read_pars(path='', parfile='final.par'):

    hasRTN = False
    isbinary = False
    hasM2 = False
    if path[-1] != '/':
        path = path +'/'
        
    fopen = open(path+parfile)
    flines = fopen.readlines()
    
    for line in flines:
        split_line = line.split()
        if split_line[0] == 'TNRedAmp':
            logA = (float(split_line[1]))
            hasRTN = True
        elif split_line[0] == 'TNRedGam':
            gamma = (float(split_line[1]))
            hasRTN = (True)
        elif split_line[0] == 'F0':
            F0 = (float(split_line[1]))
        elif split_line[0] == 'F1':
            F1 = (float(split_line[1]))
        elif split_line[0] == 'BINARY':
            BINARY = (True)
            isbinary = True
        elif split_line[0] == 'PB':
            PB = (float(split_line[1]))
        elif split_line[0] == 'M2':
            hasM2 = True
            M2 = (float(split_line[1]))
    if isbinary == False:
        BINARY = (False)
        PB = (None)
    if hasM2 == False:
        M2 = (None)
    if hasRTN == False:
        logA = (np.nan)
        gamma = (np.nan)
    
    return logA, gamma, F0, F1, BINARY, PB, M2


# In[4]:


def read_results(path='', parfile='PSR.par.results'):

    if path[-1] != '/':
        path = path +'/'
        
    fopen = open(path+parfile)
    flines = fopen.readlines()
    
    for line in flines:
        split_line = line.split()
        if 'planet1_period' in split_line[0]:
            PB_maxL = float(split_line[1])
            PB_mean = float(split_line[2])
            PB_std = float(split_line[3])

            
    return PB_maxL, PB_mean, PB_std


# In[5]:


days_per_yr = 365.25
sec_per_year = 365.25*60.*60.*24.
days_per_year = 365.25
sec_per_day = 60.*60.*24.


# In[6]:


def qrfit(dm, b):
    """Least Squares fitting useing the QR decomposition"""
    q, r = np.linalg.qr(dm,mode='reduced')
    p = np.dot(q.T, b)
    
    ri = linalg.inv(r) 
    param = np.dot(linalg.inv(r), p)
    newcvm= ri.dot(ri.T)
    model = param.dot(dm.T)
    postfit = b - model
    chisq = np.sum(np.power(postfit,2))
    newcvm *= (chisq/float(len(b)-len(param)))
    return param,newcvm,chisq


def PSD(freqs_days,fc_yr,logA,gamma):
	# logA is log10(A)
	Amp2 = np.power(np.power(10.,logA),2)
	fc_days = fc_yr/days_per_yr

	psd = 1./(12.* np.power(np.pi,2))* Amp2 *             np.power(np.power((freqs_days/fc_days),2)+1.,-gamma/2.) * (fc_yr)**(-gamma) # yr^3

	return psd

def getC(psr,times,residuals,errors,logA,gamma,fc_yr):
	# note: 'times' can be irregularly sampled; also negative
	total_time = (np.amax(times)-np.amin(times)) # days
	ndays_time = int(total_time+1)
	#fc_days = 1./total_time 
	#fc_days = 0.01/days_per_year
	fc_days = fc_yr/days_per_yr
	# create regular sequence of times..
	#.. to include 'times' after the iFFT
	npts_regtimes = 128
	delta_regtimes = 1 #day
	while npts_regtimes*delta_regtimes < (ndays_time+1)*2 or npts_regtimes*delta_regtimes < (2./fc_days):
		npts_regtimes *= 2

	# this creates an array of times to use in covFunc...
	# ... regtimes =  [0,1,...,npts_regtimes-1]*delta_regtimes
	# ... in days

	# which then gives the FT frequencies...
	# ... [0,1/npts_regtimes,...,(npts_regtimes/2)/npts_regtimes]]*(1/delta_regtimes)
	# ... in days^-1

	nu = np.fft.rfftfreq(npts_regtimes,delta_regtimes) #day^-1
	# len(nu) == npts_regtimes/2+1
	# delta_nu == 1./(npts_regtimes*delta_regtimes)
	delta_nu = (nu[-1] - nu[0])/nu.shape[0] #day^-1

	# calculate power spectral density
	psd_nu = PSD(nu, fc_yr, logA, gamma) # yr^3

	# get the covariance function:
	covFunc = np.fft.irfft(psd_nu*sec_per_year**3/(delta_regtimes*sec_per_day)) #sec^2
	# len(covFunc) == (len(nu)-1)*2 == npts_regtimes

	# need covFunc to work for abs(times[j]-times[i]), any j,i
	# ndays_time-1 <max (abs(times[j]-times[i])) =< ndays_time

	# so interpolate it for the regular times up to ndays_time:
	covFunc = interp.interp1d(np.arange(ndays_time+1), covFunc[:ndays_time+1])

	C = np.empty([len(times),len(times)])
	#alltimes = np.empty([len(times),len(times)])
	for i in range(len(times)):
		try:
			C[i] = covFunc(np.abs(times-times[i]))
			#alltimes[i] = np.abs(times-times[i])
		except ValueError as e:
			print (e)
			print (np.abs(times-times[i]))
	return C
	# have covariance matrix!

def getLikelihood(psr,C,times,residuals,errors):
	C += np.diag(errors**2)
	# Cholesky decomposition:
	L = np.linalg.cholesky(C)

	Linv = linalg.inv(L) # Not efficient, but ok for now

	w = Linv.dot(residuals) # This is the "whitened residuals"

	# First we will fit just a mean value to get the comparative log-likelihood
	# Might want to fit pulsar model here too?
	M = np.zeros((1,len(times))) # Design matrix with 1 parameter

	M[0] = 1.0 # Set all elements to 1.0

	wM = Linv.dot(M.T) # "Whitened" design matrix
	beta, cvm, chisq = qrfit(wM,w) # Solve least-squares for whitened design matrix and whitened data.
	default_ll = -0.5*chisq # Assume log-likelihood is just chi-square

	total_time = (np.amax(times)-np.amin(times)) # days
	avg_delta_times = total_time/(len(times)-1)

	#f_Nyq = 1./(2.*avg_delta_times) # per day; max freq
	#print ('MAX F = ',f_Nyq)
	f_Nyq = 2. # per day
	#delta_f = 1.0/(10.*total_time) # frequency step, per day
	delta_f = 1.0/(total_time) # frequency step, per day
	N_f = int(f_Nyq/delta_f) # number of spectral channels
	f = np.linspace(delta_f,f_Nyq, N_f, endpoint=False) # per day

	# Create the output arrays
	dpower = np.zeros(N_f) 
	loglike = np.zeros(N_f)


	angf = f*2.*np.pi # rad/day
	wt_matrix = np.outer (angf,times)
	sin_wt_matrix = np.sin(wt_matrix)
	cos_wt_matrix = np.cos(wt_matrix)
	# sin and cos have shapes = len(f),len(times)
	print ('chisq = '+ str(chisq))
	print ('redchisq = '+str(chisq/len(times)))
	print ('\n')

	for ifreq in range(len(f)):
		if ifreq%1000 == 0:
			print(ifreq, '/', len(f))
		M = np.zeros((3,len(times))) # Design Matrix of three parameters: sin, cos, mean
		#M[0] = np.sin(omega*times) 
		#M[1] = np.cos(omega*times) # Elements of DM for cos
		M[0] = sin_wt_matrix[ifreq] # Elements of DM for sin
		M[1] = cos_wt_matrix[ifreq] # Elements of DM for cos
		M[2] = 1.0            # Elements of DM for mean value.

		# Whiten and solve the least-squares problem
		wM = Linv.dot(M.T)
		beta, cvm, chisq = qrfit(wM,w)
		ll = -0.5*chisq
		delta_ll = ll - default_ll

		#Output results
		dpower[ifreq] = beta[0]**2 + beta[1]**2 # sec^2
		loglike[ifreq] = delta_ll


	#dpower = psd_f_yr3*delta_f_perday*(sec_per_year**3/sec_per_day) # sec^2
	psd_f_yr3 = dpower*sec_per_day/delta_f/(sec_per_year**3) # yr^3

	return f,psd_f_yr3,loglike


# In[ ]:





# In[7]:


PSRs = []
fopen = open('25PSRs.txt','r')
flines = fopen.readlines()
for line in flines:
    PSRs.append(line.split()[0])
        
# PBs = np.array([1473, 1417, 1102, [231, 498], 319, 350, 1013, 719, 650, 946, 117, 723, 1070, 5180, 28])
PSRs[5]


# In[8]:


flagPSRs = ['J0751+1807', 'J0900-3144', 'J1012+5307', 'J1022+1001', 'J1600-3053', 'J1640+2224', 'J1713+0747',            'J1744-1134', 'J1909-3744', 'J1918-0642']
# len(flagPSRs)


# In[9]:


nonflagPSRs = [x for x in PSRs if x not in flagPSRs]
# len(nonflagPSRs)


# In[10]:


for psr in flagPSRs:
    print(psr)

    parfile = psr+'_post_wfit.par'
    #timfile = 'final.tim'
    resfile = 'residuals_'+psr+'_tndm.txt'

    pre_times, pre_residuals, pre_errors = read_residuals(path='.', resfile=resfile)
    times_MJD = pre_times
    logA, gamma, _, _, _, _, _ = read_pars(path=psr+'/', parfile=parfile)
    fc_yr = 0.01 
    totalT = times_MJD.max() - times_MJD.min()
    fc_days = fc_yr / days_per_yr
    
    dynpath = 'planet_runs/dynesty-run_enterprise/withP/' + psr + '/'
    # _, _, _, _, _, PB, _ = read_pars(path=dynpath_short, parfile='.par')
    # PB, PBs[k], type(PBs[k]) == int
    PB_maxL, PB_mean, PB_std = read_results(path=dynpath, parfile=psr+'.par.results')
    print(PB_maxL, PB_mean, PB_std)
    
    print('Calculating C.........\n')
    C = getC(psr, pre_times, pre_residuals, pre_errors, logA, gamma, fc_yr)
    print('Calculating likelihood.........\n')
    f, psd_f_yr3, loglike = getLikelihood(psr, C, pre_times, pre_residuals, pre_errors)

    PSD_model = PSD(f, fc_days, logA, gamma)
    
    f_yr = sec_per_year*f/sec_per_day #in yr-1
    # PB_maxL in days
    # sec_per_year*f/sec_per_day in yr-1
    # psd_f_yr3
    # PSD_model in yr3
    
    np.savez(psr+'_tndm_forPSDplot.npz', PB_maxL, f_yr, psd_f_yr3, PSD_model)
    


# In[ ]:





# In[ ]:




