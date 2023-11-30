#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import emcee
import argparse
import sys
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.offsetbox import AnchoredText


# ------ PARSE ARGUMENTS:


parser=argparse.ArgumentParser(description="Find planet mass limits for a single pulsar and a single run.")

# --- SPECIFIC ARGS:
parser.add_argument('PSR', help='Name of the pulsar to analyse')
parser.add_argument('--no-enterprise',dest='ent',default=True,action='store_false', help='Disable running enterprise search')
parser.add_argument('--no-copydata',dest='copydata',default=True,action='store_false', help='Disable copying data files .par and .tim')
parser.add_argument('--machine',default='godrevy', help='Machine to copy (rsync) data from. Can also be set to \'local\'. Default is \'godrevy\'.')
parser.add_argument('--datapath', help='Path where data is on the remote machine. Default is \'/scratch/mkeith/jbo_dataset/4_ready_for_timing/+$PSR+/analysis/\'.')
parser.add_argument('--par', help='Name of the .par file. Default is \'final.par\'.')
parser.add_argument('--tim', help='Name of the .tim file. Default is \'final.tim\'.')
parser.add_argument('--plot-periodmass',default=False,action='store_true', help='Option to plot period vs mass samples.')
parser.add_argument('--plot-allmass',default=False,action='store_true', \
			help=r'Option to plot the mass distribution, also showing the 95%% limit.')
parser.add_argument('--plot-masslim',default=False,action='store_true', \
			help='Option to plot the 95%% mass limit against period.')
parser.add_argument('--folder', help='Set the folder to save results in. Default is local.')

# --- enterprise ARGS:
#parser.add_argument('par')
#parser.add_argument('tim')
parser.add_argument('--f2', type=float, default=0, help='ENTERPRISE: Range of f2 to search')
parser.add_argument('-N','--nsample', type=float, default=1e6, help='ENTERPRISE: Number of samples in MCMC')
parser.add_argument('-D','--dm',action='store_true', help='ENTERPRISE: Enable DM variation search')
parser.add_argument('--Adm-max',type=float,default=-12,help='ENTERPRISE: Max log10A_DM')
parser.add_argument('--Adm-min',type=float,default=-18,help='ENTERPRISE: Min log10A_DM')
parser.add_argument('--dm-prior-log',action='store_true',help='ENTERPRISE: Use uniform prior in log space for DM noise amplitude')
parser.add_argument('--dm-ncoeff',type=int,default=None,help='ENTERPRISE: Number of DM noise coefficients (nC)')
parser.add_argument('--no-red-noise',dest='red',default=True,action='store_false', help='ENTERPRISE: Disable Power Law Red Noise search')
parser.add_argument('--jbo','-j',action='store_true',help='ENTERPRISE: Use -be flag for splitting backends')
parser.add_argument('--be-flag','-f',help='ENTERPRISE: Use specified flag for splitting backends')
parser.add_argument('--Ared-max','-A',type=float,default=-12, help='ENTERPRISE: Max log10A_Red')
parser.add_argument('--Ared-min',type=float,default=-18, help='ENTERPRISE: Min log10A_Red')
parser.add_argument('--red-gamma-max',type=float,default=8,help='ENTERPRISE: Max gamma red')
parser.add_argument('--red-gamma-min',type=float,default=0,help='ENTERPRISE: Min gamma red')
parser.add_argument('--red-prior-log',action='store_true',help='ENTERPRISE: Use uniform prior in log space for red noise amplitude')
parser.add_argument('--red-ncoeff',type=int,default=60,help='ENTERPRISE: Number of red noise coefficients (nC)')
parser.add_argument('-n','--no-sample',dest='sample',default=True,action='store_false', help='ENTERPRISE: Disable the actual sampling...')
parser.add_argument('--no-white',dest='white',default=True,action='store_false', help='ENTERPRISE: Disable EFAC and EQUAD.')
parser.add_argument('--white-prior-log',action='store_true',help='ENTERPRISE: Use uniform prior in log space for EQUAD.')
parser.add_argument('--efac-max',type=float, default=5, help='ENTERPRISE: Max for efac prior')
parser.add_argument('--efac-min',type=float, default=0.2, help='ENTERPRISE: Min for efac prior')
parser.add_argument('--ngecorr',action='store_true', help='ENTERPRISE: Add ECORR for the nanograv backends.')
parser.add_argument('--white-corner',action='store_true', help='ENTERPRISE: Make the efac/equad corner plots.')
parser.add_argument('--all-corner',action='store_true', help='ENTERPRISE: Make corner plots with all params.')
parser.add_argument('--pm',action='store_true', help='ENTERPRISE: Fit for PMRA+PMDEC.')
parser.add_argument('--px',action='store_true', help='ENTERPRISE: Fit for parallax.')
parser.add_argument('--px-range',type=float,default=10, help='ENTERPRISE: Max parallax to search')
parser.add_argument('--px-verbiest',action='store_true', help=argparse.SUPPRESS) #help='Use Verbiest PX prior to correct for L-K bias')
parser.add_argument('--s1400',type=float,default=None, help=argparse.SUPPRESS)  #help='S1400, used for Verbiest PX prior')
parser.add_argument('--pm-angle',action='store_true', help='ENTERPRISE: Fit for PM + angle.')
parser.add_argument('--pm-range',type=float,default=10,help='ENTERPRISE: Search range for proper motion (deg/yr).')
parser.add_argument('--models','-M',nargs='+',help='ENTERPRISE: Add a model to model selection stack. Use e.g. --model "-D --f2 --no-red-noise".')
parser.add_argument('--outdir','-o',type=str,help="ENTERPRISE: Output directory for chains etc.")
parser.add_argument('--plot-chain',action='store_true', help='ENTERPRISE: Make a plot of the chains.')
parser.add_argument('--tspan-mult',type=float,default=2,help='ENTERPRISE: Multiplier for tspan.')
parser.add_argument('--pm-ecliptic',action='store_true', help='ENTERPRISE: Generate ecliptic coords from pmra/pmdec.')

parser.add_argument('--quasiperiodic','-Q', action='store_true',help='ENTERPRISE: Fit quasiperiodic (QP) model.')
parser.add_argument('--Aqp-max',type=float,default=1,help='ENTERPRISE: Max log10A_QP')
parser.add_argument('--Aqp-min',type=float,default=-4,help='ENTERPRISE: Min log10A_QP')
parser.add_argument('--qp-prior-log',action='store_true',help='ENTERPRISE: Use uniform prior in log space for QP amplitude.')
parser.add_argument('--glitch-recovery',action='store_true', help='ENTERPRISE: Fit for glitch recoveries')
parser.add_argument('--glitches',type=int, nargs='+', help='ENTERPRISE: Select glitches to fit')

parser.add_argument('--qp-f0-max',type=float,default=10.0,help='ENTERPRISE: Max QP f0')
parser.add_argument('--qp-f0-min',type=float,default=0.1,help='ENTERPRISE: Min QP f0')

parser.add_argument('--qp-sigma-max',type=float,default=10.0,help='ENTERPRISE: Max QP sigma')
parser.add_argument('--qp-sigma-min',type=float,default=0.01,help='ENTERPRISE: Min QP sigma')

parser.add_argument('--emcee',action='store_true', help='ENTERPRISE: Use emcee sampler')
parser.add_argument('--cont',action='store_true', help='ENTERPRISE: Continue existing chain (emcee)')
parser.add_argument('--nthread','-t',type=int, default=1, help="ENTERPRISE: Number of threads (emcee)")
parser.add_argument('--nwalkers',type=int, default=0, help="ENTERPRISE: Number of walkers (emcee)")
parser.add_argument('--test-threads',action='store_true', help='ENTERPRISE: Test threading options (emcee)')

parser.add_argument('--planet','-P',action='store_true', help='ENTERPRISE: Fit for planet orbit parameters')
parser.add_argument('--mass-max',type=float,default=1,help='ENTERPRISE: Max planet mass (Earth masses)')
parser.add_argument('--mass-min',type=float,default=1e-3,help='ENTERPRISE: Min planet mass (Earth masses)')
parser.add_argument('--mass-log-prior',action='store_true',help='ENTERPRISE: Use log prior for planet mass')
parser.add_argument('--period-max',type=float,default=3000,help='ENTERPRISE: Max planet period (days)')
parser.add_argument('--period-min',type=float,default=50,help='ENTERPRISE: Min planet period (days)')
parser.add_argument('--ecc-log-prior',action='store_true',help='ENTERPRISE: Use log prior for planet ecc')

args=parser.parse_args()

##### --------- (BEGIN) READ ARGS --------- #####
print('\n\n')
print(vars(args))
print('\n\n')
PSR = args.PSR

# assume run enterprise is in current dir as well
# mkdir with psr name in current dir

def create_output_folder(input_value):
	# check if it already exists and if not, mkdir
	if os.path.isdir(input_value)==False:
	
		try:
			#print ("\n--- The folder doesn't exist, creating the folder...\n")
			os.mkdir(input_value)
		
		except OSError:
			print ("\n--- The folder path doesn't exist!\n")
			sys.exit()
		
						
	input_value = str(input_value)
	if input_value[-1] != '/':
		input_value = input_value + '/'
		
	return input_value

if args.folder == None:
	currentdir = os.getcwd()
	output_folder = currentdir+'/'
else:
	output_folder = create_output_folder(args.folder)

if args.par==None:
	parfile = 'final.par'

else:	
	parfile = args.par

if args.tim==None:
	timfile = 'final.tim'

else:	
	timfile = args.tim

if args.datapath==None:
	datapath = '/scratch/mkeith/jbo_dataset/4_ready_for_timing/'+PSR+'/analysis/'
else:
	datapath = args.datapath

if datapath[-1] != '/':
	datapath = datapath+'/'

if args.copydata:
	# copy the par and tim files there
	if args.machine == 'local':
		if datapath[-1] != '/':
			loc = datapath + '/'
		else:
			loc = datapath

		#rsync_output = subprocess.run(['rsync','-avhe',loc+parfile,loc+timfile,\
		#				output_folder],capture_output=True)
		rsync_output = subprocess.run(['rsync','-avh',loc+parfile,loc+timfile,\
						output_folder], stdout=subprocess.PIPE,\
						stderr=subprocess.STDOUT)

	else:
		if datapath[-1] != '/':
			sshloc = args.machine + ':' + datapath + '/'
		else:
			sshloc = args.machine + ':' + datapath

		#rsync_output = subprocess.run(['rsync','-avhe','ssh',sshloc+parfile,sshloc+timfile,output_folder],\
#				              capture_output=True)
		rsync_output = subprocess.run(['rsync','-avhe','ssh',sshloc+parfile,sshloc+timfile,output_folder],stdout=subprocess.PIPE,\
		stderr=subprocess.STDOUT)
	#rsync_err = rsync_output.stderr
	#rsync_out = rsync_output.stdout

	print('*** rsync OUTPUT ***')
	print(rsync_output.stdout.decode('utf-8'))

	#if rsync_err.decode('utf-8'):
		#print('*** rsync ERRORS/WARNINGS ***')
		#print(rsync_err.decode('utf-8'))

# enable params & delete TN params from par file:
dumparfile = 'dum_'+parfile 
runparfile = 'run_'+parfile 


## read in the TNRedAmp:

fullparfile = parfile

fpar = open(fullparfile,'r')
flines = fpar.readlines()
for line in flines:
    split_line = line.split()
    if split_line[0] == 'TNRedAmp':
        TNRedAmp = float(split_line[1])

fpar.close()

fdum = open(output_folder+dumparfile,'w')
frun = open(output_folder+runparfile,'w')

enablepar_list = ['../../enableparam.py', fullparfile, 'RAJ', 'DECJ', 'F0', 'F1', 'F2', 'PX',\
			'PB', 'A1', 'OM','ECC','T0','PMRA','PMDEC', 'NE_SW', 'DM']

grep_list = ['grep','-v','TN',dumparfile]
print('\nenableparam.py LIST:',enablepar_list)
subprocess.run(args=enablepar_list,cwd=output_folder,stdout=fdum)
subprocess.run(args=grep_list,cwd=output_folder,stdout=frun)
frun.close()
fdum.close()
subprocess.run(['rm',dumparfile],cwd=output_folder)
print('\nUsed par file: '+output_folder+runparfile)


# run enterprise in there, so ../run_enterprise.py ....

print('\n !!!! Make sure OMP_NUM_THREADS=1.\n')

def get_optent_params(args):

	param_list=[]
	param_list.append('--f2'); param_list.append(str(args.f2))
	param_list.append('-N'); param_list.append(str(args.nsample))
	if args.dm:
		param_list.append('--dm')
		param_list.append('--Adm-max'); param_list.append(str(args.Adm_max))
		param_list.append('--Adm-min'); param_list.append(str(args.Adm_min))
		param_list.append('--dm-ncoeff'); param_list.append(str(args.dm_ncoeff))

	if args.red:
		param_list.append('--Ared-max'); param_list.append(str(args.Ared_max))
		param_list.append('--Ared-min'); param_list.append(str(args.Ared_min))
		param_list.append('--red-gamma-min'); param_list.append(str(args.red_gamma_min))
		param_list.append('--red-gamma-max'); param_list.append(str(args.red_gamma_max))
		param_list.append('--red-ncoeff'); param_list.append(str(args.red_ncoeff))
		if args.red_prior_log:
			param_list.append('--red-prior-log')

	else: 
		param_list.append('--no-red-noise')
	if args.jbo:
		param_list.append('-j')
	if args.be_flag:
		param_list.append('--be-flag'); param_list.append(str(args.be_flag))
	if args.sample:
		if args.models==None:
			if args.emcee:
				param_list.append('--emcee')
				param_list.append('--nwalkers'); param_list.append(str(args.nwalkers))
				param_list.append('-t'); param_list.append(str(args.nthread))
				if args.cont:
					param_list.append('--cont')

				if args.test_threads:
					param_list.append('--test-threads')
		else:
			param_list.append('--models'); param_list.append(str(args.models))

	else:
		param_list.append('--no-sample')

	
	if args.px:
		param_list.append('--px')
		if args.px_verbiest:
			param_list('--px-verbiest')
		else:
			param_list.append('--px-range'); param_list.append(str(args.px_range))

	if args.pm:
		param_list.append('--pm')
	if args.pm_angle:
		param_list.append('--pm-angle')
	if args.pm_angle or args.pm:
		param_list.append('--pm-range'); param_list.append(str(args.px_range))
	
	if args.glitch_recovery:
		param_list.append('--glitch-recovery')
	if args.glitches:
		param_list.append('--glitches'); 
		for g in args.glitches:
			param_list.append(str(g))

	if args.white:
		param_list.append('--efac-max'); param_list.append(str(args.efac_max))
		param_list.append('--efac-min'); param_list.append(str(args.efac_min))
		if args.white_prior_log:
			param_list.append('--white-prior-log')
	else:
		param_list.append('--no-white')
	if args.ngecorr:
		param_list.append('--ngecorr')


	if args.quasiperiodic:
		param_list.append('-Q')
		param_list.append('--Aqp-max'); param_list.append(str(args.Aqp_max))
		param_list.append('--Aqp-min'); param_list.append(str(args.Aqp_min))
		param_list.append('--qp-f0-max'); param_list.append(str(args.qp_f0_max))
		param_list.append('--qp-f0-min'); param_list.append(str(args.qp_f0_min))
		param_list.append('--qp-sigma-max'); param_list.append(str(args.qp_sigma_max))
		param_list.append('--qp-sigma-min'); param_list.append(str(args.qp_sigma_min))

		if args.qp_prior_log:
			param_list.append('--qp-prior-log')

	if args.planet:
		param_list.append('-P')
		param_list.append('--mass-max'); param_list.append(str(args.mass_max))
		param_list.append('--mass-min'); param_list.append(str(args.mass_min))
		param_list.append('--period-max'); param_list.append(str(args.period_max))
		param_list.append('--period-min'); param_list.append(str(args.period_min))
		if args.mass_log_prior:
			param_list.append('--mass-log-prior')
		if args.ecc_log_prior:
			param_list.append('--ecc-log-prior')

	param_list.append('--tspan-mult'); param_list.append(str(args.tspan_mult))
	if args.outdir != None:
		param_list.append('-o'); param_list.append(str(args.outdir))

	if args.pm_ecliptic:
		param_list.append('--pm-ecliptic')
	if args.plot_chain:
		param_list.append('--plot-chain')
	if args.all_corner:
		param_list.append('--all-corner')
	if args.white_corner:
		param_list.append('--white-corner')

	return param_list


if args.ent:

	# get list of enterprise optional parameters.
	optent_params = get_optent_params(args)
	print('\n\n OUTPUT FOLDER: ',output_folder)
	print('\n\n OPTENT list: ',optent_params)
	run_ent_list = ['../../run_enterprise.py',runparfile,timfile]
	for i in range(len(optent_params)):
		run_ent_list.append(optent_params[i])
	run_ent_list.append('--outdir')
	run_ent_list.append('chains/'+PSR+"/")
	print('\n\n run ent list: ',run_ent_list)
	#ent_output = subprocess.run(args=run_ent_list, capture_output=True, cwd=output_folder)
	ent_output = subprocess.run(args=run_ent_list, cwd=output_folder, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#ent_out = ent_output.stdout
	#ent_err = ent_output.stderr
	print('\n\n*** enterprise OUTPUT ***\n\n')
	print(ent_output.stdout.decode('utf-8'))
	#if ent_err.decode('utf-8'):
		#print('\n\n*** enterprise ERRORS/WARNINGS ***\n\n')
		#print(ent_err.decode('utf-8'))

##### --------- (END) READ ARGS --------- #####


##### --------- (BEGIN) FUNCTIONS --------- #####

def get_chain_h5(outdir='PSR/chains/PSR/',filename='chain.h5'):
    
    filename = os.path.join(outdir,filename)
    scls,offs = np.loadtxt(os.path.join(outdir,"scloff"),unpack=True)
    reader = emcee.backends.HDFBackend(filename)

    tau = reader.get_autocorr_time(tol=0)
    #burnin = int(2 * np.max(tau))
    burnin = 0
    #thin = int(0.5 * np.min(tau))
    thin = 10
    #print("tau = {} burn = {} thin = {}".format(tau, burnin, thin))
    samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
    samples *= scls
    samples += offs
    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)

    log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)
    #pars = np.loadtxt(outdir + '/pars.txt', dtype=np.unicode_)

    N=len(samples)
    chain = np.concatenate( (samples, log_prob_samples[:, None], log_prior_samples[:, None],\
                             np.zeros(N)[:,None],np.zeros(N)[:,None]), axis=1)
    # cut extra 25% from chain
    burn=int(0.25*chain.shape[0])

    return chain[burn:]
    
def read_outpars(parsfile='PSR/runparfile.results'):
        
    fopen = open(parsfile)
    flines = fopen.readlines()
    pars = []
    firstline=True
    for line in flines:
        split_line = line.split()
        if firstline:
            firstline=False
        else:
            pars.append(split_line[0])
    return pars            

def read_chain_h5(chaindir='PSR/chains/PSR/',chainfile='chain.h5',\
                         parsfile='PSR/runparfile.results'):
    # usually only use this to get mass(m) and period(p)
    # so just_mp = True
    full_chain = get_chain_h5(outdir=chaindir,filename=chainfile)
    outpars = read_outpars(parsfile=parsfile)
    outpars = np.array(outpars)
    contains_mass = [s for s in outpars if 'mass' in s]
    mass_index = np.where(outpars==contains_mass)
    contains_period = [s for s in outpars if 'period' in s]
    period_index = np.where(outpars==contains_period)
    
    mass = full_chain[:,mass_index]
    period = full_chain[:,period_index]
    return mass[:,0,0], period[:,0,0]   


def sort_all(xsort,y):
	# sorts xsort, y by increasing xsort
	xsort,y = zip(*sorted(zip(xsort,y)))
	xsort = np.array(xsort)
	y = np.array(y)
	return xsort,y

def containInt_linp(m, w, c):
    # sort (m,w) by m, from lowest m to highest:
    m,w = zip(*sorted(zip(m,w)))
    m = np.array(m)
    w = np.array(w)
    
    totalN = np.sum(w)
    
    sump= 0 
    maxind = 0
    
    while sump <= c*totalN:
        #print (maxind,m[maxind],sump)
        sump += w[maxind]
        maxind += 1
        
    lim = m[maxind]
    if maxind+1<len(m):
        errlim = m[maxind+1]-m[maxind] 
    else:
        errlim = np.nan
    return lim,errlim

def mass_w(mass):

    prior_uniflog = 1./(5*np.log(10)*mass) # prior (m|log)
    prior_uniflin = 1./(10-1e-4) # prior (m|lin)

    w = prior_uniflin/prior_uniflog
    
    return w

def f_limmass_linp(mass,w,c=0.95):
    
    # note post is not normalised
    if len(mass) < 10:
        limmass = np.nan
        errlimmass = np.nan
    else:
        limmass,errlimmass = containInt_linp(mass,w,c)

    return limmass,errlimmass

def get_avg_w(x,w):
    return np.sum(x*w)/np.sum(w)

def get_stdev_w(x,w):
    avgx = get_avg_w(x,w)
    stdev2 = np.sum(x*x*w)/np.sum(w) - avgx*avgx
    return np.sqrt(stdev2)


##### --------- (END) FUNCTIONS --------- #####

##### --------- (BEGIN) MAIN --------- #####

if args.test_threads==False:
    mass,period = read_chain_h5(chaindir=output_folder+'chains/'+PSR+'/',chainfile='chain.h5',\
		         parsfile=output_folder+runparfile+'.results')
    if args.mass_log_prior:
        mass = 10**mass

    period,mass = sort_all(period,mass)
    dperiod = (args.period_max-args.period_min)/2.
    midperiod = args.period_min+dperiod

    w = mass_w(mass)
    normw = w/np.sum(w)
    plot_periodmass=True

    if plot_periodmass:

        plt.figure()
        plt.plot(mass,period,'.k')
        plt.xlabel('mass (ME)')
        plt.ylabel('period (days)')


        #plt.axhline(midperiod,color = 'red',ls='-')
        plt.axhline(midperiod-dperiod,color = 'red',ls='--')
        plt.axhline(midperiod+dperiod,color = 'red',ls='--')

        plt.savefig(output_folder+'Period_vs_mass.pdf')
        plt.close(plt.gcf())

    limmass,errlimmass = f_limmass_linp(mass,w)

    if get_avg_w(mass,w) > 3*get_stdev_w(mass,w):
        flagdet = True
    else:
        flagdet = False


    if args.plot_allmass:

        fig = plt.figure(figsize=(16,9))
        ax = fig.add_subplot(111)
        plt.title("period = "+str(np.min(period))+" to "+str(np.max(period)))

        #plt.errorbar(mass,w,np.sqrt(w),c='k',marker='.',ls='',label="mean: %.4f\nstdev: %.4f\nflag: %s" \
        #             % (get_avg_w(mass,w), get_stdev_w(mass,w), str(flagdet)))

        #plt.hist(mass,density=True,bins=100,label="log: mean: %.4f\nstdev: %.4f\nflag: %s" \
        #             % (get_avg(mass,w), get_stdev(mass,w), str(flagdet)))
        plt.hist(mass,weights=w,density=True,bins=100,label="linear (non-log) \nmean: %.4f\nstdev: %.4f\nflag: %s" \
                     % (get_avg_w(mass,w), get_stdev_w(mass,w), str(flagdet)))

        plt.axvline(limmass,c='red')
        ylims = plt.gca().get_ylim()
        plt.fill_betweenx(ylims,limmass-errlimmass,limmass+errlimmass,alpha=0.5,color='red')

        ax.legend(loc=1)
        plt.xlabel('mass (ME)')
        plt.ylabel('weight (no. occur)')
        plt.savefig(output_folder+'all_masses_linprior.pdf')
        #plt.show()
        plt.close(plt.gcf())


    np.save(output_folder+'masslim_werr_linprior',[midperiod,dperiod,limmass,errlimmass,flagdet])
    np.savetxt(output_folder+'masslim_linprior.txt', np.c_[midperiod,dperiod,limmass,errlimmass,flagdet],header='mPeriod(days) sPeriod(days) Masslimit(Me) errMasslimit(Me) Detection(3sigma)',fmt='%.10f')    
        

    if args.plot_masslim:
	#print('\n\n\n WHY: ',avg_period,limmass)
        fig = plt.figure(figsize=(16,9))
        ax = fig.add_subplot(111)
        if flagdet:
        	plt.plot(midperiod, limmass, color='r', ls='', marker='*')
        plt.errorbar(midperiod, limmass, yerr=errlimmass, color='k', ls='-', marker='')
        plt.step(np.array([midperiod-dperiod,midperiod+dperiod]),np.array([limmass,limmass]),c='k')
         
        plt.loglog()
        ax.set_ylabel('mass (ME)')
        ax.set_xlabel('period (days)')
        #ax.set_xlim((args.period_min),(args.period_max))
        ax.set_xlim(20,7000)
        ax.set_ylim((args.mass_min),args.mass_max)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        plt.savefig(output_folder+'masslim_linprior.pdf')
        plt.close(plt.gcf())
        #plt.show()

    print('\nOutput folder: '+ output_folder)
