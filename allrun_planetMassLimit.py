#!/usr/bin/env python
from timeit import default_timer as timer
start = timer()

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import emcee
import argparse
from matplotlib.backends.backend_pdf import PdfPages

# ------ PARSE ARGUMENTS:


parser=argparse.ArgumentParser(description="Find planet mass limits for multiple pulsars and multiple runs of different period ranges.")

# --- SPECIFIC ARGS:
parser.add_argument('txtfile', help='Text file with the names of pulsars to analyse.')
parser.add_argument('--add-period-range','-r', dest='pranges', action='append', nargs=2, metavar=('MIN','MAX'), help='Add a(nother) range of planet period (in days) to be investigated by enterprise. Can be called multiple times for multiple ranges at once.')


# --- run_planetMassLimit ARGS:

parser.add_argument('--no-enterprise',dest='ent',default=True,action='store_false', help='Disable running enterprise search')
parser.add_argument('--no-copydata',dest='copydata',default=True,action='store_false', help='Disable copying data files .par and .tim')

parser.add_argument('--machine',default='godrevy', help='Machine to copy (rsync) data from. Can also be set to \'local\'. Default is \'godrevy\'.')
parser.add_argument('--datapath', help='Path where data is on the remote machine. Default is \'/scratch/mkeith/jbo_dataset/4_ready_for_timing/+$PSR+/analysis/\'.')
parser.add_argument('--par', default='final.par',help='Name of the .par file. Default is \'final.par\'.')
parser.add_argument('--tim', default='final.tim',help='Name of the .tim file. Default is \'final.tim\'.')
parser.add_argument('--plot-periodmass',default=False,action='store_true', help='Option to plot period vs mass samples.')
parser.add_argument('--plot-allmass',default=False,action='store_true', \
			help=r'Option to plot the mass distribution, also showing the 95%% limit.')
parser.add_argument('--plot-masslim',default=False,action='store_true', \
			help='Option to plot the 95%% mass limit against period.')
#parser.add_argument('--folder', help='Set the folder to save results in. Default is PSR name.')

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
parser.add_argument('--Ared-max','-A',type=float,help='ENTERPRISE: Max log10A_Red')
parser.add_argument('--Ared-min',type=float,help='ENTERPRISE: Min log10A_Red')
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
parser.add_argument('--mass-log-prior',action='store_true',help='EMTERPRISE: Use log prior for planet mass')
#parser.add_argument('--period-max',type=float,default=3000,help='ENTERPRISE: Max planet period (days)')
#parser.add_argument('--period-min',type=float,default=50,help='ENTERPRISE: Min planet period (days)')
parser.add_argument('--ecc-log-prior',action='store_true',help='EMTERPRISE: Use log prior for planet ecc')

args=parser.parse_args()

##### --------- (BEGIN) READ ARGS --------- #####
print('\n\n')
print(vars(args))
print('\n\n')

# read pulsars from textfile:
PSRs = []
psrfile = open(args.txtfile,'r')
for line in psrfile:
	if line.split()[0][0] != '#':
		PSRs.append(line.split()[0])
print('\nPSRs: ',PSRs)

pranges = np.array(args.pranges)
pranges = pranges.astype(np.float)
pranges = np.sort(pranges,axis=1) # makes sure MIN<MAX
print('\nPRANGES:')
print(pranges)

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


def get_optent_params2(args):

	param_list=[]
	param_list.append('--f2'); param_list.append(str(args.f2))
	param_list.append('-N'); param_list.append(str(args.nsample))
	if args.dm:
		param_list.append('--dm')
		param_list.append('--Adm-max'); param_list.append(str(args.Adm_max))
		param_list.append('--Adm-min'); param_list.append(str(args.Adm_min))
		param_list.append('--dm-ncoeff'); param_list.append(str(args.dm_ncoeff))

	if args.red:
		if args.Ared_max:	
			param_list.append('--Ared-max'); param_list.append(str(args.Ared_max))
		if args.Ared_min:		
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

##### -------- for each pulsar:
k=0
for PSR in PSRs:
	print('\n\n ###### ('+str(k)+') PSR '+PSR+'\n\n')
	output_folder = create_output_folder(PSR)

	for i in range(len(pranges)):
		run_list = ['../../run_planetMassLimit.py',PSR,'--par',args.par,'--tim',args.tim]
		if args.ent==False:
			run_list.append('--no-enterprise')
		if args.copydata==False:
			run_list.append('--no-copydata')
		if args.plot_periodmass:
			run_list.append('--plot-periodmass')
		if args.plot_allmass:
			run_list.append('--plot-allmass')
		if args.plot_masslim:
			run_list.append('--plot-masslim')

		minP=pranges[i][0];maxP=pranges[i][1]
		print('\n\n ########## ('+str(k)+'.'+str(i)+') range: '+str(int(minP))+' - '+str(int(maxP))+' days'+'\n\n')

		ifolder = create_output_folder(output_folder+'p'+str(int(minP))+'to'+str(int(maxP))+'/')

		#fullparfile = ifolder+args.par

		#fpar = open(fullparfile,'r')
		#flines = fpar.readlines()
		#for line in flines:
		#    split_line = line.split()
		#    if split_line[0] == 'TNRedAmp':
		#        TNRedAmp = float(split_line[1])

		#fpar.close()

		#run_list.append('--folder');run_list.append(local)
		optent_params = get_optent_params2(args)
		print('OPTENT_PARAMS: ',(optent_params))	
		if args.copydata:
			run_list.append('--machine');run_list.append(args.machine)
		if args.datapath!=None:
			run_list.append('--datapath');run_list.append(args.datapath)

		for j in range(len(optent_params)):
			# the ent params, excluding period-max and min
			run_list.append(optent_params[j])

		run_list.append('--period-min');run_list.append(str(minP))
		run_list.append('--period-max');run_list.append(str(maxP))
	
		#if not args.Ared_max:
		#	run_list.append('--Ared-max'); param_list.append(str(TNRedAmp+1.))

		#if not args.Ared_min:
		#	run_list.append('--Ared-min'); param_list.append(str(TNRedAmp-1.))

		print('\n\n RUN LIST: ',run_list)
		#run_output = subprocess.run(args=run_list,capture_output=True,cwd=ifolder)
		run_output = subprocess.run(args=run_list, cwd=ifolder, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		#run_out = run_output.stdout
		#run_err = run_output.stderr
		print('\n\n*** run OUTPUT ***\n\n')
		print(run_output.stdout.decode('utf-8'))
		#if run_err.decode('utf-8'):
			#print('\n\n*** run ERRORS/WARNINGS ***\n\n')
			#print(run_err.decode('utf-8'))
	k += 1

end = timer()
print('Time: '+ str(end - start)+' sec') # Time in seconds


