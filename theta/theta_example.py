# -*- coding: utf-8 -*-

import numpy
import ROOT
from array import array

import os
import sys

sys.path.append('./')
from PostFitUncertainties import *
from PostFitCorrelations import *
from ThetaPostFitPlot import *

def run(fname_stat, outname_stat, write_report = False):
    if write_report:report.reopen_file()

    #model_stat = build_model_from_rootfile(inputpath+fname_stat, include_mc_uncertainties=True)
    model_stat = build_model_from_rootfile(fname_stat)
    model_stat.fill_histogram_zerobins()

    rate_unc = math.log(2.0)
    rate_unc= math.log(2.0)

    model_stat.add_lognormal_uncertainty('TTbar_mergedTop_Pass_rate', rate_unc , procname='TTbar_mergedTop',obsname='Mass_pass')
    model_stat.add_lognormal_uncertainty('TTbar_semimerged_Pass_rate',  rate_unc, procname='TTbar_semimerged',obsname='Mass_pass')
    model_stat.add_lognormal_uncertainty('TTbar_notmerged_Pass_rate', rate_unc , procname='TTbar_notmerged',obsname='Mass_pass')

    model_stat.add_lognormal_uncertainty('TTbar_mergedTop_Fail_rate', rate_unc, procname='TTbar_mergedTop',obsname='Mass_fail')
    model_stat.add_lognormal_uncertainty('TTbar_semimerged_Fail_rate', rate_unc, procname='TTbar_semimerged',obsname='Mass_fail')
    model_stat.add_lognormal_uncertainty('TTbar_notmerged_Fail_rate', rate_unc, procname='TTbar_notmerged',obsname='Mass_fail') 

    options = Options()
    options.set('minimizer', 'strategy', 'robust')

    n_name_stat = deepcopy(outname_stat)
    n_name_stat = n_name_stat.replace(".root","_")
   
    print '============'
    print 'stat'
    print '============'
    mle_output_stat = mle(model_stat, input='data', n=100, with_covariance=True, chi2=True, options=options, signal_process_groups ={'background_only':[]})

    PlotPostFitCorrelations(model_stat, mle_output_stat['background_only'], "fitResults/nuisance/Corr_"+n_name_stat)  
    writeOutputFile(fname_stat,  outname_stat, mle_output_stat['background_only'], model_stat)
 
    for pf_vals in mle_output_stat.itervalues():
        del pf_vals['__nll']
        del pf_vals['__cov']
     
    postfit_stat = ThetaPostFitPlot(mle_output_stat)
    postfit_stat.make_plots("fitResults/nuisance/",n_name_stat)
   

run("thetaFile_400_PUPPI_stat.root", "Hists_400_PUPPI_stat.root")
