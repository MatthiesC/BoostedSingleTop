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

    inputpath = "./rootfiles/"

    model_stat = build_model_from_rootfile(inputpath+fname_stat, include_mc_uncertainties=True)
    #model_stat = build_model_from_rootfile(inputpath+fname_stat)
    model_stat.fill_histogram_zerobins()

    rate_unc_tt    = math.log(1.5)
    rate_unc_tW    = math.log(2.0)
    rate_unc_QCD   = math.log(2.0)
    rate_unc_w     = math.log(1.1)
    rate_unc_dy    = math.log(1.1)

    #model_stat.add_lognormal_uncertainty('luminosity_unc', math.log(1.025), procname='*', obsname='*')

    model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate',  math.log(2.0) , procname='SingleTop_tWch',obsname='*')
    model_stat.add_lognormal_uncertainty('TTbar_rate',           math.log(1.5) , procname='TTbar',obsname='*')

    #model_stat.add_lognormal_uncertainty('SingleTop_stch_rate', rate_unc_other , procname='SingleTop_stch',obsname='*')
    model_stat.add_lognormal_uncertainty('WJets_rate', rate_unc_w , procname='WJets',obsname='*')
    model_stat.add_lognormal_uncertainty('DYJets_rate', rate_unc_dy , procname='DYJets',obsname='*')
    #model_stat.add_lognormal_uncertainty('Diboson_rate', rate_unc_other , procname='Diboson',obsname='*')
    #model_stat.add_lognormal_uncertainty('QCD_rate', rate_unc_QCD , procname='QCD',obsname='*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate', rate_unc_tt , procname='TTbar',obsname='*')

    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_200to300', rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_200to300*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_200to300', rate_unc_tt , procname='TTbar',obsname='BDT_200to300*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_300to400', rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_300to400*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_300to400', rate_unc_tt , procname='TTbar',obsname='BDT_300to400*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_400to600', rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_400to600*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_400to600', rate_unc_tt , procname='TTbar',obsname='BDT_400to600*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_600to1200', rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_600to1200*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_600to1200', rate_unc_tt , procname='TTbar',obsname='BDT_600to1200*')

    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_200to300',  rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_200to300*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_200to300',           rate_unc_tt , procname='TTbar',obsname='BDT_200to300*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_300to400',  rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_300to400*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_300to400',           rate_unc_tt , procname='TTbar',obsname='BDT_300to400*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_400to600',  rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_400to600*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_400to600',           rate_unc_tt , procname='TTbar',obsname='BDT_400to600*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_600to1200', rate_unc_tW , procname='SingleTop_tWch',obsname='BDT_600to1200*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_600to1200',          rate_unc_tt , procname='TTbar',obsname='BDT_600to1200*')

    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_200to300',  math.log(2.0) , procname='SingleTop_tWch',obsname='BDT_200to300*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_200to300',           math.log(1.5) , procname='TTbar',obsname='BDT_200to300*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_300to400',  math.log(2.0) , procname='SingleTop_tWch',obsname='BDT_300to400*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_300to400',           math.log(1.4) , procname='TTbar',obsname='BDT_300to400*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_400to600',  math.log(2.0) , procname='SingleTop_tWch',obsname='BDT_400to600*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_400to600',           math.log(1.3) , procname='TTbar',obsname='BDT_400to600*')
    #model_stat.add_lognormal_uncertainty('SingleTop_tWch_rate_600to1200', math.log(2.0) , procname='SingleTop_tWch',obsname='BDT_600to1200*')
    #model_stat.add_lognormal_uncertainty('TTbar_rate_600to1200',          math.log(1.2) , procname='TTbar',obsname='BDT_600to1200*')

    options = Options()
    options.set('minimizer', 'strategy', 'robust')
    #options.set('minimizer', 'strategy')

    n_name_stat = deepcopy(outname_stat)
    n_name_stat = n_name_stat.replace(".root","_")

    model_summary(model_stat)

    print '============'
    print 'stat'
    print '============'
    mle_output_stat = mle(model_stat, input='data', n=200, with_covariance=True, chi2=True, options=options, signal_process_groups ={'background_only':[]})

    #print 'Results:'
    #print mle_output_stat
    #report.write_html('htmlout')
    
    PlotPostFitCorrelations(model_stat, mle_output_stat['background_only'], "fitResults/nuisance/Corr_"+n_name_stat)  
    writeOutputFile(inputpath+fname_stat, inputpath+outname_stat, mle_output_stat['background_only'], model_stat)

    for pf_vals in mle_output_stat.itervalues():
        del pf_vals['__nll']
        del pf_vals['__cov']

    postfit_stat = ThetaPostFitPlot(mle_output_stat)
    postfit_stat.make_plots("fitResults/nuisance/",n_name_stat)

    #Use this instead of Torben's "writeOutputFile" if you want to have wrong uncertainties: ;-)
    #parameter_values = {}
    #for p in model_stat.get_parameters([]):
    #    parameter_values[p] = mle_output_stat['background_only'][p][0][0]
    #histos = evaluate_prediction(model_stat, parameter_values, include_signal = False)
    #write_histograms_to_rootfile(histos, inputpath+outname_stat)

#run("theta-input_Muon_stat.root", "theta-output_Muon_stat.root")
#run("theta-input_Electron_stat.root", "theta-output_Electron_stat.root")

#run("theta-input_stat.root", "theta-output_stat.root")
#run("theta-input_Electron_1b_200to300.root", "theta-output_Electron_1b_200to300.root")

channels = ["Electron", "Muon"]
ptbins   = ["200to300", "300to400", "400to600", "600to1200"]
nbtags   = ["1b", "2b"]

for c in channels:
    for p in ptbins:
        for nb in nbtags:
            name = c+"_"+nb+"_"+p
            run("theta-input_"+name+".root", "theta-output_"+name+".root")
