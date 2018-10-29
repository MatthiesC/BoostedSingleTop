#include <vector>
#include <iostream>

void copyAndRenameHists() //(TString channel)
{
  //TFile* final_file = TFile::Open("rootfiles/theta-input_stat.root", "RECREATE"); // my output file
 
  for(TString channel : {"Muon", "Electron"})
    {
      TString pathName = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/"+channel+"/";
      TString fileName_prefix = "uhh2.AnalysisModuleRunner.";
      std::vector<TString> fileNames_DATA {"DATA"};
      std::vector<TString> fileNames_MC {"SingleTop_tWch", "SingleTop_stch", "TTbar", "WJets", "DYJets", "Diboson"}; //, "QCD"};
      std::vector<vector<TString>> fileNames {fileNames_DATA, fileNames_MC};
      std::vector<TString> systNames {"", "syst/btagSFbc_up/", "syst/btagSFbc_down/", "syst/btagSFudsg_up/", "syst/btagSFudsg_down/", "syst/scale_up/", "syst/scale_down/", "syst/jec_up/", "syst/jec_down/", "syst/jer_up/", "syst/jer_down/", "syst/pileup_up/", "syst/pileup_down/"};
  
      for(TString pt_bin : {"200to300", "300to400", "400to600", "600to1200"})
	{
	  for(TString n_btags : {"1b", "2b"})
	    {
	      TFile* final_file = TFile::Open("rootfiles/theta-input_"+channel+"_"+n_btags+"_"+pt_bin+".root", "RECREATE"); // my output file
	      for(int k = 0; k < systNames.size(); k++)
		{
		  for(int i = (k == 0 ? 0 : 1); i < fileNames.size(); i++)
		    {
		      for(int j = 0; j < fileNames.at(i).size(); j++)
			{
			  //for(TString n_btags : {"1b", "2b"})
			  //{
			      TString input_file = pathName+systNames.at(k)+fileName_prefix+(i==0 ? "DATA." : "MC.")+(fileNames.at(i)).at(j)+".root";
			      std::cout << "Copying hists from input file:   " << input_file << std::endl;
			      TFile *input_file_ = TFile::Open(input_file, "read");
			      TH1F *input_histo = (TH1F*)input_file_->Get("MVA_1t"+n_btags+"/"+pt_bin+"_mvaD"); // _mvaD has 40 bins!
			      ////////////////////////////////////////////////////////////////////////////////
			      double bin40 = input_histo->GetBinContent(40);
			      double bin41 = input_histo->GetBinContent(41);
			      input_histo->SetBinContent(40,bin40+bin41); // include overflow into last bin (since BDT=1 went into the overflow bin)
			      double bin40_err = input_histo->GetBinError(40);
			      double bin41_err = input_histo->GetBinError(41);
			      input_histo->SetBinError(40,sqrt(bin40_err*bin40_err+bin41_err*bin41_err));
			      input_histo->SetBinContent(41,0); // set overflow bin content to zero
			      input_histo->SetBinError(41,0);
			      ////////////////////////////////////////////////////////////////////////////////
			      if (pt_bin == "200to300") input_histo->Rebin(2); // 20 bins
			      if (pt_bin == "300to400") input_histo->Rebin(4); // 10 bins
			      if (pt_bin == "400to600") input_histo->Rebin(8); // 5 bins
			      if (pt_bin == "600to1200") input_histo->Rebin(10); // 4 bins
			      input_histo->SetTitle("BDT response");
			      TString systSuffix = "";
			      if(systNames.at(k).EndsWith("_up/"))
				{
				  systSuffix += systNames.at(k);
				  systSuffix.ReplaceAll("syst/","__");
				  systSuffix.ReplaceAll("_up/","__plus");
				}
			      else if(systNames.at(k).EndsWith("_down/"))
				{
				  systSuffix += systNames.at(k);
				  systSuffix.ReplaceAll("syst/","__");
				  systSuffix.ReplaceAll("_down/","__minus");
				}
			      input_histo->SetName("BDT_"+pt_bin+"_"+n_btags+"_"+(channel == "Electron" ? "ele" : (channel == "Muon" ? "muo" : ""))+"__"+fileNames.at(i).at(j)+systSuffix);
			      final_file->cd(); // change to this file to be able to write on it
			      if (fileNames.at(i).at(j) == "SingleTop_tWch" && systSuffix.BeginsWith("__scale")) continue; // exclude muR/muF variation for single-top tW
			      input_histo->Write();
			      input_file_->Close();
			      //}
			}
		    }
		}
	      final_file->Close(); // close my output file
	    }
	}
    }

  //final_file->Close(); // close my output file
}

void writeInputHists()
{
  std::vector<TString> pt_bins {"200to300", "300to400", "400to600", "600to1200"};
  std::cout << "WARNING:  Do not forget to run SFramePlotter first!  Else, there are no input root files.  Do not forget to merge single-top s+t channels." << std::endl;
  //for(TString channel : {"Muon", "Electron"})
  //{
  copyAndRenameHists();
  //}
}
