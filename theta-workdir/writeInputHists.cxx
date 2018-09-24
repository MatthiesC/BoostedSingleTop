#include <vector>
#include <iostream>

void copyAndRenameHists(TString channel)
{
  TFile* final_file = TFile::Open("rootfiles/theta-input_"+channel+"_stat.root", "RECREATE"); // my output file
  
  TString pathName = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/"+channel+"/";
  TString fileName_prefix = "uhh2.AnalysisModuleRunner.";
  std::vector<TString> fileNames_DATA {"DATA"};
  std::vector<TString> fileNames_MC {"SingleTop_tWch", "SingleTop_stch", "TTbar", "WJets", "DYJets", "Diboson"}; //, "QCD"};
  std::vector<vector<TString>> fileNames {fileNames_DATA, fileNames_MC};
  
  for(TString pt_bin : {"200to300", "300to400", "400to600", "600to1200"})
    {
      for(int i = 0; i < fileNames.size(); i++)
	{
	  for(int j = 0; j < fileNames.at(i).size(); j++)
	    {
	      for(TString n_btags : {"1b", "2b"})
		{
		  TString input_file = pathName+fileName_prefix+(i==0 ? "DATA." : "MC.")+(fileNames.at(i)).at(j)+".root";
		  std::cout << "Copying hists from input file:   " << input_file << std::endl;
		  TFile *input_file_ = TFile::Open(input_file, "read");
		  TH1F *input_histo = (TH1F*)input_file_->Get("MVA_1t"+n_btags+"/"+pt_bin+"_mvaD"); // _mvaD hat 40 bins!
		  if (pt_bin == "200to300") input_histo->Rebin(4); // 10 bins
		  if (pt_bin == "300to400") input_histo->Rebin(5); // 8 bins
		  if (pt_bin == "400to600") input_histo->Rebin(8); // 5 bins
		  if (pt_bin == "600to1200") input_histo->Rebin(10); // 4 bins
		  input_histo->SetTitle("BDT response");
		  input_histo->SetName("BDT_"+pt_bin+"_"+n_btags+"__"+fileNames.at(i).at(j));
		  final_file->cd(); // change to this file to be able to write on it
		  input_histo->Write();
		}
	    }
	}
    }

  final_file->Close(); // close my output file
}

void writeInputHists()
{
  std::vector<TString> pt_bins {"200to300", "300to400", "400to600", "600to1200"};
  std::cout << "WARNING:  Do not forget to run SFramePlotter first!  Else, there are no input root files.  Do not forget to merge single-top s+t channels." << std::endl;
  for(TString channel : {"Muon", "Electron"})
    {
      copyAndRenameHists(channel);
    }
}
