#include <vector>
#include <iostream>

void copyAndRenameHists(TString pt_bin)
{
  std::cout << "WARNING:  Do not forget to run SFramePlotter first!  Else, there are no input root files.  Do not forget to merge single-top s+t channels." << std::endl;

  TFile* final_file = TFile::Open("theta-input_PT-"+pt_bin+".root", "RECREATE"); // my output file
  
  TString pathName = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/Muon/";
  TString fileName_prefix = "uhh2.AnalysisModuleRunner.";
  std::vector<TString> fileNames_DATA {"DATA"};
  std::vector<TString> fileNames_MC {"SingleTop_tWch", "SingleTop_stch", "TTbar", "WJets", "DYJets", "Diboson", "QCD"};
  std::vector<vector<TString>> fileNames {fileNames_DATA, fileNames_MC};
  
  for(int i = 0; i < fileNames.size(); i++)
    {
      for(int j = 0; j < fileNames.at(i).size(); j++)
	{
	  TString input_file = pathName+fileName_prefix+(i==0 ? "DATA." : "MC.")+(fileNames.at(i)).at(j)+".root";
	  std::cout << "Copying hists from input file:   " << input_file << std::endl;
	  TFile *input_file_ = TFile::Open(input_file, "read");
	  TH1F *input_histo = (TH1F*)input_file_->Get("MVA/"+pt_bin+"_mvaD_rebin01");
	  input_histo->SetName(pt_bin+"__"+fileNames.at(i).at(j)+"__for-theta");
	  final_file->cd(); // change to this file to be able to write on it
	  input_histo->Write();
	}
    }

  final_file->Close(); // close my output file
}

void writeInputHists()
{
  std::vector<TString> pt_bins {"200to300", "300to400", "400to600", "600to1200"};
  for(int bin = 0; bin < pt_bins.size(); ++bin)
    {
      copyAndRenameHists(pt_bins.at(bin));
    }
}
