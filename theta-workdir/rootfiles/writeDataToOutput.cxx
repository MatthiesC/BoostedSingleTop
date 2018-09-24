#include <vector>
#include <iostream>

void copyDataHists(TString channel)
{
  TFile* input_file = TFile::Open("theta-input_"+channel+"_stat.root", "READ");
  
  TFile* final_file = TFile::Open("theta-output_"+channel+"_stat.root", "UPDATE");

  for(TString pt_bin : {"200to300", "300to400", "400to600", "600to1200"})
    {
      TH1F *input_histo = (TH1F*)input_file->Get("BDT_"+pt_bin+"__DATA");
      final_file->cd();
      input_histo->Write();
    }
  final_file->Close();
  input_file->Close();
}

void sortHists(TString channel)
{
  TFile* input_file = TFile::Open("theta-input_"+channel+"_stat.root", "READ");
  
  TFile* theta_file = TFile::Open("theta-output_"+channel+"_stat.root", "READ");

  TFile* final_file = TFile::Open("theta-output_"+channel+"_stat_sorted.root", "RECREATE");

  std::vector<TString> names_of_histos;
  TList* list_of_histos = input_file->GetListOfKeys();
  TObject* obj;
  for(auto key : *list_of_histos)
    {
      obj = ((TKey*)key)->ReadObj();
      names_of_histos.push_back(obj->GetName());
    }
  for(TString name_of_hist : names_of_histos)
    {
      TH1F *hist = (TH1F*)theta_file->Get(name_of_hist);
      final_file->cd();
      hist->Write();
    }

  final_file->Close();
  theta_file->Close();
  input_file->Close();
}

void writeDataToOutput()
{
  for(TString channel : {"Muon", "Electron"})
    {
      copyDataHists(channel);
      sortHists(channel);
    }
}
