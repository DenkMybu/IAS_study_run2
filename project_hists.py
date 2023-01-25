import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors,     TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math



#python project_hists.py template_2018A_15mars_selection_5_bin_eta_2p1_pmin_5_pmax_45_final_HSCP_noiso_new.root


if __name__ == "__main__":

     file_trf = sys.argv[1]

     split_name = file_trf.split('_')
     name_p_range = "proj_1D_template_" + split_name[8] + split_name[9] + split_name[10] + split_name[11] + ".root"
     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", file_trf)
     else:
         print("Successfully opened root : ",file_trf)

     str_template = "Charge_Vs_Path"

     loaded_template = file0.Get(str_template)

     num_binx = loaded_template.GetNbinsX() 
     num_biny = loaded_template.GetNbinsY() 
     num_binz = loaded_template.GetNbinsZ() 

     print("x bins :", num_binx, " , y bins :", num_biny, " , z bins :", num_binz)

     print("Creating root file to store results ...")

     hfile = TFile(name_p_range, 'RECREATE' )

     name_nb_stat = "AA_Stat_per_1d_projection"
     nb_stat_per_tpt = TH1D(name_nb_stat,name_nb_stat,10000,0,1000000)
     for i in range(0,num_binx+1):
         for p in range(0,num_biny+1):
	     loaded_template.GetXaxis().SetRange(i,i)
             loaded_template.GetYaxis().SetRange(p,p)
             name = "projection_xbin_" + str(i) + "_ybin_" + str(p)
             test = TH1D(name,name,500,0,5000)
             test = loaded_template.Project3D("z")
             nb_stat_per_tpt.Fill(test.Integral())
             test.SetName(name)
             if(test.GetEntries() != 0):
                 print(name, " has an integral z entries = ", test.Integral())
                 test.Write()

     nb_stat_per_tpt.Write()
     hfile.Write()
