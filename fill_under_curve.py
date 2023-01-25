import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math




#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     file_trf = sys.argv[1]

     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", file_trf)
     else:
         print("Successfully opened root : ",file_trf)




     str_hist = "Ias_all_base_cutIH"
     hfile = TFile( str_hist+"_filled.root", 'RECREATE' )
     
     c1 = TCanvas("Ias_all_p_10_50","Ias_all_p_10_50",800,800)
 
     hist = file0.Get(str_hist) 
     hist.SetLineColor(4)
     hist.Draw()

     hist_cl = hist.Clone("color")
     hist_cl.GetXaxis().SetRange(6,49)
     hist_cl.SetFillColor(2)
     hist_cl.SetFillStyle(1001) 
     hist_cl.Draw("FC SAME")

     c1.Write()
     hfile.Write()

