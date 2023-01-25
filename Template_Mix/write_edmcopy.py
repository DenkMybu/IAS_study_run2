import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math




#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":
     
    txt_file = open("list_file.txt","w")
    for i in range(999):         
        filename = "file:Histos_" + str(i) +".root"
        if not filename:
            print("could not open " ,filename)

        txt_file.write(filename)
        txt_file.write("\n")
    print("Sucessfully terminated")     
