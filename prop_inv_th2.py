import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3
import sys
import argparse

def poissonHisto(h,p,RNG):
    for i in range(p.GetNbinsX()+1):
        for j in range(p.GetNbinsY()+1):
            for k in range(p.GetNbinsZ()+1):
                h.SetBinContent(i,j,k,RNG.Poisson(p.GetBinContent(i,j,k)))

#python3 propdraw.py template_2018B_15mars_selection_5_bin_eta_2p1_pmin_10.root





if __name__ == "__main__":
     hfile = TFile( 'inversed.root', 'RECREATE' )
     file0 = ROOT.TFile.Open("2d_test.root")

     for x in range(file0.GetNbinsX()+1):
         for y in range(file0.GetNbinsY()+1):
             file0.GetBinContent(x,y)
             hfile.Fill(y,x)
     '''

     test = TH2D("a_vs_b","a_vs_b",100,0,100,100,0,100)
     inversed = TH2D("b_vs_a","b_vs_a",100,0,100,100,0,100)

     for u in range(100):
             test.Fill(u,20)
     '''
     
     hfile.Write()
