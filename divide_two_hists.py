import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse





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

     hfile = TFile( 'Ratio_plots_num_denom.root', 'RECREATE' )

     hist_num_raw = file0.Get("HSCP_Mass_Single_ihcut_ptsupp50_qtl_80_90_ias") 
     hist_denom_raw = file0.Get("HSCP_Mass_Triple_ihcut_ptsupp50_qtl_80_90_ias")

     low_bin = 200
     high_bin = 1000
     low_val_num_raw = hist_num_raw.GetXaxis().FindBin(low_bin)
     up_val_num_raw = hist_num_raw.GetXaxis().FindBin(high_bin)

     low_val_denom_raw = hist_denom_raw.GetXaxis().FindBin(low_bin)
     up_val_denom_raw = hist_denom_raw.GetXaxis().FindBin(high_bin)

    
     itg_num_raw = hist_num_raw.Integral(low_val_num_raw,up_val_num_raw)
     itg_denom_raw = hist_denom_raw.Integral(low_val_denom_raw,up_val_denom_raw)
     print("BEFORE SCALING")
     print("Integral ",low_bin," to ", high_bin, " from single  : ", itg_num_raw, " and triple : ", itg_denom_raw)
     
     print("Difference ITG(simple) - ITG(triple) = ", itg_num_raw - itg_denom_raw)

     hist_num = hist_num_raw.Clone("HSCP_Mass_Single_cutIH_PT50_quantile_80_90")
     hist_denom = hist_denom_raw.Clone("HSCP_Mass_Multiple_cutIH_PT50_quantile_80_90")

     hname = hist_num.GetName()
     hnameden = hist_denom.GetName()

     hist_num.Scale(1./hist_num.Integral())
     hist_denom.Scale(1./hist_denom.Integral())


     
     raw_divide = TCanvas("single/multiple","single/multiple",800,800)
     #all_iso_diff = TH1D("ias_var_noiso_write_tpt","aha",50,0,1)
     hist_num.SetMarkerStyle(2)
     hist_num.SetMarkerColor(4)
     #hist_num.SetLineColor(2)
     hist_num.SetLineWidth(1)
     hist_num.Draw()
     hist_denom.SetMarkerStyle(2)
     hist_denom.SetMarkerColor(2)
     #hist_denom.SetLineColor(4)
     hist_denom.SetLineWidth(1)
     hist_denom.Draw("same")
     lej_denom = TLegend(0.1,0.7,0.5,0.9)
     lej_denom.AddEntry(hist_num,"Single","l")
     lej_denom.AddEntry(hist_denom,"Multiple","l")
     lej_denom.Draw("same")
     raw_divide.Write()


     #all_iso_diff.Write()
     #diff_iso_tpt.Write()
     #all_iso_diff2 = TH1D("ias_var_iso_write_tpt","aha",50,0,1)




     low_val_num = hist_num.GetXaxis().FindBin(low_bin)
     up_val_num = hist_num.GetXaxis().FindBin(high_bin)

     low_val_denom = hist_denom.GetXaxis().FindBin(low_bin)
     up_val_denom = hist_denom.GetXaxis().FindBin(high_bin)

    
     itg_num = hist_num.Integral(low_val_num,up_val_num)
     itg_denom = hist_denom.Integral(low_val_denom,up_val_denom)
     
     print("Integral ",low_bin," to ", high_bin, " from single  : ", itg_num, " and triple : ", itg_denom)
     
     print("Difference ITG(simple) - ITG(triple) = ", itg_num - itg_denom)


     hist_num.GetYaxis().SetTitle("#events/tot")
     hist_num.GetYaxis().SetTitleOffset(2.8)
     hist_num.GetYaxis().SetTitleFont(43)
     hist_num.GetYaxis().SetTitleSize(20)
     
     hist_num.SetLineColor(4)
     hist_num.SetLineWidth(2)
     hist_num.GetXaxis().SetTitle("MASS")
     hist_num.GetXaxis().SetTitleSize(15)
     hist_num.GetXaxis().SetTitleFont(43)
     hist_num.GetXaxis().SetTitleOffset(2)
     hist_num.Write()
  

     hist_denom.GetYaxis().SetTitle("#events/tot")
     hist_denom.GetYaxis().SetTitleOffset(2.8)
     hist_denom.GetYaxis().SetTitleFont(43)
     hist_denom.GetYaxis().SetTitleSize(20)

     hist_denom.GetXaxis().SetTitle("MASS")
     hist_denom.GetXaxis().SetTitleSize(15)
     hist_denom.GetXaxis().SetTitleFont(43)
     hist_denom.GetXaxis().SetTitleOffset(1)

     hist_denom.SetLineColor(2)
     hist_denom.SetLineWidth(2)
     hist_denom.Write()
     

     c =  TCanvas("c","canvas",800,800)
     pad1 =  TPad("pad1","pad1",0,0.5,1,0.95)
     pad1.SetLogy()

     pad1.SetGridx()
     pad1.Draw()
     pad1.cd()

     labellumi = TPaveLabel(0.6,0.95,0.5,1,"SingleMu 2018A, 13.6 TeV","brNDC")
     labellumi.SetFillColor(0)
     labellumi.SetTextColor(1)
     labellumi.SetFillStyle(0)
     labellumi.SetBorderSize(0)
     labellumi.SetTextSize(0.25)
     labellumi.SetTextAlign(12)
     labellumi.Draw("same")

     hist_num.SetStats(0)
     hist_num.Draw()
     hist_denom.Draw("same")


     lej_brut = TLegend(0.5,0.8,0.6,0.9)
     lej_brut.AddEntry(hist_num,"Single","l")
     lej_brut.AddEntry(hist_denom,"Multiple","l")
     lej_brut.Draw("same")

     axis = hist_num.GetYaxis()
     axis.ChangeLabel(1,-1,-1,-1,-1,-1,"")
     axis.SetLabelFont(43)
     axis.SetLabelSize(15)
     axis.SetTitleOffset(1.5)

     c.cd()
     pad2 = TPad("pad2","pad2",0,0.3,1,0.5)
     #pad2.SetTopMargin(0)
     pad2.SetBottomMargin(0.1)
     pad2.SetLogy()
     pad2.SetGridx()
     pad2.Draw()
     pad2.cd()
 
     ratio_plot = hist_num.Clone("h3")
     
     #ratio_plot.SetLineColor(1)
     ratio_plot.SetTitle("")


          
     ratio_plot.GetYaxis().SetTitle(" blue / red ")
     ratio_plot.GetYaxis().SetNdivisions(505)
     ratio_plot.GetYaxis().SetTitleSize(15)
     ratio_plot.GetYaxis().SetTitleFont(43)
     ratio_plot.GetYaxis().SetTitleOffset(2)
     ratio_plot.GetYaxis().SetLabelFont(43)
     ratio_plot.GetYaxis().SetLabelSize(15)
     ratio_plot.GetXaxis().SetTitle("")
     ratio_plot.GetXaxis().SetTitleSize(15)
     ratio_plot.GetXaxis().SetTitleFont(43)
     ratio_plot.GetXaxis().SetTitleOffset(2)
     ratio_plot.GetXaxis().SetLabelFont(43)
     ratio_plot.GetXaxis().SetLabelSize(15)
     ratio_plot.SetMinimum(0.3)
     ratio_plot.SetMaximum(3)
     ratio_plot.Sumw2()
     ratio_plot.SetStats(0)
     ratio_plot.Divide(hist_denom)
     ratio_plot.SetMarkerStyle(2)
     ratio_plot.SetMarkerSize(1)
     ratio_plot.SetLineWidth(1)
     ratio_plot.SetLineColor(15)
     ratio_plot.SetMarkerColor(15)
    
     ratio_plot.Draw("ep")

     line = TLine(0,1,1000,1)
     line.SetLineStyle(4)

     
     line_q1 = TLine(200,0.3,200,3)
     line_q1.SetLineStyle(3)
     line_q1.SetLineColor(2)
     line_q1.SetLineWidth(2)
     
     line_q2 = TLine(1000,0.3,1000,3)
     line_q2.SetLineStyle(3)
     line_q2.SetLineColor(2)
     line_q2.SetLineWidth(2)
     
     line_q1.Draw("same")
     line_q2.Draw("same")
     line.Draw("same")
     
     c.cd()
     pad3 = TPad("pad3","pad3",0,0.07,1,0.27)
     #pad3.SetTopMargin(0)
     #pad3.SetLogy()
     #pad3.SetBottomMargin(0.1)
     
     pad3.SetGridx()
     pad3.Draw()
     pad3.cd()

     ratio_itg = TH1D("ratio_itg","ratio_itg",200,0,1000)
     ratio_itg.Sumw2()
     ratio_itg.SetTitle("")
     ratio_itg.GetYaxis().SetTitle("Right Integrals blue / red ")
     ratio_itg.GetYaxis().SetNdivisions(505)
     ratio_itg.GetYaxis().SetTitleSize(15)
     ratio_itg.GetYaxis().SetTitleFont(43)
     ratio_itg.GetYaxis().SetTitleOffset(2)
     ratio_itg.GetYaxis().SetLabelFont(43)
     ratio_itg.GetYaxis().SetLabelSize(15)

     ratio_itg.SetMinimum(0.8)
     ratio_itg.SetMaximum(1.3)
     ratio_itg.SetMarkerStyle(2)
     ratio_itg.SetMarkerSize(1)
     
     print("there are ", ratio_itg.GetNbinsX(), " x bins")
     print("ITG num bins 200-200 = ", hist_num.Integral(200,200))
     print("ITG den bins 200-200 = ", hist_denom.Integral(200,200))
     print("ITG num bins 0-0 = ", hist_num.Integral(0,0))
     print("ITG den bins 0-0 = ", hist_denom.Integral(0,0))
     for i in range(0,ratio_itg.GetNbinsX()+1):
         it = hist_num.Integral(i,50)
         #print("it num",i, " = ", it)
         it2 = hist_denom.Integral(i,50)
         #print("it den",i, " = ", it2)
         if(it2 !=0):
             ratio_itg.SetBinContent(i,it/it2)
             print("bin #",i, " gets ratio = ",it, " / ", it2, " = ", it/it2)
             #SetBinError -> propagate errors each number


     ratio_itg.GetYaxis().SetNdivisions(505)
     
     ratio_itg.GetYaxis().SetLabelFont(43)
     ratio_itg.GetYaxis().SetLabelSize(15)
     ratio_itg.GetYaxis().SetTitleOffset(1.2)
     ratio_itg.GetYaxis().SetTitleSize(20)

     ratio_itg.GetXaxis().SetLabelFont(43)
     ratio_itg.GetXaxis().SetLabelSize(15)
     #ratio_itg.GetXaxis().SetTitleOffset(1)
     ratio_itg.GetXaxis().SetTitle("IAS")
     ratio_itg.GetXaxis().SetTitleSize(20)

     ratio_itg.SetStats(0)

     ratio_itg.Draw("ep")
     line2 = TLine(0,1,1000,1)
     line2.SetLineStyle(4)
     line2.Draw("same")

     line3 = TLine(200,0.8,200,1.3)
     line3.SetLineStyle(3)
     line3.SetLineColor(2)
     line3.SetLineWidth(2)
     line3.Draw("same")
     
     line4 = TLine(1000,0.8,1000,1.3)
     line4.SetLineStyle(3)
     line4.SetLineColor(2)
     line4.SetLineWidth(2)
     line4.Draw("same")
     c.Write()


     c.SaveAs(hname+"_over_"+hnameden+".pdf")
     hfile.Write()

