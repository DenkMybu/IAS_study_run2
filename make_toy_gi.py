import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH1D,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3
import sys
import argparse



def pull(ntoys,n_tot_entries,nbins,do_chi2):
    method_prefix = TString("Log-Likelihood ")
    if (do_chi2):
        method_prefix="#chi^{2} "

    #Create histo
    h4 = TH1F(method_prefix+"h4",method_prefix+" Random Gauss",nbins,-4,4);
    h4.SetMarkerStyle(21)
    h4.SetMarkerSize(0.8)
    h4.SetMarkerColor(kRed)


    #Histogram for sigma and pull
    sigma = TH1F(method_prefix+"sigma",method_prefix+"sigma from gaus fit",50,0.5,1.5)
    pull = TH1F(method_prefix+"pull",method_prefix+"pull from gaus fit",50,-4.,4.)

    c0 = TCanvas(method_prefix+"Gauss",ethod_prefix+"Gauss",0,0,320,240)
    c0.SetGrid()
    c0.cd()

    sig,mean = 0,0
    for i in range (ntoys):
       h4.Reset()
       for u in range(n_tot_entries):
           h4.Fill(gRandom.Gaus())
       if(do_chi2):
           h4.Fit("gaus","q")
       else:
           h4.Fit("gaus","lq")

       if i%100 == 0:
           h4.Draw("ep")
           c0.Update()


       myfit = h4.GetFunction("gaus")

       sig = myfit.GetPArameter(2)
       mean = myfit.GetParameter(1)
       sigma.Fill(sig)
       pull.Fill(mean/sig * sqrt(n_tot_entries))
       
    c1.cd()
    pull.DrawClone()

#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 propdraw.py template_2018B_15mars_selection_5_bin_eta_2p1_pmin_10.root

def poissonHisto(h,p,RNG,binx,biny):
    for k in range(p.GetNbinsZ()+1):
        h.SetBinContent(binx,biny,k,RNG.Poisson(p.GetBinContent(binx,biny,k)))


def toy_charge(p,RNG,binx,biny,binz):      
    return RNG.Poisson(p.GetBinContent(binx,biny,binz))

def toy_charge_1d(h,p,RNG):
    for i in range(p.GetNbinsX()+1):  
        h.SetBinContent(i,RNG.Poisson(p.GetBinContent(i)))



def get_contents(p):
    for i in range(p.GetNbinsX()+1):
        print("Bin X #",i, "\n")
        for j in range(p.GetNbinsY()+1):
            print("Bin Y #",j,"\n")
            for k in range(p.GetNbinsZ()+1):
                if p.GetBinContent(i,j,k) != 0:
                    print("X : ", i , " , Y : ", j , ", Z :", k, " has bin content = ", p.GetBinContent(i,j,k), "\n")

P_Min = 1
P_Max = 16 # 1 + 14 + 1; final one is for pixel!
P_NBins = 15 #15th bin = pixel; 0 is underflow
Path_Min = 0.2
Path_Max = 1.6
Path_NBins = 42
Charge_Min = 0
Charge_Max = 5000
Charge_NBins = 500 


if __name__ == "__main__":

     file_trf = sys.argv[1]
     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", sys.argv[1])
     else:
         print("Successfully opened root file containing templates")


     loaded_template = file0.Get("Charge_Vs_Path_PU_between_0_20")
     if not loaded_template:
         print("could not open histogram Charge_Vs_Path_PU_between_0_20")
     else:
         print("Successfully opened histogram")


     #get_contents(loaded_template)
     num_binx = loaded_template.GetNbinsX()
     num_biny = loaded_template.GetNbinsY()
     num_binz = loaded_template.GetNbinsZ()

     n_toy = 100

     hfile = TFile( 'Toy_templates_2018A.root', 'RECREATE' )

     toy_templates = [0] * n_toy
     charge_projection = [0] * n_toy
   
     charge_toy = [0] * n_toy 


     x_study = 13
     y_study = 11

     RNG = TRandom3()

             

     # k = sys 1 (0.1 etc..) 
     for i in range(n_toy):


         name = "toy_binx_4_biny_8_prod_" + str(i)
         name_charge = "charge_given_x_y_bins" + str(i)
         toy_templates[i] = TH3D (name,name,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         charge_projection[i] = TH1D (name_charge,name_charge, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(toy_templates[i],loaded_template,RNG,x_study,y_study)
         print("Done generating Gi toys #",i, " with set x bin : 4, and set y bin :8 ")         
         


         for k in range(num_binz+1):
             tmp = toy_charge(loaded_template,RNG,x_study,y_study,k)
             if tmp != 0:
                 print("random poisson drawn from charge distribution = ", tmp)

             charge_projection[i].Fill(tmp)
         


     for i in range(n_toy):
         for j in range(0,num_binx+1):
             for k in range(0,num_biny+1):
                 toy_templates[i].GetXaxis().SetRange(j,j)
                 toy_templates[i].GetYaxis().SetRange(k,k)
                  
     ''' 
     #Test projection for bin x = 13 and biny = 11
     clone = loaded_template.Clone("Test_projection_charge")
     for i in range(0,num_binx+1):
         for p in range(0,num_biny+1):
             clone.GetXaxis().SetRange(i,i)
             clone.GetYaxis().SetRange(p,p)
             name = "projection_xbin_" + str(i) + "_ybin_" + str(p)
             test = TH1D(name,name,500,0,5000)
             test = clone.Project3D("z")
             test.SetName(name)
             if(test.GetEntries() != 0):
                 print(name, " has an integral z entries = ", test.Integral())
                 test.Write()
                 for u in range(n_toy):
                     name_charge_toy = "poisson_charge_binx_" + str(i) + "_biny_" + str(p) + "_toy" + str(u)
                     charge_toy[u] = TH1D(name_charge_toy,name_charge_toy,500,0,5000)
                     toy_charge_1d(charge_toy[u],test,RNG)
                        
     '''

     hfile.Write()
