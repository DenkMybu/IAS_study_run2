#define run2analysis_cxx
#include "run2analysis.h"
#include <TMatrix.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <TEfficiency.h>
#include "TRandom3.h"
#include <TRatioPlot.h>
#include <THStack.h>

enum TrackQuality {
      undefQuality = -1,
      loose = 0,
      tight = 1,
      highPurity = 2,
      confirmed = 3,      // means found by more than one iteration
      goodIterative = 4,  // meaningless
      looseSetWithPV = 5,
      highPuritySetWithPV = 6,
      discarded = 7,  // because a better track found. kept in the collection for reference....
      qualitySize = 8
    };


bool writeTemplateOnDisk = false;
bool writeTptHSCP = false;

bool computeSpecial= false;
bool boolDeDxTemp= false;


bool UseTemplatesForPUReweighting = false; //When running over all tracks, but no iso available
bool get_list_of_bins = false; //computes first all events and their NPV distributions to get equal statistic bins

bool compute_PE = false; //compute all 100 pseudo experiments (nPE can be changed, but need to use python script propdraw.py
bool UsePURwtHSCP = (!writeTptHSCP);//same but when running on HSCPs

bool TemplateIso = false;
bool TemplateNoIso = (!TemplateIso);
 
bool StudyIso = true;
bool StudyNoIso = (!StudyIso);

string low_bound = "20";
string high_bound = "48";
// Modification of the code in September 2021 to 
// align it with the use of xtalk inversion (only for cluster cleaning) and saturation 
//
void run2analysis::Loop(int year, TString Letter, bool dataFlag=true)
{
//   In a ROOT session, you can do:
//      root> .L run2analysis.C
//      root> run2analysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   Long64_t nbytes_pu = 0, nbpu = 0;

   double pi=acos(-1);
   double dEdxSF [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF
   
   if (writeTptHSCP){
       std::cout << "Generating templates...." << endl;
       std::cout << "P range [" <<low_bound << "-"<< high_bound << "], with 8th september preselection (except iso)" << endl;
   }
   std::string a = "iso tk 15";
   std::string b = a;
   if(!StudyIso) a = "no iso";
   if(!TemplateIso) b = "no iso";
   std::cout << "Study will use template with " << b << " and study with " << a;


   // values for 2017 UL data & MC
   if (dataFlag) {
      if (year==2017) {
      dEdxSF[0]=1.;
      dEdxSF[1]=1.0325;
      }
      else if (year==2018) { 
      dEdxSF[0]=1.;
      dEdxSF[1]=1.0817;
      }
      else {
         cout << " AIE AIE AIE : Year not well specified for dEdXSF in data !!!! " << endl;
         return;
      }
   }
   else {
      if (year==2017) {
        dEdxSF[0]=1.0079;
        dEdxSF[1]=1.0875;
      } 
      else if (year==2018) {
        dEdxSF[0]=1.0047;
        dEdxSF[1]=1.1429;
      }
      else {
         cout << " AIE AIE AIE : Year not well specified for dEdXSF in MC !!!! " << endl;
         return;
      }
   }
   cout << "  dEdXSF : " << dEdxSF[0] << "   &   " << dEdxSF[1] << endl;

   // update all on feb 25 (on jan21 rootfile)
   // with Cval extracted on 3-5 
   // MC 2017
   float Kval_ld = 2.33;
   float Cval_ld =3.37 ;
   float Kval_ldstrip = 2.45;
   float Cval_ldstrip = 3.41;
   float Kval_all = 2.19;
   float Cval_all = 3.18;

   float Cval_nol1 = 3.31; // fit C in 5-25
   float Kval_nol1 = 2.20; 
   float Cval_nol1_2 = 3.22; // fit C in 3-5
   float Kval_nol1_2 = 2.26;
   float Cval_nol1_3 = 2.88; // with log term
   float Kval_nol1_3 = 2.57;
   float Nval_nol1_3 = 0.099;

   float Kval_strip = 2.43;
   float Cval_strip = 3.23;
   float Kval_hdnol1 = 2.20; 
   float Cval_hdnol1 = 3.07; 

   float Kval_pix = 1.75;
   float Cval_pix = 3.03;
   float Kval_pixnol1 = 1.62;
   float Cval_pixnol1 = 3.13;

   if (dataFlag) {
    if (year==2017) {
     // extracted on analysis_ul_2017_21jan.root
     Kval_ld = 2.36 ;
     Cval_ld = 3.37 ;
     Kval_ldstrip = 2.58 ;
     Cval_ldstrip = 3.37 ;
     Kval_all = 2.19 ;
     Cval_all = 3.10 ;

     Cval_nol1 = 3.19 ; // fit C in 5-20
     Kval_nol1 = 2.28 ;
     Cval_nol1_2 = 3.17 ; // fit C in 3-5
     Kval_nol1_2 = 2.30 ;
     Cval_nol1_3 = 2.80 ;
     Kval_nol1_3 = 2.83 ;
     Nval_nol1_3 = 0.090 ;

     Kval_strip = 2.50 ;
     Cval_strip = 3.19 ;
     Kval_hdnol1 = 2.23 ;
     Cval_hdnol1 = 3.03 ;

/* not possible to extract K and C values on the pixel only samples  ==> do not update and keep the MC values
     Kval_pix = 2.173 ;
     Cval_pix = 3.242 ;
     Kval_pixnol1 = 2.173 ;
     Cval_pixnol1 = 3.242 ;
*/
    }
    else if (year==2018) {
     //extracted on analysis_ul_2018_28feb.root
     Kval_ld = 2.35 ;
     Cval_ld = 3.36 ;
     Kval_ldstrip = 2.60 ;
     Cval_ldstrip = 3.38 ;
     Kval_all = 2.17 ;
     Cval_all = 3.09 ;

     Cval_nol1 = 3.18 ; // fit C in 5-20
     Kval_nol1 = 2.25 ;
     Cval_nol1_2 = 3.16 ; // fit C in 3-5
     Kval_nol1_2 = 2.27 ;
     Cval_nol1_3 = 2.90 ;
     Kval_nol1_3 = 2.52 ;
     Nval_nol1_3 = 0.066 ;

     Kval_strip = 2.50 ;
     Cval_strip = 3.19 ;
     Kval_hdnol1 = 2.21 ;
     Cval_hdnol1 = 3.02 ;

// not possible to extract K and C values on the pixel only samples  ==> do not update and keep the MC values
     Kval_pix = 1.72 ;
     Cval_pix = 2.96 ;
// MC 2018 not available 
/*
     Kval_pixnol1 = .. ;
     Cval_pixnol1 = .. ;
*/
    }
   }
   else if (year==2018) {
     // MC : extracted on analysis_ul_2018MC_w18_MC_28feb.root
     Kval_ld = 2.34 ;
     Cval_ld = 3.35 ;
     Kval_ldstrip = 2.48 ;
     Cval_ldstrip = 3.38 ;
     Kval_all = 2.20 ;
     Cval_all = 3.14 ;

     Cval_nol1 = 3.29 ; // fit C in 5-20
     Kval_nol1 = 2.22 ;
     Cval_nol1_2 = 3.16 ; // fit C in 3-5
     Kval_nol1_2 = 2.27 ;
     Cval_nol1_3 = 2.86 ;
     Kval_nol1_3 = 2.61 ;
     Nval_nol1_3 = 0.103 ;

     Kval_strip = 2.40 ;
     Cval_strip = 3.21 ;
     Kval_hdnol1 = 2.23 ;
     Cval_hdnol1 = 3.05 ;

// MC values : no convergence... oscillation between 1.64 and 1.72 for Kval_pix :(
     Kval_pix = 1.72 ;
     Cval_pix = 2.96 ;
// MC values : not possible to fit --> K =1 at limit :(
/*
     Kval_pixnol1 = .. ;
     Cval_pixnol1 = .. ;
*/
   }
   bool blind_data= false;
   if (dataFlag)  blind_data=true;


// histograms
   TH1D* NHSCP = new TH1D("N_HSCP_POST_PRESEL", "N_HSCP_POST_PRESEL", 3, 0,3);
   TH1D* NTRK = new TH1D("N_TRK_POST_PRESEL", "N_TRK_POST_PRESEL", 3, 0,3);
   TH1D* N_CLU_HSCP = new TH1D("N_CLU_HSCP_POST_PRESEL", "N_CLU_HSCP_POST_PRESEL", 100, 0,100);
   TH1D* N_CLU_TRK = new TH1D("N_CLU_TRK_POST_PRESEL", "N_CLU_TRK_POST_PRESEL", 100, 0,100);

   TH1D* HNtracks = new TH1D("HNtracks", "HNtracks", 50, -0.5,49.5);
   TH1D* HNtracks1 = new TH1D("HNtracks1", "Ntracks with pT>1", 40, -0.5,39.5);
   TH1D* HNtracks20 = new TH1D("HNtracks20", "Ntracks with pT>20", 30, -0.5,29.5);

   TH1D* Htrackpt = new TH1D("Htrackpt", "track pT", 500, 0.,500);
   TH1D* Htracketa = new TH1D("Htracketa", "track eta", 24, -3.,3.);
   TH1D* Htracketa_lowp = new TH1D("Htracketa_lowp", "track eta", 24, -3.,3.);
   TH1D* Htrackphi = new TH1D("Htrackphi", "track phi", 24, -1.*pi,pi);
   TH1D* Htracknhit = new TH1D("Htracknhit", "track nhit", 50, -0.5,49.5);

   TH1D* Htrackih_reco = new TH1D("Htrackih_reco", "track ih (reco)", 50, 0.,20);
   TH1D* Htrackih_pix = new TH1D("Htrackih_pix", "track ih in pix", 50, 0.,20);
   TH1D* Htrackih_strip = new TH1D("Htrackih_strip", "track ih in strip", 50, 0.,20);
   TH1D* Htrackdedx_pix = new TH1D("Htrackdedx_pix", "track dedx in pix", 100, 0.,500000);
   TH1D* Htrackdedx_strip = new TH1D("Htrackdedx_strip", "track dedx in strip", 100, 0.,2000);
   TH1D* Htrackias = new TH1D("Htrackias", "track ias (noL1)", 80, 0.,1);
   TH1D* Htrackiasall = new TH1D("Htrackiasall", "track ias (all)", 80, 0.,1);

   TH1D* Htrackih_lowp = new TH1D("Htrackih_lowp", "track ih", 50, 0.,20);
   TH1D* Htrackih_pix_lowp = new TH1D("Htrackih_pix_lowp", "track ih in pix", 50, 0.,20);
   TH1D* Htrackih_strip_lowp = new TH1D("Htrackih_strip_lowp", "track ih in strip", 50, 0.,20);
   TH1D* Htrackih0_lowp = new TH1D("Htrackih0_lowp", "track ih0 ", 50, 0.,20);
   TH1D* Htrackih0noL1_lowp = new TH1D("Htrackih0noL1_lowp", "track ih0 noL1", 50, 0.,20);
   TH1D* Htrackias_lowp = new TH1D("Htrackias_lowp", "track ias (noL1)", 80, 0.,1);
   TH1D* Htrackiasall_lowp = new TH1D("Htrackiasall_lowp", "track ias (all)", 80, 0.,1);

   TH1D* Htrackdedx_pix_lowp = new TH1D("Htrackdedx_pix_lowp", "track dedx in pix", 100, 0.,500000);
   TH1D* Htrackdedx_strip_lowp = new TH1D("Htrackdedx_strip_lowp", "track dedx in strip", 100, 0.,2000);
   TH1D* Htrackdedx_strip_lowp1 = new TH1D("Htrackdedx_strip_lowp1", "track dedx in strip", 100, 0.,20);
   TH1D* Htrackdedx_strip_lowp2 = new TH1D("Htrackdedx_strip_lowp2", "track dedx in strip", 100, 0.,20);

   TH1D* Nsat = new TH1D("Nsat", "Nsat", 20, -0.5,19.5);
   TH1D* NPix = new TH1D("NPix", "NPix", 10, -0.5,9.5);
   TH1D* NStrip = new TH1D("NStrip", "NStrip", 30, -0.5,29.5);
   TH2D* dEdXVsP = new TH2D("dEdXVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXVsP_lowp = new TH2D("dEdXVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXVsP_lowp2 = new TH2D("dEdXVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdXpixVsP = new TH2D("dEdXpixVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXpixVsP_lowp = new TH2D("dEdXpixVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXpixVsP_lowp2 = new TH2D("dEdXpixVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdXstripVsP = new TH2D("dEdXstripVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXstripVsP_lowp = new TH2D("dEdXstripVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXstripVsP_lowp2 = new TH2D("dEdXstripVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0VsP_lowp = new TH2D("dEdX0VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0VsP_lowp2 = new TH2D("dEdX0VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0stripVsP_lowp = new TH2D("dEdX0stripVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0stripVsP_lowp2 = new TH2D("dEdX0stripVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_lowp = new TH2D("dEdX0noL1VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_lowp2 = new TH2D("dEdX0noL1VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta1_lowp = new TH2D("dEdX0noL1VsP_eta1_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta1_lowp2 = new TH2D("dEdX0noL1VsP_eta1_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta2_lowp = new TH2D("dEdX0noL1VsP_eta2_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta2_lowp2 = new TH2D("dEdX0noL1VsP_eta2_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta3_lowp = new TH2D("dEdX0noL1VsP_eta3_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta3_lowp2 = new TH2D("dEdX0noL1VsP_eta3_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu1_lowp = new TH2D("dEdX0noL1VsP_pu1_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu1_lowp2 = new TH2D("dEdX0noL1VsP_pu1_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu2_lowp = new TH2D("dEdX0noL1VsP_pu2_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu2_lowp2 = new TH2D("dEdX0noL1VsP_pu2_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu3_lowp = new TH2D("dEdX0noL1VsP_pu3_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu3_lowp2 = new TH2D("dEdX0noL1VsP_pu3_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu4_lowp = new TH2D("dEdX0noL1VsP_pu4_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu4_lowp2 = new TH2D("dEdX0noL1VsP_pu4_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdXHDnoL1VsP_lowp = new TH2D("dEdXHDnoL1VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXHDnoL1VsP_lowp2 = new TH2D("dEdXHDnoL1VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0pixnoL1VsP_lowp = new TH2D("dEdX0pixnoL1VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0pixnoL1VsP_lowp2 = new TH2D("dEdX0pixnoL1VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEstrVsdE_lowp = new TH2D("dEstrVsdE_lowp", "dEdX(all):dEdX(strip)", 50,2.,7., 80, 2.,10.);
   TH2D* dEdXstripVsEta_lowp = new TH2D("dEdXstripVsEta_lowp", "dEdX:P", 12, -3.,3., 80, 2.,10.);
   TH2D* EtaVsPhi_nhit = new TH2D("EtaVsPhi_nhit", "Eta:Phi", 12, -1.*pi,pi,12, -3.,3.);


   TH2D* dEdXstripVsNhit_lowp  = new TH2D("dEdXstripVsNhit_lowp", "dEdX(strip):Nhit", 20, -0.5,19.5,80, 2.,10.);
   TH2D* dEdXstripVsNhittrunc_lowp  = new TH2D("dEdXstripVsNhittrunc_lowp", "dEdX(strip):Nhit", 20, -0.5,19.5,80, 2.,10.);
   TH2D* dEdXstripVsCharge_lowp = new TH2D("dEdXstripVsCharge_lowp", "dEdX(strip):Charge", 100, 0.,20, 80, 2.,10.);

   TH1D* Charge_pixl1 = new TH1D("Charge_pixl1","Charge per layer",400, 0.,20);
   TH1D* Charge_pixl2 = new TH1D("Charge_pixl2","Charge per layer",400, 0.,20);
   TH1D* Charge_pixl3 = new TH1D("Charge_pixl3","Charge per layer",400, 0.,20);
   TH1D* Charge_pixl4 = new TH1D("Charge_pixl4","Charge per layer",400, 0.,20);
   TH1D* Charge_pixd1 = new TH1D("Charge_pixd1","Charge per layer",400, 0.,20);
   TH1D* Charge_pixd2 = new TH1D("Charge_pixd2","Charge per layer",400, 0.,20);
   TH1D* Charge_pixd3 = new TH1D("Charge_pixd3","Charge per layer",400, 0.,20);
   TH1D* Charge_pixr1 = new TH1D("Charge_pixr1","Charge per layer",400, 0.,20);
   TH1D* Charge_pixr2 = new TH1D("Charge_pixr2","Charge per layer",400, 0.,20);
   TH1D* Charge_tib1 = new TH1D("Charge_tib1","Charge per layer",400, 0.,20);
   TH1D* Charge_tib2 = new TH1D("Charge_tib2","Charge per layer",400, 0.,20);
   TH1D* Charge_tib3 = new TH1D("Charge_tib3","Charge per layer",400, 0.,20);
   TH1D* Charge_tib4 = new TH1D("Charge_tib4","Charge per layer",400, 0.,20);
   TH1D* Charge_tob1 = new TH1D("Charge_tob1","Charge per layer",400, 0.,20);
   TH1D* Charge_tob2 = new TH1D("Charge_tob2","Charge per layer",400, 0.,20);
   TH1D* Charge_tob3 = new TH1D("Charge_tob3","Charge per layer",400, 0.,20);
   TH1D* Charge_tob4 = new TH1D("Charge_tob4","Charge per layer",400, 0.,20);
   TH1D* Charge_tob5 = new TH1D("Charge_tob5","Charge per layer",400, 0.,20);
   TH1D* Charge_tob6 = new TH1D("Charge_tob6","Charge per layer",400, 0.,20);
   TH1D* Charge_tid1 = new TH1D("Charge_tid1","Charge per layer",400, 0.,20);
   TH1D* Charge_tid2 = new TH1D("Charge_tid2","Charge per layer",400, 0.,20);
   TH1D* Charge_tid3 = new TH1D("Charge_tid3","Charge per layer",400, 0.,20);
   TH1D* Charge_tec1 = new TH1D("Charge_tec1","Charge per layer",400, 0.,20);
   TH1D* Charge_tec2 = new TH1D("Charge_tec2","Charge per layer",400, 0.,20);
   TH1D* Charge_tec3 = new TH1D("Charge_tec3","Charge per layer",400, 0.,20);
   TH1D* Charge_tec4 = new TH1D("Charge_tec4","Charge per layer",400, 0.,20);
   TH1D* Charge_tec5 = new TH1D("Charge_tec5","Charge per layer",400, 0.,20);
   TH1D* Charge_tec6 = new TH1D("Charge_tec6","Charge per layer",400, 0.,20);
   TH1D* Charge_tec7 = new TH1D("Charge_tec7","Charge per layer",400, 0.,20);
   TH1D* Charge_tec8 = new TH1D("Charge_tec8","Charge per layer",400, 0.,20);
   TH1D* Charge_tec9 = new TH1D("Charge_tec9","Charge per layer",400, 0.,20);

   TH1D* LowpCharge_tib1 = new TH1D("LowpCharge_tib1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tib2 = new TH1D("LowpCharge_tib2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tib3 = new TH1D("LowpCharge_tib3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tib4 = new TH1D("LowpCharge_tib4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob1 = new TH1D("LowpCharge_tob1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob2 = new TH1D("LowpCharge_tob2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob3 = new TH1D("LowpCharge_tob3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob4 = new TH1D("LowpCharge_tob4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob5 = new TH1D("LowpCharge_tob5","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob6 = new TH1D("LowpCharge_tob6","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tid1 = new TH1D("LowpCharge_tid1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tid2 = new TH1D("LowpCharge_tid2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tid3 = new TH1D("LowpCharge_tid3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec1 = new TH1D("LowpCharge_tec1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec2 = new TH1D("LowpCharge_tec2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec3 = new TH1D("LowpCharge_tec3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec4 = new TH1D("LowpCharge_tec4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec5 = new TH1D("LowpCharge_tec5","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec6 = new TH1D("LowpCharge_tec6","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec7 = new TH1D("LowpCharge_tec7","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec8 = new TH1D("LowpCharge_tec8","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec9 = new TH1D("LowpCharge_tec9","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl1 = new TH1D("LowpCharge_pixl1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl2 = new TH1D("LowpCharge_pixl2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl3 = new TH1D("LowpCharge_pixl3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl4 = new TH1D("LowpCharge_pixl4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixd1 = new TH1D("LowpCharge_pixd1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixd2 = new TH1D("LowpCharge_pixd2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixd3 = new TH1D("LowpCharge_pixd3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixr1 = new TH1D("LowpCharge_pixr1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixr2 = new TH1D("LowpCharge_pixr2","Charge per layer",400, 0.,20);
   TH2D* LowpCharge_Eta_pix = new TH2D("LowpCharge_Eta_pix", "Charge/pathlength:eta", 12, -3.,3., 200, 0.,10);
   TH2D* LowpCharge_Eta_strip = new TH2D("LowpCharge_Eta_strip", "Charge/pathlength:eta", 12, -3.,3., 200, 0.,10);

   TH2D* ChargeVsRun_pixl1 = new TH2D("ChargeVsRun_pixl1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixl2 = new TH2D("ChargeVsRun_pixl2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixl3 = new TH2D("ChargeVsRun_pixl3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixl4 = new TH2D("ChargeVsRun_pixl4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixd1 = new TH2D("ChargeVsRun_pixd1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixd2 = new TH2D("ChargeVsRun_pixd2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixd3 = new TH2D("ChargeVsRun_pixd3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixr1 = new TH2D("ChargeVsRun_pixr1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixr2 = new TH2D("ChargeVsRun_pixr2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib1 = new TH2D("ChargeVsRun_tib1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib2 = new TH2D("ChargeVsRun_tib2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib3 = new TH2D("ChargeVsRun_tib3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib4 = new TH2D("ChargeVsRun_tib4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob1 = new TH2D("ChargeVsRun_tob1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob2 = new TH2D("ChargeVsRun_tob2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob3 = new TH2D("ChargeVsRun_tob3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob4 = new TH2D("ChargeVsRun_tob4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob5 = new TH2D("ChargeVsRun_tob5","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob6 = new TH2D("ChargeVsRun_tob6","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tid1 = new TH2D("ChargeVsRun_tid1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tid2 = new TH2D("ChargeVsRun_tid2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tid3 = new TH2D("ChargeVsRun_tid3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec1 = new TH2D("ChargeVsRun_tec1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec2 = new TH2D("ChargeVsRun_tec2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec3 = new TH2D("ChargeVsRun_tec3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec4 = new TH2D("ChargeVsRun_tec4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec5 = new TH2D("ChargeVsRun_tec5","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec6 = new TH2D("ChargeVsRun_tec6","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec7 = new TH2D("ChargeVsRun_tec7","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec8 = new TH2D("ChargeVsRun_tec8","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec9 = new TH2D("ChargeVsRun_tec9","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);

   TH2D* ZooChargeVsRun_pixl1 = new TH2D("ZooChargeVsRun_pixl1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl2 = new TH2D("ZooChargeVsRun_pixl2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl3 = new TH2D("ZooChargeVsRun_pixl3","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl4 = new TH2D("ZooChargeVsRun_pixl4","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd1 = new TH2D("ZooChargeVsRun_pixd1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd2 = new TH2D("ZooChargeVsRun_pixd2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd3 = new TH2D("ZooChargeVsRun_pixd3","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixr1 = new TH2D("ZooChargeVsRun_pixr1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixr2 = new TH2D("ZooChargeVsRun_pixr2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_tib1 = new TH2D("ZooChargeVsRun_tib1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tib2 = new TH2D("ZooChargeVsRun_tib2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tib3 = new TH2D("ZooChargeVsRun_tib3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tib4 = new TH2D("ZooChargeVsRun_tib4","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob1 = new TH2D("ZooChargeVsRun_tob1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob2 = new TH2D("ZooChargeVsRun_tob2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob3 = new TH2D("ZooChargeVsRun_tob3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob4 = new TH2D("ZooChargeVsRun_tob4","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob5 = new TH2D("ZooChargeVsRun_tob5","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob6 = new TH2D("ZooChargeVsRun_tob6","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tid1 = new TH2D("ZooChargeVsRun_tid1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tid2 = new TH2D("ZooChargeVsRun_tid2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tid3 = new TH2D("ZooChargeVsRun_tid3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec1 = new TH2D("ZooChargeVsRun_tec1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec2 = new TH2D("ZooChargeVsRun_tec2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec3 = new TH2D("ZooChargeVsRun_tec3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec4 = new TH2D("ZooChargeVsRun_tec4","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec5 = new TH2D("ZooChargeVsRun_tec5","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec6 = new TH2D("ZooChargeVsRun_tec6","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec7 = new TH2D("ZooChargeVsRun_tec7","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec8 = new TH2D("ZooChargeVsRun_tec8","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec9 = new TH2D("ZooChargeVsRun_tec9","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);

   TH1D* HSCP_dEdX = new TH1D("HSCP_dEdX", "dEdX", 60, 0.,15.);
   TH1D* HSCP_dEdXpix = new TH1D("HSCP_dEdXpix", "dEdX(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdXstrip = new TH1D("HSCP_dEdXstrip", "dEdX(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdX0 = new TH1D("HSCP_dEdX0", "dEdX0", 60, 0.,15.);
   TH1D* HSCP_dEdX0pix = new TH1D("HSCP_dEdX0pix", "dEdX0(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdX0strip = new TH1D("HSCP_dEdX0strip", "dEdX0(strip)", 60, 0.,15.);

   TH1D* HSCP_MassIh = new TH1D("HSCP_MassIh", "Mass via Ih", 200, 0.,4000.);
   TH1D* HSCP_MassIh0 = new TH1D("HSCP_MassIh0", "Mass via Ih(no drop)", 200, 0.,4000.);
   TH1D* HSCP_MassIhstrip = new TH1D("HSCP_MassIhstrip", "Mass via Ih(strip)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1 = new TH1D("HSCP_MassIh0noL1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2 = new TH1D("HSCP_MassIh0noL1_2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3 = new TH1D("HSCP_MassIh0noL1_3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_11 = new TH1D("HSCP_MassIh0noL1_11", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_12 = new TH1D("HSCP_MassIh0noL1_12", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_13 = new TH1D("HSCP_MassIh0noL1_13", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_1s1 = new TH1D("HSCP_MassIh0noL1_1s1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_1s2 = new TH1D("HSCP_MassIh0noL1_1s2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_1s3 = new TH1D("HSCP_MassIh0noL1_1s3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2s1 = new TH1D("HSCP_MassIh0noL1_2s1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2s2 = new TH1D("HSCP_MassIh0noL1_2s2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2s3 = new TH1D("HSCP_MassIh0noL1_2s3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3s1 = new TH1D("HSCP_MassIh0noL1_3s1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3s2 = new TH1D("HSCP_MassIh0noL1_3s2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3s3 = new TH1D("HSCP_MassIh0noL1_3s3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0strip = new TH1D("HSCP_MassIh0strip", "Mass via Ih(no drop strip)", 200, 0.,4000.);
   TH1D* HSCP_MassIhHDnoL1 = new TH1D("HSCP_MassIhHDnoL1", "Mass via Ih(High drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassTOF = new TH1D("HSCP_MassTOF", "Mass via 1/beta", 200, 0.,4000.);
   TH2D* HSCP2d_MassTOFvsIh = new TH2D("HSCP2d_MassTOFvsIh", "Mass Ih vs 1/beta",50, 0.,4000., 50, 0.,4000.);
   TH2D* HSCP2d_MassIh = new TH2D("HSCP2d_MassIh", "Mass via Ih",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIh0 = new TH2D("HSCP2d_MassIh0", "Mass via Ih(no drop)",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIhstrip = new TH2D("HSCP2d_MassIhstrip", "Mass via Ih(strip)",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIh0noL1 = new TH2D("HSCP2d_MassIh0noL1", "Mass via Ih(no drop noL1)",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIh0strip = new TH2D("HSCP2d_MassIh0strip", "Mass via Ih(no drop strip)",30, 0, 3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIhHDnoL1 = new TH2D("HSCP2d_MassIhHDnoL1", "Mass via Ih(High drop noL1)",30,0,3000, 100, 0.,4000.);

   TH2D* HSCP2d_Mass_pix_strip15 = new TH2D("HSCP2d_Mass_pix_strip15", "Mass Pixel vs Mass Strip (15p drop)",100,0,4000, 100, 0.,4000.);
   TH2D* HSCP2d_Mass_pix_strip0 = new TH2D("HSCP2d_Mass_pix_strip0", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,4000, 100, 0.,4000.);
   TH2D* HSCP2d_Mass_pix_strip = new TH2D("HSCP2d_Mass_pix_strip", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,500, 100, 0.,500.);
   TH1D* HSCP_MassDiff_pix_strip0 = new TH1D("HSCP_MassDiff_pix_strip0", "MassDiff Pixel-Strip (no drop noL1)",200,-4000,4000);
   TH1D* HSCP_MassDiff_pix_strip15 = new TH1D("HSCP_MassDiff_pix_strip15", "MassDiff Pixel-Strip (15p drop)",200,-4000,4000);
   TH1D* HSCP_MassResol_pix_strip0 = new TH1D("HSCP_MassResol_pix_strip0", "MassResol Pixel-Strip (no drop noL1)",200,-2,2);
   TH1D* HSCP_MassResol_pix_strip15 = new TH1D("HSCP_MassResol_pix_strip15", "MassResol Pixel-Strip (15p drop)",200,-2,2);

   TH1D* lowp_MassIh = new TH1D("lowp_MassIh", "Mass via Ih", 100, 0.,5.);
   TH1D* lowp_MassIh0 = new TH1D("lowp_MassIh0", "Mass via Ih(no drop)", 100, 0.,5.);
   TH1D* lowp_MassIhstrip = new TH1D("lowp_MassIhstrip", "Mass via Ih(strip)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1 = new TH1D("lowp_MassIh0noL1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2 = new TH1D("lowp_MassIh0noL1_2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3 = new TH1D("lowp_MassIh0noL1_3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_11 = new TH1D("lowp_MassIh0noL1_11", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_12 = new TH1D("lowp_MassIh0noL1_12", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_13 = new TH1D("lowp_MassIh0noL1_13", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_1s1 = new TH1D("lowp_MassIh0noL1_1s1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_1s2 = new TH1D("lowp_MassIh0noL1_1s2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_1s3 = new TH1D("lowp_MassIh0noL1_1s3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2s1 = new TH1D("lowp_MassIh0noL1_2s1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2s2 = new TH1D("lowp_MassIh0noL1_2s2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2s3 = new TH1D("lowp_MassIh0noL1_2s3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3s1 = new TH1D("lowp_MassIh0noL1_3s1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3s2 = new TH1D("lowp_MassIh0noL1_3s2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3s3 = new TH1D("lowp_MassIh0noL1_3s3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0strip = new TH1D("lowp_MassIh0strip", "Mass via Ih(no drop strip)", 100, 0.,5.);
   TH1D* lowp_MassIhHDnoL1 = new TH1D("lowp_MassIhHDnoL1", "Mass via Ih(High drop noL1)", 100, 0.,5.);
   TH2D* lowp2d_MassIh = new TH2D("lowp2d_MassIh", "Mass via Ih", 30,0,3, 100, 0.,2.);
   TH2D* lowp2d_MassIh0 = new TH2D("lowp2d_MassIh0", "Mass via Ih(no drop)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIhstrip = new TH2D("lowp2d_MassIhstrip", "Mass via Ih(strip)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIh0noL1 = new TH2D("lowp2d_MassIh0noL1", "Mass via Ih(no drop noL1)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIh0strip = new TH2D("lowp2d_MassIh0strip", "Mass via Ih(no drop strip)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIhHDnoL1 = new TH2D("lowp2d_MassIhHDnoL1", "Mass via Ih(High drop noL1)",30,0,3, 100, 0.,5.);
   TH2D* lowp_dEdXpixVsstrip = new TH2D("lowp_dEdXpixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* lowp_dEdX0pixVsstrip = new TH2D("lowp_dEdX0pixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* lowp2d_Mass_pix_strip15 = new TH2D("lowp2d_Mass_pix_strip15", "Mass Pixel vs Mass Strip (15p drop)",100,0,5, 100, 0.,5.);
   TH2D* lowp2d_Mass_pix_strip0 = new TH2D("lowp2d_Mass_pix_strip0", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,5, 100, 0.,5.);
   TH2D* bg_lowp2d_Mass_pix_strip0 = new TH2D("bg_lowp2d_Mass_pix_strip0", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,100, 100, 0.,100.);
   TH2D* bg_transf_Mass = new TH2D("bg_transf_Mass", "Diff Mass Pixel-Strip  vs Mean (no drop noL1)",100,0,100, 100, -2.,2.);
   TH1D* lowp_MassDiff_pix_strip0 = new TH1D("lowp_MassDiff_pix_strip0", "MassDiff Pixel-Strip (no drop noL1)",200,-5,5);
   TH1D* lowp_MassDiff_pix_strip15 = new TH1D("lowp_MassDiff_pix_strip15", "MassDiff Pixel-Strip (15p drop)",200,-5,5);
   TH2D* bg_dEdXVsIas = new TH2D("bg_dEdXVsIas","Ih:Ias ", 100, 0.,1., 100, 0.,10.);
   TH1D* bg_test_event_Mass = new TH1D("bg_test_event_Mass", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* bg_test_prefiring_Mass = new TH1D("bg_test_prefiring_Mass", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* bg_test_nohem_Mass = new TH1D("bg_test_nohem_Mass", "Mass via Ih(no drop noL1)", 200, 0.,4000.);

   TH1D* HSCP_dEdXHiDrop  = new TH1D("HSCP_dEdXHiDrop", "dEdX HighDrop", 60, 0.,15.);
   TH1D* HSCP_dEdXstripHiDrop  = new TH1D("HSCP_dEdXstripHiDrop", "dEdX HighDrop(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdXpixHiDrop = new TH1D("HSCP_dEdXpixHiDrop", "dEdX HighDrop(pix)", 60, 0.,15.);
   TH1D* HSCP_dEdXHiDropNoL1  = new TH1D("HSCP_dEdXHiDropNoL1", "dEdX HighDrop No L1", 60, 0.,15.);
   TH1D* HSCP_dEdX0NoL1  = new TH1D("HSCP_dEdX0NoL1", "dEdX No L1", 60, 0.,15.);
   //CUSTOM HISTOGRAMS

   TH1D* NB_PV = new TH1D("PU_PV", "PrimaryVertices", 80, 0,80); 



   //END custom histograms
   TH1D* HSCP_FMIP4 = new TH1D("HSCP_FMIP4", "FMIP(4)", 50, 0.,1.);
   TH1D* HSCP_FMIP3p5 = new TH1D("HSCP_FMIP3p5", "FMIP(3.5)", 50, 0.,1.);
   TH1D* HSCP_FMIP3p2 = new TH1D("HSCP_FMIP3p2", "FMIP(3.2)", 50, 0.,1.);

   TH1D* HSCP_iasnol1 = new TH1D("HSCP_iasnol1", "Ias (noL1)", 80, 0.,1.);
   TH1D* HSCP_iasall = new TH1D("HSCP_iasall", "Ias (all)", 80, 0.,1.);
   TH1D* HSCP_iasstrip = new TH1D("HSCP_iasstrip", "Ias (strip)", 80, 0.,1.);
   TH1D* HSCP_iaspix = new TH1D("HSCP_iaspix", "Ias (pix)", 80, 0.,1.);
   TH1D* HSCP_probQ = new TH1D("HSCP_probQ", "ProbQ", 80, 0.,1.);
   TH1D* HSCP_probQNoL1 = new TH1D("HSCP_probQNoL1", "ProbQNoL1", 80, 0.,1.);
   TH1D* HSCP_probXY = new TH1D("HSCP_probXY", "ProbXY", 80, 0.,1.);
   TH1D* HSCP_probXYNoL1 = new TH1D("HSCP_probXYNoL1", "ProbXYNoL1", 80, 0.,1.);

   TH1D* HSCP_pt = new TH1D("HSCP_pt", "pT",  50, 55.,1550);
   TH1D* HSCP_eta = new TH1D("HSCP_eta", "eta",  24, -3.,3.);
   TH1D* HSCP_iso_eop = new TH1D("HSCP_iso_eop", "Isolation (ECAL+HCAL)/p",  50, 0., 5.);
   TH1D* nPV = new TH1D("nPV", "nPV",  40,0,80);
   TH1D* HSCP_invB = new TH1D("HSCP_invB", "invBeta",  300,-1,2);
   TH1D* HSCP_errinvB = new TH1D("HSCP_errinvB", "err_invBeta",  50,0,0.5);
   TH1D* HSCP_invBDT = new TH1D("HSCP_invBDT", "invBeta(DT)",  90,-1,2);
   TH1D* HSCP_invBCSC = new TH1D("HSCP_invBCSC", "invBeta(CSC)",  90,-1,2);
   TH1D* HSCP_time = new TH1D("HSCP_time", "VertexTiming",  200,-100,100);
   TH1D* HSCP_npix = new TH1D("HSCP_npix", "#(pix)", 15, -0.5,14.5);
   TH1D* HSCP_nstrip = new TH1D("HSCP_nstrip", "#(strip)", 40, -0.5,39.5);
   TH1D* HSCP_nmpix = new TH1D("HSCP_nmpix", "#m(pix)", 15, -0.5,14.5);
   TH1D* HSCP_nmstrip = new TH1D("HSCP_nmstrip", "#m(strip)", 40, -0.5,39.5);
   TH1D* HSCP_nratio = new TH1D("HSCP_nratio", "#(ratio)", 100, 0.,1.);
   TH1D* HSCP_nmratio = new TH1D("HSCP_nmratio", "#m(ratio)", 100, 0.,1.);

   TH2D* HSCP_dEdXpixVsstrip = new TH2D("HSCP_dEdXpixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* HSCP_dEdX0pixVsstrip = new TH2D("HSCP_dEdX0pixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* HSCP_dEdXstripVsall  = new TH2D("HSCP_dEdXstripVsall","Strip dE/dX: Ih(all)", 100, 0.,10., 100, 0.,10.);
   TH2D* HSCP_dEdXpixVsall   = new TH2D("HSCP_dEdXpixVsall","Pix dE/dX: Ih(all)", 100, 0.,10., 100, 0.,10.);

   TH2D* dEdXVsRun = new TH2D("dEdXVsRun", "dEdX:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXpixVsRun = new TH2D("dEdXpixVsRun", "dEdX(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXstripVsRun = new TH2D("dEdXstripVsRun", "dEdX(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXNoL1VsRun = new TH2D("dEdXNoL1VsRun", "dEdX(noL1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXNoL1pixVsRun = new TH2D("dEdXNoL1pixVsRun", "dEdX(Pix only with no L1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXHiDropVsRun = new TH2D("dEdXHiDropVsRun","dEdX HighDrop:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXpixHiDropVsRun = new TH2D("dEdXpixHiDropVsRun","dEdX HighDrop(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXstripHiDropVsRun = new TH2D("dEdXstripHiDropVsRun","dEdX HighDrop(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXHiDropNoL1VsRun = new TH2D("dEdXHiDropNoL1VsRun","dEdX HighDrop:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0VsRun = new TH2D("dEdX0VsRun", "dEdX0:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0pixVsRun = new TH2D("dEdX0pixVsRun", "dEdX0(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0stripVsRun = new TH2D("dEdX0stripVsRun", "dEdX0(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* MassStripVsRun = new TH2D("MassStripVsRun", "Mass:Run", 545, 271000,325500, 80, 0.,4000.);
   TH2D* MassNoL1VsRun = new TH2D("MassNoL1VsRun", "Mass:Run", 545, 271000,325500, 80, 0.,4000.);
   TH2D* dEdX0NoL1pixVsRun = new TH2D("dEdX0NoL1pixVsRun", "dEdX(Pix only with no L1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0NoL1VsRun = new TH2D("dEdX0NoL1VsRun", "dEdX(noL1Pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4VsRun = new TH2D("dEdXV4sRun", "dEdX:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4pixVsRun = new TH2D("dEdX4pixVsRun", "dEdX(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4stripVsRun = new TH2D("dEdX4stripVsRun", "dEdX(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40VsRun = new TH2D("dEdX40VsRun", "dEdX0:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40pixVsRun = new TH2D("dEdX40pixVsRun", "dEdX0(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40stripVsRun = new TH2D("dEdX40stripVsRun", "dEdX0(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* bg_dEdX0NoL1VsRun = new TH2D("bg_dEdX0NoL1VsRun", "dEdX(noL1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* iasNoL1VsRun = new TH2D("iasNoL1VsRun", "Ias(noL1):Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* iasAllVsRun = new TH2D("iasAllVsRun", "Ias(all):Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probQVsRun = new TH2D("probQVsRun", "ProbQ:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probQNoL1VsRun = new TH2D("probQNoL1VsRun", "ProbQ:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probXYVsRun = new TH2D("probXYVsRun", "ProbXY:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probXYNoL1VsRun = new TH2D("probXYNoL1VsRun", "ProbXY:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probQVsIas = new TH2D("probQVsIas", "ProbQ:Ias", 80, 0.,1., 80, 0.,1.);

   TH2D* FMIP4VsRun = new TH2D("FMIP4VsRun", "FMIP(4):Run", 545, 271000,325500, 50, 0.,1.);
   TH2D* FMIP3p5VsRun = new TH2D("FMIP3p5VsRun", "FMIP(3.5):Run", 545, 271000,325500, 50, 0.,1.);
   TH2D* FMIP3p2VsRun = new TH2D("FMIP3p2VsRun", "FMIP(3.2):Run", 545, 271000,325500, 50, 0.,1.);
   TH2D* FMIP4VsEta = new TH2D("FMIP4VsEta", "FMIP(4):P", 12, -3.,3., 50, 0.,1.);

   TH2D* NmeasVsRun = new TH2D("NmeasVsRun", "#meas:Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasPixVsRun = new TH2D("NmeasPixVsRun", "#meas(pix):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasStrVsRun = new TH2D("NmeasStrVsRun", "#meas(strip):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* Nmeas0VsRun = new TH2D("Nmeas0VsRun", "#meas:Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasPix0VsRun = new TH2D("NmeasPix0VsRun", "#meas(pix):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasStr0VsRun = new TH2D("NmeasStr0VsRun", "#meas(strip):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NsatVsRun = new TH2D("NsatVsRun", "#sat:Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatPixVsRun = new TH2D("NsatPixVsRun", "#sat(pix):Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatStrVsRun = new TH2D("NsatStrVsRun", "#sat(strip):Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* Nsat0VsRun = new TH2D("Nsat0VsRun", "#sat:Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatPix0VsRun = new TH2D("NsatPix0VsRun", "#sat(pix):Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatStr0VsRun = new TH2D("NsatStr0VsRun", "#sat(strip):Run", 545, 271000,325500, 10, -0.5,9.5);
/*
   TH2D* dEdXVsIL = new TH2D("dEdXVsIL", "dEdX:InstLumi", 28, 0,14000, 60, 0.,15.);
   TH2D* dEdXpixVsIL = new TH2D("dEdXpixVsIL", "dEdX(strip):InstLumi", 28, 0,14000, 60, 0.,15.);
   TH2D* dEdXstripVsIL = new TH2D("dEdXstripVsIL", "dEdX(strip):InstLumi", 28, 0,14000, 60, 0.,15.);
*/

   TH2D* ptVsRun = new TH2D("ptVsRun", "pT:Run", 545, 271000,325500, 50, 55.,1550);
   TH2D* nPVVsRun = new TH2D("nPVVsRun", "nPV:Run", 545, 271000,325500, 40,0,80);
   TH2D* invBVsRun = new TH2D("invBVsRun", "invBeta:Run", 545, 271000,325500, 90,-1,2);
   TH2D* errinvBVsRun = new TH2D("errinvBVsRun", "err_invBeta:Run", 545, 271000,325500, 50,0,0.5);
   TH2D* invBDTVsRun = new TH2D("invBDTVsRun", "invBeta(DT):Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBCSCVsRun = new TH2D("invBCSCVsRun", "invBeta(CSC):Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBnewVsRun = new TH2D("invnewBVsRun", "invBeta:Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBnewDTVsRun = new TH2D("invnewBDTVsRun", "invBeta(DT):Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBnewCSCVsRun = new TH2D("invnewBCSCVsRun", "invBeta(CSC):Run", 545, 271000,325500, 90,-1,2);
   TH2D* timeVsRun = new TH2D("timeVsRun", "VertexTimung:Run", 545, 271000,325500, 100,-100,100);
   TH2D* lumiVsRun = new TH2D("lumiVsRun", "Lumi:Run", 545, 271000,325500, 56, 0,14000);


   TH2D* R1_StdEdXVsEvent = new TH2D("R1_StdEdXVsEvent", "dEdX(strip):Event",500,0.,2000000000,60, 0.,15.); 
   TH2D* R1_StdEdXVsLumi = new TH2D("R1_StdEdXVsLumi", "dEdX(strip):Lumi",80,0.,16000,60, 0.,15.); 
   TH2D* R1_LumiVsEvent = new TH2D("R1_LumiVsEvent", "Lumi:Event",500,0.,2000000000,80,0.,16000); 
   TH2D* R1_nPVVsEvent = new TH2D("R1_nPVVsEvent", "nPV:Event",500,0.,2000000000,40, 0.,80.); 
   TH1D* R1_CandVsEvent = new TH1D("R1_CandVsEvent", "Cand per 4M events vs Event", 500,0.,2000000000);

   TH2D* R2_StdEdXVsEvent = new TH2D("R2_StdEdXVsEvent", "dEdX(strip):Event",500,0.,2000000000,60, 0.,15.); 
   TH2D* R2_StdEdXVsLumi = new TH2D("R2_StdEdXVsLumi", "dEdX(strip):Lumi",80,0.,16000,60, 0.,15.); 
   TH2D* R2_LumiVsEvent = new TH2D("R2_LumiVsEvent", "Lumi:Event",500,0.,2000000000,80,0.,16000); 
   TH2D* R2_nPVVsEvent = new TH2D("R2_nPVVsEvent", "nPV:Event",500,0.,2000000000,40, 0.,80.); 
   TH1D* R2_CandVsEvent = new TH1D("R2_CandVsEvent", "Cand per 4M events vs Event", 500,0.,2000000000);

   TH2D* R3_StdEdXVsEvent = new TH2D("R3_StdEdXVsEvent", "dEdX(strip):Event",500,0.,2000000000,60, 0.,15.); 
   TH2D* R3_StdEdXVsLumi = new TH2D("R3_StdEdXVsLumi", "dEdX(strip):Lumi",80,0.,16000,60, 0.,15.); 
   TH2D* R3_LumiVsEvent = new TH2D("R3_LumiVsEvent", "Lumi:Event",500,0.,2000000000,80,0.,16000); 
   TH2D* R3_nPVVsEvent = new TH2D("R3_nPVVsEvent", "nPV:Event",500,0.,2000000000,40, 0.,80.); 
   TH1D* R3_CandVsEvent = new TH1D("R3_CandVsEvent", "Cand per 4M events vs Event", 500,0.,2000000000);


   TH1D* HHitPix = new TH1D( "HHitPix", "HHitPix", 200, 0, 20);
   TProfile* HHitProfilePix = new TProfile( "HHitProfilePix", "HHitProfilePix", 50, 0, 100);
   TH2D* HHit2DPix = new TH2D( "HHit2DPix", "HHit2DPix", 50, 0, 100,200, 0, 20);
   TH2D* HHit2DPix_NoM = new TH2D( "HHit2DPix_NoM", "HHit2DPix", 50, 0, 100,200, 0, 20);
   TH1D* HHitStrip = new TH1D( "HHitStrip", "HHitStrip", 200, 0, 20);
   TProfile* HHitProfileStrip = new TProfile( "HHitProfileStrip", "HHitProfileStrip", 50, 0, 100);
   TH2D* HHit2DStrip = new TH2D( "HHit2DStrip", "HHit2DStrip", 50, 0, 100,200, 0, 20);
   TH2D* HHit2DStrip_NoM = new TH2D( "HHit2DStrip_NoM", "HHit2DStrip", 50, 0, 100,200, 0, 20);

   TH2D* IhStripVsP = new TH2D("IhStrip_VsP_EtaCentral", "IhStripVsP_EtaCentral",100,0,10.0,50,0,100 );
   TH2D* IhStripVsPtight = new TH2D("IhStrip_VsP_EtaTight", "IhStripVsP_EtaTight",100,0,10.0,50,0,100 );

   TH2D* IhBestStripVsP = new TH2D("Ih0_noL1_VsP_EtaCentral", "Ih0_noL1StripVsP_EtaCentral",100,0,10.0,50,0,100 );
   TH2D* IhBestStripVsPtight = new TH2D("Ih0_noL1_VsP_EtaTight", "Ih0_noL1StripVsP_EtaTight",100,0,10.0,50,0,100 );


   TH2D* IhStripVsP_presk = new TH2D("IhStrip_VsP_EtaCentral_presk", "IhStripVsP_EtaCentral_presk",100,0,10.0,50,0,100 );
   TH2D* IhStripVsPtight_presk = new TH2D("IhStrip_VsP_EtaTight_presk", "IhStripVsP_EtaTight_presk",100,0,10.0,50,0,100 );

   TH2D* IhBestStripVsP_presk = new TH2D("Ih0_noL1_VsP_EtaCentral_presk", "Ih0_noL1StripVsP_EtaCentral_presk",100,0,10.0,50,0,100 );
   TH2D* IhBestStripVsPtight_presk = new TH2D("Ih0_noL1_VsP_EtaTight_presk", "Ih0_noL1StripVsP_EtaTight_presk",100,0,10.0,50,0,100 );

   TH2D* IhBestvsIas_p_5_100 = new TH2D("Ih0_noL1_vs_Ias_strip_only_p_5_100", "Ih0_noL1_vs_Ias_strip_p_5_100",100,0,10.0,50,0,1 );
   TH2D* IhBestvsIas_p_10_45 = new TH2D("Ih0_noL1_vs_Ias_strip_only_p_10_45", "Ih0_noL1_vs_Ias_strip_p_10_45",100,0,10.0,50,0,1 );


   TH2D* IasStripVsIh0noL1_p_10_45 = new TH2D("Ias_strip_VS_Ih0_noL1_p_10_45", "Ias_strip_VS_Ih0_noL1_p_10_45",50,0,1.0,100,0,10.0 );

   TH2D* IhBestvsP_10_45_central = new TH2D("Ih0_noL1_vs_P_10_45_central", "Ih0_noL1_vs_P_10_45_central",100,0,10.0,200,0,1000 );
   TH2D* IhBestvsP_10_45_tight = new TH2D("Ih0_noL1_vs_P_10_45_tight", "Ih0_noL1_vs_P_10_45_tight",100,0,10.0,200,0,1000 );

   TH1D* P_after_presel = new TH1D("P_after_presel","P_after_presel",50,0,100);
   TH1D* NPV_all = new TH1D("NPV_all","NPV all events",100,0,100);
   TH1D* NPV_presel = new TH1D("NPV_presel","NPV after presel",50,0,100);

   TH1D* MASS_SINGLE_ihcut = new TH1D("HSCP_Mass_Single_ihcut", "Mass via Ih0_noL1, cut IAS (single) > 0.115, cut IH 3.47", 200, 0.,200.);
   TH1D* MASS_TRIPLE_ihcut = new TH1D("HSCP_Mass_Triple_ihcut", "Mass via Ih0_noL1, cut IAS (triple) > 0.115, cut IH 3.47", 200, 0.,200.);


   TH1D* MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias = new TH1D("HSCP_Mass_Single_ihcut_ptsupp50_qtl_80_90_ias", "Mass via Ih0_noL1,IAS IN 80-90 QTL, cut IH 3.47, PT > 50", 200, 0.,1000.);

   TH1D* MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias = new TH1D("HSCP_Mass_Triple_ihcut_ptsupp50_qtl_80_90_ias", "Mass via Ih0_noL1, IAS IN 80-90 QTL, cut IH 3.47, PT > 50", 200, 0.,1000.);


   TH1D* IAS_TRIPLE_PT_50_60_selection = new TH1D("IAS_TRIPLE_PT_50_60_selection", "IAS triple, cut IH 3.47, PT [50-60]", 50, 0.,1.0);
   TH1D* IAS_SINGLE_PT_50_60_selection = new TH1D("IAS_SINGLE_PT_50_60_selection", "IAS single, cut IH 3.47, PT [50-60]", 50, 0.,1.0);

   

   TH1D* MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias = new TH1D("HSCP_Mass_Single_ihcut_ptsupp100_qtl_80_90_ias", "Mass via Ih0_noL1,IAS IN 80-90 QTL, cut IH 3.47, PT > 100", 200, 0.,1000.);
   TH1D* MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias = new TH1D("HSCP_Mass_Triple_ihcut_ptsupp100_qtl_80_90_ias", "Mass via Ih0_noL1, IAS IN 80-90 QTL, cut IH 3.47, PT > 100", 200, 0.,1000.);

   TH1D* MASS_SINGLE_ihcut_base_notriple_ias_0p1_ptsupp50_qtl_80_90_ias = new TH1D("HSCP_Mass_Single_ihcut_base_above_triple_below_ias_0p115_ptsupp50_qtl_80_90_ias", "Mass via Ih0_noL1, 0.07 < IAS < 0.1 on base but triple has IAS < 0.08, cut IH 3.47, PT > 50", 200, 0.,1000.);

   TH1D* MASS_TRIPLE_ihcut_triple_nobase_ias_0p1_ptsupp50_qtl_80_90_ias = new TH1D("HSCP_Mass_Triple_ihcut_triple_above_base_below_ias_0p115_ptsupp50_qtl_80_90_ias", "Mass via Ih0_noL1, 0.07 < IAS < 0.1 on triple but base has IAS < 0.08, PT > 50, cut IH 3.47", 200, 0.,1000.);



   TH1D* MASS_SINGLE_ihcut_base_notriple_ias_0p1 = new TH1D("HSCP_Mass_Single_ihcut_base_above_triple_below_ias_0p115", "Mass via Ih0_noL1, cut IAS (single) > 0.1 on base but not on triple, cut IH 3.47, P [10-45]", 200, 0.,200.);
   TH1D* MASS_TRIPLE_ihcut_triple_nobase_ias_0p1 = new TH1D("HSCP_Mass_Single_ihcut_triple_above_base_below_ias_0p115", "Mass via Ih0_noL1, cut IAS (single) > 0.1 on base but not on triple, cut IH 3.47, P [ 10-45]", 200, 0.,200.);


   TH1D* MASS_SINGLE_nocut = new TH1D("HSCP_Mass_Single_nocut", "Mass via Ih0_noL1, cut IAS (single) > 0.115, no cut IH", 200, 0.,200.);

   TH1D* MASS_TRIPLE_nocut = new TH1D("HSCP_Mass_Triple_nocut", "Mass via Ih0_noL1, cut IAS (triple) > 0.115, no cut IH", 200, 0.,200.);



   TH1D* IAS_single_pt50 = new TH1D("IAS_single_pt50", "IAS single pt > 50, IH > 3.47 and track preselection", 100, 0,1);
   TH1D* IAS_triple_pt50 = new TH1D("IAS_triple_pt50", "IAS triple pt > 50", 100, 0,1);

   TH1D* IAS_triple_when_simple_above_0p115 = new TH1D("IAS_triple_when_simple_above_0p115", "IAS triple when IAS simple > 0.115", 200, 0.08,0.4);
   TH1D* IAS_simple_when_triple_above_0p115 = new TH1D("IAS_simple_when_triple_above_0p115", "IAS simple when IAS triple > 0.115", 200, 0.08,0.4);

   TH2D* Ias_simple_Vs_Ih0noL1_p_10_45 = new TH2D("Ias_strip_simple_VS_Ih0_noL1_p_10_45_nocutIH", "Ias_strip_simple_VS_Ih0_noL1_p_10_45_nocutIH",50,0,1.0,100,0,10.0 );
   TH2D* Ias_triple_Vs_Ih0noL1_p_10_45 = new TH2D("Ias_strip_triple_VS_Ih0_noL1_p_10_45_nocutIH", "Ias_strip_triple_VS_Ih0_noL1_p_10_45_nocutIH",50,0,1.0,100,0,10.0 );

   TH2D* NPV_ias_triple_cutih = new TH2D("NPV_vs_ias_triple_cutih", "NPV vs ias triple, IH > 3.47",100,0,100,50,0,1);
   TH2D* NPV_ias_single_cutih = new TH2D("NPV_vs_ias_single_cutih", "NPV vs ias single, IH > 3.47",100,0,100,50,0,1);
 
   TH2D* NPV_ih_cutih = new TH2D("NPV_vs_ih0_noL1_cutih", "NPV vs ih0_noL1, IH > 3.47",100,0,100,100,0,10);
   TH2D* NPV_ih_cut3p29 = new TH2D("NPV_vs_ih0_noL1_cutih_3p29", "NPV vs ih0_noL1, IH > 3.29",100,0,100,100,0,10);
   TH2D* NPV_ih_nocut = new TH2D("NPV_vs_ih0_noL1_nocut", "NPV vs ih0_noL1, no cut IH",100,0,100,100,0,10);
   TH2D* NPV_ih_cut1 = new TH2D("NPV_vs_ih0_noL1_cut_1", "NPV vs ih0_noL1, IH > 1",100,0,100,100,0,10);

   TH2D* NPV_ih_strip_cutih = new TH2D("NPV_vs_ih_strip_cutih", "NPV vs ih strip, IH > 3.47",100,0,100,100,0,10);
   TH2D* NPV_ih_strip_cut3p29 = new TH2D("NPV_vs_ih_strip_cutih_3p29", "NPV vs ih strip, IH > 3.29",100,0,100,100,0,10);
   TH2D* NPV_ih_strip_nocut = new TH2D("NPV_vs_ih_strip_nocut", "NPV vs ih strip, no cut IH",100,0,100,100,0,10);
   TH2D* NPV_ih_strip_cut1 = new TH2D("NPV_vs_ih_strip_cut_1", "NPV vs ih strip, IH > 1",100,0,100,100,0,10);



   TH1D* IH0_noL1 = new TH1D("IH0_noL1_selection", "IH0_noL1 for selection and pt > 50", 100, 0,10);

   TH1D* P_LOWPU_pt50 = new TH1D("P_LOWPU_pt50_cutih", "P for NPV < 10, Ih > 3.47 and pt > 50", 100, 0,1000);
   TH1D* P_HIGHPU_pt50 = new TH1D("P_HIGHPU_pt50_cutih", "P for NPV > 60, Ih > 3.47 and pt > 50", 100, 0,1000);


   TH1D* P_PU_28_29_pt50 = new TH1D("P_PU_28_29_pt50_cutih", "P for NPV in [28-29], Ih > 3.47 and pt > 50", 100, 0,1000);
   TH1D* P_PU_30_31_pt50 = new TH1D("P_PU_30_31_pt50_cutih", "P for NPV in [30-31], Ih > 3.47 and pt > 50", 100, 0,1000);


   TH1D* P_PU_28_30_pt50 = new TH1D("P_PU_28_30_pt50_cutih", "P for NPV in [28-30], Ih > 3.47 and pt > 50", 100, 0,1000);
   TH1D* P_PU_32_34_pt50 = new TH1D("P_PU_32_34_pt50_cutih", "P for NPV in [32-34], Ih > 3.47 and pt > 50", 100, 0,1000);


   TH1D* P_PU_58_60_pt50 = new TH1D("P_PU_58_60_pt50_cutih", "P for NPV in [58-60], Ih > 3.47 and pt > 50", 100, 0,1000);
   TH1D* P_PU_08_10_pt50 = new TH1D("P_PU_08_10_pt50_cutih", "P for NPV in [08-10], Ih > 3.47 and pt > 50", 100, 0,1000);


   TH1D* P_single_pt50_ias_qtl_80_90 = new TH1D("P_single_pt50_ias_qtl_80_90", "P single for IAS single in 80-90% qtl [0.07-1] and pt > 50", 100, 0,1000);
   TH1D* P_triple_pt50_ias_qtl_80_90 = new TH1D("P_triple_pt50_ias_qtl_80_90", "P triple for IAS triple in 80-90% qtl [0.07-1] and pt > 50", 100, 0,1000);

   TH1D* PU_ias_single_0p3 = new TH1D("PU_when_ias_single_supp_0p3", "PU for events with IAS single > 0.3, Ih > 3.47 and pt > 50", 100, 0,100);
   TH1D* PU_ias_triple_0p3 = new TH1D("PU_when_ias_triple_supp_0p3", "PU for events with IAS triple > 0.3, Ih > 3.47 and pt > 50", 100, 0,100);


   TH2D* NPV_mass_cutih = new TH2D("NPV_vs_mass_cutih_ias_single", "NPV vs MASS, Ih0_noL1 > 3.47, IAS single in 80-90% qtl [0.07-1]",100,0,100,200,0,1000);
   TH2D* NPV_mass_cut3p29 = new TH2D("NPV_vs_mass_cutih_3p29_ias_single", "NPV vs MASS, Ih0_noL1 > 3.29, IAS single in 80-90% qtl [0.07-1]",100,0,100,200,0,1000);
   

   const double P_Min               = 1   ;
   const double P_Max               = 16  ; // 1 + 14 + 1; final one is for pixel!
   const int    P_NBins             = 15  ; // 15th bin = pixel; 0 is underflow
   const double Path_Min            = 0.2 ;
   const double Path_Max            = 1.6 ;
   const int    Path_NBins          = 42  ;
   const double Charge_Min          = 0   ;
   const double Charge_Max          = 5000;
   const int Charge_NBins = 500 ;
   TH3D* Charge_Vs_Path = new TH3D( "Charge_Vs_Path", "Charge_Vs_Path", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH1D* CHARGE_OVER_PATHLENGTH_GIVEN_MODULE = new TH1D("Charge_over_Pathlength", "Charge for pathlength 0.5-0.6 for mod[geom] == 1", Charge_NBins, Charge_Min,Charge_Max);

   TH1D* CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP = new TH1D("Charge_over_Pathlength_STRIP", "STRIP Charge for pathlength 0.5-0.6 for mod[geom] == 1", Charge_NBins, Charge_Min,Charge_Max);

   TH1D* CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_LOWPU = new TH1D("Charge_over_Pathlength_STRIP_LOWPU", "STRIP Charge for pathlength 0.5-0.6 for mod[geom] == 1, NPV < 20", Charge_NBins, Charge_Min,Charge_Max);
 
   TH1D* CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_HIGHPU = new TH1D("Charge_over_Pathlength_STRIP_HIGHPU", "STRIP Charge for pathlength 0.5-0.6 for mod[geom] == 1, NPV > 50", Charge_NBins, Charge_Min,Charge_Max);

   TH1D* CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_NOSATURATED = new TH1D("Charge_over_Pathlength_STRIP_nosaturated_only", "STRIP Charge for pathlength 0.5-0.6 for mod[geom] == 1, clusters not saturated", Charge_NBins, Charge_Min,Charge_Max);
   TH1D* CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_PIXEL = new TH1D("Charge_over_Pathlength_PIXEL", "PXL Charge for pathlength 0.5-0.6 for mod[geom] == 1", Charge_NBins, Charge_Min,Charge_Max);
   TH1D* PATHLENGTH_GIVEN_MODULE = new TH1D("Pathlength", "Pathlength 0.5-0.6 for mod[geom] == 1", Path_NBins, Path_Min,Path_Max);

   // PU DEPENDENCIES HISTOGRAMS 
   TH2D* TptGivenPathlenght1 = new TH2D("TptPathlenght1","Template for given pathlenght = 1 ", P_NBins, P_Min ,P_Max, Charge_NBins, Charge_Min , Charge_Max);


   TH3D* Charge_Vs_Path_PU_corr = new TH3D( "Charge_Vs_Path_PU_corr", "Charge_Vs_Path_PU_corr", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   
   TH3D* Charge_Vs_Path_PU_below_20 = new TH3D( "Charge_Vs_Path_Low_PU", "Charge_Vs_Path_low_pu", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);

   TH3D* Charge_Vs_Path_PU_between = new TH3D( "Charge_Vs_Path_Middle_PU", "Charge_Vs_Path_medium_pu", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);

   TH3D* Charge_Vs_Path_PU_above_40 = new TH3D( "Charge_Vs_Path_High_PU", "Charge_Vs_Path_high_pu", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);


   int what_bins[ias_intervals+1]={ 0 }; 
   int nb_ias_0p1_triple = 0,nb_ias_0p1_single = 0;
   int nb_ias_0p1_triple_nocut = 0,nb_ias_0p1_single_nocut = 0;
   int nb_ias_0p1_triple_cut3p29 = 0, nb_ias_0p1_single_cut3p29 = 0;

   int nb_ias_0p2_triple = 0,nb_ias_0p2_single = 0;
   int nb_ias_0p2_triple_nocut = 0,nb_ias_0p2_single_nocut = 0;
   int nb_ias_0p2_triple_cut3p29 = 0, nb_ias_0p2_single_cut3p29 = 0;


   int nb_ias_0p3_triple = 0,nb_ias_0p3_single = 0;
   int nb_ias_0p3_triple_nocut = 0,nb_ias_0p3_single_nocut = 0;
   int nb_ias_0p3_triple_cut3p29 = 0, nb_ias_0p3_single_cut3p29 = 0;


   int nb_ias_0p4_triple = 0,nb_ias_0p4_single = 0;
   int nb_ias_0p4_triple_nocut = 0,nb_ias_0p4_single_nocut = 0;
   int nb_ias_0p4_triple_cut3p29 = 0, nb_ias_0p4_single_cut3p29 = 0;

   TRandom3* RNG = new TRandom3();

   if ( get_list_of_bins){
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
          Long64_t pentry = LoadTree(jentry);
          if (pentry < 0) break;
          nbpu = b_npv->GetEntry(pentry); nbytes_pu += nbpu;
          //cout << "npv for this event:" << npv << endl;
          NPV_all->Fill(npv);
       }
      
       int all_entries = NPV_all->GetEntries();
       int nb_off = (all_entries%ias_intervals);
    
       int per_bin = (all_entries / ias_intervals);
    
       int bins_pu[ias_intervals];    
       int real_events[ias_intervals] = { 0 };
       std::fill_n(bins_pu, ias_intervals, per_bin);
       for (int i = 0; i < nb_off ; i++){
           bins_pu[i]+=1;
       }
       cout << "Finnished splitting into equal stats, there will be " << ias_intervals << " Pile-Up bins that should have : " <<endl;
    
       for(int i = 0; i< ias_intervals ; i++){
           cout << "Bin " << i << " : " << bins_pu[i] << "events " << endl;
       }
       int nbinx = NPV_all->GetXaxis()->GetNbins();
       cout << " this hist has " << nbinx << " bins" << endl;
       what_bins[ias_intervals] = nbinx-1;
       int mean_bin = NPV_all->GetXaxis()->FindBin(NPV_all->GetMean()); // 30 around the mean of NPV
       int whole_integal = NPV_all->Integral();
       int all_test = 0;
       for (int j = 0 ; j < ias_intervals-1; j++ ){
           bool bin_filled = false;
           for (int i= what_bins[j]+1; i< nbinx ; i++){
               int itgtr = NPV_all->Integral(what_bins[j]+1,i);
               //cout << "integrating between bins" << what_bins[j]+1 << " and " << i << " = " << itgtr << endl;
               //cout << "is result above " << bins_pu[j] << endl;
               if ( (abs( itgtr - bins_pu[j])) < (per_bin*0.2)){
                   //cout << "YES so superior bin is : " << i << " and nb events in this PU bin = " << NPV_all->Integral(what_bins[j],i)<<endl; 
                   what_bins[j+1] = i;
                   real_events[j] = NPV_all->Integral(what_bins[j]+1,i);
                   bin_filled = true;
                   all_test+=real_events[j];
               }
               else if(itgtr > bins_pu[j]){
                   cout << "warning nb of events in this bin > than 120 percent of expected amound" << endl;
                   what_bins[j+1] = i;
                   real_events[j] = NPV_all->Integral(what_bins[j]+1,i);
                   bin_filled = true;
               }
               if (bin_filled) break;
           }
       }
    
       real_events[ias_intervals-1] = NPV_all->Integral(what_bins[ias_intervals-1]+1,what_bins[ias_intervals]);
       all_test+=real_events[ias_intervals-1];
       cout <<"Last bin integral : " << NPV_all->Integral(nbinx,nbinx) << " so the sum of both is : " << NPV_all->Integral(nbinx,nbinx) + all_test <<endl;
       cout << "Total number of entries : " << nentries <<endl; 
       cout << "********** since we have a finite number of bins, we can not provide an exact splitting, the best one can do is : " << endl;
       for (int i = 0; i < ias_intervals ; i++){
           double valbinmin = NPV_all->GetXaxis()->GetBinCenter(what_bins[i]);
           double valbin = NPV_all->GetXaxis()->GetBinCenter(what_bins[i+1]);
           cout << "Bin " << i+1 << " NPV :  " << what_bins[i] << " - " << what_bins[i+1] << " with " << real_events[i] << " events" << endl;
       }
       cout << "Sum of all integrated events : " << all_test << " difference (in last bin probably)" << endl; 
   }
   else{
       std::ifstream fin("list_5_bins_pu.txt");
       
       std::vector<int> results;
       int result, pe = 0;
       if(fin.is_open()){
           while (fin >> result){
               if (pe > ias_intervals) cout << "warning, too many numbers in txt file" <<endl;
               what_bins[pe] = int(result);
               pe++;
           }
       
           std::cout << "We read values from txt files, and the borns are : " << endl;
           for (int i = 0; i< ias_intervals ; i++){
               cout << "BIN " << i+1 << " - >   " << what_bins[i] << "-" << what_bins[i+1] <<endl;
           }
       }
   }

   vector <TH1D*> PU_distrib(ias_intervals);
 
   vector < vector <TH1D*> > PE_Ias, PE_Ias_cutih;
   PE_Ias.resize(ias_intervals, vector<TH1D*>(nPE));
   PE_Ias_cutih.resize(ias_intervals,vector<TH1D*>(nPE));

   dEdxTemplatesPU.resize(ias_intervals, NULL);
   dEdxTemplatesPUBis.resize(ias_intervals, NULL);
 
   PE_dEdxTemplatesPU.resize(ias_intervals, vector<TH3F*>(nPE));

   //to change hardcoded ias_intervals * 4 for TGraphError 

   double means_ias_pe[ias_intervals][nPE] = { 0 };
   double err_means_ias_pe[ias_intervals][nPE] = { 0 };
   double means_ias_pe_cutih[ias_intervals][nPE] = { 0 };
   double err_means_ias_pe_cutih[ias_intervals][nPE] = { 0 };

   double std_dev_means_ias_pe[ias_intervals][nPE] = { 0 };
   double std_dev_means_ias_pe_cutih[ias_intervals][nPE] = { 0 };
   
   



   for (int j = 0; j < ias_intervals; j++){ 
      for (int i = 0 ; i < nPE ; i++){
          string name_ias_pe = "IAS_strip_PE_" + to_string(i+1) + "_PU_" + to_string(what_bins[j]) + "_" + to_string(what_bins[j+1]);
          PE_Ias[j][i] = new TH1D(name_ias_pe.c_str(),name_ias_pe.c_str(),50,0,1);
          string name_ias_pe_cut = "IAS_strip_PE_" + to_string(i+1) +"_PU_" + to_string(what_bins[j]) + "_" + to_string(what_bins[j+1]) + "_cut_Ih";
          PE_Ias_cutih[j][i]  = new TH1D(name_ias_pe_cut.c_str(),name_ias_pe_cut.c_str(),50,0,1);
       }
   }

   vector <TH3D*> Charge_Vs_Path_PU;
   Charge_Vs_Path_PU.resize(ias_intervals);
   vector <string> names_templates;names_templates.resize(ias_intervals);
   cout << "after declaring all templates and for PE" << endl;

   double ias_top_born[ias_intervals+1];
   for (int k = 0 ; k <= ias_intervals ; k++){
       ias_top_born[k] = what_bins[k];
   }

   for (int i = 0; i < ias_intervals; i++){
       string pren = "Charge_Vs_Path_PU_between_", pudistrib = "PU_distrib_between_";
       string uds = "_";
       double lowerlim = what_bins[i];
       double upperlim = what_bins[i+1];
       int precisionVal = 1; 
       string bg = std::to_string(lowerlim).substr(0, std::to_string(lowerlim).find(".") + precisionVal + 1);
       string end = std::to_string(upperlim).substr(0, std::to_string(upperlim).find(".") + precisionVal + 1);
       int trunc_bg = stoi(bg);
       int trunc_end = stoi(end);
       string name = pren + to_string(trunc_bg) + uds + to_string(trunc_end);
       string pu_name = pudistrib + to_string(trunc_bg) + uds + to_string(trunc_end);
       names_templates[i] = name;
       Charge_Vs_Path_PU[i] = new TH3D(name.c_str(),name.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
       PU_distrib[i] = new TH1D(pu_name.c_str(),pu_name.c_str(),100,0,100);
   }

   //STATISTICAL TESTS

   int count_tenth = 0;
   int choice_third = 3; //1 and 2 for first two thirds to produce templates, 3 for last reading third
   bool tenth_stat_prod = false;
   bool tenth_stat_study = true;

   TH1D* Ias_when_PU_in_16_18_base  = new TH1D("Ias_when_PU_between_16_18_base","Ias_when_PU_small_interval_base",50,0,1);
   TH1D* Ias_below_5GeV_PU_below_20_base  = new TH1D("Ias_for_p_below_5GeV_when_low_PU_base","Ias_for_p_below_5GeV_when_low_PU_base",50,0,1);


   TH1D* Ias_when_PU_between_base  = new TH1D("Ias_when_Mid_PU_base","Ias_when_Mid_PU_base",50,0,1);
   TH1D* Ias_when_PU_in_30_32_base  = new TH1D("Ias_when_PU_between_30_32_base","Ias_when_PU_medium_interval_base",50,0,1);
   TH1D* Ias_below_5GeV_PU_between_base  = new TH1D("Ias_for_p_below_5GeV_when_Mid_PU_base","Ias_for_p_below_5GeV_when_Mid_PU_base",50,0,1);


   TH1D* Ias_when_PU_above_40_base  = new TH1D("Ias_when_High_PU_base","Ias_when_High_PU_base",50,0,1);
   TH1D* Ias_below_5GeV_PU_above_40_base  = new TH1D("Ias_for_p_below_5GeV_when_High_PU_base","Ias_for_p_below_5GeV_when_High_PU_base",50,0,1);

   
   std::string name_delta_ias = "Ias(stat2) vs [Ias(stat2) - Ias(stat1)]" 
   if(tenth_
   TH2D* Ias_1_VS_Delta_Ias_1_2  = new TH2D("IAS_1_VS_Delta_Ias"," Ias(stat1) vs Ias(stat1) - Ias(stat2)",50,0,1,500,-0.5,0.5);

   TH1D* Ias_toy  = new TH1D("Ias_toy_modgeom_1","Ias_toy_modgeom1",50,0,1);

   TEfficiency* eff_ih_PU = new TEfficiency("Eff_ih_PU","Eff of IH cut;PU;#epsilon",80,0,80);



   TH1D* Is_all_base  = new TH1D("Is_strip_all_inclusive","Is_strip_all_inclusive",50,0,1);

   vector <TH1D*> Ias_when_PU,Ias_when_PU_triple,Ias_when_PU_ih_cut,Ias_when_PU_ih_cut_triple;

   Ias_when_PU.resize(ias_intervals);
   Ias_when_PU_triple.resize(ias_intervals);
   Ias_when_PU_ih_cut.resize(ias_intervals);
   Ias_when_PU_ih_cut_triple.resize(ias_intervals);

   vector <TH1D*> Is_when_PU,Is_when_PU_triple,Is_when_PU_ih_cut,Is_when_PU_ih_cut_triple;
   Is_when_PU.resize(ias_intervals);
   Is_when_PU_triple.resize(ias_intervals);
   Is_when_PU_ih_cut.resize(ias_intervals);
   Is_when_PU_ih_cut_triple.resize(ias_intervals);
   


   double sums_ias[ias_intervals] = { 0 };
   double means_ias[ias_intervals] = { 0 };


   double sums_ias_nocut[ias_intervals] = { 0 };
   double means_ias_nocut[ias_intervals] = { 0 };
 



   double sums_ias_triple[ias_intervals] = { 0 };
   double means_ias_triple[ias_intervals] = { 0 };

   double sums_ias_triple_nocut[ias_intervals] = { 0 };
   double means_ias_triple_nocut[ias_intervals] = { 0 };

   double err_mean_x_base[ias_intervals] = { 0 } ;  
   double err_mean_x_triple[ias_intervals] = { 0 } ;

   double error_mean_base[ias_intervals] = { 0 } ;
   double error_mean_triple[ias_intervals] = { 0 } ;

   double error_mean_base_cutih[ias_intervals] = { 0 } ;
   double error_mean_triple_cutih[ias_intervals] = { 0 } ;

   int ih_eff_num[ias_intervals] = { 0 };
   int ih_eff_denom[ias_intervals] = { 0 };   

   double ratio_above_qtl[ias_intervals] = {0};
   double ratio_above_qtl_triple[ias_intervals] = {0};
   double ratio_above_qtl_cut[ias_intervals] = {0};
   double ratio_above_qtl_cut_triple[ias_intervals] = {0};





   int nb_above_99_npv = 0;
   double quantiles_PU_80[ias_intervals] = { 0 };
   double quantiles_PU_triple_80[ias_intervals] = { 0 };

   double quantiles_PU_90[ias_intervals] = { 0 };
   double quantiles_PU_triple_90[ias_intervals] = { 0 };

   double quantiles_PU_cut_80[ias_intervals] = { 0 };
   double quantiles_PU_triple_cut_80[ias_intervals] = { 0 };

   double quantiles_PU_cut_90[ias_intervals] = { 0 };
   double quantiles_PU_triple_cut_90[ias_intervals] = { 0 };


   double quantiles_PU_99[ias_intervals] = { 0 };
   double quantiles_PU_triple_99[ias_intervals] = { 0 };

   double quantiles_PU_cut_99[ias_intervals] = { 0 };
   double quantiles_PU_triple_cut_99[ias_intervals] = { 0 };

   int precisionVal = 1;
   string ias_axis_names[ias_intervals];
   double ias_diff_2by2[ias_intervals-1] = { 0 };

   string names_pu[ias_intervals] = {}; 

   for (int i = 0; i < ias_intervals; i++){
       string is = "Is_when_PU_between_";
       string pren = "Ias_when_PU_between_";
       string uds = "_", base = "_base", triple = "_triple", cutIH = "_cutIH";
       string puname = "PU_";
       double lowerlim = what_bins[i];
       double upperlim = what_bins[i+1];
        
       string bg = std::to_string(lowerlim).substr(0, std::to_string(lowerlim).find(".") + precisionVal + 1);
       string end = std::to_string(upperlim).substr(0, std::to_string(upperlim).find(".") + precisionVal + 1);
       int trunc_bg = stoi(bg);
       int trunc_end = stoi(end);

       string fnlnamepu = puname + to_string(trunc_bg) + uds + to_string(trunc_end);
       names_pu[i] = fnlnamepu;

      
       string name_is = is + to_string(trunc_bg) + uds + to_string(trunc_end) + base;
       string name_is_triple = is + to_string(trunc_bg) + uds + to_string(trunc_end) + triple;
       string name_is_cutIH = is + to_string(trunc_bg) + uds + to_string(trunc_end) + base + cutIH;
       string name_is_cutIH_triple = is + to_string(trunc_bg) + uds + to_string(trunc_end) + triple + cutIH;

       string name = pren + to_string(trunc_bg) + uds + to_string(trunc_end) + base;
       string nametriple = pren + to_string(trunc_bg) + uds + to_string(trunc_end) + triple;
       string namecutIH = pren + to_string(trunc_bg) + uds + to_string(trunc_end) + base + cutIH;
       string namecutIHtriple  = pren + to_string(trunc_bg) + uds + to_string(trunc_end) + triple + cutIH;

       Ias_when_PU[i] = new TH1D(name.c_str(),name.c_str(),50,0,1);
       Ias_when_PU_ih_cut[i] = new TH1D(namecutIH.c_str(),namecutIH.c_str(),50,0,1);
       Ias_when_PU_ih_cut_triple[i] = new TH1D(namecutIHtriple.c_str(),namecutIHtriple.c_str(),50,0,1);
       Ias_when_PU_triple[i] = new TH1D(nametriple.c_str(),nametriple.c_str(),50,0,1);

       Is_when_PU[i] = new TH1D(name_is.c_str(),name_is.c_str(),50,0,1);
       Is_when_PU_ih_cut[i] = new TH1D(name_is_cutIH.c_str(),name_is_cutIH.c_str(),50,0,1);
       Is_when_PU_ih_cut_triple[i] = new TH1D(name_is_cutIH_triple.c_str(),name_is_cutIH_triple.c_str(),50,0,1);
       Is_when_PU_triple[i] = new TH1D(name_is_triple.c_str(),name_is_triple.c_str(),50,0,1);
       
   }




   TH1D* Ias_when_PU_in_16_18_triple  = new TH1D("Ias_when_PU_between_16_18_triple","Ias_when_PU_small_interval_triple",50,0,1);
   TH1D* Ias_below_5GeV_PU_below_20_triple  = new TH1D("Ias_for_p_below_5GeV_when_low_PU_triple","Ias_for_p_below_5GeV_when_low_PU_triple",50,0,1);


   TH1D* Ias_when_PU_in_30_32_triple  = new TH1D("Ias_when_PU_between_30_32_triple","Ias_when_PU_medium_interval_triple",50,0,1);
   TH1D* Ias_below_5GeV_PU_between_triple  = new TH1D("Ias_for_p_below_5GeV_when_Mid_PU_triple","Ias_for_p_below_5GeV_when_Mid_PU_triple",50,0,1);

   TH2D* Ias_vs_PU_above_40_triple  = new TH2D("Ias_vs_High_PU_triple","Ias_vs_High_PU_triple",50,0,1,50,0,100);

   TH1D* Ias_below_5GeV_PU_above_40_triple  = new TH1D("Ias_for_p_below_5GeV_when_High_PU_triple","Ias_for_p_below_5GeV_when_High_PU_triple",50,0,1);

   TH1D* Ias_all_triple_cutIH  = new TH1D("Ias_all_triple_cutIH","Ias_all_triple_cutIH",50,0,1);
   TH1D* Ias_all_base_cutIH  = new TH1D("Ias_all_base_cutIH","Ias_all_base_cutIH",50,0,1);
   
   TH1D* Ias_all_triple_nocut  = new TH1D("Ias_all_triple_nocut","Ias_all_triple_nocut",50,0,1);
   TH1D* Ias_all_base_nocut  = new TH1D("Ias_all_base_nocut","Ias_all_base_nocut",50,0,1);

   TH1D* PU_distrib_cut_ias_01  = new TH1D("Ias_all_base_cut01","Ias_all_base_cut01",50,0,1);
   

   Ias_vs_PU_above_40_triple->GetXaxis()->SetTitle("Ias");
   Ias_vs_PU_above_40_triple->GetYaxis()->SetTitle("NPV");

   // END PU DEPENDENCIES HISTOGRAMS

   TH3D* Charge_Vs_Path_noL1 = new TH3D( "Charge_Vs_Path_noL1", "Charge_Vs_Path_noL1", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH3D* Charge_Vs_Path_NoM = new TH3D( "Charge_Vs_Path_NoM", "Charge_Vs_Path", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH3D* Charge_Vs_Path_noL1_NoM = new TH3D( "Charge_Vs_Path_noL1_NoM", "Charge_Vs_Path_noL1", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);


   TH1D* PT_compute_ias = new TH1D("PT_for_tracks_ias","PT_for_tracks_ias",1000,0,1000);
   TH1D* P_compute_ias = new TH1D("P_for_tracks_ias","P_for_tracks_ias",1000,0,1000);

   TH1D* PT_hscp_compute_ias = new TH1D("HSCP_PT_for_tracks_ias_preselection","HSCP_PT_for_tracks_ias_post_presel",1000,0,1000);
   TH1D* P_hscp_compute_ias = new TH1D("HSCP_P_for_tracks_ias_preselection","HSCP_P_for_tracks_ias_post_presel",1000,0,1000);

   TH2D* PU_VS_NPV = new TH2D( "PU_VS_NPV","PU_VS_NPV",100,0,100,100,0,100);

   // P and PT distrib for IAS tracks choosen
   TH2D* IAS_VS_PT = new TH2D( "IAS_VS_PT","IAS_VS_PT",50,0,1,200,0,1000);
   TH1D* PATHLENGHT_BIN_ETA_0_01 = new TH1D( "PATHLENGTH_BIN_ETA_0_01","Pathlength_for_eta_bin_0_01",Path_NBins, Path_Min, Path_Max);

   vector <TH1F*> ias_bin_pt;
   ias_bin_pt.resize(10);

   vector <TH1F*> ias_bin_pt_PUcorr;
   ias_bin_pt_PUcorr.resize(10);


   double itg_norm[10] = { 0 };
   double itg_norm_PU[10] = { 0 };

   
   for (int m = 0 ; m < ias_intervals ; m++){
      double lowerlim = what_bins[m];
      double upperlim = what_bins[m+1];
      string bg = std::to_string(lowerlim).substr(0, std::to_string(lowerlim).find(".") + precisionVal + 1);
      string end = std::to_string(upperlim).substr(0, std::to_string(upperlim).find(".") + precisionVal + 1);
      string trf = " < PU < ";
      string nameax = bg + trf + end; 
      //ias_axis_names[m] = nameax;
   }

   for(int i = 0; i < 10;i++){
       int p = (i+1) * 20;
       string trf = "IAS_BIN_PT_"+ to_string(p+40) + "_to_" + to_string(p+60);
       string trfPU = "IAS_BIN_PT_"+ to_string(p+40) + "_to_" + to_string(p+60) + "_PU_corr";
       ias_bin_pt[i] = new TH1F(trf.c_str(),trf.c_str(),50,0,1);
       ias_bin_pt[i]->Sumw2();
   
       ias_bin_pt_PUcorr[i] = new TH1F(trfPU.c_str(),trfPU.c_str(),50,0,1);
       ias_bin_pt_PUcorr[i]->Sumw2();
   } 

   if(UseTemplatesForPUReweighting) {
       Charge_Vs_Path_PU_corr->Sumw2();
   }
   HHitPix->Sumw2();
   HHitProfilePix->Sumw2();
   HHit2DPix->Sumw2();
   HHit2DPix_NoM->Sumw2();
   HHitStrip->Sumw2();
   HHitProfileStrip->Sumw2();
   HHit2DStrip->Sumw2();
   HHit2DStrip_NoM->Sumw2();
   Charge_Vs_Path->Sumw2();

   Charge_Vs_Path_PU_below_20->Sumw2();
   Charge_Vs_Path_PU_between->Sumw2();
   Charge_Vs_Path_PU_above_40->Sumw2();

   Charge_Vs_Path_noL1->Sumw2();
   Charge_Vs_Path_NoM->Sumw2();
   Charge_Vs_Path_noL1_NoM->Sumw2();
   NHSCP->Sumw2();
   NTRK->Sumw2();
   N_CLU_HSCP->Sumw2();
   N_CLU_TRK->Sumw2();
   IAS_triple_pt50->Sumw2();
   IAS_single_pt50->Sumw2();


   Ias_all_triple_cutIH->Sumw2();
   Ias_all_base_cutIH->Sumw2();
   Ias_all_triple_nocut->Sumw2();
   Ias_all_base_nocut->Sumw2();

   Ias_toy->Sumw2();
   for(int j = 0; j< ias_intervals; j++){
       for (int i = 0; i < nPE ; i++){
           PE_Ias[j][i]->Sumw2();
           PE_Ias_cutih[j][i]->Sumw2();
       }
   }

   IhStripVsP->Sumw2();
   IhStripVsPtight->Sumw2();

   IhBestStripVsP->Sumw2();
   IhBestStripVsPtight->Sumw2();

   IhStripVsP_presk->Sumw2();
   IhStripVsPtight_presk->Sumw2();

   IhBestStripVsP_presk->Sumw2();
   IhBestStripVsPtight_presk->Sumw2();

   IhBestvsP_10_45_central->Sumw2();
   IhBestvsP_10_45_tight->Sumw2();
 
   IhBestvsIas_p_5_100->Sumw2();
   IhBestvsIas_p_10_45->Sumw2();
   IasStripVsIh0noL1_p_10_45->Sumw2();
   cout << "after cut ihvsIas" << endl;
   PT_compute_ias->Sumw2();
   P_compute_ias->Sumw2();
   P_hscp_compute_ias->Sumw2();
   PT_hscp_compute_ias->Sumw2();

   NPV_all->Sumw2();
   NPV_presel->Sumw2();



   MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->Sumw2();
   MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Sumw2();
   IAS_TRIPLE_PT_50_60_selection->Sumw2();
   IAS_SINGLE_PT_50_60_selection->Sumw2();


   MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->GetXaxis()->SetTitle("MASS");
   MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->GetXaxis()->SetTitle("MASS");

   MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->Sumw2();
   MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->Sumw2();
   MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->GetXaxis()->SetTitle("MASS");
   MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->GetXaxis()->SetTitle("MASS");




   Ias_1_VS_Delta_Ias_1_2->Sumw2();


   Ias_vs_PU_above_40_triple->Sumw2();





   Is_all_base->Sumw2();

   for (int i = 0 ; i < ias_intervals ; i++){
      Ias_when_PU[i]->Sumw2();
      Ias_when_PU_triple[i]->Sumw2();
      Ias_when_PU_ih_cut[i]->Sumw2();
      Ias_when_PU_ih_cut_triple[i]->Sumw2();
      Is_when_PU[i]->Sumw2();
      Is_when_PU_triple[i]->Sumw2();
      Is_when_PU_ih_cut[i]->Sumw2();
      Is_when_PU_ih_cut_triple[i]->Sumw2();
      Charge_Vs_Path_PU[i]->Sumw2();   
      PU_distrib[i]->Sumw2();
   }


   
   Ias_below_5GeV_PU_below_20_base->Sumw2();
   Ias_below_5GeV_PU_between_base->Sumw2();
   Ias_below_5GeV_PU_above_40_base->Sumw2();


   Ias_below_5GeV_PU_below_20_triple->Sumw2();
   Ias_below_5GeV_PU_between_triple->Sumw2();
   Ias_below_5GeV_PU_above_40_triple->Sumw2();

   Ias_when_PU_between_base->Sumw2();
   Ias_when_PU_above_40_base->Sumw2();
   
   Ias_when_PU_in_16_18_base->Sumw2();
   Ias_when_PU_in_16_18_triple->Sumw2();

   Ias_when_PU_in_30_32_base->Sumw2();
   Ias_when_PU_in_30_32_triple->Sumw2();
   

   HNtracks->Sumw2();
   HNtracks1->Sumw2();
   HNtracks20->Sumw2();

   Htrackpt->Sumw2();
   Htracketa->Sumw2();
   Htracketa_lowp->Sumw2();
   Htrackphi->Sumw2();
   Htracknhit->Sumw2();

   Htrackih_reco->Sumw2();
   Htrackih_pix->Sumw2();
   Htrackih_strip->Sumw2();
   Htrackdedx_pix->Sumw2();
   Htrackdedx_strip->Sumw2();
   Htrackias->Sumw2();
   Htrackiasall->Sumw2();

   Htrackih_lowp->Sumw2();
   Htrackih0noL1_lowp->Sumw2();
   Htrackih0_lowp->Sumw2();
   Htrackih_pix_lowp->Sumw2();
   Htrackih_strip_lowp->Sumw2();
   Htrackdedx_pix_lowp->Sumw2();
   Htrackdedx_strip_lowp->Sumw2();
   Htrackdedx_strip_lowp1->Sumw2();
   Htrackdedx_strip_lowp2->Sumw2();
   Htrackias_lowp->Sumw2();
   Htrackiasall_lowp->Sumw2();


   PU_VS_NPV->Sumw2();
   PATHLENGHT_BIN_ETA_0_01->Sumw2();   
   TptGivenPathlenght1->Sumw2();

   Nsat->Sumw2();
   NPix->Sumw2();
   NStrip->Sumw2();
   dEdXVsP->Sumw2();
   dEdXVsP_lowp->Sumw2();
   dEdXVsP_lowp2->Sumw2();
   dEdXpixVsP->Sumw2();
   dEdXpixVsP_lowp->Sumw2();
   dEdXpixVsP_lowp2->Sumw2();
   dEdXstripVsP->Sumw2();
   dEdXstripVsP_lowp->Sumw2();
   dEdXstripVsP_lowp2->Sumw2();
   dEdX0VsP_lowp->Sumw2();
   dEdX0VsP_lowp2->Sumw2();
   dEdX0stripVsP_lowp->Sumw2();
   dEdX0stripVsP_lowp2->Sumw2();
   dEdX0noL1VsP_lowp->Sumw2();
   dEdX0noL1VsP_lowp2->Sumw2();
   dEdX0noL1VsP_eta1_lowp->Sumw2();
   dEdX0noL1VsP_eta1_lowp2->Sumw2();
   dEdX0noL1VsP_eta2_lowp->Sumw2();
   dEdX0noL1VsP_eta2_lowp->Sumw2();
   dEdX0noL1VsP_eta3_lowp->Sumw2();
   dEdX0noL1VsP_eta3_lowp2->Sumw2();
   dEdX0noL1VsP_pu1_lowp->Sumw2();
   dEdX0noL1VsP_pu1_lowp2->Sumw2();
   dEdX0noL1VsP_pu2_lowp->Sumw2();
   dEdX0noL1VsP_pu2_lowp2->Sumw2();
   dEdX0noL1VsP_pu3_lowp->Sumw2();
   dEdX0noL1VsP_pu3_lowp2->Sumw2();
   dEdX0noL1VsP_pu4_lowp->Sumw2();
   dEdX0noL1VsP_pu4_lowp2->Sumw2();
   dEdXHDnoL1VsP_lowp->Sumw2();
   dEdXHDnoL1VsP_lowp2->Sumw2();
   dEdX0pixnoL1VsP_lowp->Sumw2();
   dEdX0pixnoL1VsP_lowp2->Sumw2();
   dEdXstripVsEta_lowp->Sumw2();
   dEstrVsdE_lowp->Sumw2();
   dEdXstripVsNhit_lowp->Sumw2();
   dEdXstripVsNhittrunc_lowp->Sumw2();
   dEdXstripVsCharge_lowp->Sumw2();
   EtaVsPhi_nhit->Sumw2();
   HSCP_dEdXpixVsstrip->Sumw2();
   HSCP_dEdX0pixVsstrip->Sumw2();
   HSCP_dEdXstripVsall->Sumw2();
   HSCP_dEdXpixVsall->Sumw2();

   Charge_pixl1->Sumw2();
   Charge_pixl2->Sumw2();
   Charge_pixl3->Sumw2();
   Charge_pixl4->Sumw2();
   Charge_pixd1->Sumw2();
   Charge_pixd2->Sumw2();
   Charge_pixd3->Sumw2();
   Charge_pixr1->Sumw2();
   Charge_pixr2->Sumw2();
   Charge_tib1->Sumw2();
   Charge_tib2->Sumw2();
   Charge_tib3->Sumw2();
   Charge_tib4->Sumw2();
   Charge_tob1->Sumw2();
   Charge_tob2->Sumw2();
   Charge_tob3->Sumw2();
   Charge_tob4->Sumw2();
   Charge_tob5->Sumw2();
   Charge_tob6->Sumw2();
   Charge_tid1->Sumw2();
   Charge_tid2->Sumw2();
   Charge_tid3->Sumw2();
   Charge_tec1->Sumw2();
   Charge_tec2->Sumw2();
   Charge_tec3->Sumw2();
   Charge_tec4->Sumw2();
   Charge_tec5->Sumw2();
   Charge_tec6->Sumw2();
   Charge_tec7->Sumw2();
   Charge_tec8->Sumw2();
   Charge_tec9->Sumw2();

   LowpCharge_tib1->Sumw2();
   LowpCharge_tib2->Sumw2();
   LowpCharge_tib3->Sumw2();
   LowpCharge_tib4->Sumw2();
   LowpCharge_tob1->Sumw2();
   LowpCharge_tob2->Sumw2();
   LowpCharge_tob3->Sumw2();
   LowpCharge_tob4->Sumw2();
   LowpCharge_tob5->Sumw2();
   LowpCharge_tob6->Sumw2();
   LowpCharge_tid1->Sumw2();
   LowpCharge_tid2->Sumw2();
   LowpCharge_tid3->Sumw2();
   LowpCharge_tec1->Sumw2();
   LowpCharge_tec2->Sumw2();
   LowpCharge_tec3->Sumw2();
   LowpCharge_tec4->Sumw2();
   LowpCharge_tec5->Sumw2();
   LowpCharge_tec6->Sumw2();
   LowpCharge_tec7->Sumw2();
   LowpCharge_tec8->Sumw2();
   LowpCharge_tec9->Sumw2();
   LowpCharge_pixl1->Sumw2();
   LowpCharge_pixl2->Sumw2();
   LowpCharge_pixl3->Sumw2();
   LowpCharge_pixl4->Sumw2();
   LowpCharge_pixd1->Sumw2();
   LowpCharge_pixd2->Sumw2();
   LowpCharge_pixd3->Sumw2();
   LowpCharge_pixr1->Sumw2();
   LowpCharge_pixr2->Sumw2();
   LowpCharge_Eta_pix->Sumw2();
   LowpCharge_Eta_strip->Sumw2();

   ChargeVsRun_pixl1->Sumw2();
   ChargeVsRun_pixl2->Sumw2();
   ChargeVsRun_pixl3->Sumw2();
   ChargeVsRun_pixl4->Sumw2();
   ChargeVsRun_pixd1->Sumw2();
   ChargeVsRun_pixd2->Sumw2();
   ChargeVsRun_pixd3->Sumw2();
   ChargeVsRun_pixr1->Sumw2();
   ChargeVsRun_pixr2->Sumw2();
   ChargeVsRun_tib1->Sumw2();
   ChargeVsRun_tib2->Sumw2();
   ChargeVsRun_tib3->Sumw2();
   ChargeVsRun_tib4->Sumw2();
   ChargeVsRun_tob1->Sumw2();
   ChargeVsRun_tob2->Sumw2();
   ChargeVsRun_tob3->Sumw2();
   ChargeVsRun_tob4->Sumw2();
   ChargeVsRun_tob5->Sumw2();
   ChargeVsRun_tob6->Sumw2();
   ChargeVsRun_tid1->Sumw2();
   ChargeVsRun_tid2->Sumw2();
   ChargeVsRun_tid3->Sumw2();
   ChargeVsRun_tec1->Sumw2();
   ChargeVsRun_tec2->Sumw2();
   ChargeVsRun_tec3->Sumw2();
   ChargeVsRun_tec4->Sumw2();
   ChargeVsRun_tec5->Sumw2();
   ChargeVsRun_tec6->Sumw2();
   ChargeVsRun_tec7->Sumw2();
   ChargeVsRun_tec8->Sumw2();
   ChargeVsRun_tec9->Sumw2();

   ZooChargeVsRun_pixl1->Sumw2();
   ZooChargeVsRun_pixl2->Sumw2();
   ZooChargeVsRun_pixl3->Sumw2();
   ZooChargeVsRun_pixl4->Sumw2();
   ZooChargeVsRun_pixd1->Sumw2();
   ZooChargeVsRun_pixd2->Sumw2();
   ZooChargeVsRun_pixd3->Sumw2();
   ZooChargeVsRun_pixr1->Sumw2();
   ZooChargeVsRun_pixr2->Sumw2();
   ZooChargeVsRun_tib1->Sumw2();
   ZooChargeVsRun_tib2->Sumw2();
   ZooChargeVsRun_tib3->Sumw2();
   ZooChargeVsRun_tib4->Sumw2();
   ZooChargeVsRun_tob1->Sumw2();
   ZooChargeVsRun_tob2->Sumw2();
   ZooChargeVsRun_tob3->Sumw2();
   ZooChargeVsRun_tob4->Sumw2();
   ZooChargeVsRun_tob5->Sumw2();
   ZooChargeVsRun_tob6->Sumw2();
   ZooChargeVsRun_tid1->Sumw2();
   ZooChargeVsRun_tid2->Sumw2();
   ZooChargeVsRun_tid3->Sumw2();
   ZooChargeVsRun_tec1->Sumw2();
   ZooChargeVsRun_tec2->Sumw2();
   ZooChargeVsRun_tec3->Sumw2();
   ZooChargeVsRun_tec4->Sumw2();
   ZooChargeVsRun_tec5->Sumw2();
   ZooChargeVsRun_tec6->Sumw2();
   ZooChargeVsRun_tec7->Sumw2();
   ZooChargeVsRun_tec8->Sumw2();
   ZooChargeVsRun_tec9->Sumw2();


   dEdXVsRun->Sumw2();
   dEdXpixVsRun->Sumw2();
   dEdXstripVsRun->Sumw2();
   dEdXNoL1pixVsRun->Sumw2();
   dEdXNoL1VsRun->Sumw2();
   dEdXHiDropVsRun->Sumw2();
   dEdXpixHiDropVsRun->Sumw2();
   dEdXstripHiDropVsRun->Sumw2();
   dEdXHiDropNoL1VsRun->Sumw2();
   dEdX0VsRun->Sumw2();
   dEdX0pixVsRun->Sumw2();
   dEdX0stripVsRun->Sumw2();
   dEdX0NoL1pixVsRun->Sumw2();
   dEdX0NoL1VsRun->Sumw2();
   MassStripVsRun->Sumw2();
   MassNoL1VsRun->Sumw2();
   dEdX4VsRun->Sumw2();
   dEdX4pixVsRun->Sumw2();
   dEdX4stripVsRun->Sumw2();
   dEdX40VsRun->Sumw2();
   dEdX40pixVsRun->Sumw2();
   dEdX40stripVsRun->Sumw2();
   bg_dEdX0NoL1VsRun->Sumw2();
   iasNoL1VsRun->Sumw2();
   iasAllVsRun->Sumw2();

   FMIP4VsRun->Sumw2();
   FMIP3p5VsRun->Sumw2();
   FMIP3p2VsRun->Sumw2();
   FMIP4VsEta->Sumw2();

   NmeasVsRun->Sumw2();
   NmeasPixVsRun->Sumw2();
   NmeasStrVsRun->Sumw2();
   Nmeas0VsRun->Sumw2();
   NmeasPix0VsRun->Sumw2();
   NmeasStr0VsRun->Sumw2();
   NsatVsRun->Sumw2();
   NsatPixVsRun->Sumw2();
   NsatStrVsRun->Sumw2();
   Nsat0VsRun->Sumw2();
   NsatPix0VsRun->Sumw2();
   NsatStr0VsRun->Sumw2();

   ptVsRun->Sumw2();
   nPVVsRun->Sumw2();
   invBVsRun->Sumw2();
   timeVsRun->Sumw2();
   lumiVsRun->Sumw2();
   errinvBVsRun->Sumw2();
   invBDTVsRun->Sumw2();
   invBCSCVsRun->Sumw2();
   invBnewVsRun->Sumw2();
   invBnewDTVsRun->Sumw2();
   invBnewCSCVsRun->Sumw2();

   HSCP_dEdX->Sumw2(); 
   HSCP_dEdXpix->Sumw2(); 
   HSCP_dEdXstrip->Sumw2(); 
   HSCP_dEdX0->Sumw2(); 
   HSCP_dEdX0pix->Sumw2(); 
   HSCP_dEdX0strip->Sumw2(); 

   NPV_ias_single_cutih->Sumw2();
   NPV_ias_triple_cutih->Sumw2();

   PU_ias_single_0p3->Sumw2();
   PU_ias_triple_0p3->Sumw2();

   P_LOWPU_pt50->Sumw2();
   P_HIGHPU_pt50->Sumw2();

   P_PU_28_29_pt50->Sumw2();
   P_PU_30_31_pt50->Sumw2();
 
   P_PU_28_29_pt50->GetXaxis()->SetTitle("P [GeV/c]");
   P_PU_30_31_pt50->GetXaxis()->SetTitle("P [GeV/c]");


   P_PU_28_30_pt50->Sumw2();
   P_PU_32_34_pt50->Sumw2();

   P_PU_28_30_pt50->GetXaxis()->SetTitle("P [GeV/c]");
   P_PU_32_34_pt50->GetXaxis()->SetTitle("P [GeV/c]");



   P_PU_58_60_pt50->Sumw2();
   P_PU_08_10_pt50->Sumw2();

   P_PU_58_60_pt50->GetXaxis()->SetTitle("P [GeV/c]");
   P_PU_08_10_pt50->GetXaxis()->SetTitle("P [GeV/c]");

   NPV_mass_cutih->Sumw2();
   NPV_mass_cut3p29->Sumw2();

   P_triple_pt50_ias_qtl_80_90->Sumw2();

   NPV_ih_cutih->Sumw2();
   NPV_ih_cut3p29->Sumw2();
   NPV_ih_nocut->Sumw2();
   NPV_ih_cut1->Sumw2();
   NPV_ih_strip_cutih->Sumw2();
   NPV_ih_strip_cut3p29->Sumw2();
   NPV_ih_strip_nocut->Sumw2();
   NPV_ih_strip_cut1->Sumw2();
   IH0_noL1->Sumw2();

   MASS_SINGLE_ihcut->Sumw2();
   MASS_TRIPLE_ihcut->Sumw2();
   MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Sumw2();
   MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Sumw2();
   MASS_SINGLE_ihcut_base_notriple_ias_0p1_ptsupp50_qtl_80_90_ias->Sumw2();
   MASS_TRIPLE_ihcut_triple_nobase_ias_0p1_ptsupp50_qtl_80_90_ias->Sumw2();
   MASS_SINGLE_ihcut_base_notriple_ias_0p1->Sumw2();
   MASS_TRIPLE_ihcut_triple_nobase_ias_0p1->Sumw2();
   MASS_SINGLE_nocut->Sumw2();
   MASS_TRIPLE_nocut->Sumw2();
   Ias_simple_Vs_Ih0noL1_p_10_45->Sumw2();
   Ias_triple_Vs_Ih0noL1_p_10_45->Sumw2();
   IAS_simple_when_triple_above_0p115->Sumw2();
   IAS_triple_when_simple_above_0p115->Sumw2();

   CHARGE_OVER_PATHLENGTH_GIVEN_MODULE->Sumw2();
   CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP->Sumw2();
   CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_LOWPU->Sumw2();
   CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_HIGHPU->Sumw2();
   CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_NOSATURATED->Sumw2();
   CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_PIXEL->Sumw2();
   HSCP_MassIh->Sumw2(); 
   HSCP_MassIh0->Sumw2(); 
   HSCP_MassIhstrip->Sumw2(); 
   HSCP_MassIh0noL1->Sumw2(); 
   HSCP_MassIh0noL1_2->Sumw2(); 
   HSCP_MassIh0noL1_3->Sumw2(); 
   HSCP_MassIh0noL1_11->Sumw2(); 
   HSCP_MassIh0noL1_12->Sumw2(); 
   HSCP_MassIh0noL1_13->Sumw2(); 
   HSCP_MassIh0noL1_1s1->Sumw2(); 
   HSCP_MassIh0noL1_1s2->Sumw2(); 
   HSCP_MassIh0noL1_1s3->Sumw2(); 
   HSCP_MassIh0noL1_2s1->Sumw2(); 
   HSCP_MassIh0noL1_2s2->Sumw2(); 
   HSCP_MassIh0noL1_2s3->Sumw2(); 
   HSCP_MassIh0noL1_3s1->Sumw2(); 
   HSCP_MassIh0noL1_3s2->Sumw2(); 
   HSCP_MassIh0noL1_3s3->Sumw2(); 
   HSCP_MassIh0strip->Sumw2(); 
   HSCP_MassIhHDnoL1->Sumw2(); 
   HSCP_MassTOF->Sumw2(); 
   HSCP2d_MassTOFvsIh->Sumw2(); 
   HSCP2d_MassIh->Sumw2(); 
   HSCP2d_MassIh0->Sumw2(); 
   HSCP2d_MassIhstrip->Sumw2(); 
   HSCP2d_MassIh0noL1->Sumw2(); 
   HSCP2d_MassIh0strip->Sumw2(); 
   HSCP2d_MassIhHDnoL1->Sumw2(); 
   HSCP2d_Mass_pix_strip15->Sumw2(); 
   HSCP2d_Mass_pix_strip0->Sumw2(); 
   HSCP2d_Mass_pix_strip->Sumw2(); 
   HSCP_MassDiff_pix_strip0->Sumw2(); 
   HSCP_MassDiff_pix_strip15->Sumw2(); 
   HSCP_MassResol_pix_strip0->Sumw2(); 
   HSCP_MassResol_pix_strip15->Sumw2(); 
   HSCP_dEdXHiDrop->Sumw2();
   HSCP_dEdXstripHiDrop->Sumw2();
   HSCP_dEdXpixHiDrop->Sumw2();
   HSCP_dEdXHiDropNoL1->Sumw2();
   HSCP_dEdX0NoL1->Sumw2();
   NB_PV->Sumw2();
   lowp_MassIh->Sumw2(); 
   lowp_MassIh0->Sumw2(); 
   lowp_MassIhstrip->Sumw2(); 
   lowp_MassIh0noL1->Sumw2(); 
   lowp_MassIh0noL1_2->Sumw2(); 
   lowp_MassIh0noL1_3->Sumw2(); 
   lowp_MassIh0noL1_11->Sumw2(); 
   lowp_MassIh0noL1_12->Sumw2(); 
   lowp_MassIh0noL1_13->Sumw2(); 
   lowp_MassIh0noL1_1s1->Sumw2(); 
   lowp_MassIh0noL1_1s2->Sumw2(); 
   lowp_MassIh0noL1_1s3->Sumw2(); 
   lowp_MassIh0noL1_2s1->Sumw2(); 
   lowp_MassIh0noL1_2s2->Sumw2(); 
   lowp_MassIh0noL1_2s3->Sumw2(); 
   lowp_MassIh0noL1_3s1->Sumw2(); 
   lowp_MassIh0noL1_3s2->Sumw2(); 
   lowp_MassIh0noL1_3s3->Sumw2(); 
   lowp_MassIh0strip->Sumw2(); 
   lowp_MassIhHDnoL1->Sumw2(); 
   lowp2d_MassIh->Sumw2(); 
   lowp2d_MassIh0->Sumw2(); 
   lowp2d_MassIhstrip->Sumw2(); 
   lowp2d_MassIh0noL1->Sumw2(); 
   lowp2d_MassIh0strip->Sumw2(); 
   lowp2d_MassIhHDnoL1->Sumw2(); 
   lowp_dEdXpixVsstrip->Sumw2();
   lowp_dEdX0pixVsstrip->Sumw2();
   lowp2d_Mass_pix_strip15->Sumw2(); 
   lowp2d_Mass_pix_strip0->Sumw2(); 
   bg_lowp2d_Mass_pix_strip0->Sumw2(); 
   bg_transf_Mass->Sumw2(); 
   lowp_MassDiff_pix_strip0->Sumw2(); 
   lowp_MassDiff_pix_strip15->Sumw2(); 
   bg_dEdXVsIas->Sumw2(); 
   bg_test_event_Mass->Sumw2(); 
   bg_test_prefiring_Mass->Sumw2(); 
   bg_test_nohem_Mass->Sumw2(); 

   HSCP_FMIP4->Sumw2(); 
   HSCP_FMIP3p5->Sumw2(); 
   HSCP_FMIP3p2->Sumw2(); 
   HSCP_iasnol1->Sumw2(); 
   HSCP_iasall->Sumw2(); 
   HSCP_iasstrip->Sumw2(); 
   HSCP_iaspix->Sumw2(); 
   HSCP_probQ->Sumw2();
   HSCP_probQNoL1->Sumw2();
   HSCP_probXY->Sumw2();
   HSCP_probXYNoL1->Sumw2();
   probQVsRun->Sumw2();
   probQNoL1VsRun->Sumw2();
   probXYVsRun->Sumw2();
   probXYNoL1VsRun->Sumw2();
   probQVsIas->Sumw2();

   IAS_VS_PT->Sumw2();


   HSCP_pt->Sumw2(); 
   HSCP_eta->Sumw2(); 
   HSCP_iso_eop->Sumw2(); 
   nPV->Sumw2(); 
   HSCP_invB->Sumw2(); 
   HSCP_errinvB->Sumw2(); 
   HSCP_invBDT->Sumw2(); 
   HSCP_invBCSC->Sumw2(); 
   HSCP_time->Sumw2(); 
   HSCP_npix->Sumw2(); 
   HSCP_nstrip->Sumw2(); 
   HSCP_nmpix->Sumw2(); 
   HSCP_nmstrip->Sumw2(); 
   HSCP_nratio->Sumw2(); 
   HSCP_nmratio->Sumw2(); 

   R1_StdEdXVsEvent->Sumw2();
   R1_StdEdXVsLumi->Sumw2();
   R1_LumiVsEvent->Sumw2();
   R1_nPVVsEvent->Sumw2();
   R1_CandVsEvent->Sumw2();
   R2_StdEdXVsEvent->Sumw2();
   R2_StdEdXVsLumi->Sumw2();
   R2_LumiVsEvent->Sumw2();
   R2_nPVVsEvent->Sumw2();
   R2_CandVsEvent->Sumw2();
   R3_StdEdXVsEvent->Sumw2();
   R3_StdEdXVsLumi->Sumw2();
   R3_LumiVsEvent->Sumw2();
   R3_nPVVsEvent->Sumw2();
   R3_CandVsEvent->Sumw2();


   TString outputfilename="analysis_ul";
   if (year==2016) outputfilename+="_2016";
   else if (year==2017) outputfilename+="_2017";
   else if (year==2018) outputfilename+="_2018";
   outputfilename+=Letter;
   if (!dataFlag) outputfilename+="_MC";
   outputfilename+="_15mars.root";
   TFile* OutputHisto = new TFile(outputfilename,"RECREATE");
   TString templateFileName="template";
   if (year==2016) templateFileName+="_2016";
   else if (year==2017) templateFileName+="_2017";
   else if (year==2018) templateFileName+="_2018";
   templateFileName+=Letter;
   if (!dataFlag) templateFileName+="_MC";

   if(TemplateIso) templateFileName+="_15mars_selection_5_bin_eta_2p1_pmin_" + low_bound +"_pmax_" + high_bound + "_final_HSCP_iso15.root";
   else templateFileName+="_15mars_selection_5_bin_eta_2p1_pmin_" + low_bound +"_pmax_" + high_bound + "_final_HSCP_noiso_new_stat1_tenth.root";

   TFile* OutputTemplate;
   if (writeTemplateOnDisk || writeTptHSCP) OutputTemplate = new TFile(templateFileName,"RECREATE");
   if(!(writeTemplateOnDisk || writeTptHSCP)) OutputTemplate = new TFile(templateFileName);

   TString outputIAS;   
   if(!(writeTemplateOnDisk || writeTptHSCP)){
       outputIAS ="ias_study_2018";
       outputIAS+=Letter;
       outputIAS+="_template_";

       if(TemplateIso) outputIAS+="iso15_";
       else outputIAS+="noiso_";
      
       outputIAS+="2018A_15mars_eta_1_pubins_";
       string tr2 = std::to_string(ias_intervals);
       outputIAS+=tr2;
       outputIAS+="_minp_";
       outputIAS+=low_bound;
       outputIAS+="_pmax_";
       outputIAS+=high_bound;
       outputIAS+="_final_HSCP_";

       if(StudyIso) outputIAS+="iso15_v2_tpt_stat1_tenth.root"; 
       else outputIAS+="noiso_v2.root";
       
   }
   else{
       outputIAS = "not_filled_ias_study_writing_template_HSCP.root";
   }
   TFile* OutIasStudy = new TFile(outputIAS,"RECREATE");
//   loadSFPixelCalib();
   loadSFPixelTamas();
   bool PixelCorr2Apply = true;
   if (!dataFlag) PixelCorr2Apply = false;  // No Pixel correction for MC !


//   loadDeDxTemplates

    if (boolDeDxTemp) {
     if (!dataFlag) {
//     dEdxTemplatesNoL1 = loadDeDxTemplate("templateMC_w17_MC_21jan.root", "Charge_Vs_Path_noL1_NoM", true);
      if (year==2017) {
       dEdxTemplatesNoL1 = loadDeDxTemplate("templateMC_w17_MC_21jan.root", "Charge_Vs_Path_noL1", true);
       dEdxTemplatesAll = loadDeDxTemplate("templateMC_w17_MC_21jan.root", "Charge_Vs_Path", true);
      }
      else if (year==2018){
       dEdxTemplatesNoL1 = loadDeDxTemplate("template_2018MC_w18_MC_28feb.root", "Charge_Vs_Path_noL1", true);
       dEdxTemplatesAll = loadDeDxTemplate("template_2018MC_w18_MC_28feb.root", "Charge_Vs_Path", true);
      }
     }
     else {
     // note : for data, we should write something smart to read template per era
      std::string TemplateDataName = "template";
      if (year>2016) { 
//      if (year==2016) TemplateDataName+="_2016";
       if (year==2017) TemplateDataName+="_2017";
       else if (year==2018) TemplateDataName+="_2018";
       TemplateDataName+=Letter;
       if (year==2017) TemplateDataName+="_21jan.root";
       else if (year==2018) TemplateDataName+="_28feb.root";
      }
      else { 
        TemplateDataName+="_2017B_21jan.root";
      }

//      dEdxTemplatesNoL1 = loadDeDxTemplate(TemplateDataName, "Charge_Vs_Path_noL1_NoM", true);
      dEdxTemplatesNoL1 = loadDeDxTemplate(TemplateDataName, "Charge_Vs_Path_noL1", true);
      dEdxTemplatesAll = loadDeDxTemplate(TemplateDataName, "Charge_Vs_Path", true);
     }
    }

   //Added for tests 
   if (UsePURwtHSCP || UseTemplatesForPUReweighting){ 
       // --------- LOADING ROOT FILE FOR TEMPLATES -------------
       std::string TemplateStudyIasBis;
       std::string TemplateStudyIas = "template_2018";
      
       TemplateStudyIas+=Letter; // chose for cross template, need to change by hand (lazy :s)
       if(UseTemplatesForPUReweighting){ 
           TemplateStudyIas+="_15mars_selection_5_bin_eta_2p1_pmin_20_pmax_50_final";
       }
       else{
           TemplateStudyIas+="_15mars_selection_";
           string tr1 = std::to_string(ias_intervals);
           TemplateStudyIas+= tr1;
           TemplateStudyIas+= "_bin_eta_2p1_pmin_";
           TemplateStudyIas+=low_bound;
           TemplateStudyIas+="_pmax_";
           TemplateStudyIas+=high_bound;
           TemplateStudyIas+="_final_HSCP";
       }
      
       TemplateStudyIasBis = TemplateStudyIas;

       if (TemplateIso) TemplateStudyIas+="_iso50_new.root";
       else{
           TemplateStudyIas+="_noiso_new_stat1_tenth.root";
           TemplateStudyIasBis+="_noiso_new_stat2_tenth.root";
       }
       std::string PE_TemplateStudyIas = "PE_templates_2018";
       PE_TemplateStudyIas+=Letter;
       if (TemplateIso) PE_TemplateStudyIas+="_iso50.root";
       else PE_TemplateStudyIas+="_noiso.root"; 


       dEdxTemplatesAll = loadDeDxTemplate(TemplateStudyIas,"Charge_Vs_Path",true);

       for (int i = 0; i < ias_intervals ; i++){
           dEdxTemplatesPU[i] = loadDeDxTemplate(TemplateStudyIas,names_templates[i].c_str(),true); 
           dEdxTemplatesPUBis[i] = loadDeDxTemplate(TemplateStudyIasBis,names_templates[i].c_str(),true); 
       }  

       if (compute_PE){ 
           for (int j = 0; j< ias_intervals; j++){ 
               for (int i = 0; i < nPE ; i++){
                   std::string name_pe_tpt = "PE_" + to_string(i) + "_PU_" + to_string(what_bins[j]) + "_" + to_string(what_bins[j+1]); 
                   PE_dEdxTemplatesPU[j][i] = loadDeDxTemplate(PE_TemplateStudyIas,name_pe_tpt.c_str(),true);
               }
           }
       }

       dEdxTemplatesPuLow = loadDeDxTemplate(TemplateStudyIas,"Charge_Vs_Path_Low_PU",true);
       dEdxTemplatesPuMedium = loadDeDxTemplate(TemplateStudyIas,"Charge_Vs_Path_Middle_PU",true);
       dEdxTemplatesPuHigh = loadDeDxTemplate(TemplateStudyIas,"Charge_Vs_Path_High_PU",true);
   }

//   nentries = 200000;
  // cout << "run on  " << nentries << " entries " << endl;


   //cout << "run on first third entries, nb = " << nentries/3 << endl;
   //cout << "run on second third entries, between " << (2*nentries/3) << " and " << nentries << endl;


   for (Long64_t jentry= 0 ; jentry < nentries ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%100000 ==0 && jentry!=0) cout << " number of processed events is " << jentry <<  " = " << (100.*jentry)/(1.*nentries) << "%" <<endl;
//      if(jentry%10000 ==0 && jentry!=0) cout << " number of processed events is " << jentry <<  " = " << (100.*jentry)/(1.*nentries) << "%" <<endl;

      if (!dataFlag) runNumber = 300000;

//      if (runNumber<305150 || runNumber>305300) continue;
//      if (runNumber>323232) continue;  // because 2018 is rereco at the moment --> remove this if when running on UL
      
      if ( ! (jentry+(choice_third-1)%3 == 0 && jentry != 0) ){          
          continue;
      }
      
      if(tenth_stat_prod){
          count_tenth+=1; 
          if( !(count_tenth%10 == 0)){
              continue;
          }
      }
      // jentry%3 is stat1
      // (jentry+1)%3 is stat2
      // (jentry+2)%3 is stat3

      // only Muon Trigger at the moment (because I am running on SingleMu data)
//      if (!hlt_mu50 && !hlt_tkmu100 && !hlt_oldmu100) continue;
//     Mix SingleMu & PFMET_MHT
      

      //commented to test because all events have P way too high 

      //if (!hlt_mu50 && !hlt_tkmu100 && !hlt_oldmu100 && !hlt_pfmet_mht) continue;

      // test PU

      /* 
      unsigned int pileup_fromLumi = 0;
      const Handle<LumiScalersCollection> lumiScalers = iEvent.getHandle(lumiScalersToken_);
      if (lumiScalers.isValid() &&  !lumiScalers->empty()){
          LumiScalersCollection::const_iterator scalit = lumiScalers->begin();
          pileup_fromLumi = scalit->pileup();
      }
      cout << "PU from token = " << pileup_fromLumi << " and NPV = : " << nPV << endl;
      PU_VS_NPV->Fill(pileup_fromLumi,nPV);
      */
      

      nPVVsRun->Fill(runNumber,npv);
      nPV->Fill(npv);
      NPV_all->Fill(npv);
      if (npv<1) continue;
      if(npv > 99) nb_above_99_npv+=1;
      //if (jentry > 100000) break;
      float HighestMass=-1;
      // Loop on HSCP candidate

      for (int ihs=0; ihs<nhscp; ihs++) {
          int index_of_the_track=hscp_track_idx[ihs];
          //no presk before, below 10 we dont care right ?
          int presk_hscp = 1;
          if(year!=2016) presk_hscp=track_prescale[index_of_the_track];

          if (index_of_the_track>-1) { 
             bool selection=true;

             //if (track_pt[index_of_the_track]<55) selection=false;
             if(writeTptHSCP){
                 if (abs(track_eta[index_of_the_track])>2.1) selection=false;
             }
             else{
                 if (abs(track_eta[index_of_the_track])>1) selection=false;
             }
             //if (track_nvalidhits[index_of_the_track]<10) selection=false;
             if (track_npixhits[index_of_the_track]<2) selection=false;
             if (track_validfraction[index_of_the_track]<0.8) selection=false;
             if (track_missing[index_of_the_track]>99999) selection=false;
             if (track_validlast[index_of_the_track]<-99999) selection=false;

             
             bool is_high_qual =  (track_qual[index_of_the_track] & (1 << TrackQuality::highPurity)) >> TrackQuality::highPurity ;
             if (!is_high_qual) selection=false;
             if (track_chi2[index_of_the_track]>5) selection=false;

             // which selection for dz and dyx w/r to the primary vertex ? here 1st PV to pass the selection
             if (abs(track_dxy[index_of_the_track])>0.02) selection=false;
             if (abs(track_dz[index_of_the_track])>0.1) selection=false;

             // which DR used for iso ? here 0.3
             
             if (hscp_iso2_tk[ihs]>15) selection=false;


             float eop=(hscp_iso2_ecal[ihs] + hscp_iso2_hcal[ihs])/track_p[index_of_the_track];
             //if (eop>0.3) selection=false;

             // no cut on relative iso : false if hscpIso.Get_TK_SumEt() / track->pt() > 9999999

             // no cut on sigma pT/pT for signal yet :
             //if (track_pterr[index_of_the_track]/track_pt[index_of_the_track]>0.25) selection=false; // commented on September 22, 2021
             
             // no cut in TOF yet : false if tof->nDof() < 8
             //     &&  (dttof->nDof() < 6 || csctof->nDof() < 6)
             // no cut in TOF yet : false if tof->inverseBetaErr() > GlobalMaxTOFErr
            
             

             // no cut yet :  false if  dedxSObj->numberOfMeasurements() < 6 ;
             
             if (selection) {
                PT_hscp_compute_ias->Fill(track_pt[index_of_the_track]);
                P_hscp_compute_ias->Fill(track_p[index_of_the_track]);
        
                NHSCP->Fill(2);
                ptVsRun->Fill(runNumber,track_pt[index_of_the_track]);

                std::vector <float> charge_corr;
                std::vector <float> pathlength;
                std::vector <int> subdetId;
                std::vector <int> moduleGeometry;
                std::vector <bool> bool_cleaning;
                std::vector <bool> mustBeInside;
             
                std::vector <float> charge_corr1;
                std::vector <float> pathlength1;
                std::vector <int> subdetId1;
                std::vector <UInt_t> detId1;
                std::vector <int> moduleGeometry1;
                std::vector <bool> bool_cleaning1;
                std::vector <bool> mustBeInside1;

                std::vector <float> charge_corr2;
                std::vector <float> pathlength2;
                std::vector <int> subdetId2;
                std::vector <int> moduleGeometry2;
                std::vector <bool> bool_cleaning2;
                std::vector <bool> mustBeInside2;

                std::vector <float> charge_corr3;
                std::vector <float> pathlength3;
                std::vector <int> subdetId3;
                std::vector <UInt_t> detId3;
                std::vector <int> moduleGeometry3;
                std::vector <bool> bool_cleaning3;
                std::vector <bool> mustBeInside3;

                std::vector <float> charge_corr5;
                std::vector <float> pathlength5;
                std::vector <int> subdetId5;
                std::vector <UInt_t> detId5;
                std::vector <int> moduleGeometry5;
                std::vector <bool> bool_cleaning5;
                std::vector <bool> mustBeInside5;

                int nstip_=0;
                int npix_=0;
                N_CLU_HSCP->Fill(track_index_hit[hscp_track_idx[ihs]]+track_nhits[hscp_track_idx[ihs]]);
                for (int iclu=track_index_hit[hscp_track_idx[ihs]]; iclu<track_index_hit[hscp_track_idx[ihs]]+track_nhits[hscp_track_idx[ihs]]; iclu++) {
                     float ch1=dedx_charge[iclu];
                     bool clean1=true;
                     if (dedx_subdetid[iclu]>=3) {
                        // strip
                        // without any correction :
//                        ch1=sclus_charge[iclu];  // charge without any correction
//                        clean1=sclus_clusclean[iclu];  
//                      // with Saturation only (but no Xtalk inversion) :
                          nstip_++;
                          float check_charge=0;
                          vector<int> Quncor;
                          for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++) {
                            check_charge+=strip_ampl[istrip];
                            Quncor.push_back(strip_ampl[istrip]);
                          }
                          float deltaDif=check_charge-ch1;
                          if (deltaDif<0) deltaDif*=-1;
                          if (deltaDif>0.001) std::cout << "difference dans le calcul de la charge " << ch1 << " " << check_charge << " --> probleme acces ampl ???? " << std::endl; 
                          vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
                          float newcharge =0;
                          for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
                          ch1=newcharge;
                          clean1=sclus_clusclean2[iclu]; // clusterCleaning with Xtalk inversion and Saturation (September 22, 2021)
                          // for the record : there is a possibility to recompute the cluster cleaning with the Saturation only : clusterCleaning(Qcor, 1, &exitCode)
                          // but it is not what we want here
                     }
                     else {
                        // pixel
                          // float scaling =GetSFPixel(dedx_subdetid[iclu], dedx_detid[iclu], year, abs(track_eta[hscp_track_idx[ihs]]), runNumber);
                          npix_++;
                          if (PixelCorr2Apply) {
                            float scaling =GetSFPixelTamas(dedx_subdetid[iclu], dedx_detid[iclu], year, runNumber);
                            ch1*=scaling;
                          }
                     }
                     if (clean1 && dedx_insideTkMod[iclu]) {
                       // fill both Strip and Pixel
                       charge_corr.push_back(ch1);
                       pathlength.push_back(dedx_pathlength[iclu]);
                       subdetId.push_back(dedx_subdetid[iclu]);
                       moduleGeometry.push_back(dedx_modulgeom[iclu]);
                       mustBeInside.push_back(dedx_insideTkMod[iclu]);
                       bool_cleaning.push_back(clean1);

                       if (dedx_isstrip[iclu]) {
                        // fill only Strip
                        charge_corr1.push_back(ch1);
                        pathlength1.push_back(dedx_pathlength[iclu]);
                        subdetId1.push_back(dedx_subdetid[iclu]);
                        detId1.push_back(dedx_detid[iclu]);
                        moduleGeometry1.push_back(dedx_modulgeom[iclu]);
                        mustBeInside1.push_back(dedx_insideTkMod[iclu]);
                        bool_cleaning1.push_back(clean1);
                       }
                       else {
                        // fill only Pixel
                        charge_corr2.push_back(ch1);
                        pathlength2.push_back(dedx_pathlength[iclu]);
                        subdetId2.push_back(dedx_subdetid[iclu]);
                        moduleGeometry2.push_back(dedx_modulgeom[iclu]);
                        mustBeInside2.push_back(dedx_insideTkMod[iclu]);
                        bool_cleaning2.push_back(clean1);
                       }
                     
                       bool no_in_L1_pixel =true;
                       int info_layr=GetLayerLabel(dedx_subdetid[iclu], dedx_detid[iclu],year);
                       if (info_layr==23) no_in_L1_pixel=false;
                       if (no_in_L1_pixel) {
                         // fill both Strip and Pixel but without L1 Pixel
                         charge_corr3.push_back(ch1);
                         pathlength3.push_back(dedx_pathlength[iclu]);
                         subdetId3.push_back(dedx_subdetid[iclu]);
                         detId3.push_back(dedx_detid[iclu]);
                         moduleGeometry3.push_back(dedx_modulgeom[iclu]);
                         mustBeInside3.push_back(dedx_insideTkMod[iclu]);
                         bool_cleaning3.push_back(clean1);
                         if (dedx_subdetid[iclu]<3) {  // pixel only
                            // fill only Pixel, without L1 Pixel
                            charge_corr5.push_back(ch1);
                            pathlength5.push_back(dedx_pathlength[iclu]);
                            subdetId5.push_back(dedx_subdetid[iclu]);
                            detId5.push_back(dedx_detid[iclu]);
                            moduleGeometry5.push_back(dedx_modulgeom[iclu]);
                            mustBeInside5.push_back(dedx_insideTkMod[iclu]);
                            bool_cleaning5.push_back(clean1);
                         }
                       }
                       double SFtest = dEdxSF[0];
                       float norm_mult_test = 265;
                       //cout << "modul geom is strip : "<< dedx_isstrip[iclu] << " , geom = " << dedx_modulgeom[iclu] << " , track p = " << track_p[index_of_the_track] << " , pathlength * 10 = " << dedx_pathlength[iclu]*10 << " charge*SF / pathlength*10 = " << SFtest*ch1/(dedx_pathlength[iclu]*10)<<endl; 


                       double low_bin_center = PATHLENGTH_GIVEN_MODULE->GetXaxis()->GetBinCenter(6); 
                       double high_bin_center = PATHLENGTH_GIVEN_MODULE->GetXaxis()->GetBinCenter(8); 

                     
                       double bin_center = PATHLENGTH_GIVEN_MODULE->GetXaxis()->GetBinCenter(7); 
                       double widthbinnb = PATHLENGTH_GIVEN_MODULE->GetBinWidth(7);
                       //cout << "bin 7 has center : " << bin_center << " and that bin has width = " << widthbinnb <<endl; 
                       if(dedx_modulgeom[iclu]==1){
                           if (track_p[index_of_the_track]>10 && track_p[index_of_the_track]<45) {
                               if(dedx_isstrip[iclu]){
                                   if( (dedx_pathlength[iclu]*10) > (low_bin_center-(widthbinnb/2)) && (dedx_pathlength[iclu]*10) < (high_bin_center+(widthbinnb/2))){
                                       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10));
                                       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10));
                                       if(sclus_charge[iclu]<254){
                                           CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_NOSATURATED->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10));
                                       }
                                       if(npv < 20){
                                           CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_LOWPU->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10));
                                       }
                                       if(npv > 50){
                                           CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_HIGHPU->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10));
                                       }
                                   }
                               }
                               else{
                                   SFtest *=dEdxSF[1];
                                   if( (dedx_pathlength[iclu]*10) > (low_bin_center-(widthbinnb/2)) && (dedx_pathlength[iclu]*10) < (high_bin_center+(widthbinnb/2))){
                                       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10*norm_mult_test));
                                       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_PIXEL->Fill(SFtest*ch1/(dedx_pathlength[iclu]*10*norm_mult_test));
                                   }
                               }
                           }
                       }
                       if(writeTptHSCP){
                           float norm_mult = 265; // 247 or 265?
                           double Norm = 3.61e-06*norm_mult;
                           double scaleFactor = dEdxSF[0];
                           int lower_bd = stoi(low_bound);
                           int high_bd = stoi(high_bound);
                           if(dedx_isstrip[iclu]){
                               if (track_p[index_of_the_track]>=lower_bd && track_p[index_of_the_track]<=high_bd && track_pt[index_of_the_track] > 10) {
                                   Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk_hscp);
                                   for (int l = 0 ; l < ias_intervals; l++){
                                       if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                           Charge_Vs_Path_PU[l]->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk_hscp);
                                       }
                                   }
                               }
                           }
                           else{
                               scaleFactor *=dEdxSF[1];
                               if (track_p[index_of_the_track]>=lower_bd && track_p[index_of_the_track]<=high_bd && track_pt[index_of_the_track] > 10) {
                                   Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk_hscp);
                                   for (int l = 0 ; l < ias_intervals; l++){
                                       if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                           Charge_Vs_Path_PU[l]->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk_hscp);
                                       }
                                   }
                               }
                           }
                       }//END TPT
                    }
                    //HERE TO CHANGE TRACK TO HSCP

        
                } // end loop iclu

                //TEST nb measurements for production template
                /*
                for (int iclu=track_index_hit[hscp_track_idx[ihs]]; iclu<track_index_hit[hscp_track_idx[ihs]]+track_nhits[hscp_track_idx[ihs]]; iclu++){
                    if (dedx_subdetid[iclu]>=3) {

                    }
                    else{

                    } 
                    float ch1=dedx_charge[iclu];
                    bool clean1=true;

                    vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
                    float newcharge =0;
                    for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
                    ch1=newcharge;
                    clean1=sclus_clusclean2[iclu]; // clusterCleaning with Xtalk inversion and Saturation (September 22, 2021)
                }
                */

                if(UsePURwtHSCP){

                    // APPLY Selection on the Number of Measurements :
                    if (charge_corr3.size()>9) {
    
                        int nval20_0=0;
                        int nsat20_0=0;
                        int anval20_0=0;
                        int ansat20_0=0;
                        int nval4_0=0;
                        int nsatv4_0=0;
                        int nv = 0;
                        int ns = 0;
                        double scaleFactor = dEdxSF[0];
                        float norm_mult = 265; // 247 or 265?
                        double ias_all = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, dEdxTemplatesAll,2, 0., anval20_0, ansat20_0);
                        double ias_strip_pu_base = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
                        double is_strip_pu_base = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
          
                        double ih0_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0., nv, ns);
                        double ih_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF,  NULL,2, 0.15,  nv, ns); 
    
                        if (track_p[index_of_the_track] > 5 && track_p[index_of_the_track] < 100){
                            if(abs(track_eta[index_of_the_track])<1){
                                IhStripVsP_presk->Fill(ih_strip,track_p[index_of_the_track]);
                                IhBestStripVsP_presk->Fill(ih0_noL1,track_p[index_of_the_track]);
                                IhBestvsIas_p_5_100->Fill(ih0_noL1,ias_strip_pu_base);
                            }
                            if(abs(track_eta[index_of_the_track])<0.2){
                                IhStripVsPtight_presk->Fill(ih_strip,track_p[index_of_the_track]);
                                IhBestStripVsPtight_presk->Fill(ih0_noL1,track_p[index_of_the_track]);
                            }
           
                        }
                        if (track_p[index_of_the_track] > 10 && track_p[index_of_the_track] < 45){
                            if(abs(track_eta[index_of_the_track])<1){
                               IhBestvsIas_p_10_45->Fill(ih0_noL1,ias_strip_pu_base);
                               IasStripVsIh0noL1_p_10_45->Fill(ias_strip_pu_base,ih0_noL1);
                               IhBestvsP_10_45_central->Fill(ih0_noL1,track_p[index_of_the_track]);
                            }
                            if(abs(track_eta[index_of_the_track])<0.2){ 
                                IhBestvsP_10_45_tight->Fill(ih0_noL1,track_p[index_of_the_track]);
                            }
                        }
           
                        if(track_p[index_of_the_track] > 5 && track_p[index_of_the_track] < 100){
                            if(abs(track_eta[index_of_the_track])<1){
                                //add presk for prescale in case of p < 5 GeV
                                IhStripVsP->Fill(ih_strip,track_p[index_of_the_track]);
                                IhBestStripVsP->Fill(ih0_noL1,track_p[index_of_the_track]);
                            }
                            if(abs(track_eta[index_of_the_track])<0.2){
                                IhStripVsPtight->Fill(ih_strip,track_p[index_of_the_track]);
                                IhBestStripVsPtight->Fill(ih0_noL1,track_p[index_of_the_track]);
                            }
                            
                        }
                    
    
                        PT_compute_ias->Fill(track_pt[index_of_the_track]);
                        P_compute_ias->Fill(track_p[index_of_the_track]);
                        P_after_presel->Fill(track_p[index_of_the_track]);
                         
                        //Cut on track_p pour avoir la meme statistique que celle pour les templates
                        int lower_bd = stoi(low_bound);
                        int high_bd = stoi(high_bound);
                        if (track_p[index_of_the_track]>=lower_bd && track_p[index_of_the_track]<=high_bd) {
                            if(compute_PE){
                                //cout << "EVENT #" << jentry << " ,HSCP #" << ihs << " has track_p = " << track_p[index_of_the_track]<<endl;
                                for (int i = 0; i < ias_intervals; i++){
                                    if ( npv > ias_top_born[i] && npv <= ias_top_born[i+1]){
                                        //cout << " NPV for this event : " << npv << endl;
                                        for (int k = 0 ; k < nPE ; k++){
                                            int nval20_0 = 0,nsat20_0=0;
                                            double ias_poisson_err = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, PE_dEdxTemplatesPU[i][k],2, 0., nval20_0, nsat20_0);
                                            PE_Ias[i][k]->Fill(ias_poisson_err);
                                            //if (k%50 == 0) cout << "below IH cut : PE #"<<k<< " from PU bin #" << i+1 << " has ias = " << ias_poisson_err <<endl;
                                            
                                            if (ih0_noL1 > 3.47){
                                                PE_Ias_cutih[i][k]->Fill(ias_poisson_err);
                                                //if(k%50 == 0) cout << "above IH cut : PE #" <<k<< " from PU bin #" << i+1 << " has ias = " << ias_poisson_err <<endl;
                                            }
                                        }
                                    }
                                }
                            }
           
           
           
                            if (npv >= 16 && npv <= 18){
                               double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuLow,2, 0., nval20_0, nsat20_0);
                               Ias_when_PU_in_16_18_base->Fill(ias_strip_pu_base);    
                               Ias_when_PU_in_16_18_triple->Fill(ias_strip_pu_triple);    
                            }
                            if ( npv >= 30 && npv <= 32){
                               double ias_strip_pu_triple =getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuMedium,2, 0., nval20_0, nsat20_0);
                               Ias_when_PU_in_30_32_base->Fill(ias_strip_pu_base);
                               Ias_when_PU_in_30_32_triple->Fill(ias_strip_pu_triple);
                 
                            }
                            bool ihcut = false;
                            if (ih0_noL1 > 3.47 ) ihcut = true;
           
                            eff_ih_PU->Fill(ihcut,npv);
                            for (int l = 0 ; l < ias_intervals; l++){
                                if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                    ih_eff_denom[l] += 1;
                                    if (ih0_noL1 > 3.47) ih_eff_num[l] +=1;
           
                                }
                            }
                            double trf = 0;
                            for (int l = 0 ; l < ias_intervals; l++){
                                if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                    double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);

                                    double ias_strip_pu_triple_bis = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPUBis[l],2, 0., nval20_0, nsat20_0);
                                    Ias_1_VS_Delta_Ias_1_2->Fill(ias_strip_pu_triple_bis,(ias_strip_pu_triple_bis - ias_strip_pu_triple));
                                    double is_strip_pu_triple = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                                    //cut sur IAS 
                                    // 0.115
                                    float mass_single_nocut=-1;
                                    if(ias_strip_pu_base>0.115){
                                        if (ih0_noL1 - Cval_nol1>0) mass_single_nocut= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                        MASS_SINGLE_nocut->Fill(mass_single_nocut);
                                    }
                                    float mass_triple_nocut=-1;
                                    if(ias_strip_pu_triple>0.115){
                                        if (ih0_noL1 - Cval_nol1>0) mass_triple_nocut= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                        MASS_TRIPLE_nocut->Fill(mass_triple_nocut);
                                    }
    
                                    Ias_triple_Vs_Ih0noL1_p_10_45->Fill(ias_strip_pu_triple,ih0_noL1);
                                    Ias_simple_Vs_Ih0noL1_p_10_45->Fill(ias_strip_pu_base,ih0_noL1);
                                    //check pT > 50, et IAS entre 80-90 quantile 
                                    Ias_when_PU_triple[l]->Fill(ias_strip_pu_triple);
                                    Ias_when_PU[l]->Fill(ias_strip_pu_base);
           
                                    Is_when_PU_triple[l]->Fill(is_strip_pu_triple);
                                    Is_when_PU[l]->Fill(is_strip_pu_base);
           
                                    PU_distrib[l]->Fill(npv);
                                    sums_ias_nocut[l] += ias_strip_pu_base;
                                    sums_ias_triple_nocut[l] += ias_strip_pu_triple;
                                    Ias_all_triple_nocut->Fill(ias_strip_pu_triple);
           
                                }
                            }
                            if (npv > ias_top_born[ias_intervals]) trf = -1;
           
                            Ias_all_base_nocut->Fill(ias_strip_pu_base);
           
                            for (int l = 0 ; l < ias_intervals; l++){
                                if (ih0_noL1 > 3.47){
                                    if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                        double ias_strip_pu_triple_l = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                                        double is_strip_pu_triple_l = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
           
                                        float mass_single=-1;
                                        if(ias_strip_pu_base>0.1){
                                            IAS_triple_when_simple_above_0p115->Fill(ias_strip_pu_triple_l);
                                            if (ih0_noL1 - Cval_nol1>0) mass_single= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                            MASS_SINGLE_ihcut->Fill(mass_single);
                                            if(ias_strip_pu_triple_l < 0.1){
                                                MASS_SINGLE_ihcut_base_notriple_ias_0p1->Fill(mass_single);
                                            }
                                        }
                                        //BLINDED, WE STAY BETWEEN 80 - 90 QTL, X - 0.115
                                        

                                        float mass_triple=-1;
                                        if(ias_strip_pu_triple_l>0.1){
                                            IAS_simple_when_triple_above_0p115->Fill(ias_strip_pu_base);
                                            if (ih0_noL1 - Cval_nol1>0) mass_triple= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                            MASS_TRIPLE_ihcut->Fill(mass_triple);
                                            if(ias_strip_pu_base < 0.1){
                                                MASS_TRIPLE_ihcut_triple_nobase_ias_0p1->Fill(mass_triple);
                                            }
                                        }
                                        Ias_all_triple_cutIH->Fill(ias_strip_pu_triple_l);
                                        Ias_all_base_cutIH->Fill(ias_strip_pu_base);
                                        sums_ias_triple[l]+=ias_strip_pu_triple_l;
                                        sums_ias[l]+=ias_strip_pu_base;
           
                                        Ias_when_PU_ih_cut[l]->Fill(ias_strip_pu_base);
                                        Ias_when_PU_ih_cut_triple[l]->Fill(ias_strip_pu_triple_l);
           
                                        Is_when_PU_ih_cut[l]->Fill(is_strip_pu_base);
                                        Is_when_PU_ih_cut_triple[l]->Fill(is_strip_pu_triple_l);
                                    }
                                }
                            }
           
    
                        }//end track p cut
                         
    
                        else{
                            if (npv <= 20){
                                double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuHigh,2, 0., nval20_0, nsat20_0);
                                Ias_below_5GeV_PU_below_20_base->Fill(ias_strip_pu_base);
                                Ias_below_5GeV_PU_below_20_triple->Fill(ias_strip_pu_triple);
           
           
                            }
                            else if (npv > 20 && npv <=40){
                                double ias_strip_pu_triple =getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuMedium,2, 0., nval20_0, nsat20_0);
                                Ias_below_5GeV_PU_between_base->Fill(ias_strip_pu_base);
                                Ias_below_5GeV_PU_between_triple->Fill(ias_strip_pu_triple);
                                
           
                            }
                            else if (npv > 40){
                                double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuHigh,2, 0., nval20_0, nsat20_0); 
                                Ias_below_5GeV_PU_above_40_base->Fill(ias_strip_pu_base);
                                Ias_below_5GeV_PU_above_40_triple->Fill(ias_strip_pu_triple);
           
                            }              
           
                        }
                        if(track_pt[index_of_the_track]>50){
                            IH0_noL1->Fill(ih0_noL1); 
                            for (int l = 0 ; l < ias_intervals; l++){
                                if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                    double ias_strip_pu_triple_l = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                                    double is_strip_pu_triple_l = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                                    NPV_ih_nocut->Fill(npv,ih0_noL1);
                                    NPV_ih_strip_nocut->Fill(npv,ih_strip);

                                    if(ih_strip > 1){
                                        NPV_ih_strip_cut1->Fill(npv,ih_strip);
                                        if(ih_strip > 3.29){
                                            NPV_ih_strip_cut3p29->Fill(npv,ih_strip);
                                            if(ih_strip > 3.47){
                                                NPV_ih_strip_cutih->Fill(npv,ih_strip);
                                            }
                                        }
                                    }

                                    if(ih0_noL1 > 1){
                                        NPV_ih_cut1->Fill(npv,ih0_noL1);
                                    }
                                    if (ih0_noL1 > 3.47){
                                        float mass_single_pt_sup50=-1;
                                        NPV_ias_single_cutih->Fill(npv,ias_strip_pu_base); 
                                        NPV_ih_cutih->Fill(npv,ih0_noL1); 
                                        NPV_ias_triple_cutih->Fill(npv,ias_strip_pu_triple_l); 
    

                                        if( npv <= 10) P_LOWPU_pt50->Fill(track_p[hscp_track_idx[ihs]]);
                                        if( npv >= 60) P_HIGHPU_pt50->Fill(track_p[hscp_track_idx[ihs]]);


                                        if(npv >= 28 && npv <=29 ) P_PU_28_29_pt50->Fill(track_p[hscp_track_idx[ihs]]);
                                        if(npv >= 30 && npv <=31 ) P_PU_30_31_pt50->Fill(track_p[hscp_track_idx[ihs]]);

                                        if(npv >= 28 && npv <=30 ) P_PU_28_30_pt50->Fill(track_p[hscp_track_idx[ihs]]);
                                        if(npv >= 32 && npv <=34 ) P_PU_32_34_pt50->Fill(track_p[hscp_track_idx[ihs]]);
        
                                        if(npv >= 8 && npv <=10 ) P_PU_08_10_pt50->Fill(track_p[hscp_track_idx[ihs]]);
                                        if(npv >= 58 && npv <=60 ) P_PU_58_60_pt50->Fill(track_p[hscp_track_idx[ihs]]);

                                        IAS_single_pt50->Fill(ias_strip_pu_base);
                                        IAS_triple_pt50->Fill(ias_strip_pu_triple_l);
                                        if(ias_strip_pu_triple_l > 0.1) nb_ias_0p1_triple+=1;
                                        if(ias_strip_pu_base > 0.1) nb_ias_0p1_single+=1;
    
                                        if(ias_strip_pu_triple_l > 0.2) nb_ias_0p2_triple+=1;
                                        if(ias_strip_pu_base > 0.2) nb_ias_0p2_single+=1;
    
    
                                        if(ias_strip_pu_triple_l > 0.3){
                                            nb_ias_0p3_triple+=1;
                                            PU_ias_triple_0p3->Fill(npv);
                                        }
                                        if(ias_strip_pu_base > 0.3){
                                            nb_ias_0p3_single+=1;
                                            PU_ias_single_0p3->Fill(npv);
                                        }

                                        if(ias_strip_pu_base > 0.4){
                                            nb_ias_0p4_single+=1;
                                        }
        
                                        if(ias_strip_pu_triple_l > 0.4){
                                            nb_ias_0p4_triple+=1;
                                        }

                                        if (ias_strip_pu_base > 0.07 && ias_strip_pu_base < 0.1){
                                            P_single_pt50_ias_qtl_80_90->Fill(track_p[hscp_track_idx[ihs]]);
                                            if (ih0_noL1 - Cval_nol1>0) mass_single_pt_sup50= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                            MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->Fill(mass_single_pt_sup50);
                                            NPV_mass_cutih->Fill(npv,mass_single_pt_sup50);
                                            if(track_pt[index_of_the_track]>100) MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->Fill(mass_single_pt_sup50);
                                            if(ias_strip_pu_triple_l < 0.07){
                                                MASS_SINGLE_ihcut_base_notriple_ias_0p1_ptsupp50_qtl_80_90_ias->Fill(mass_single_pt_sup50);
                                            } 
                                        }
    
                                        float mass_triple_pt_sup50=-1;
                                        if( ias_strip_pu_triple_l >0.07 && ias_strip_pu_triple_l < 0.1){
                                            P_triple_pt50_ias_qtl_80_90->Fill(track_p[hscp_track_idx[ihs]]);
                                            if (ih0_noL1 - Cval_nol1>0) mass_triple_pt_sup50= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                            MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Fill(mass_triple_pt_sup50);
                                            if(track_pt[index_of_the_track]>100) MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->Fill(mass_triple_pt_sup50);
                                            if(ias_strip_pu_base < 0.07){
                                                MASS_TRIPLE_ihcut_triple_nobase_ias_0p1_ptsupp50_qtl_80_90_ias->Fill(mass_triple_pt_sup50);
                                            }
                                        }
                                    }
                                    if(ih0_noL1 - Cval_nol1>0){
                                        float mass_single_pt_sup50_cut3p29=-1;
                                        if (ias_strip_pu_base > 0.07 && ias_strip_pu_base < 0.1){
                                            mass_single_pt_sup50_cut3p29 = sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                                            NPV_mass_cut3p29->Fill(npv,mass_single_pt_sup50_cut3p29);
                                        }

                                        NPV_ih_cut3p29->Fill(npv,ih0_noL1);
                                        if(ias_strip_pu_triple_l > 0.1) nb_ias_0p1_triple_cut3p29+=1;
                                        if(ias_strip_pu_base > 0.1) nb_ias_0p1_single_cut3p29+=1;
    
                                        if(ias_strip_pu_triple_l > 0.2) nb_ias_0p2_triple_cut3p29+=1;
                                        if(ias_strip_pu_base > 0.2) nb_ias_0p2_single_cut3p29+=1;
    
                                        if(ias_strip_pu_triple_l > 0.3) nb_ias_0p3_triple_cut3p29+=1;
                                        if(ias_strip_pu_base > 0.3) nb_ias_0p3_single_cut3p29+=1;

                                        if(ias_strip_pu_triple_l > 0.4) nb_ias_0p4_triple_cut3p29+=1;
                                        if(ias_strip_pu_base > 0.4) nb_ias_0p4_single_cut3p29+=1;

                                    }
                                     
                                    if(ias_strip_pu_triple_l > 0.1) nb_ias_0p1_triple_nocut+=1;
                                    if(ias_strip_pu_base > 0.1) nb_ias_0p1_single_nocut+=1;
                                    if(ias_strip_pu_triple_l > 0.2) nb_ias_0p2_triple_nocut+=1;
                                    if(ias_strip_pu_base > 0.2) nb_ias_0p2_single_nocut+=1;
    
                                    if(ias_strip_pu_triple_l > 0.3) nb_ias_0p3_triple_nocut+=1;
                                    if(ias_strip_pu_base > 0.3) nb_ias_0p3_single_nocut+=1;

                                    if(ias_strip_pu_triple_l > 0.4) nb_ias_0p4_triple_nocut+=1;
                                    if(ias_strip_pu_base > 0.4) nb_ias_0p4_single_nocut+=1;


                                    if (track_pt[index_of_the_track] <= 60){
                                        if(ih0_noL1 > 3.47){
                                            IAS_TRIPLE_PT_50_60_selection->Fill(ias_strip_pu_triple_l);
                                            IAS_SINGLE_PT_50_60_selection->Fill(ias_strip_pu_base); 
                                        }
                                    }
       

                                }

                            }
                        } //END IF PT > 50
                        IAS_VS_PT->Fill(ias_all,track_pt[index_of_the_track]);
                        //cout << "filling ias_vs_pt with ias : " << ias_all << " for a pt = " << track_pt[index_of_the_track] << endl;
                        
                        if(track_pt[index_of_the_track]>60 && track_pt[index_of_the_track]<=80){
                           ias_bin_pt[0]->Fill(ias_all);
                           //cout << "filling ias for pt bin 60-80 with ias : " << ias_all << " and pt = " << track_pt[index_of_the_track] << endl;
                        }
                        if(track_pt[index_of_the_track] > 80 && track_pt[index_of_the_track] <= 100){
                           ias_bin_pt[1]->Fill(ias_all);
                           //cout << "filling ias for pt bin 80-100 with ias : " << ias_all << " and pt = " << track_pt[index_of_the_track] << endl;
                        }
                        if(track_pt[index_of_the_track] > 100 && track_pt[index_of_the_track] <= 120){
                           ias_bin_pt[2]->Fill(ias_all);
                        }
                        if(track_pt[index_of_the_track] > 120 && track_pt[index_of_the_track] <= 140){
                           ias_bin_pt[3]->Fill(ias_all); 
                        }
                        if(track_pt[index_of_the_track] > 140 && track_pt[index_of_the_track] <= 160){
                           ias_bin_pt[4]->Fill(ias_all);
                        }
                        if(track_pt[index_of_the_track] > 160 && track_pt[index_of_the_track] <= 180){
                           ias_bin_pt[5]->Fill(ias_all);
                        }
                        if(track_pt[index_of_the_track] > 180 && track_pt[index_of_the_track] <= 200){
                           ias_bin_pt[6]->Fill(ias_all);
                        }
                        if(track_pt[index_of_the_track] > 200 && track_pt[index_of_the_track] <= 220){
                           ias_bin_pt[7]->Fill(ias_all);
                        }
                        if(track_pt[index_of_the_track] > 220 && track_pt[index_of_the_track] <= 240){
                           ias_bin_pt[8]->Fill(ias_all);
                        }
                        if(track_pt[index_of_the_track] > 240 && track_pt[index_of_the_track] <= 260){
                           ias_bin_pt[9]->Fill(ias_all);
                        }
                    }
                } // end PU reweighting

                if (charge_corr.size()>6) {  // the one for which we don't remove any hits
                  int nval1=0;
                  int nsatv1=0;
                  int nval1_0=0;
                  int nsatv1_0=0;
                  // Ih 15% drop of lowest values
                  double ih_LD        = getdEdX(charge_corr,    pathlength,  subdetId,  moduleGeometry,  bool_cleaning,  mustBeInside,  dEdxSF, NULL,2, 0.15, nval1, nsatv1);
                  dEdXVsRun->Fill(       runNumber,ih_LD);
                  NmeasVsRun->Fill(runNumber,nval1);
                  if (nsatv1>9) nsatv1=9;
                  // Ih no drop
                  double ih0_cor     = getdEdX(charge_corr,  pathlength,  subdetId,  moduleGeometry,  bool_cleaning,  mustBeInside,  dEdxSF, NULL,2, 0., nval1_0, nsatv1_0);
                  dEdX0VsRun->Fill(runNumber,ih0_cor);
                  int nval1HD=0;
                  int nsatv1HD=0;
                  // Ih  15% drop of highest values
                  double ih_corHiDrop = getdEdX(charge_corr,  pathlength,  subdetId,  moduleGeometry,  bool_cleaning,  mustBeInside,  dEdxSF, NULL,2, 0., 0.15, nval1HD, nsatv1HD);
                  dEdXHiDropVsRun->Fill(     runNumber,ih_corHiDrop);

                  int nval2=0;
                  int nval2_0=0;
                  int nsatv2=0;
                  int nsatv2_0=0;
                  int nval2HD=0;
                  int nsatv2HD=0;
                  // Ih 15% drop of lowest values
                  double ih_LDstrip   = getdEdX(charge_corr1,   pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0.15, nval2, nsatv2);
                  dEdXstripVsRun->Fill(  runNumber,ih_LDstrip);
                  NmeasStrVsRun->Fill(runNumber,nval2);
                  if (nsatv2>9) nsatv2=9;
                  NsatVsRun->Fill(runNumber,nsatv1);
                  NsatStrVsRun->Fill(runNumber,nsatv2);
                  // Ih no drop
                  double ih0_strip   = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0., nval2_0, nsatv2_0);
                  dEdX0stripVsRun->Fill(runNumber,ih0_strip);


                  // Ih  15% drop of highest values
                  double ih_stripHiDrop = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0., 0.15, nval2HD, nsatv2HD);
                  dEdXstripHiDropVsRun->Fill(runNumber,ih_stripHiDrop );

                  int nval3=0;
                  int nval3_0=0;
                  int nsatv3=0;
                  int nsatv3_0=0;
                  int nval3HD=0;
                  int nsatv3HD=0;
                  // Ih 15% drop of lowest values
                  double ih_LDpix     = getdEdX(charge_corr2,   pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0.15, nval3, nsatv3);
                  dEdXpixVsRun->Fill(    runNumber,ih_LDpix);
                  NmeasPixVsRun->Fill(runNumber,nval3);
                  if (nsatv3>9) nsatv3=9;
                  NsatPixVsRun->Fill(runNumber,nsatv3);
                  // Ih no drop
                   double ih0_pix     = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0., nval3_0, nsatv3_0);
                  dEdX0pixVsRun->Fill(runNumber,ih0_pix);



                  // Ih  15% drop of highest values
                  double ih_pixHiDrop   = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0., 0.15, nval3HD, nsatv3HD);
                  dEdXpixHiDropVsRun->Fill(  runNumber,ih_pixHiDrop);


                  int nval6=0;
                  int nval6_0=0;
                  int nsatv6=0;
                  int nsatv6_0=0;
                  // Ih 15% drop of lowest values
                  double ih_LDnoL1pix = getdEdX(charge_corr5, pathlength5, subdetId5, moduleGeometry5, bool_cleaning5, mustBeInside5, dEdxSF, NULL,2, 0.15, nval6, nsatv6);
                  dEdXNoL1pixVsRun->Fill(runNumber,ih_LDnoL1pix);
                  // Ih no drop
                  double ih0_noL1pix = getdEdX(charge_corr5, pathlength5, subdetId5, moduleGeometry5, bool_cleaning5, mustBeInside5, dEdxSF, NULL,2, 0., nval6_0, nsatv6_0);
                  dEdX0NoL1pixVsRun->Fill(runNumber,ih0_noL1pix);
                 

                  float mass_ih0noL1=-1;
                  float mass_ihHDnoL1=-1;
                  float mass_ih0noL1_2=-1;
                  float mass_ih0noL1_3=-1;


                  // Ias
                  int nval20_0=0;
                  int nsat20_0=0;
                  double ias_all =-1;
                  if (boolDeDxTemp)  ias_all=getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
                  iasAllVsRun->Fill(   runNumber,ias_all);

                  double ias_strip =-1;
                  if (boolDeDxTemp)  ias_strip=getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
                  double ias_pix =-1;
                  if (boolDeDxTemp)  ias_pix=getdEdX(charge_corr5, pathlength5, subdetId5, moduleGeometry5, bool_cleaning5, mustBeInside5, dEdxSF, dEdxTemplatesNoL1,2, 0., nval20_0, nsat20_0);


                  nval20_0=0;
                  nsat20_0=0;
                  double ias_noL1 = -1;
                  if (boolDeDxTemp) ias_noL1= getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, dEdxTemplatesNoL1,2, 0., nval20_0, nsat20_0);
                  iasNoL1VsRun->Fill(   runNumber,ias_noL1);

                  float mass_strip=-1;
                  if (ih_LDstrip - Cval_ldstrip>0) mass_strip=sqrt((ih_LDstrip - Cval_ldstrip)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_ldstrip);

                  float mass_ih=-1;
                  if (ih_LD - Cval_ld>0) mass_ih=sqrt((ih_LD - Cval_ld)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_ld);

                  float mass_ih0=-1;
                  if (ih0_cor - Cval_all>0) mass_ih0= sqrt((ih0_cor - Cval_all)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_all);


                  float mass_0strip=-1;
                  if (ih0_strip - Cval_strip>0) mass_0strip= sqrt((ih0_strip - Cval_strip)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_strip);


                  float mass_pixel=-1;
                  if (ih_LDpix - Cval_pix>0) mass_pixel= sqrt((ih_LDpix - Cval_pix)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_pix);

                  float mass_pixnL1=-1;
                  if (ih0_noL1pix - Cval_pixnol1>0) mass_pixnL1= sqrt((ih0_noL1pix - Cval_pixnol1)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_pixnol1);


                  // additional cut on the number of measurement points
                  if (charge_corr3.size() >9) {

                     int nval4=0;
                     int nval4_0=0;
                     int nsatv4=0;
                     int nsatv4_0=0;

                     /*
                     for (int iclu=track_index_hit[hscp_track_idx[ihs]]; iclu<track_index_hit[hscp_track_idx[ihs]]+track_nhits[hscp_track_idx[ihs]]; iclu++) {
                         float norm_mult = 265; // 247 or 265?
                         double Norm = 3.61e-06*norm_mult;
                         double scaleFactor = dEdxSF[0];
                         if(dedx_isstrip[iclu]){
                             if (track_p[index_of_the_track]>10 && track_p[index_of_the_track]<45) {
                                 Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10));
                                 for (int l = 0 ; l < ias_intervals; l++){
                                     if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                         Charge_Vs_Path_PU[l]->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10));
                                     }
                                 }
                             }
                         }
                         else{
                             scaleFactor *=dEdxSF[1];
                             if (track_p[index_of_the_track]>10 && track_p[index_of_the_track]<45) {
                                 for (int l = 0 ; l < ias_intervals; l++){
                                     if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                                         Charge_Vs_Path_PU[l]->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult));
                                     }
                                 }
                             }
                         }
                     }
                     */

                     // Ih 15% drop of lowest values
                     double ih_LDnoL1    = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0.15, nval4, nsatv4);
                     dEdXNoL1VsRun->Fill(   runNumber,ih_LDnoL1);
 
                     // Ih no drop
                     double ih0_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0., nval4_0, nsatv4_0);

                     dEdX0NoL1VsRun->Fill(runNumber,ih0_noL1);  //  <==== our best estimate for Ih

                     Nmeas0VsRun->Fill(runNumber,nval4_0);    // no L1
                     NmeasPix0VsRun->Fill(runNumber,nval6_0); // no L1
                     NmeasStr0VsRun->Fill(runNumber,nval2_0);  

                     // Ih  15% drop of highest values
                     nval1HD=0;
                     nsatv1HD=0;
                     double ih_corHiDropNoL1 = getdEdX(charge_corr3,  pathlength3,  subdetId3,  moduleGeometry3,  bool_cleaning3,  mustBeInside3,  dEdxSF, NULL,2, 0., 0.15, nval1HD, nsatv1HD);
                     dEdXHiDropNoL1VsRun->Fill(     runNumber,ih_corHiDropNoL1);

                     if (ih0_noL1 - Cval_nol1>0) mass_ih0noL1= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                     if (ih0_noL1 - Cval_nol1_2>0) mass_ih0noL1_2= sqrt((ih0_noL1 - Cval_nol1_2)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1_2);

                     if (mass_ih0noL1_2>HighestMass) HighestMass=mass_ih0noL1_2;

                     float val_nol1_3 = Cval_nol1_3 + (0.92*0.92*Kval_nol1_3)/(6.5*6.5) +Nval_nol1_3*log(6.5/0.92);
                     if (ih0_noL1 - val_nol1_3 >0 && computeSpecial) 
                                          mass_ih0noL1_3= getMassSpecial(ih0_noL1,track_p[hscp_track_idx[ihs]], Kval_nol1_3, Cval_nol1_3,Nval_nol1_3);
                     if (ih_corHiDropNoL1 - Cval_hdnol1>0) mass_ihHDnoL1= sqrt((ih_corHiDropNoL1 - Cval_hdnol1)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_hdnol1);

                     float fmip_strip4 =   FMIP(charge_corr1, pathlength1,dEdxSF[0],4);
                     float fmip_strip3p5 = FMIP(charge_corr1, pathlength1,dEdxSF[0], 3.5);
                     float fmip_strip3p2 = FMIP(charge_corr1, pathlength1,dEdxSF[0], 3.2);

                     FMIP4VsRun->Fill(   runNumber,fmip_strip4);
                     FMIP3p5VsRun->Fill( runNumber,fmip_strip3p5);
                     FMIP3p2VsRun->Fill( runNumber,fmip_strip3p2);

                     if (nsatv4_0>9) nsatv4_0=9;
                     if (nsatv2_0>9) nsatv2_0=9;
                     if (nsatv6_0>9) nsatv6_0=9;
                     Nsat0VsRun->Fill(runNumber,nsatv4_0);
                     NsatPix0VsRun->Fill(runNumber,nsatv6_0);
                     NsatStr0VsRun->Fill(runNumber,nsatv2_0);


                     if ((!blind_data) || (mass_ih0noL1_2<500)) {

                       HSCP_pt->Fill(track_pt[index_of_the_track]);
                       HSCP_eta->Fill(track_eta[index_of_the_track]);
                       HSCP_iso_eop->Fill(eop);
                       if (boolILumi) lumiVsRun->Fill(runNumber,InstLumi);
                       HSCP_dEdX->Fill(ih_LD);
                       HSCP_dEdX0->Fill(ih0_cor);
                       HSCP_dEdXHiDrop->Fill(ih_corHiDrop);
                       HSCP_dEdXstrip->Fill(ih_LDstrip);
                       HSCP_dEdX0strip->Fill(ih0_strip);
                       HSCP_dEdXHiDropNoL1->Fill(ih_corHiDropNoL1);
                       // comparison
                       HSCP_dEdXstripHiDrop->Fill(ih_stripHiDrop);
                       HSCP_dEdX0pix->Fill(ih0_pix);
                       HSCP_dEdXpix->Fill(ih_LDpix);
                       HSCP_dEdXpixHiDrop->Fill(ih_pixHiDrop);
                       HSCP_dEdXpixVsstrip->Fill(ih_LDstrip,ih_LDpix);
                       HSCP_dEdXstripVsall->Fill(ih_LD,ih_LDstrip);
                       HSCP_dEdXpixVsall->Fill(ih_LD,ih_LDpix);
  
                       // our best Ih estimate :
                       HSCP_dEdX0NoL1->Fill(ih0_noL1);
                       HSCP_dEdX0pixVsstrip->Fill(ih0_strip,ih0_noL1pix);




                       HSCP_iasall->Fill(ias_all);
                       HSCP_iasnol1->Fill(ias_noL1);
                       HSCP_iasstrip->Fill(ias_strip);
                       HSCP_iaspix->Fill(ias_pix);
                       if (boolProbQ) {
                        HSCP_probQ->Fill(track_probQ[index_of_the_track]);
                        HSCP_probQNoL1->Fill(track_probQNoL1[index_of_the_track]);
                        HSCP_probXY->Fill(track_probXY[index_of_the_track]);
                        HSCP_probXYNoL1->Fill(track_probXYNoL1[index_of_the_track]);
                        probQVsRun->Fill(   runNumber,track_probQ[index_of_the_track]);
                        probQNoL1VsRun->Fill(   runNumber,track_probQNoL1[index_of_the_track]);
                        probXYVsRun->Fill(   runNumber,track_probXY[index_of_the_track]);
                        probXYNoL1VsRun->Fill(   runNumber,track_probXYNoL1[index_of_the_track]);
                        probQVsIas->Fill(ias_all,track_probQ[index_of_the_track]);
                       }

                       HSCP_nstrip->Fill(nstip_);
                       HSCP_npix->Fill(npix_);
                       if (nstip_>0) HSCP_nratio->Fill(npix_*1./nstip_);
                       HSCP_nmstrip->Fill(charge_corr1.size());
                       HSCP_nmpix->Fill(charge_corr5.size());
                       if (charge_corr1.size()>0) HSCP_nmratio->Fill(charge_corr5.size()*1./charge_corr1.size());

                       HSCP_MassIh0noL1->Fill(mass_ih0noL1);
                       HSCP_MassIh0noL1_2->Fill(mass_ih0noL1_2);
                       HSCP_MassIh0noL1_3->Fill(mass_ih0noL1_3);

                       if (ih0_noL1> Cval_nol1) {
                       HSCP_MassIh0noL1_11->Fill(mass_ih0noL1);
                       HSCP_MassIh0noL1_12->Fill(mass_ih0noL1_2);
                       HSCP_MassIh0noL1_13->Fill(mass_ih0noL1_3);
                       MassNoL1VsRun->Fill(runNumber,mass_ih0noL1_2);
                       }
                       if (ih0_noL1>3.27 + 0.21) {
                       HSCP_MassIh0noL1_1s1->Fill(mass_ih0noL1);
                       HSCP_MassIh0noL1_1s2->Fill(mass_ih0noL1_2);
                       HSCP_MassIh0noL1_1s3->Fill(mass_ih0noL1_3);
                         if (ih0_noL1>3.27 + 2*0.21) {
                         HSCP_MassIh0noL1_2s1->Fill(mass_ih0noL1);
                         HSCP_MassIh0noL1_2s2->Fill(mass_ih0noL1_2);
                         HSCP_MassIh0noL1_2s3->Fill(mass_ih0noL1_3);
                           if (ih0_noL1>3.27 + 3*0.21) {
                           HSCP_MassIh0noL1_3s1->Fill(mass_ih0noL1);
                           HSCP_MassIh0noL1_3s2->Fill(mass_ih0noL1_2);
                           HSCP_MassIh0noL1_3s3->Fill(mass_ih0noL1_3);
                           }
                         }
                       }

                       HSCP_MassIhHDnoL1->Fill(mass_ihHDnoL1);
                       HSCP2d_MassIh0noL1->Fill(track_p[index_of_the_track],mass_ih0noL1);
                       HSCP2d_MassIhHDnoL1->Fill(track_p[index_of_the_track],mass_ihHDnoL1);
                       HSCP_MassIhstrip->Fill(mass_strip);
                       HSCP_MassIh->Fill(mass_ih);
                       HSCP_MassIh0->Fill(mass_ih0);
                       HSCP_MassIh0strip->Fill(mass_0strip);
                       HSCP2d_MassIhstrip->Fill(track_p[index_of_the_track],mass_strip);
                       HSCP2d_MassIh->Fill(track_p[index_of_the_track],mass_ih);
                       HSCP2d_MassIh0->Fill(track_p[index_of_the_track],mass_ih0);
                       HSCP2d_MassIh0strip->Fill(track_p[index_of_the_track],mass_0strip);

                       HSCP2d_Mass_pix_strip15->Fill(mass_strip,mass_pixel);
                       HSCP2d_Mass_pix_strip0->Fill(mass_0strip,mass_pixnL1);
                       HSCP2d_Mass_pix_strip->Fill(mass_0strip,mass_pixnL1);
                       if (mass_0strip>0 && mass_pixnL1>0 ) {
                        HSCP_MassDiff_pix_strip0->Fill(mass_0strip - mass_pixnL1);
                        HSCP_MassResol_pix_strip0->Fill((mass_0strip - mass_pixnL1)/mass_0strip);
                       }
                       if (mass_strip>0 && mass_pixel>0) {
                        HSCP_MassDiff_pix_strip15->Fill(mass_strip - mass_pixel);
                        HSCP_MassResol_pix_strip15->Fill((mass_strip - mass_pixel)/mass_strip);
                       }

                       MassStripVsRun->Fill(runNumber,mass_strip);

                       FMIP4VsEta->Fill(track_eta[index_of_the_track],fmip_strip4);

                       HSCP_FMIP4->Fill(fmip_strip4);
                       HSCP_FMIP3p5->Fill(fmip_strip3p5);
                       HSCP_FMIP3p2->Fill(fmip_strip3p2);
                     }




/*
                     if (boolILumi) {
                      dEdXVsIL->Fill(InstLumi,ih_LD);
                      dEdXpixVsIL->Fill(InstLumi,ih_pix);
                      dEdXstripVsIL->Fill(InstLumi,ih_strip);
                     }
*/

                  } // end cut on #meas charge_corr3.size


                  if (runNumber==305186) {
                        R1_StdEdXVsEvent->Fill(event,ih_LDstrip);
                        if (boolILumi) {
                         R1_StdEdXVsLumi->Fill(InstLumi,ih_LDstrip);
                         R1_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R1_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R1_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R1_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
                  }
                  else if (runNumber==305188) {
                        R2_StdEdXVsEvent->Fill(event,ih_LDstrip);
                        if (boolILumi) {
                         R2_StdEdXVsLumi->Fill(InstLumi,ih_LDstrip);
                         R2_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R2_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R2_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R2_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
                  }
                  else if (runNumber==305204) {
                        R3_StdEdXVsEvent->Fill(event,ih_LDstrip);
                        if (boolILumi) {
                         R3_StdEdXVsLumi->Fill(InstLumi,ih_LDstrip);
                         R3_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R3_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R3_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R3_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
		  }
                    
                  // Test with a power 4 in place of 2 in the formula for Ih, 15% drop of the lowest values
                  double ih4_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,4, 0.15,  nval1, nsatv1);
                  double ih4_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,4, 0.15, nval2, nsatv2);
                  double ih4_pix = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,4, 0.15, nval3, nsatv3);

                  double ih40_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,4, 0., nval1_0, nsatv1_0);
                  double ih40_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,4, 0., nval2_0, nsatv2_0);
                  double ih40_pix = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,4, 0.,nval3_0, nsatv3_0);

                  dEdX4VsRun->Fill(runNumber,ih4_cor);
                  dEdX4pixVsRun->Fill(runNumber,ih4_pix);
                  dEdX4stripVsRun->Fill(runNumber,ih4_strip);
                  dEdX40VsRun->Fill(runNumber,ih40_cor);
                  dEdX40pixVsRun->Fill(runNumber,ih40_pix);
                  dEdX40stripVsRun->Fill(runNumber,ih40_strip);

                  if (hscp_muon_idx[ihs]>-1 && 
                       charge_corr3.size() >=6) {
                       // un muon existe pour le candidat HSCP
             // no cut in TOF yet : false if tof->nDof() < 8
             //     &&  (dttof->nDof() < 6 || csctof->nDof() < 6)
             // no cut in TOF yet : false if tof->inverseBetaErr() > GlobalMaxTOFErr
                       bool tof_sel= true;
                       if (muon_comb_tofndof[hscp_muon_idx[ihs]] <8 &&
                          (muon_dt_tofndof[hscp_muon_idx[ihs]]<6 || muon_csc_tofndof[hscp_muon_idx[ihs]]<6)) tof_sel = false;
                       if (muon_comb_inversebetaerr[hscp_muon_idx[ihs]]>0.15 ) tof_sel = false;
                       if (fabs(muon_comb_inversebeta[hscp_muon_idx[ihs]]-1)>=50) tof_sel = false;  // not in the official selection in 
                       if (tof_sel) {
                            // https://github.com/enibigir/SUSYBSMAnalysis-HSCP/blob/dev/Analyzer/plugins/Analyzer.cc  ???
/*
 // what I did before
                       if (muon_comb_tofndof[hscp_muon_idx[ihs]]>=8
                             && (muon_dt_tofndof[hscp_muon_idx[ihs]]>=6 || muon_csc_tofndof[hscp_muon_idx[ihs]]>=6) 
                             && muon_comb_inversebetaerr[hscp_muon_idx[ihs]]<=0.15 
                             && fabs(muon_comb_inversebeta[hscp_muon_idx[ihs]]-1)<50) {
*/

                         if (muon_comb_inversebeta[hscp_muon_idx[ihs]]!=0.) {
                          float beta_for_mass = 1./muon_comb_inversebeta[hscp_muon_idx[ihs]];
                          float gamma_for_mass = 1/sqrt(1-beta_for_mass*beta_for_mass);
                          float mass_tof = track_p[index_of_the_track]/(beta_for_mass*gamma_for_mass) ;
//                          if (!blind_data) {
                          if ((!blind_data) || (mass_ih0noL1_2<500)) {
                              HSCP_MassTOF->Fill(mass_tof); 
                              HSCP2d_MassTOFvsIh->Fill(mass_tof,mass_ih0noL1); 
                          }
                         }

                         invBVsRun->Fill(runNumber,muon_comb_inversebeta[hscp_muon_idx[ihs]]);
                         errinvBVsRun->Fill(runNumber,muon_comb_inversebetaerr[hscp_muon_idx[ihs]]);
                         timeVsRun->Fill(runNumber,muon_comb_vertextime[hscp_muon_idx[ihs]]);
                         if ((!blind_data) || (mass_ih0noL1_2<500)) {
                           HSCP_invB->Fill(muon_comb_inversebeta[hscp_muon_idx[ihs]]);
                           HSCP_errinvB->Fill(muon_comb_inversebetaerr[hscp_muon_idx[ihs]]);
                           HSCP_time->Fill(muon_comb_vertextime[hscp_muon_idx[ihs]]);
                           if (muon_dt_tofndof[hscp_muon_idx[ihs]]>=6) {
                            invBDTVsRun->Fill(runNumber,muon_dt_inversebeta[hscp_muon_idx[ihs]]);
                            HSCP_invBDT->Fill(muon_dt_inversebeta[hscp_muon_idx[ihs]]);
                           }
                           if (muon_csc_tofndof[hscp_muon_idx[ihs]]>=6) {
                            invBCSCVsRun->Fill(runNumber,muon_csc_inversebeta[hscp_muon_idx[ihs]]);
                            HSCP_invBCSC->Fill(muon_csc_inversebeta[hscp_muon_idx[ihs]]);
                           }
                         }
                       }
                       if (year==2016) {
                        bool tof_sel= true;
                        if (muon_newcomb_tofndof[hscp_muon_idx[ihs]] <8 &&
                          (muon_newdt_tofndof[hscp_muon_idx[ihs]]<6 || muon_newcsc_tofndof[hscp_muon_idx[ihs]]<6)) tof_sel = false;
                        if (muon_newcomb_inversebetaerr[hscp_muon_idx[ihs]]>0.15 ) tof_sel = false;
                        if (fabs(muon_newcomb_inversebeta[hscp_muon_idx[ihs]]-1)>=50) tof_sel = false;  // not in the official selection in 
                        if (tof_sel) {

                         if ((!blind_data) || (mass_ih0noL1_2<500)) {
                          invBnewVsRun->Fill(runNumber,muon_newcomb_inversebeta[hscp_muon_idx[ihs]]);
                          if (muon_newdt_tofndof[hscp_muon_idx[ihs]]>=6) {
                           invBnewDTVsRun->Fill(runNumber,muon_newdt_inversebeta[hscp_muon_idx[ihs]]);
                          }
                          if (muon_newcsc_tofndof[hscp_muon_idx[ihs]]>=6) {
                           invBnewCSCVsRun->Fill(runNumber,muon_newcsc_inversebeta[hscp_muon_idx[ihs]]);
                          }
                         }
                        }
                       }
                         
                 } // muon
               } // #measurements
             } // selection
          } // index
      } // hscp
      if (HighestMass>0) {
        bg_test_event_Mass->Fill(HighestMass);   // fill the highest mass per event (in case there are multiple candidates)
        bool noprefiring = true;
        bool nohem = true;
        for (int ijet=0; ijet<njets; ijet++) {
           if (jet_pt[ijet]>100 && fabs(jet_eta[ijet])>2.25 && fabs(jet_eta[ijet])<3) noprefiring=false;
           if (jet_pt[ijet]>20 &&  -3<jet_eta[ijet]  && jet_eta[ijet]<-1.3 && -1.57<jet_phi[ijet] && jet_phi[ijet]<-0.87) nohem=false;
        }
        if (noprefiring) bg_test_prefiring_Mass->Fill(HighestMass);  
        if (nohem) bg_test_nohem_Mass->Fill(HighestMass);  
      }

      int ntracks1=0;
      int ntracks20=0;
      //std::cout << "There are " << ntracks << " in this event, looping on them" << std::endl;
      bool npv_presel = true;
      for (int itr=0; itr<ntracks; itr++) {
         //cout << " debug loop tracks "  << itr << endl;
         int presk= 1;
         if (year!=2016) presk=track_prescale[itr];
         int index_of_the_track=itr;
         bool selection=true;
 
         //PT cut a comparer 

         //if (track_pt[index_of_the_track]<55) selection=false;

         if (abs(track_eta[index_of_the_track])>1) selection=false; // Cutting tighter on eta means our templates will not be filled between bins 5-14


         //Garder les coupures de qualite

         if (track_nvalidhits[index_of_the_track]<10) selection=false;
         if (track_npixhits[index_of_the_track]<2) selection=false;
         if (track_validfraction[index_of_the_track]<0.8) selection=false;
         if (track_missing[index_of_the_track]>99999) selection=false;
         if (track_validlast[index_of_the_track]<-99999) selection=false;

             
         bool is_high_qual =  (track_qual[index_of_the_track] & (1 << TrackQuality::highPurity)) >> TrackQuality::highPurity ;
         if (!is_high_qual) selection=false;
         if (track_chi2[index_of_the_track]>5) selection=false;

         // which selection for dz and dyx w/r to the primary vertex ? here 1st PV to pass the selection
         if (abs(track_dxy[index_of_the_track])>0.02) selection=false;
         if (abs(track_dz[index_of_the_track])>0.1) selection=false;

         // Iso cut autour du cone -> track_only sur la relecture mais pas lors de l'Ãcriture
         // refaire lors de l'Ãcriture + lecture 
         //
         //norm + superposer les Ias_when_PU_xx_xx
         //echelle normale entre 0 - 0.1 
         //distrib 0.1 a 1 et first bin doit etre egal dans les 3 plots 
         // ajouter cut IH loose et ressortir les meme plots 
         // distributions IAS en fonction de eta slides 11 meme plots 

         // no cut on the isolation

         // cut on sigma pT/pT for signal  :

         //if (track_pterr[index_of_the_track]/track_pt[index_of_the_track]>0.25) selection=false;

         if (!selection) continue;
         NTRK->Fill(2);
         //Cut same as Analysis
         //p cut -> [5,45] post MIP but not in SR (change cut to 30 or 40 GeV jeu entre B et 
         //
         if (npv_presel){
             NPV_presel->Fill(npv);
             npv_presel = false;
         }
         Htrackpt->Fill(track_pt[itr],presk);

         if (track_pt[itr]>1) ntracks1+=presk;
         if (track_pt[itr]>20) ntracks20+=presk;

      
         std::vector <float> charge_corr;
         std::vector <float> pathlength;
         std::vector <int> subdetId;
         std::vector <UInt_t> detId;
         std::vector <int> moduleGeometry;
         std::vector <bool> bool_cleaning;
         std::vector <bool> mustBeInside;

         std::vector <float> charge_corr1;
         std::vector <float> pathlength1;
         std::vector <int> subdetId1;
         std::vector <UInt_t> detId1;
         std::vector <int> moduleGeometry1;
         std::vector <bool> bool_cleaning1;
         std::vector <bool> mustBeInside1;

         std::vector <float> charge_corr2;
         std::vector <float> pathlength2;
         std::vector <int> subdetId2;
         std::vector <int> moduleGeometry2;
         std::vector <bool> bool_cleaning2;
         std::vector <bool> mustBeInside2;

         std::vector <float> charge_corr3;
         std::vector <float> pathlength3;
         std::vector <int> subdetId3;
         std::vector <UInt_t> detId3;
         std::vector <int> moduleGeometry3;
         std::vector <bool> bool_cleaning3;
         std::vector <bool> mustBeInside3;

         std::vector <float> charge_corr4;
         std::vector <float> pathlength4;
         std::vector <int> subdetId4;
         std::vector <UInt_t> detId4;
         std::vector <int> moduleGeometry4;
         std::vector <bool> bool_cleaning4;
         std::vector <bool> mustBeInside4;


         int nsatclu=0;
         N_CLU_TRK->Fill(track_index_hit[index_of_the_track]+track_nhits[itr]);
         for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
               //cout << "running over all hits from given track" << std::endl

               if(abs(track_eta[itr]) > 0 && abs(track_eta[itr]) < 0.1){
                   PATHLENGHT_BIN_ETA_0_01->Fill(dedx_pathlength[iclu]*10); // *10 from mm to cm
               }
               float ch1=dedx_charge[iclu];
               bool clean1=true;
               if (dedx_subdetid[iclu]>=3) {
                        // strip
                        // without any correction :
//                        ch1=sclus_charge[iclu];  // charge without any correction
//                        clean1=sclus_clusclean[iclu];  
//                      // with Saturation only (but no Xtalk inversion) :
                          float check_charge=0;
                          vector<int> Quncor;
                          for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++) {
                            check_charge+=strip_ampl[istrip];
                            Quncor.push_back(strip_ampl[istrip]);
                          }
                          float deltaDif=check_charge-ch1;
                          if (deltaDif<0) deltaDif*=-1;
                          if (deltaDif>0.001) std::cout << "difference dans le calcul de la charge " << ch1 << " " << check_charge << " --> probleme acces ampl ???? " << std::endl; 
                          vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
                          float newcharge =0;
                          for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
                          ch1=newcharge;
                          clean1=sclus_clusclean2[iclu]; // clusterCleaning with Xtalk inversion and Saturation (September 22, 2021) 

                          
                          // for the record : there is a possibility to recompute the cluster cleaning with the Saturation only : clusterCleaning(Qcor, 1, &exitCode)
                          // but it is not what we want here
                          if (dedx_isstrip[iclu] && clean1 && dedx_insideTkMod[iclu] && (sclus_sat254[iclu] || sclus_sat255[iclu])) nsatclu++;
               }
               else {
                        // pixel
                          // float scaling =GetSFPixel(dedx_subdetid[iclu], dedx_detid[iclu], year, abs(track_eta[hscp_track_idx[ihs]]), runNumber);
                          if (PixelCorr2Apply) {
                            float scaling =GetSFPixelTamas(dedx_subdetid[iclu], dedx_detid[iclu], year, runNumber);
                            ch1*=scaling;
                          }
               }
               if (clean1 && dedx_insideTkMod[iclu]) {
                 // fill both Strip and Pixel
                 charge_corr.push_back(ch1);
                 pathlength.push_back(dedx_pathlength[iclu]);
                 subdetId.push_back(dedx_subdetid[iclu]);
                 detId.push_back(dedx_detid[iclu]);
                 moduleGeometry.push_back(dedx_modulgeom[iclu]);
                 mustBeInside.push_back(dedx_insideTkMod[iclu]);
                 bool_cleaning.push_back(clean1);

                 if (dedx_isstrip[iclu]) {
                 // fill only Strip
                     charge_corr1.push_back(ch1);
                     pathlength1.push_back(dedx_pathlength[iclu]);
                     subdetId1.push_back(dedx_subdetid[iclu]);
                     detId1.push_back(dedx_detid[iclu]);
                     moduleGeometry1.push_back(dedx_modulgeom[iclu]);
                     mustBeInside1.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning1.push_back(clean1);
                     if (track_p[itr]<5) Htrackdedx_strip_lowp->Fill(ch1,presk);
                     if (track_pt[itr]>55) Htrackdedx_strip->Fill(ch1, presk);
                 }
                 else {
                 // fill only Pixel
                     charge_corr2.push_back(ch1);
                     pathlength2.push_back(dedx_pathlength[iclu]);
                     subdetId2.push_back(dedx_subdetid[iclu]);
                     moduleGeometry2.push_back(dedx_modulgeom[iclu]);
                     mustBeInside2.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning2.push_back(clean1);
                     if (track_p[itr]<5) Htrackdedx_pix_lowp->Fill(ch1,presk);
                     if (track_pt[itr]>55) Htrackdedx_pix->Fill(ch1,presk);
                 }
                 bool no_in_L1_pixel =true;
                 int info_layr=GetLayerLabel(dedx_subdetid[iclu], dedx_detid[iclu],year);
                 if (info_layr==23) no_in_L1_pixel=false;
                 if (no_in_L1_pixel) {
                       // fill both Strip and Pixel but without L1 Pixel
                       charge_corr3.push_back(ch1);
                       pathlength3.push_back(dedx_pathlength[iclu]);
                       subdetId3.push_back(dedx_subdetid[iclu]);
                       detId3.push_back(dedx_detid[iclu]);
                       moduleGeometry3.push_back(dedx_modulgeom[iclu]);
                       mustBeInside3.push_back(dedx_insideTkMod[iclu]);
                       bool_cleaning3.push_back(clean1);
                       if (!dedx_isstrip[iclu]) {
                       // fill only Pixel without L1
                         charge_corr4.push_back(ch1);
                         pathlength4.push_back(dedx_pathlength[iclu]);
                         subdetId4.push_back(dedx_subdetid[iclu]);
                         detId4.push_back(dedx_detid[iclu]);
                         moduleGeometry4.push_back(dedx_modulgeom[iclu]);
                         mustBeInside4.push_back(dedx_insideTkMod[iclu]);
                         bool_cleaning4.push_back(clean1);
                       }
                 }


              float norm_mult = 265; // 247 or 265?
              double Norm = 3.61e-06*norm_mult;
              double scaleFactor = dEdxSF[0];
              //std::cout << "right before dedx_ispixel = " << endl; 
              if (dedx_ispixel[iclu]) {
                  //std::cout << "Bool IS pixel, track p = " << track_p[itr] << endl;
                  Norm = 3.61e-06;
                  scaleFactor *=dEdxSF[1];
                  double ChargeOverPathlength = scaleFactor*Norm*ch1;
                  if (dedx_pathlength[iclu]>0) ChargeOverPathlength/=dedx_pathlength[iclu];
                  else ChargeOverPathlength=0;
                  HHitPix->Fill(ChargeOverPathlength, presk);
                  if (no_in_L1_pixel) {
                   if (fabs(track_eta[itr])<0.4){
                    HHitProfilePix->Fill(track_p[itr], ChargeOverPathlength, presk);
                    HHit2DPix->Fill(track_p[itr], ChargeOverPathlength, presk);
		   }
                   if (track_p[itr]>5 && track_p[itr]<45) {
                      Charge_Vs_Path_noL1->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                   }
                  }

                  /*
                  if (track_p[itr]>10 && track_p[itr]<45) {
                      Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                      for (int l = 0 ; l < ias_intervals; l++){
                          if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                              Charge_Vs_Path_PU[l]->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                          }
                      }

                      if(npv <= 20){
                          Charge_Vs_Path_PU_below_20->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                      }
                      else if (npv > 20 && npv <= 40){
                          Charge_Vs_Path_PU_between->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                      }
                      else if (npv > 40){
                          Charge_Vs_Path_PU_above_40->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                      }
                      
                  }
                  */
              }
              else {
                  //std::cout << "bool is NOT pixel, track p = " << track_p[itr] << std::endl;
                  double ChargeOverPathlength = scaleFactor*Norm*ch1;
		  if (dedx_pathlength[iclu]>0) ChargeOverPathlength/=dedx_pathlength[iclu];
	  	  else ChargeOverPathlength=0;
		  HHitStrip->Fill(ChargeOverPathlength, presk);
		   if (fabs(track_eta[itr])<0.4 && sclus_clusclean2[iclu] && dedx_insideTkMod[iclu] ){
			    HHitProfileStrip->Fill(track_p[itr], ChargeOverPathlength, presk);
			    HHit2DStrip->Fill(track_p[itr], ChargeOverPathlength, presk);
		  }
                  //cout << " NO PIXL Before track_p cut Charge_Vs_Path with modulgeom " << dedx_modulgeom[iclu] << endl;
                  
                  if (track_p[itr]>10 && track_p[itr]<45) {
                 //Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk);
                   Charge_Vs_Path_noL1->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10), presk);

                   /*
                   for (int l = 0 ; l < ias_intervals; l++){
                       if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                           Charge_Vs_Path_PU[l]->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk);
                       }
                   }
                   */

                   if(npv <= 20){
                       Charge_Vs_Path_PU_below_20->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk);
                   }
                   else if (npv > 20 && npv <= 40){
                       Charge_Vs_Path_PU_between->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk);
                   }
                   else if (npv > 40){
                       Charge_Vs_Path_PU_above_40->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk);
                   }
                  }
	      }

            }  // end if cleaning & inside

			     
         }
  
         // APPLY Selection on the Number of Measurements :
         if (charge_corr3.size()>9) {
            // Pixel no L1;
            float norm_mult = 265; // 247 or 265?
            for (unsigned int jch=0;jch<charge_corr4.size();jch++) {
		   float Norm = 3.61e-06;
		   float scaleFactor =dEdxSF[0]*dEdxSF[1];
		   double ChargeOverPathlength = scaleFactor*Norm*charge_corr4[jch];
		   if (pathlength4[jch]>0) ChargeOverPathlength/=pathlength4[jch];
		   else ChargeOverPathlength=0;
                   if (fabs(track_eta[itr])<0.4) HHit2DPix_NoM->Fill(track_p[itr], ChargeOverPathlength, presk);
                   if (track_p[itr]>5 && track_p[itr]<45) {
                      Charge_Vs_Path_noL1_NoM->Fill (moduleGeometry4[jch], pathlength4[jch]*10, scaleFactor*charge_corr4[jch]/(pathlength4[jch]*10*norm_mult),presk);
                      Charge_Vs_Path_NoM->Fill (moduleGeometry4[jch], pathlength4[jch]*10, scaleFactor*charge_corr4[jch]/(pathlength4[jch]*10*norm_mult),presk);
                   }
            }
            // Strip
            //std::cout << "before loop strip , npv = " << npv << std::endl;
            for (unsigned int jch=0;jch<charge_corr1.size();jch++) {
                  double Norm = 3.61e-06*norm_mult;
                  double scaleFactor = dEdxSF[0];
                  double ChargeOverPathlength = scaleFactor*Norm*charge_corr1[jch];
		  if (pathlength1[jch]>0) ChargeOverPathlength/=pathlength1[jch];
	  	  else ChargeOverPathlength=0;
		  HHitStrip->Fill(ChargeOverPathlength, presk);

                  if (fabs(track_eta[itr])<0.4 && bool_cleaning1[jch] && mustBeInside1[jch] )
                             HHit2DStrip_NoM->Fill(track_p[itr], ChargeOverPathlength, presk);
                  if (track_p[itr]>5 && track_p[itr]<45) {
                    Charge_Vs_Path_NoM->Fill (moduleGeometry1[jch], pathlength1[jch]*10, scaleFactor*charge_corr1[jch]/(pathlength1[jch]*10),presk);
                    Charge_Vs_Path_noL1_NoM->Fill (moduleGeometry1[jch], pathlength1[jch]*10, scaleFactor*charge_corr1[jch]/(pathlength1[jch]*10),presk);
                  }
             }
         }
         // endAPPLY
                     

         int nv=0;
         int ns=0;
         double ih_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF,  NULL,2, 0.15,  nv, ns);
         double ih_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF,  NULL,2, 0.15,  nv, ns);
         double ih_pix = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF,  NULL,2, 0.15,  nv, ns);

         double ih0_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF,  NULL,2, 0.,  nv, ns);
         double ih0_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF,  NULL,2, 0.,  nv, ns);
         double ih0_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0., nv, ns);  // <======= important one
         double ih0_pixnoL1 = getdEdX(charge_corr4, pathlength4, subdetId4, moduleGeometry4, bool_cleaning4, mustBeInside4, dEdxSF, NULL,2, 0., nv, ns);

         double ih_corHiDropNoL1 = getdEdX(charge_corr3,  pathlength3,  subdetId3,  moduleGeometry3,  bool_cleaning3,  mustBeInside3,  dEdxSF, NULL,2, 0., 0.15, nv, ns);
         // Ias
         int nval20_0=0;
         int nsat20_0=0;
         double ias_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, dEdxTemplatesNoL1,2, 0., nval20_0, nsat20_0);
         double ias_all = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
          
         double ias_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);

         //double ias_strip_toy_test = getdEdX(//GetRandom()random from template GetRandom->(dEdxTemplatesAll[module][pathlebnght]->Find), pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
         //Reweighting of templates depending on the PU ( # npv )
  
         if(UseTemplatesForPUReweighting){
             double scaleFactor = dEdxSF[0];
             float norm_mult = 265; // 247 or 265?
             double ias_strip_pu_base = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
             double is_strip_pu_base = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);

            
             if (track_p[itr] > 5 && track_p[itr] < 100){
                 if(abs(track_eta[itr])<1){
                     if(track_p[itr] < 10){
                         IhStripVsP_presk->Fill(ih_strip,track_p[itr],presk);
                         IhBestStripVsP_presk->Fill(ih0_noL1,track_p[itr],presk);
                     }
                     else{
                         IhStripVsP_presk->Fill(ih_strip,track_p[itr]);
                         IhBestStripVsP_presk->Fill(ih0_noL1,track_p[itr]);
                         IhBestvsIas_p_5_100->Fill(ih0_noL1,ias_strip_pu_base);
                     }
                 }
                 if(abs(track_eta[itr])<0.2){
                     if(track_p[itr] < 10){
                         IhStripVsPtight_presk->Fill(ih_strip,track_p[itr],presk);
                         IhBestStripVsPtight_presk->Fill(ih0_noL1,track_p[itr],presk);
                     }
                     else{
                         IhStripVsPtight_presk->Fill(ih_strip,track_p[itr]);
                         IhBestStripVsPtight_presk->Fill(ih0_noL1,track_p[itr]);
                          
                     }
                 }

             }
             if (track_p[itr] > 10 && track_p[itr] < 45){
                 if(abs(track_eta[itr])<1){
                    IhBestvsIas_p_10_45->Fill(ih0_noL1,ias_strip_pu_base);
                    //IhBestvsP_10_45->Fill(ih0_noL1,);
                 }
             }

             if(track_p[itr] > 5 && track_p[itr] < 100){
                 if(abs(track_eta[itr])<1){
                     //add presk for prescale in case of p < 5 GeV
                     IhStripVsP->Fill(ih_strip,track_p[itr]);
                     IhBestStripVsP->Fill(ih0_noL1,track_p[itr]);
                 }
                 if(abs(track_eta[itr])<0.2){
                     IhStripVsPtight->Fill(ih_strip,track_p[itr]);
                     IhBestStripVsPtight->Fill(ih0_noL1,track_p[itr]);
                 }
                 
             }
         
             //TOY FOR IAS MC 
             //
             // Get a given modul geom, draw from TH2D random values of charge vs Pathlength 


             //   TOY MC

             

             PT_compute_ias->Fill(track_pt[itr]);
             P_compute_ias->Fill(track_p[itr]);


             P_after_presel->Fill(track_p[itr]);


             //Cut on track_p pour avoir la meme statistique que celle pour les templates
             if (track_p[itr]>10 && track_p[itr]<45) {
                 if(compute_PE){
                     for (int i = 0; i < ias_intervals; i++){
                         if ( npv > ias_top_born[i] && npv <= ias_top_born[i+1]){
                             for (int k = 0 ; k < nPE ; k++){
                                 double ias_poisson_err = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, PE_dEdxTemplatesPU[i][k],2, 0., nval20_0, nsat20_0);
                                 PE_Ias[i][k]->Fill(ias_poisson_err);
                                 if (ih0_noL1 > 3.47){
                                     PE_Ias_cutih[i][k]->Fill(ias_poisson_err);
                                 }
                             }
                         }
                     }
                 }



                 if (npv >= 16 && npv <= 18){
                    double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuLow,2, 0., nval20_0, nsat20_0);
                    Ias_when_PU_in_16_18_base->Fill(ias_strip_pu_base);    
                    Ias_when_PU_in_16_18_triple->Fill(ias_strip_pu_triple);    
                 }
                 if ( npv >= 30 && npv <= 32){
                    double ias_strip_pu_triple =getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuMedium,2, 0., nval20_0, nsat20_0);
                    Ias_when_PU_in_30_32_base->Fill(ias_strip_pu_base);
                    Ias_when_PU_in_30_32_triple->Fill(ias_strip_pu_triple);
      
                 }
                 bool ihcut = false;
                 if (ih0_noL1 > 3.47 ) ihcut = true;

                 eff_ih_PU->Fill(ihcut,npv);
                 for (int l = 0 ; l < ias_intervals; l++){
                     if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                         ih_eff_denom[l] += 1;
                         if (ih0_noL1 > 3.47) ih_eff_num[l] +=1;

                     }
                 }
                 double trf = 0;
                 for (int l = 0 ; l < ias_intervals; l++){
                     if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                         double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                         double is_strip_pu_triple = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                         
                         Ias_when_PU_triple[l]->Fill(ias_strip_pu_triple);
                         Ias_when_PU[l]->Fill(ias_strip_pu_base);

                         Is_when_PU_triple[l]->Fill(is_strip_pu_triple);
                         Is_when_PU[l]->Fill(is_strip_pu_base);

                         PU_distrib[l]->Fill(npv);
                         sums_ias_nocut[l] += ias_strip_pu_base;
                         sums_ias_triple_nocut[l] += ias_strip_pu_triple;

                     }
                 }
                 if (npv > ias_top_born[ias_intervals]) trf = -1;

                 Ias_all_triple_nocut->Fill(trf);
                 Ias_all_base_nocut->Fill(ias_strip_pu_base);

                 for (int l = 0 ; l < ias_intervals; l++){
                     if (ih0_noL1 > 3.47){
                         if ( npv > ias_top_born[l] && npv <= ias_top_born[l+1] ){
                             double ias_strip_pu_triple_l = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);
                             double is_strip_pu_triple_l = getdEdXIs(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPU[l],2, 0., nval20_0, nsat20_0);

                             Ias_all_triple_cutIH->Fill(ias_strip_pu_triple_l);
                             Ias_all_base_cutIH->Fill(ias_strip_pu_base);
                             sums_ias_triple[l]+=ias_strip_pu_triple_l;
                             sums_ias[l]+=ias_strip_pu_base;

                             Ias_when_PU_ih_cut[l]->Fill(ias_strip_pu_base);
                             Ias_when_PU_ih_cut_triple[l]->Fill(ias_strip_pu_triple_l);

                             Is_when_PU_ih_cut[l]->Fill(is_strip_pu_base);
                             Is_when_PU_ih_cut_triple[l]->Fill(is_strip_pu_triple_l);
                         }
                     }
                 }

             }

             else{
                 if (npv <= 20){
                     double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuHigh,2, 0., nval20_0, nsat20_0);
                     Ias_below_5GeV_PU_below_20_base->Fill(ias_strip_pu_base);
                     Ias_below_5GeV_PU_below_20_triple->Fill(ias_strip_pu_triple);


                 }
                 else if (npv > 20 && npv <=40){
                     double ias_strip_pu_triple =getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuMedium,2, 0., nval20_0, nsat20_0);
                     Ias_below_5GeV_PU_between_base->Fill(ias_strip_pu_base);
                     Ias_below_5GeV_PU_between_triple->Fill(ias_strip_pu_triple);
                     

                 }
                 else if (npv > 40){
                     double ias_strip_pu_triple = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, dEdxTemplatesPuHigh,2, 0., nval20_0, nsat20_0); 
                     Ias_below_5GeV_PU_above_40_base->Fill(ias_strip_pu_base);
                     Ias_below_5GeV_PU_above_40_triple->Fill(ias_strip_pu_triple);

                 }              

             }

             IAS_VS_PT->Fill(ias_all,track_pt[itr]);
             //cout << "filling ias_vs_pt with ias : " << ias_all << " for a pt = " << track_pt[itr] << endl;
             
             if(track_pt[itr]>60 && track_pt[itr]<=80){
                ias_bin_pt[0]->Fill(ias_all);
                //cout << "filling ias for pt bin 60-80 with ias : " << ias_all << " and pt = " << track_pt[itr] << endl;
             }
             if(track_pt[itr] > 80 && track_pt[itr] <= 100){
                ias_bin_pt[1]->Fill(ias_all);
                //cout << "filling ias for pt bin 80-100 with ias : " << ias_all << " and pt = " << track_pt[itr] << endl;
             }
             if(track_pt[itr] > 100 && track_pt[itr] <= 120){
                ias_bin_pt[2]->Fill(ias_all);
             }
             if(track_pt[itr] > 120 && track_pt[itr] <= 140){
                ias_bin_pt[3]->Fill(ias_all); 
             }
             if(track_pt[itr] > 140 && track_pt[itr] <= 160){
                ias_bin_pt[4]->Fill(ias_all);
             }
             if(track_pt[itr] > 160 && track_pt[itr] <= 180){
                ias_bin_pt[5]->Fill(ias_all);
             }
             if(track_pt[itr] > 180 && track_pt[itr] <= 200){
                ias_bin_pt[6]->Fill(ias_all);
             }
             if(track_pt[itr] > 200 && track_pt[itr] <= 220){
                ias_bin_pt[7]->Fill(ias_all);
             }
             if(track_pt[itr] > 220 && track_pt[itr] <= 240){
                ias_bin_pt[8]->Fill(ias_all);
             }
             if(track_pt[itr] > 240 && track_pt[itr] <= 260){
                ias_bin_pt[9]->Fill(ias_all);
             }
         } // end PU reweighting


         float mass_strip=-1;
         if (ih_strip - Cval_ldstrip>0.) mass_strip=sqrt((ih_strip - Cval_ldstrip)*track_p[itr]*track_p[itr]/Kval_ldstrip);

         float mass_ih=-1;
         if (ih_cor - Cval_ld>0.) mass_ih=sqrt((ih_cor - Cval_ld)*track_p[itr]*track_p[itr]/Kval_ld);

         float mass_ih0=-1;
         if (ih0_cor - Cval_all>0.) mass_ih0= sqrt((ih0_cor - Cval_all)*track_p[itr]*track_p[itr]/Kval_all);

         float mass_ih0noL1=-1;
         float mass_ih0noL1_2=-1;
         float mass_ih0noL1_3=-1;
         if (ih0_noL1 - Cval_nol1>0.) mass_ih0noL1= sqrt((ih0_noL1 - Cval_nol1)*track_p[itr]*track_p[itr]/Kval_nol1);
         if (ih0_noL1 - Cval_nol1_2>0.) mass_ih0noL1_2= sqrt((ih0_noL1 - Cval_nol1_2)*track_p[itr]*track_p[itr]/Kval_nol1_2);
         float val_nol1_3 = Cval_nol1_3 + (0.92*0.92*Kval_nol1_3)/(6.5*6.5) +Nval_nol1_3*log(6.5/0.92);
         if (track_p[itr]<3) {
           if (ih0_noL1 - val_nol1_3  >0. && computeSpecial) mass_ih0noL1_3= getMassSpecial(ih0_noL1,track_p[itr], Kval_nol1_3, Cval_nol1_3,Nval_nol1_3);
         }

        
         float mass_0strip=-1;
         if (ih0_strip - Cval_strip>0.) mass_0strip= sqrt((ih0_strip - Cval_strip)*track_p[itr]*track_p[itr]/Kval_strip);

         float mass_ihHDnoL1=-1;
         if (ih_corHiDropNoL1 - Cval_hdnol1>0.) mass_ihHDnoL1= sqrt((ih_corHiDropNoL1 - Cval_hdnol1)*track_p[itr]*track_p[itr]/Kval_hdnol1);

         float mass_pixel=-1;
         if (ih_pix - Cval_pix>0) mass_pixel= sqrt((ih_pix - Cval_pix)*track_p[itr]*track_p[itr]/Kval_pix);

         float mass_pixnL1=-1;
         if (ih0_pixnoL1 - Cval_pixnol1>0) mass_pixnL1= sqrt((ih0_pixnoL1 - Cval_pixnol1)*track_p[itr]*track_p[itr]/Kval_pixnol1);


         if (track_p[itr]<20) {
            dEdXVsP_lowp2->Fill(track_p[itr],ih_cor,presk);

            dEdX0VsP_lowp2->Fill(track_p[itr],ih0_cor,presk);
            dEdX0noL1VsP_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            // eta bins motivated by Dylan's presentation on March 12, 2021
            if (fabs(track_eta[itr])<0.91) {
                 dEdX0noL1VsP_eta1_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (fabs(track_eta[itr])<1.74) {
                 dEdX0noL1VsP_eta2_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_eta3_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            if (npv<20.) {
                 dEdX0noL1VsP_pu1_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<30.) {
                 dEdX0noL1VsP_pu2_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<40.) {
                 dEdX0noL1VsP_pu3_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_pu4_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
         
            dEdXHDnoL1VsP_lowp2->Fill(track_p[itr],ih_corHiDropNoL1,presk);
            dEdXstripVsP_lowp2->Fill(track_p[itr],ih_strip,presk);
            dEdXpixVsP_lowp2->Fill(track_p[itr],ih_pix,presk);
            dEdX0stripVsP_lowp2->Fill(track_p[itr],ih0_strip,presk);
            dEdX0pixnoL1VsP_lowp2->Fill(track_p[itr],ih0_pixnoL1,presk);
         }
         if (track_p[itr]<5) {
            Htrackih_lowp->Fill(ih_cor, presk);
            Htrackih_pix_lowp->Fill(ih_pix,presk);
            Htrackih_strip_lowp->Fill(ih_strip,presk);
            Htrackih0_lowp->Fill(ih0_cor,presk);
            Htrackih0noL1_lowp->Fill(ih0_noL1,presk);
            Htracketa_lowp->Fill(track_eta[itr],presk);
            Htrackias_lowp->Fill(ias_noL1,presk);
            Htrackiasall_lowp->Fill(ias_all,presk);

            dEdXVsP_lowp->Fill(track_p[itr],ih_cor,presk);
            dEdXpixVsP_lowp->Fill(track_p[itr],ih_pix,presk);
            dEdXstripVsP_lowp->Fill(track_p[itr],ih_strip,presk);
            dEdX0stripVsP_lowp->Fill(track_p[itr],ih0_strip,presk);
            dEdX0VsP_lowp->Fill(track_p[itr],ih0_cor,presk);
            dEdX0noL1VsP_lowp->Fill(track_p[itr],ih0_noL1,presk);
            // eta bins motivated by Dylan's presentation on March 12, 2021
            if (fabs(track_eta[itr])<0.91) {
                 dEdX0noL1VsP_eta1_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (fabs(track_eta[itr])<1.74) {
                 dEdX0noL1VsP_eta2_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_eta3_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            if (npv<20.) {
                 dEdX0noL1VsP_pu1_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<30.) {
                 dEdX0noL1VsP_pu2_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<40.) {
                 dEdX0noL1VsP_pu3_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_pu4_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            dEdXHDnoL1VsP_lowp->Fill(track_p[itr],ih_corHiDropNoL1,presk);
            dEdX0pixnoL1VsP_lowp->Fill(track_p[itr],ih0_pixnoL1,presk);

            dEstrVsdE_lowp->Fill(ih_cor,ih_strip,presk);
            dEdXstripVsEta_lowp->Fill(track_eta[itr],ih_strip,presk);

            // test proton Mass 

            if (track_p[itr]<3) {
                    lowp_MassIhstrip->Fill(mass_strip,presk);
                    lowp_MassIh->Fill(mass_ih,presk);
                    lowp_MassIh0->Fill(mass_ih0,presk);
                    lowp_MassIh0noL1->Fill(mass_ih0noL1,presk);
                    lowp_MassIh0noL1_2->Fill(mass_ih0noL1_2,presk);
                    lowp_MassIh0noL1_3->Fill(mass_ih0noL1_3,presk);
                       if (ih0_noL1> Cval_nol1) {
                       lowp_MassIh0noL1_11->Fill(mass_ih0noL1);
                       lowp_MassIh0noL1_12->Fill(mass_ih0noL1_2);
                       lowp_MassIh0noL1_13->Fill(mass_ih0noL1_3);
                       }
                       if (ih0_noL1>3.27 + 0.21) {
                       lowp_MassIh0noL1_1s1->Fill(mass_ih0noL1);
                       lowp_MassIh0noL1_1s2->Fill(mass_ih0noL1_2);
                       lowp_MassIh0noL1_1s3->Fill(mass_ih0noL1_3);
                         if (ih0_noL1>3.27 + 2*0.21) {
                         lowp_MassIh0noL1_2s1->Fill(mass_ih0noL1);
                         lowp_MassIh0noL1_2s2->Fill(mass_ih0noL1_2);
                         lowp_MassIh0noL1_2s3->Fill(mass_ih0noL1_3);
                           if (ih0_noL1>3.27 + 3*0.21) {
                           lowp_MassIh0noL1_3s1->Fill(mass_ih0noL1);
                           lowp_MassIh0noL1_3s2->Fill(mass_ih0noL1_2);
                           lowp_MassIh0noL1_3s3->Fill(mass_ih0noL1_3);
                           }
                         }
                       }
                    lowp_MassIh0strip->Fill(mass_0strip,presk);
                    lowp_MassIhHDnoL1->Fill(mass_ihHDnoL1,presk);
                    lowp2d_MassIhstrip->Fill(track_p[itr],mass_strip,presk);
                    lowp2d_MassIh->Fill(track_p[itr],mass_ih,presk);
                    lowp2d_MassIh0->Fill(track_p[itr],mass_ih0,presk);
                    lowp2d_MassIh0noL1->Fill(track_p[itr],mass_ih0noL1,presk);
                    lowp2d_MassIh0strip->Fill(track_p[itr],mass_0strip,presk);
                    lowp2d_MassIhHDnoL1->Fill(track_p[itr],mass_ihHDnoL1,presk);
                    lowp_dEdXpixVsstrip->Fill(ih_strip,ih_pix,presk);
                    lowp_dEdX0pixVsstrip->Fill(ih0_strip,ih0_pixnoL1,presk);
                    lowp2d_Mass_pix_strip15->Fill(mass_strip,mass_pixel,presk);
                    lowp2d_Mass_pix_strip0->Fill(mass_0strip,mass_pixnL1,presk);
                    if (mass_0strip>0 && mass_pixnL1>0 ) {
                      lowp_MassDiff_pix_strip0->Fill(mass_0strip - mass_pixnL1,presk);
                    }
                    if (mass_strip>0 && mass_pixel>0) {
                      lowp_MassDiff_pix_strip15->Fill(mass_strip - mass_pixel,presk);
                    }
              }

         } 
         else {
          // track_p>5
          //
          // test for the background
                bg_lowp2d_Mass_pix_strip0->Fill(mass_0strip,mass_pixnL1,presk);
                if (mass_0strip>0 && mass_pixnL1>0) {
                     float meanval=(mass_0strip+mass_pixnL1)/2;
                     float alpha= (mass_0strip-mass_pixnL1)/meanval;
                     bg_transf_Mass->Fill(meanval,alpha,presk);
                }
          // stability
                bg_dEdX0NoL1VsRun->Fill(runNumber,ih0_noL1,presk);  //  <==== our best estimate for Ih
                bg_dEdXVsIas->Fill(ias_noL1,ih0_noL1,presk);

          std::vector<double> vect;
          std::vector <int> subdetvec;
          std::vector <UInt_t> detvec;
//          for (int h=0; h<charge_corr1.size(); h++) {   //only strip
          for (unsigned int h=0; h<charge_corr3.size(); h++) {      //all cluster
            if(!bool_cleaning[h])continue;
            if(!mustBeInside[h])continue;
            double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
            double scaleFactor = dEdxSF[0];
            if (subdetId[h]<3) scaleFactor*=dEdxSF[1];
            double ChargeOverPathlength = 0;
            if (pathlength[h]>0) ChargeOverPathlength = scaleFactor*Norm*charge_corr[h]/pathlength[h];
            vect.push_back(ChargeOverPathlength); //save charge
            subdetvec.push_back(subdetId[h]);
            detvec.push_back(detId[h]);
            if (10<track_pt[itr] && track_pt[itr]<30) {
              if (subdetId[h]<3) {
                 LowpCharge_Eta_pix->Fill(track_eta[itr],ChargeOverPathlength,presk);
              }
              else { 
                 LowpCharge_Eta_strip->Fill(track_eta[itr],ChargeOverPathlength,presk);
              }
            }
          }
          std::vector <double> tmp (vect.size());
          std::copy (vect.begin(), vect.end(), tmp.begin());
//          std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
//          int nTrunc = tmp.size()*0.15;
//          dEdXstripVsNhit_lowp->Fill(vect.size(),ih_strip,presk);
//          for(unsigned int t=0;t+nTrunc<tmp.size();t++){
          for(unsigned int t=0;t<tmp.size();t++){
/*
              if (subdetvec[t]>2) {
                dEdXstripVsCharge_lowp->Fill(tmp[t],ih_strip,presk);
                if (ih_strip<3.4) Htrackdedx_strip_lowp1->Fill(tmp[t],presk);
                else Htrackdedx_strip_lowp2->Fill(tmp[t],presk);
              }
*/
              int lab=0;
/*
              int indexx=-1;
              for (int h=0; h<vect.size(); h++) {
                float dif = tmp[t]-vect[h];
                if (dif<0) dif*=-1;
                if (dif<0.000001) indexx=h;
              }
//              if (indexx>-1) lab=GetLayerLabel(subdetvec[t], detvec[t]);
              if (indexx>-1) lab=GetLayerLabel(subdetvec[indexx], detvec[indexx]);
*/

              lab=GetLayerLabel(subdetvec[t], detvec[t],year);
              if (lab==1)      Charge_tib1->Fill(tmp[t],presk);
              else if (lab==2) Charge_tib2->Fill(tmp[t],presk);
              else if (lab==3) Charge_tib3->Fill(tmp[t],presk);
              else if (lab==4) Charge_tib4->Fill(tmp[t],presk);
              else if (lab==5) Charge_tob1->Fill(tmp[t],presk);
              else if (lab==6) Charge_tob2->Fill(tmp[t],presk);
              else if (lab==7) Charge_tob3->Fill(tmp[t],presk);
              else if (lab==8) Charge_tob4->Fill(tmp[t],presk);
              else if (lab==9) Charge_tob5->Fill(tmp[t],presk);
              else if (lab==10) Charge_tob6->Fill(tmp[t],presk);
              else if (lab==11) Charge_tid1->Fill(tmp[t],presk);
              else if (lab==12) Charge_tid2->Fill(tmp[t],presk);
              else if (lab==13) Charge_tid3->Fill(tmp[t],presk);
              else if (lab==14) Charge_tec1->Fill(tmp[t],presk);
              else if (lab==15) Charge_tec2->Fill(tmp[t],presk);
              else if (lab==16) Charge_tec3->Fill(tmp[t],presk);
              else if (lab==17) Charge_tec4->Fill(tmp[t],presk);
              else if (lab==18) Charge_tec5->Fill(tmp[t],presk);
              else if (lab==19) Charge_tec6->Fill(tmp[t],presk);
              else if (lab==20) Charge_tec7->Fill(tmp[t],presk);
              else if (lab==21) Charge_tec8->Fill(tmp[t],presk);
              else if (lab==22) Charge_tec9->Fill(tmp[t],presk);
              else if (lab==23) Charge_pixl1->Fill(tmp[t],presk);
              else if (lab==24) Charge_pixl2->Fill(tmp[t],presk);
              else if (lab==25) Charge_pixl3->Fill(tmp[t],presk);
              else if (lab==26) Charge_pixl4->Fill(tmp[t],presk);
              else if (lab==27) Charge_pixd1->Fill(tmp[t],presk);
              else if (lab==28) Charge_pixd2->Fill(tmp[t],presk);
              else if (lab==29) Charge_pixd3->Fill(tmp[t],presk);

   
              if (subdetvec[t]==2) {
                 int layerorside = int((detvec[t]>>23)&0x3); // 1=FPIX- 2=FPIX+
                 if (layerorside==1) {  
                     Charge_pixr1->Fill(tmp[t],presk);
                     ChargeVsRun_pixr1->Fill(runNumber,tmp[t],presk);
                     ZooChargeVsRun_pixr1->Fill(runNumber,tmp[t],presk);
                 }
                 else if (layerorside==2) {
                     Charge_pixr2->Fill(tmp[t],presk);
                     ChargeVsRun_pixr2->Fill(runNumber,tmp[t],presk);
                     ZooChargeVsRun_pixr2->Fill(runNumber,tmp[t],presk);
                 }
              }
/*
              else {
                cout << " lab " << subdetvec[t] << "  "  << detvec[t] ;
                cout << "   "  << lab << endl;
              }
*/

              if (lab==1)      ChargeVsRun_tib1->Fill(runNumber,tmp[t],presk);
              else if (lab==2) ChargeVsRun_tib2->Fill(runNumber,tmp[t],presk);
              else if (lab==3) ChargeVsRun_tib3->Fill(runNumber,tmp[t],presk);
              else if (lab==4) ChargeVsRun_tib4->Fill(runNumber,tmp[t],presk);
              else if (lab==5) ChargeVsRun_tob1->Fill(runNumber,tmp[t],presk);
              else if (lab==6) ChargeVsRun_tob2->Fill(runNumber,tmp[t],presk);
              else if (lab==7) ChargeVsRun_tob3->Fill(runNumber,tmp[t],presk);
              else if (lab==8) ChargeVsRun_tob4->Fill(runNumber,tmp[t],presk);
              else if (lab==9) ChargeVsRun_tob5->Fill(runNumber,tmp[t],presk);
              else if (lab==10) ChargeVsRun_tob6->Fill(runNumber,tmp[t],presk);
              else if (lab==11) ChargeVsRun_tid1->Fill(runNumber,tmp[t],presk);
              else if (lab==12) ChargeVsRun_tid2->Fill(runNumber,tmp[t],presk);
              else if (lab==13) ChargeVsRun_tid3->Fill(runNumber,tmp[t],presk);
              else if (lab==14) ChargeVsRun_tec1->Fill(runNumber,tmp[t],presk);
              else if (lab==15) ChargeVsRun_tec2->Fill(runNumber,tmp[t],presk);
              else if (lab==16) ChargeVsRun_tec3->Fill(runNumber,tmp[t],presk);
              else if (lab==17) ChargeVsRun_tec4->Fill(runNumber,tmp[t],presk);
              else if (lab==18) ChargeVsRun_tec5->Fill(runNumber,tmp[t],presk);
              else if (lab==19) ChargeVsRun_tec6->Fill(runNumber,tmp[t],presk);
              else if (lab==20) ChargeVsRun_tec7->Fill(runNumber,tmp[t],presk);
              else if (lab==21) ChargeVsRun_tec8->Fill(runNumber,tmp[t],presk);
              else if (lab==22) ChargeVsRun_tec9->Fill(runNumber,tmp[t],presk);
              else if (lab==23) ChargeVsRun_pixl1->Fill(runNumber,tmp[t],presk);
              else if (lab==24) ChargeVsRun_pixl2->Fill(runNumber,tmp[t],presk);
              else if (lab==25) ChargeVsRun_pixl3->Fill(runNumber,tmp[t],presk);
              else if (lab==26) ChargeVsRun_pixl4->Fill(runNumber,tmp[t],presk);
              else if (lab==27) ChargeVsRun_pixd1->Fill(runNumber,tmp[t],presk);
              else if (lab==28) ChargeVsRun_pixd2->Fill(runNumber,tmp[t],presk);
              else if (lab==29) ChargeVsRun_pixd3->Fill(runNumber,tmp[t],presk);

              if (lab==1)      ZooChargeVsRun_tib1->Fill(runNumber,tmp[t],presk);
              else if (lab==2) ZooChargeVsRun_tib2->Fill(runNumber,tmp[t],presk);
              else if (lab==3) ZooChargeVsRun_tib3->Fill(runNumber,tmp[t],presk);
              else if (lab==4) ZooChargeVsRun_tib4->Fill(runNumber,tmp[t],presk);
              else if (lab==5) ZooChargeVsRun_tob1->Fill(runNumber,tmp[t],presk);
              else if (lab==6) ZooChargeVsRun_tob2->Fill(runNumber,tmp[t],presk);
              else if (lab==7) ZooChargeVsRun_tob3->Fill(runNumber,tmp[t],presk);
              else if (lab==8) ZooChargeVsRun_tob4->Fill(runNumber,tmp[t],presk);
              else if (lab==9) ZooChargeVsRun_tob5->Fill(runNumber,tmp[t],presk);
              else if (lab==10) ZooChargeVsRun_tob6->Fill(runNumber,tmp[t],presk);
              else if (lab==11) ZooChargeVsRun_tid1->Fill(runNumber,tmp[t],presk);
              else if (lab==12) ZooChargeVsRun_tid2->Fill(runNumber,tmp[t],presk);
              else if (lab==13) ZooChargeVsRun_tid3->Fill(runNumber,tmp[t],presk);
              else if (lab==14) ZooChargeVsRun_tec1->Fill(runNumber,tmp[t],presk);
              else if (lab==15) ZooChargeVsRun_tec2->Fill(runNumber,tmp[t],presk);
              else if (lab==16) ZooChargeVsRun_tec3->Fill(runNumber,tmp[t],presk);
              else if (lab==17) ZooChargeVsRun_tec4->Fill(runNumber,tmp[t],presk);
              else if (lab==18) ZooChargeVsRun_tec5->Fill(runNumber,tmp[t],presk);
              else if (lab==19) ZooChargeVsRun_tec6->Fill(runNumber,tmp[t],presk);
              else if (lab==20) ZooChargeVsRun_tec7->Fill(runNumber,tmp[t],presk);
              else if (lab==21) ZooChargeVsRun_tec8->Fill(runNumber,tmp[t],presk);
              else if (lab==22) ZooChargeVsRun_tec9->Fill(runNumber,tmp[t],presk);
              else if (lab==23) ZooChargeVsRun_pixl1->Fill(runNumber,tmp[t],presk);
              else if (lab==24) ZooChargeVsRun_pixl2->Fill(runNumber,tmp[t],presk);
              else if (lab==25) ZooChargeVsRun_pixl3->Fill(runNumber,tmp[t],presk);
              else if (lab==26) ZooChargeVsRun_pixl4->Fill(runNumber,tmp[t],presk);
              else if (lab==27) ZooChargeVsRun_pixd1->Fill(runNumber,tmp[t],presk);
              else if (lab==28) ZooChargeVsRun_pixd2->Fill(runNumber,tmp[t],presk);
              else if (lab==29) ZooChargeVsRun_pixd3->Fill(runNumber,tmp[t],presk);

            if (10<track_pt[itr] && track_pt[itr]<30) {
              if (lab==1)      LowpCharge_tib1->Fill(tmp[t],presk);
              else if (lab==2) LowpCharge_tib2->Fill(tmp[t],presk);
              else if (lab==3) LowpCharge_tib3->Fill(tmp[t],presk);
              else if (lab==4) LowpCharge_tib4->Fill(tmp[t],presk);
              else if (lab==5) LowpCharge_tob1->Fill(tmp[t],presk);
              else if (lab==6) LowpCharge_tob2->Fill(tmp[t],presk);
              else if (lab==7) LowpCharge_tob3->Fill(tmp[t],presk);
              else if (lab==8) LowpCharge_tob4->Fill(tmp[t],presk);
              else if (lab==9) LowpCharge_tob5->Fill(tmp[t],presk);
              else if (lab==10) LowpCharge_tob6->Fill(tmp[t],presk);
              else if (lab==11) LowpCharge_tid1->Fill(tmp[t],presk);
              else if (lab==12) LowpCharge_tid2->Fill(tmp[t],presk);
              else if (lab==13) LowpCharge_tid3->Fill(tmp[t],presk);
              else if (lab==14) LowpCharge_tec1->Fill(tmp[t],presk);
              else if (lab==15) LowpCharge_tec2->Fill(tmp[t],presk);
              else if (lab==16) LowpCharge_tec3->Fill(tmp[t],presk);
              else if (lab==17) LowpCharge_tec4->Fill(tmp[t],presk);
              else if (lab==18) LowpCharge_tec5->Fill(tmp[t],presk);
              else if (lab==19) LowpCharge_tec6->Fill(tmp[t],presk);
              else if (lab==20) LowpCharge_tec7->Fill(tmp[t],presk);
              else if (lab==21) LowpCharge_tec8->Fill(tmp[t],presk);
              else if (lab==22) LowpCharge_tec9->Fill(tmp[t],presk);
              else if (lab==23) LowpCharge_pixl1->Fill(tmp[t],presk);
              else if (lab==24) LowpCharge_pixl2->Fill(tmp[t],presk);
              else if (lab==25) LowpCharge_pixl3->Fill(tmp[t],presk);
              else if (lab==26) LowpCharge_pixl4->Fill(tmp[t],presk);
              else if (lab==27) LowpCharge_pixd1->Fill(tmp[t],presk);
              else if (lab==28) LowpCharge_pixd2->Fill(tmp[t],presk);
              else if (lab==29) LowpCharge_pixd3->Fill(tmp[t],presk);

   
              if (subdetvec[t]==2) {
                 int layerorside = int((detvec[t]>>23)&0x3); // 1=FPIX- 2=FPIX+
                 if (layerorside==1) {  
                     LowpCharge_pixr1->Fill(tmp[t],presk);
                 }
                 else if (layerorside==2) {
                     LowpCharge_pixr2->Fill(tmp[t],presk);
                 }
              }
            }



          } // end loop on tmp
//          dEdXstripVsNhittrunc_lowp->Fill(vect.size()-nTrunc,ih_strip,presk);
         } // end cut on p



         Htracketa->Fill(track_eta[itr],presk);
         Htrackphi->Fill(track_phi[itr],presk);
         dEdXVsP->Fill(track_p[itr],ih0_noL1,presk);
         dEdXpixVsP->Fill(track_p[itr],ih0_pixnoL1,presk);
         dEdXstripVsP->Fill(track_p[itr],ih0_strip,presk);
         EtaVsPhi_nhit->Fill(track_phi[itr],track_eta[itr],presk);

         Nsat->Fill(nsatclu,presk);
         NPix->Fill(charge_corr4.size(),presk);
         NStrip->Fill(charge_corr1.size(),presk);

         Htrackih_reco->Fill(ih0_noL1,presk);
         Htrackih_pix->Fill(ih0_pixnoL1,presk);
         Htrackih_strip->Fill(ih0_strip,presk);
         Htrackias->Fill(ias_noL1,presk);
         Htrackiasall->Fill(ias_all,presk);
      }
      HNtracks->Fill(ntracks);
      HNtracks1->Fill(ntracks1);
      HNtracks20->Fill(ntracks20);

   }
   OutputHisto->cd();
   HNtracks->Write();
   HNtracks1->Write();
   HNtracks20->Write();

   Htrackpt->Write();
   Htracketa->Write();
   Htracketa_lowp->Write();
   Htrackphi->Write();
   Htracknhit->Write();

   Htrackih_reco->Write();
   Htrackih_pix->Write();
   Htrackih_strip->Write();
   Htrackdedx_pix->Write();
   Htrackdedx_strip->Write();
   Htrackias->Write();
   Htrackiasall->Write();

   Htrackih_lowp->Write();
   Htrackih_pix_lowp->Write();
   Htrackih_strip_lowp->Write();
   Htrackih0_lowp->Write();
   Htrackih0noL1_lowp->Write();
   Htrackdedx_pix_lowp->Write();
   Htrackdedx_strip_lowp->Write();
   Htrackdedx_strip_lowp1->Write();
   Htrackdedx_strip_lowp2->Write();
   Htrackias_lowp->Write();
   Htrackiasall_lowp->Write();


   Nsat->Write();
   NPix->Write();
   NStrip->Write();
   dEdXVsP->Write();
   dEdXpixVsP->Write();
   dEdXstripVsP->Write();

   dEdXVsP_lowp->Write();
   dEdXVsP_lowp2->Write();

   dEdX0VsP_lowp->Write();
   dEdX0VsP_lowp2->Write();
   dEdX0noL1VsP_lowp->Write();
   dEdX0noL1VsP_lowp2->Write();
   dEdX0noL1VsP_eta1_lowp->Write();
   dEdX0noL1VsP_eta1_lowp2->Write();
   dEdX0noL1VsP_eta2_lowp->Write();
   dEdX0noL1VsP_eta2_lowp2->Write();
   dEdX0noL1VsP_eta3_lowp->Write();
   dEdX0noL1VsP_eta3_lowp2->Write();
   dEdX0noL1VsP_pu1_lowp->Write();
   dEdX0noL1VsP_pu1_lowp2->Write();
   dEdX0noL1VsP_pu2_lowp->Write();
   dEdX0noL1VsP_pu2_lowp2->Write();
   dEdX0noL1VsP_pu3_lowp->Write();
   dEdX0noL1VsP_pu3_lowp2->Write();
   dEdX0noL1VsP_pu4_lowp->Write();
   dEdX0noL1VsP_pu4_lowp2->Write();
   dEdXHDnoL1VsP_lowp->Write();
   dEdXHDnoL1VsP_lowp2->Write();
   dEdX0pixnoL1VsP_lowp->Write();
   dEdX0pixnoL1VsP_lowp2->Write();

   dEdXpixVsP_lowp->Write();
   dEdXpixVsP_lowp2->Write();
   dEdXstripVsP_lowp->Write();
   dEdXstripVsP_lowp2->Write();
   dEdX0stripVsP_lowp->Write();
   dEdX0stripVsP_lowp2->Write();

   dEdXstripVsEta_lowp->Write();
   dEstrVsdE_lowp->Write();
   dEdXstripVsNhit_lowp->Write();
   dEdXstripVsNhittrunc_lowp->Write();
   dEdXstripVsCharge_lowp->Write();
   EtaVsPhi_nhit->Write();

   Charge_pixl1->Write();
   Charge_pixl2->Write();
   Charge_pixl3->Write();
   Charge_pixl4->Write();
   Charge_pixd1->Write();
   Charge_pixd2->Write();
   Charge_pixd3->Write();
   Charge_pixr1->Write();
   Charge_pixr2->Write();
   Charge_tib1->Write();
   Charge_tib2->Write();
   Charge_tib3->Write();
   Charge_tib4->Write();
   Charge_tob1->Write();
   Charge_tob2->Write();
   Charge_tob3->Write();
   Charge_tob4->Write();
   Charge_tob5->Write();
   Charge_tob6->Write();
   Charge_tid1->Write();
   Charge_tid2->Write();
   Charge_tid3->Write();
   Charge_tec1->Write();
   Charge_tec2->Write();
   Charge_tec3->Write();
   Charge_tec4->Write();
   Charge_tec5->Write();
   Charge_tec6->Write();
   Charge_tec7->Write();
   Charge_tec8->Write();
   Charge_tec9->Write();

   LowpCharge_tib1->Write();
   LowpCharge_tib2->Write();
   LowpCharge_tib3->Write();
   LowpCharge_tib4->Write();
   LowpCharge_tob1->Write();
   LowpCharge_tob2->Write();
   LowpCharge_tob3->Write();
   LowpCharge_tob4->Write();
   LowpCharge_tob5->Write();
   LowpCharge_tob6->Write();
   LowpCharge_tid1->Write();
   LowpCharge_tid2->Write();
   LowpCharge_tid3->Write();
   LowpCharge_tec1->Write();
   LowpCharge_tec2->Write();
   LowpCharge_tec3->Write();
   LowpCharge_tec4->Write();
   LowpCharge_tec5->Write();
   LowpCharge_tec6->Write();
   LowpCharge_tec7->Write();
   LowpCharge_tec8->Write();
   LowpCharge_tec9->Write();
   LowpCharge_pixl1->Write();
   LowpCharge_pixl2->Write();
   LowpCharge_pixl3->Write();
   LowpCharge_pixl4->Write();
   LowpCharge_pixd1->Write();
   LowpCharge_pixd2->Write();
   LowpCharge_pixd3->Write();
   LowpCharge_pixr1->Write();
   LowpCharge_pixr2->Write();
   LowpCharge_Eta_pix->Write();
   LowpCharge_Eta_strip->Write();

   ChargeVsRun_pixl1->Write();
   ChargeVsRun_pixl2->Write();
   ChargeVsRun_pixl3->Write();
   ChargeVsRun_pixl4->Write();
   ChargeVsRun_pixd1->Write();
   ChargeVsRun_pixd2->Write();
   ChargeVsRun_pixd3->Write();
   ChargeVsRun_pixr1->Write();
   ChargeVsRun_pixr2->Write();
   ChargeVsRun_tib1->Write();
   ChargeVsRun_tib2->Write();
   ChargeVsRun_tib3->Write();
   ChargeVsRun_tib4->Write();
   ChargeVsRun_tob1->Write();
   ChargeVsRun_tob2->Write();
   ChargeVsRun_tob3->Write();
   ChargeVsRun_tob4->Write();
   ChargeVsRun_tob5->Write();
   ChargeVsRun_tob6->Write();
   ChargeVsRun_tid1->Write();
   ChargeVsRun_tid2->Write();
   ChargeVsRun_tid3->Write();
   ChargeVsRun_tec1->Write();
   ChargeVsRun_tec2->Write();
   ChargeVsRun_tec3->Write();
   ChargeVsRun_tec4->Write();
   ChargeVsRun_tec5->Write();
   ChargeVsRun_tec6->Write();
   ChargeVsRun_tec7->Write();
   ChargeVsRun_tec8->Write();
   ChargeVsRun_tec9->Write();

   ZooChargeVsRun_pixl1->Write();
   ZooChargeVsRun_pixl2->Write();
   ZooChargeVsRun_pixl3->Write();
   ZooChargeVsRun_pixl4->Write();
   ZooChargeVsRun_pixd1->Write();
   ZooChargeVsRun_pixd2->Write();
   ZooChargeVsRun_pixd3->Write();
   ZooChargeVsRun_pixr1->Write();
   ZooChargeVsRun_pixr2->Write();
   ZooChargeVsRun_tib1->Write();
   ZooChargeVsRun_tib2->Write();
   ZooChargeVsRun_tib3->Write();
   ZooChargeVsRun_tib4->Write();
   ZooChargeVsRun_tob1->Write();
   ZooChargeVsRun_tob2->Write();
   ZooChargeVsRun_tob3->Write();
   ZooChargeVsRun_tob4->Write();
   ZooChargeVsRun_tob5->Write();
   ZooChargeVsRun_tob6->Write();
   ZooChargeVsRun_tid1->Write();
   ZooChargeVsRun_tid2->Write();
   ZooChargeVsRun_tid3->Write();
   ZooChargeVsRun_tec1->Write();
   ZooChargeVsRun_tec2->Write();
   ZooChargeVsRun_tec3->Write();
   ZooChargeVsRun_tec4->Write();
   ZooChargeVsRun_tec5->Write();
   ZooChargeVsRun_tec6->Write();
   ZooChargeVsRun_tec7->Write();
   ZooChargeVsRun_tec8->Write();
   ZooChargeVsRun_tec9->Write();


   dEdXVsRun->Write();
   dEdXpixVsRun->Write();
   dEdXstripVsRun->Write();
   dEdXNoL1pixVsRun->Write();
   dEdXNoL1VsRun->Write();
   dEdXHiDropVsRun->Write();
   dEdXpixHiDropVsRun->Write();
   dEdXstripHiDropVsRun->Write();
   dEdXHiDropNoL1VsRun->Write();
/*
   dEdXVsIL->Write();
   dEdXpixVsIL->Write();
   dEdXstripVsIL->Write();
*/
   dEdX0VsRun->Write();
   dEdX0pixVsRun->Write();
   dEdX0stripVsRun->Write();
   MassStripVsRun->Write();
   MassNoL1VsRun->Write();
   dEdX0NoL1pixVsRun->Write();
   dEdX0NoL1VsRun->Write();
   dEdX4VsRun->Write();
   dEdX4pixVsRun->Write();
   dEdX4stripVsRun->Write();
   dEdX40VsRun->Write();
   dEdX40pixVsRun->Write();
   dEdX40stripVsRun->Write();
   bg_dEdX0NoL1VsRun->Write();
   iasNoL1VsRun->Write();
   iasAllVsRun->Write();
   NmeasVsRun->Write();
   NmeasPixVsRun->Write();
   NmeasStrVsRun->Write();
   Nmeas0VsRun->Write();
   NmeasPix0VsRun->Write();
   NmeasStr0VsRun->Write();
   NsatVsRun->Write();
   NsatPixVsRun->Write();
   NsatStrVsRun->Write();
   Nsat0VsRun->Write();
   NsatPix0VsRun->Write();
   NsatStr0VsRun->Write();
   ptVsRun->Write();
   nPVVsRun->Write();
   invBVsRun->Write();
   errinvBVsRun->Write();
   invBDTVsRun->Write();
   invBCSCVsRun->Write();
   invBnewVsRun->Write();
   invBnewDTVsRun->Write();
   invBnewCSCVsRun->Write();
   timeVsRun->Write();
   lumiVsRun->Write();
   HSCP_dEdX->Write(); 
   HSCP_dEdXpix->Write(); 
   HSCP_dEdXstrip->Write(); 
   HSCP_dEdX0->Write(); 
   HSCP_dEdX0pix->Write(); 
   HSCP_dEdX0strip->Write(); 
   HSCP_MassIh->Write(); 
   HSCP_MassIh0->Write(); 
   HSCP_MassIhstrip->Write(); 
   HSCP_MassIh0noL1->Write(); 
   HSCP_MassIh0noL1_2->Write(); 
   HSCP_MassIh0noL1_3->Write(); 
   HSCP_MassIh0noL1_11->Write(); 
   HSCP_MassIh0noL1_12->Write(); 
   HSCP_MassIh0noL1_13->Write(); 
   HSCP_MassIh0noL1_1s1->Write(); 
   HSCP_MassIh0noL1_1s2->Write(); 
   HSCP_MassIh0noL1_1s3->Write(); 
   HSCP_MassIh0noL1_2s1->Write(); 
   HSCP_MassIh0noL1_2s2->Write(); 
   HSCP_MassIh0noL1_2s3->Write(); 
   HSCP_MassIh0noL1_3s1->Write(); 
   HSCP_MassIh0noL1_3s2->Write(); 
   HSCP_MassIh0noL1_3s3->Write(); 
   HSCP_MassIhHDnoL1->Write(); 
   HSCP_MassTOF->Write(); 
   HSCP2d_MassTOFvsIh->Write(); 
   HSCP_MassIh0strip->Write(); 
   HSCP2d_MassIh->Write(); 
   HSCP2d_MassIh0->Write(); 
   HSCP2d_MassIhstrip->Write(); 
   HSCP2d_MassIh0noL1->Write(); 
   HSCP2d_MassIhHDnoL1->Write(); 
   HSCP2d_MassIh0strip->Write(); 
   HSCP2d_Mass_pix_strip15->Write(); 
   HSCP2d_Mass_pix_strip0->Write(); 
   HSCP2d_Mass_pix_strip->Write(); 
   HSCP_MassDiff_pix_strip0->Write(); 
   HSCP_MassDiff_pix_strip15->Write(); 
   HSCP_MassResol_pix_strip0->Write(); 
   HSCP_MassResol_pix_strip15->Write(); 
   lowp_MassIh->Write(); 
   lowp_MassIh0->Write(); 
   lowp_MassIhstrip->Write(); 
   lowp_MassIh0noL1->Write(); 
   lowp_MassIh0noL1_2->Write(); 
   lowp_MassIh0noL1_3->Write(); 
   lowp_MassIh0noL1_11->Write(); 
   lowp_MassIh0noL1_12->Write(); 
   lowp_MassIh0noL1_13->Write(); 
   lowp_MassIh0noL1_1s1->Write(); 
   lowp_MassIh0noL1_1s2->Write(); 
   lowp_MassIh0noL1_1s3->Write(); 
   lowp_MassIh0noL1_2s1->Write(); 
   lowp_MassIh0noL1_2s2->Write(); 
   lowp_MassIh0noL1_2s3->Write(); 
   lowp_MassIh0noL1_3s1->Write(); 
   lowp_MassIh0noL1_3s2->Write(); 
   lowp_MassIh0noL1_3s3->Write(); 
   lowp_MassIhHDnoL1->Write(); 
   lowp_MassIh0strip->Write(); 
   lowp2d_MassIh->Write(); 
   lowp2d_MassIh0->Write(); 
   lowp2d_MassIhstrip->Write(); 
   lowp2d_MassIh0noL1->Write(); 
   lowp2d_MassIhHDnoL1->Write(); 
   lowp2d_MassIh0strip->Write(); 
   lowp_dEdXpixVsstrip->Write();
   lowp_dEdX0pixVsstrip->Write();
   lowp2d_Mass_pix_strip15->Write(); 
   lowp2d_Mass_pix_strip0->Write(); 
   bg_lowp2d_Mass_pix_strip0->Write(); 
   bg_transf_Mass->Write(); 
   lowp_MassDiff_pix_strip0->Write(); 
   lowp_MassDiff_pix_strip15->Write(); 
   bg_dEdXVsIas->Write(); 
   bg_test_event_Mass->Write(); 
   bg_test_prefiring_Mass->Write(); 
   bg_test_nohem_Mass->Write(); 

   HSCP_dEdXpixVsstrip->Write();
   HSCP_dEdX0pixVsstrip->Write();
   HSCP_dEdXstripVsall->Write();
   HSCP_dEdXpixVsall->Write();
   HSCP_dEdXHiDrop->Write();
   HSCP_dEdXstripHiDrop->Write();
   HSCP_dEdXpixHiDrop->Write();
   HSCP_dEdXHiDropNoL1->Write();
   HSCP_dEdX0NoL1->Write();
   NB_PV->Write();

   HSCP_FMIP4->Write(); 
   HSCP_FMIP3p5->Write(); 
   HSCP_FMIP3p2->Write(); 
   FMIP4VsRun->Write();
   FMIP3p5VsRun->Write();
   FMIP3p5VsRun->Write();
   FMIP4VsEta->Write();
   HSCP_iasnol1->Write(); 
   HSCP_iasall->Write(); 
   HSCP_iasstrip->Write(); 
   HSCP_iaspix->Write(); 
   HSCP_probQ->Write();
   HSCP_probQNoL1->Write();
   HSCP_probXY->Write();
   HSCP_probXYNoL1->Write();
   probQVsRun->Write();
   probQNoL1VsRun->Write();
   probXYVsRun->Write();
   probXYNoL1VsRun->Write();
   probQVsIas->Write();

   HSCP_pt->Write(); 
   HSCP_eta->Write(); 
   HSCP_iso_eop->Write(); 
   nPV->Write(); 
   HSCP_invB->Write(); 
   HSCP_errinvB->Write(); 
   HSCP_invBDT->Write(); 
   HSCP_invBCSC->Write(); 
   HSCP_time->Write(); 
   HSCP_npix->Write(); 
   HSCP_nstrip->Write(); 
   HSCP_nmpix->Write(); 
   HSCP_nmstrip->Write(); 
   HSCP_nratio->Write(); 
   HSCP_nmratio->Write(); 

   R1_StdEdXVsEvent->Write();
   R1_StdEdXVsLumi->Write();
   R1_LumiVsEvent->Write();
   R1_nPVVsEvent->Write();
   R1_CandVsEvent->Write();
   R2_StdEdXVsEvent->Write();
   R2_StdEdXVsLumi->Write();
   R2_LumiVsEvent->Write();
   R2_nPVVsEvent->Write();
   R2_CandVsEvent->Write();
   R3_StdEdXVsEvent->Write();
   R3_StdEdXVsLumi->Write();
   R3_LumiVsEvent->Write();
   R3_nPVVsEvent->Write();
   R3_CandVsEvent->Write();

   HHitPix->Write();
   HHitProfilePix->Write();
   HHit2DPix->Write();
   HHit2DPix_NoM->Write();
   HHitStrip->Write();
   HHitProfileStrip->Write();
   HHit2DStrip->Write();
   HHit2DStrip_NoM->Write();

   OutputHisto->Close();


   

   if (writeTemplateOnDisk || writeTptHSCP) {
       OutputTemplate->cd();
  
       Charge_Vs_Path->Write();
       Charge_Vs_Path_PU_above_40->Write();
       Charge_Vs_Path_PU_between->Write();
       Charge_Vs_Path_PU_below_20->Write();
       for (int i = 0; i< ias_intervals; i++){
          Charge_Vs_Path_PU[i]->Write();
       }
       Charge_Vs_Path_noL1->Write();
       Charge_Vs_Path_NoM->Write();
       Charge_Vs_Path_noL1_NoM->Write();
       if(UseTemplatesForPUReweighting) Charge_Vs_Path_PU_corr->Write();

       OutputTemplate->Close();     
   }
   else{
       OutIasStudy->cd();
       for (int i = 0 ; i<ias_intervals ; i++){
          //normal are wit cut on IH
          means_ias[i] = Ias_when_PU_ih_cut[i]->GetMean();
          means_ias_triple[i] = Ias_when_PU_ih_cut_triple[i]->GetMean();


          //Following are without cut on IH
          
          means_ias_nocut[i] = Ias_when_PU[i]->GetMean();

          means_ias_triple_nocut[i] = Ias_when_PU_triple[i]->GetMean();
          
          
          error_mean_triple_cutih[i] = Ias_when_PU_ih_cut_triple[i]->GetMeanError(); 
          error_mean_base_cutih[i] = Ias_when_PU_ih_cut[i]->GetMeanError(); 
   
          error_mean_base[i] = Ias_when_PU[i]->GetMeanError();
          error_mean_triple[i] = Ias_when_PU_triple[i]->GetMeanError();

       }

       double ias_bins[ias_intervals];   
       for (int j = 0 ; j < ias_intervals ; j++){
            ias_bins[j] = PU_distrib[j]->GetMean();
            ias_axis_names[j] = std::to_string(PU_distrib[j]->GetMean()).substr(0, std::to_string(PU_distrib[j]->GetMean()).find(".") + precisionVal + 1);
       }

       double diff_mean_pe_pu_cutih[ias_intervals] = {0};
       double diff_mean_pe_pu[ias_intervals] = {0};
       double diff_qtl_pe_pu[ias_intervals] = {0};
       double diff_qtl_pe_pu_cutih[ias_intervals] = {0};
       
       double sigmas_pe_pu_cutih[ias_intervals] = {0};
       double sigmas_pe_pu[ias_intervals] = {0};
       double sigmas_qtl_pu[ias_intervals] = {0};
       double sigmas_qtl_pu_cutih[ias_intervals] = {0};

       double mean_of_means[ias_intervals] = {0};
       double mean_of_means_cutih[ias_intervals] = {0};
       if (compute_PE){
           for(int j = 0; j < ias_intervals; j++){
               cout << "pu bin #" << j <<endl;
               for (int i = 0; i<nPE; i++){
                   mean_of_means[j] += PE_Ias[j][i]->GetMean();
                   cout << "mean of means nb " << i << " = " << PE_Ias[j][i]->GetMean() << endl;
                   mean_of_means_cutih[j] += PE_Ias_cutih[j][i]->GetMean();
                   cout << "mean of means cut ih nb " << i << " = " << PE_Ias_cutih[j][i]->GetMean() << endl;
               }
               cout << "mean of means #" << j << " before division : " << mean_of_means[j] << endl;
               mean_of_means[j] = (mean_of_means[j]*1.0/nPE);
               cout << "mean of means cut ih #" << j << " before division : " << mean_of_means_cutih[j] << endl;
               mean_of_means_cutih[j] = (mean_of_means_cutih[j]*1.0/nPE);
           }
           for(int j = 0; j < ias_intervals; j++){ 
               for (int i = 0; i<nPE; i++){
                   means_ias_pe[j][i] = PE_Ias[j][i]->GetMean();
                   err_means_ias_pe[j][i] = PE_Ias[j][i]->GetMeanError();
                   std_dev_means_ias_pe[j][i] = PE_Ias[j][i]->GetStdDev();
                   diff_mean_pe_pu[j]+= pow((PE_Ias[j][i]->GetMean() - mean_of_means[j]),2);
    
                   means_ias_pe_cutih[j][i] = PE_Ias_cutih[j][i]->GetMean();
                   err_means_ias_pe_cutih[j][i] = PE_Ias_cutih[j][i]->GetMeanError();
                   std_dev_means_ias_pe_cutih[j][i] = PE_Ias_cutih[j][i]->GetStdDev();
                   diff_mean_pe_pu_cutih[j]+= pow((PE_Ias_cutih[j][i]->GetMean() - mean_of_means_cutih[j]),2);
               }
           }
           for (int i = 0 ; i < ias_intervals ; i++){
               sigmas_pe_pu_cutih[i] = sqrt(diff_mean_pe_pu_cutih[i]*1.0/nPE);   
               sigmas_pe_pu[i] = sqrt(diff_mean_pe_pu[i]*1.0/nPE);
               std::cout << "PU bin #"<<i+1<<" ,has mean of means : " << mean_of_means_cutih[i] << " and PE std dev = "<< sigmas_pe_pu_cutih[i] <<endl; 
           }
       }
       


       
       IAS_triple_when_simple_above_0p115->Write();
       IAS_simple_when_triple_above_0p115->Write();
       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE->Write();
       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP->Write();
       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_LOWPU->Write();
       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_HIGHPU->Write();
       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_STRIP_NOSATURATED->Write();
       CHARGE_OVER_PATHLENGTH_GIVEN_MODULE_PIXEL->Write();

       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->Write();
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Write();
      
       IAS_TRIPLE_PT_50_60_selection->Write();
       IAS_SINGLE_PT_50_60_selection->Write();


 
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->Write();
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->Write();

       TCanvas *masses_pu_pt50_qtl80_90 = new TCanvas("Masses_single_multiple_ihcut_pt_supp50_qtlias_80_90","Masses for pt > 50 GeV, with Ias in [80-90]qtl (0.08,0.115) single/multiple, and ih cut 3.47");
       masses_pu_pt50_qtl80_90->SetLogy();
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->SetMarkerColor(2);
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->SetMarkerStyle(2);
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->SetLineColor(2);
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->SetLineWidth(2);
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Draw("HIST");
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Scale(1./MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->Integral());
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->SetMarkerColor(4);
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->SetMarkerStyle(2);
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->SetLineColor(4);
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->SetLineWidth(2);
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->Scale(1./MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->Integral());
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->Draw("same HIST");
       auto leg_masses_pt50_qtl80_90 = new TLegend(0.2, 0.6, 0.35, 0.8);
       leg_masses_pt50_qtl80_90->AddEntry(MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias,"Triple","l");
       leg_masses_pt50_qtl80_90->AddEntry(MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias,"Single","l");
       leg_masses_pt50_qtl80_90->Draw("same");
       masses_pu_pt50_qtl80_90->Write();



       TCanvas *masses_pu_pt100_qtl80_90 = new TCanvas("Masses_single_multiple_ihcut_pt_sup100_qtlias_80_90","Masses for pt > 100 GeV, with Ias in [80-90]qtl (0.08,0.115) single/multiple, and ih cut 3.47");
       masses_pu_pt100_qtl80_90->SetLogy();
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->SetMarkerColor(2);
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->SetMarkerStyle(2);
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->SetLineColor(2);
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->SetLineWidth(2);
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->Draw("HIST");
       MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->Scale(1./MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias->Integral());
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->SetMarkerColor(4);
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->SetMarkerStyle(2);
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->SetLineColor(4);
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->SetLineWidth(2);
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->Scale(1./MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->Integral());
       MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias->Draw("same HIST");
       auto leg_masses_pt100_qtl80_90 = new TLegend(0.2, 0.6, 0.35, 0.8);
       leg_masses_pt100_qtl80_90->AddEntry(MASS_TRIPLE_ihcut_ptsupp100_qtl_80_90_ias,"Triple","l");
       leg_masses_pt100_qtl80_90->AddEntry(MASS_SINGLE_ihcut_ptsupp100_qtl_80_90_ias,"Single","l");
       leg_masses_pt100_qtl80_90->Draw("same");
       masses_pu_pt100_qtl80_90->Write();


       TCanvas *masses_pu = new TCanvas("Masses_single_multiple_ihcut","Masses for p [10-45] GeV, after cut on Ias > 0.115 (ias single/multiple) and ih cut 3.47");
       masses_pu->SetLogy();
       MASS_TRIPLE_ihcut->SetMarkerColor(2);
       MASS_TRIPLE_ihcut->SetMarkerStyle(2);
       MASS_TRIPLE_ihcut->SetLineColor(2);
       MASS_TRIPLE_ihcut->Draw();
       MASS_SINGLE_ihcut->SetMarkerColor(4);
       MASS_SINGLE_ihcut->SetMarkerStyle(2);
       MASS_SINGLE_ihcut->SetLineColor(4);
       MASS_SINGLE_ihcut->Draw("SAME");


       auto leg_masses = new TLegend(0.2, 0.6, 0.35, 0.8);
       leg_masses->AddEntry(MASS_TRIPLE_ihcut,"Triple","l");
       leg_masses->AddEntry(MASS_SINGLE_ihcut,"Single","l");
       leg_masses->Draw("same");
       masses_pu->Write();
       // ***** plot all normalized IAS distribution ******
      
     
       NPV_ias_triple_cutih->Write();
       NPV_ias_single_cutih->Write();
       NPV_ih_cutih->Write();
       NPV_ih_cut3p29->Write();
       NPV_ih_nocut->Write();
       NPV_ih_cut1->Write();   

       NPV_ih_strip_cutih->Write();
       NPV_ih_strip_cut3p29->Write();
       NPV_ih_strip_nocut->Write();
       NPV_ih_strip_cut1->Write();   
       IH0_noL1->Write();

       PU_ias_single_0p3->Write();
       PU_ias_triple_0p3->Write();
       P_LOWPU_pt50->Write();
       P_HIGHPU_pt50->Write();
       
       P_PU_28_29_pt50->Write();
       P_PU_30_31_pt50->Write();


       P_PU_28_30_pt50->Write();
       P_PU_32_34_pt50->Write();

       P_PU_08_10_pt50->Write();
       P_PU_58_60_pt50->Write();

       NPV_mass_cutih->Write();
       NPV_mass_cut3p29->Write();

       P_triple_pt50_ias_qtl_80_90->Write();
       P_single_pt50_ias_qtl_80_90->Write();
       IAS_single_pt50->Write();
       IAS_triple_pt50->Write();


       MASS_TRIPLE_ihcut->Write();
       MASS_SINGLE_ihcut->Write();

       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->GetXaxis()->SetTitle("MASS");
       MASS_SINGLE_ihcut_ptsupp50_qtl_80_90_ias->GetYaxis()->SetTitle("events");

       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->GetXaxis()->SetTitle("MASS");
       MASS_TRIPLE_ihcut_ptsupp50_qtl_80_90_ias->GetYaxis()->SetTitle("events");



       MASS_SINGLE_ihcut_base_notriple_ias_0p1_ptsupp50_qtl_80_90_ias->Write();
       MASS_TRIPLE_ihcut_triple_nobase_ias_0p1_ptsupp50_qtl_80_90_ias->Write();
         
       MASS_SINGLE_ihcut_base_notriple_ias_0p1->Write();
       MASS_TRIPLE_ihcut_triple_nobase_ias_0p1->Write();

       MASS_SINGLE_nocut->Write();
       MASS_TRIPLE_nocut->Write();
   

  

        

       TCanvas *ias_pu = new TCanvas("mean_IH_cut","<Ias> VS PU scenarios from templates, cut Ih > 3.47"); 
       auto ias_graph = new TGraphErrors(ias_intervals,ias_bins,means_ias,err_mean_x_base,error_mean_base_cutih);
       ias_graph->SetTitle("<Ias> single/multiple templates, with IH > 3.47");
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph->GetXaxis()->SetBinLabel(ias_graph->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
          if (u != (ias_intervals-1)){
             ias_diff_2by2[u]= means_ias[u+1]/means_ias[u];
             std::cout << "Ratio ( Ias_high / Ias_low ) for each consecutive bins in % :" << ias_diff_2by2[u] * 100 << " %" << std::endl; 

          }
       }
       for (int k = 0 ; k < ias_intervals ; k++){
           cout << "Eff of IH cut for PU between " << int(ias_top_born[k]) << " and " << int(ias_top_born[k+1]) << " is " << ih_eff_num[k] << " / " <<ih_eff_denom[k] << " = " << ih_eff_num[k]*1.0/ih_eff_denom[k] << endl;
       }
       
       auto ias_graph_triple = new TGraphErrors(ias_intervals,ias_bins,means_ias_triple,err_mean_x_triple,sigmas_pe_pu_cutih);
       
       ias_graph_triple->SetLineColor(2);
       ias_graph->SetLineColor(4);
       ias_graph->Draw("AC*");
       ias_graph_triple->Draw("SAME C*");

       ias_pu->Write(); 

       TCanvas *ias_pu_nocut = new TCanvas("mean_no_cut","<Ias> VS PU scenarios from templates, no cut IH");

       auto ias_graph_nocutih = new TGraphErrors(ias_intervals,ias_bins,means_ias_nocut,err_mean_x_base,error_mean_base);
       ias_graph_nocutih->SetTitle("<Ias> single/multiple templates, no cut IH"); 
       auto ias_graph_triple_nocutih = new TGraphErrors(ias_intervals,ias_bins,means_ias_triple_nocut,err_mean_x_base,sigmas_pe_pu);

       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_nocutih->GetXaxis()->SetBinLabel(ias_graph->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }
      

       ias_graph_nocutih->SetLineColor(4);
       ias_graph_triple_nocutih->SetLineColor(2);
       ias_graph_nocutih->Draw("AC*");
       ias_graph_triple_nocutih->Draw("SAME C*");

       ias_pu_nocut->Write();

       double bin_pe[nPE] = { 0 };
       for (int i = 0 ; i < nPE; i++){
           bin_pe[i] = i+1;
       }
       cout << "after bin PE 1 to 101"<<endl;
       double ias_mean_0[nPE] = {0}, ias_mean_0_cutih[nPE] = {0};
       double err_ias_mean_0[nPE] = {0}, err_ias_mean_0_cutih[nPE] = {0};
       double ias_stddev_0[nPE] = {0}, ias_stddev_0_cutih[nPE] = {0};


       double ias_mean_1[nPE] = {0}, ias_mean_1_cutih[nPE] = {0};
       double err_ias_mean_1[nPE] = {0}, err_ias_mean_1_cutih[nPE] = {0};
       double ias_stddev_1[nPE] = {0}, ias_stddev_1_cutih[nPE] = {0};

       double ias_mean_2[nPE] = {0}, ias_mean_2_cutih[nPE] = {0};
       double err_ias_mean_2[nPE] = {0}, err_ias_mean_2_cutih[nPE] = {0};
       double ias_stddev_2[nPE] = {0}, ias_stddev_2_cutih[nPE] = {0};

       double ias_mean_3[nPE] = {0}, ias_mean_3_cutih[nPE] = {0};
       double err_ias_mean_3[nPE] = {0}, err_ias_mean_3_cutih[nPE] = {0};
       double ias_stddev_3[nPE] = {0}, ias_stddev_3_cutih[nPE] = {0};

       double ias_mean_4[nPE] = {0}, ias_mean_4_cutih[nPE] = {0};
       double err_ias_mean_4[nPE] = {0}, err_ias_mean_4_cutih[nPE] = {0};
       double ias_stddev_4[nPE] = {0}, ias_stddev_4_cutih[nPE] = {0};

       if(compute_PE){
           for (int i = 0; i < nPE; i++){
               ias_mean_0[i] = means_ias_pe[0][i];
               ias_mean_0_cutih[i] = means_ias_pe_cutih[0][i];
               err_ias_mean_0[i] = err_means_ias_pe[0][i];
               err_ias_mean_0_cutih[i] = err_means_ias_pe_cutih[0][i];
               ias_stddev_0[i] = std_dev_means_ias_pe[0][i];
               ias_stddev_0_cutih[i] = std_dev_means_ias_pe_cutih[0][i];
    
    
               ias_mean_1[i] = means_ias_pe[1][i];
               ias_mean_1_cutih[i] = means_ias_pe_cutih[1][i];
               err_ias_mean_1[i] = err_means_ias_pe[1][i];
               err_ias_mean_1_cutih[i] = err_means_ias_pe_cutih[1][i];
               ias_stddev_1[i] = std_dev_means_ias_pe[1][i];
               ias_stddev_1_cutih[i] = std_dev_means_ias_pe_cutih[1][i];
    
    
               ias_mean_2[i] = means_ias_pe[2][i];
               ias_mean_2_cutih[i] = means_ias_pe_cutih[2][i];
               err_ias_mean_2[i] = err_means_ias_pe[2][i];
               err_ias_mean_2_cutih[i] = err_means_ias_pe_cutih[2][i];
               ias_stddev_2[i] = std_dev_means_ias_pe[2][i];
               ias_stddev_2_cutih[i] = std_dev_means_ias_pe_cutih[2][i];
    
    
               ias_mean_3[i] = means_ias_pe[3][i];
               ias_mean_3_cutih[i] = means_ias_pe_cutih[3][i];
               err_ias_mean_3[i] = err_means_ias_pe[3][i];
               err_ias_mean_3_cutih[i] = err_means_ias_pe_cutih[3][i];
               ias_stddev_3[i] = std_dev_means_ias_pe[3][i];
               ias_stddev_3_cutih[i] = std_dev_means_ias_pe_cutih[3][i];
    
    
               ias_mean_4[i] = means_ias_pe[4][i];
               ias_mean_4_cutih[i] = means_ias_pe_cutih[4][i];
               err_ias_mean_4[i] = err_means_ias_pe[4][i];
               err_ias_mean_4_cutih[i] = err_means_ias_pe_cutih[4][i];
               ias_stddev_4[i] = std_dev_means_ias_pe[4][i];
               ias_stddev_4_cutih[i] = std_dev_means_ias_pe_cutih[4][i];

           }
       

       
           TCanvas *pe_ias_mean_0 = new TCanvas("pe_ias_mean_PU_bin_1","pe_ias_mean_PU_bin_1");
           auto pe_ias_graph_0 = new TGraphErrors(nPE,bin_pe,ias_mean_0,err_mean_x_base,err_ias_mean_0);
           pe_ias_graph_0->SetMarkerStyle(21);
           pe_ias_graph_0->Draw("APL");
           pe_ias_mean_0->Write();
    
           TCanvas *pe_ias_mean_1 = new TCanvas("pe_ias_mean_PU_bin_2","pe_ias_mean_PU_bin_2");
           auto pe_ias_graph_1 = new TGraphErrors(nPE,bin_pe,ias_mean_1,err_mean_x_base,err_ias_mean_1);
           pe_ias_graph_1->SetMarkerStyle(21);
           pe_ias_graph_1->Draw("APL");
           pe_ias_mean_1->Write();
           
           TCanvas *pe_ias_mean_2 = new TCanvas("pe_ias_mean_PU_bin_3","pe_ias_mean_PU_bin_3");
           auto pe_ias_graph_2 = new TGraphErrors(nPE,bin_pe,ias_mean_2,err_mean_x_base,err_ias_mean_2);
           pe_ias_graph_2->SetMarkerStyle(21);
           pe_ias_graph_2->Draw("APL");
           pe_ias_mean_2->Write();
    
           TCanvas *pe_ias_mean_3 = new TCanvas("pe_ias_mean_PU_bin_4","pe_ias_mean_PU_bin_4");
           auto pe_ias_graph_3 = new TGraphErrors(nPE,bin_pe,ias_mean_3,err_mean_x_base,err_ias_mean_3);
           pe_ias_graph_3->SetMarkerStyle(21);
           pe_ias_graph_3->Draw("APL");
           pe_ias_mean_3->Write();
    
    
           TCanvas *pe_ias_mean_4 = new TCanvas("pe_ias_mean_PU_bin_5","pe_ias_mean_PU_bin_5");
           auto pe_ias_graph_4 = new TGraphErrors(nPE,bin_pe,ias_mean_4,err_mean_x_base,err_ias_mean_4);
           pe_ias_graph_4->SetMarkerStyle(21);
           pe_ias_graph_4->Draw("APL");
           pe_ias_mean_4->Write();
    
           //SAME WITH IH CUT 
    
           TCanvas *pe_ias_mean_0_cutih = new TCanvas("pe_ias_mean_PU_bin_1_cutih","pe_ias_mean_PU_bin_1_cutih");
           auto pe_ias_graph_0_cutih = new TGraphErrors(nPE,bin_pe,ias_mean_0_cutih,err_mean_x_base,err_ias_mean_0_cutih);
           pe_ias_graph_0_cutih->SetMarkerStyle(21);
           pe_ias_graph_0_cutih->Draw("APL");
           pe_ias_mean_0_cutih->Write();
    
           TCanvas *pe_ias_mean_1_cutih = new TCanvas("pe_ias_mean_PU_bin_2_cutih","pe_ias_mean_PU_bin_2_cutih");
           auto pe_ias_graph_1_cutih = new TGraphErrors(nPE,bin_pe,ias_mean_1_cutih,err_mean_x_base,err_ias_mean_1_cutih);
           pe_ias_graph_1_cutih->SetMarkerStyle(21);
           pe_ias_graph_1_cutih->Draw("APL");
           pe_ias_mean_1_cutih->Write();
           
           TCanvas *pe_ias_mean_2_cutih = new TCanvas("pe_ias_mean_PU_bin_3_cutih","pe_ias_mean_PU_bin_3_cutih");
           auto pe_ias_graph_2_cutih = new TGraphErrors(nPE,bin_pe,ias_mean_2_cutih,err_mean_x_base,err_ias_mean_2_cutih);
           pe_ias_graph_2_cutih->SetMarkerStyle(21);
           pe_ias_graph_2_cutih->Draw("APL");
           pe_ias_mean_2_cutih->Write();
    
           TCanvas *pe_ias_mean_3_cutih = new TCanvas("pe_ias_mean_PU_bin_4_cutih","pe_ias_mean_PU_bin_4_cutih");
           auto pe_ias_graph_3_cutih = new TGraphErrors(nPE,bin_pe,ias_mean_3_cutih,err_mean_x_base,err_ias_mean_3_cutih);
           pe_ias_graph_3_cutih->SetMarkerStyle(21);
           pe_ias_graph_3_cutih->Draw("APL");
           pe_ias_mean_3_cutih->Write();
    
    
           TCanvas *pe_ias_mean_4_cutih = new TCanvas("pe_ias_mean_PU_bin_5_cutih","pe_ias_mean_PU_bin_5_cutih");
           auto pe_ias_graph_4_cutih = new TGraphErrors(nPE,bin_pe,ias_mean_4_cutih,err_mean_x_base,err_ias_mean_4_cutih);
           pe_ias_graph_4_cutih->SetMarkerStyle(21);
           pe_ias_graph_4_cutih->Draw("APL");
           pe_ias_mean_4_cutih->Write();
       }

       NHSCP->Write();
       NTRK->Write();
       N_CLU_HSCP->Write();
       N_CLU_TRK->Write();
       std::cout << "after all PE canvases" << endl;

       IhBestvsIas_p_5_100->Write();
       IhBestvsIas_p_10_45->Write();
       IasStripVsIh0noL1_p_10_45->Write();

       IhBestvsP_10_45_central->Write();
       IhBestvsP_10_45_tight->Write();

       TCanvas *cstack = new TCanvas("c1","stacked + normalized distribution"); 

       THStack *hs = new THStack("hs","");

       ias_bin_pt[0]->SetFillColor(2);
       ias_bin_pt[1]->SetFillColor(41);
       ias_bin_pt[2]->SetFillColor(43);
       ias_bin_pt[3]->SetFillColor(45);
       ias_bin_pt[4]->SetFillColor(46);
       ias_bin_pt[5]->SetFillColor(30);
       ias_bin_pt[6]->SetFillColor(33);
       ias_bin_pt[7]->SetFillColor(40);
       ias_bin_pt[8]->SetFillColor(38);
       ias_bin_pt[9]->SetFillColor(4);
    
    
    
       for(int i=0;i<10;i++){
           hs->Add(ias_bin_pt[i]);
       }
       hs->Draw("pfc nostack");
       cstack->Write();
    
       TCanvas *c1 = new TCanvas("stacked+hist","normalised distribution"); 
       // ***** plot all normalized IAS distribution ******
       
       
    
       ias_bin_pt[0]->SetLineColor(2);
       ias_bin_pt[0]->SetLineWidth(2);
       ias_bin_pt[0]->SetMarkerStyle(32);
    
       ias_bin_pt[0]->Draw("HIST");
    
    
       ias_bin_pt[1]->SetLineColor(41);
       ias_bin_pt[1]->SetLineWidth(2);
       ias_bin_pt[1]->SetMarkerStyle(9);
    
       ias_bin_pt[2]->SetLineColor(43);
       ias_bin_pt[2]->SetLineWidth(2);
       ias_bin_pt[2]->SetMarkerStyle(2);
    
       ias_bin_pt[3]->SetLineColor(45);
       ias_bin_pt[3]->SetLineWidth(2);
       ias_bin_pt[3]->SetMarkerStyle(4);
    
       ias_bin_pt[4]->SetLineColor(46);
       ias_bin_pt[4]->SetLineWidth(2);
       ias_bin_pt[4]->SetMarkerStyle(5);
    
       ias_bin_pt[5]->SetLineColor(30);
       ias_bin_pt[5]->SetLineWidth(2);
       ias_bin_pt[5]->SetMarkerStyle(34);
    
       ias_bin_pt[6]->SetLineColor(33);
       ias_bin_pt[6]->SetLineWidth(2);
       ias_bin_pt[6]->SetMarkerStyle(49);
    
       ias_bin_pt[7]->SetLineColor(40);
       ias_bin_pt[7]->SetLineWidth(2);
       ias_bin_pt[7]->SetMarkerStyle(27);
    
       ias_bin_pt[8]->SetLineColor(38);
       ias_bin_pt[8]->SetLineWidth(2);
       ias_bin_pt[8]->SetMarkerStyle(22);
    
       ias_bin_pt[9]->SetLineColor(4);
       ias_bin_pt[9]->SetLineWidth(2);
       ias_bin_pt[9]->SetMarkerStyle(23);
    
       if(UsePURwtHSCP || UseTemplatesForPUReweighting){
           for(int i=0;i<10;i++){ 
               itg_norm[i] = (1./ias_bin_pt[i]->Integral());
               ias_bin_pt[i]->Write(); //write raw hist and then normalize it
               ias_bin_pt[i]->Scale(itg_norm[i]);//normalization
        
               itg_norm_PU[i] = (1./ias_bin_pt_PUcorr[i]->Integral());
               ias_bin_pt_PUcorr[i]->Write();
               ias_bin_pt_PUcorr[i]->Scale(itg_norm_PU[i]);//normalization 
               if(i != 0) ias_bin_pt[i]->Draw("SAME HIST");
           }
       } 
    
    
       auto leg = new TLegend(0.2, 0.6, 0.35, 0.8);
       leg->AddEntry(ias_bin_pt[0],"PT_60_to_80","l");
       leg->AddEntry(ias_bin_pt[1],"PT_80_to_100","l");
       leg->AddEntry(ias_bin_pt[2],"PT_100_to_120","l");
       leg->AddEntry(ias_bin_pt[3],"PT_120_to_140","l");
       leg->AddEntry(ias_bin_pt[4],"PT_140_to_160","l");
       leg->AddEntry(ias_bin_pt[5],"PT_160_to_180","l");
       leg->AddEntry(ias_bin_pt[6],"PT_180_to_200","l");
       leg->AddEntry(ias_bin_pt[7],"PT_200_to_220","l");
       leg->AddEntry(ias_bin_pt[8],"PT_220_to_240","l");
       leg->AddEntry(ias_bin_pt[9],"PT_240_to_260","l");
       leg->Draw();
    
       c1->Write();
       IAS_VS_PT->Write();
       eff_ih_PU->Write();
    
       /*
       Ias_when_PU_ih_cut_triple[0]->SetMarkerStyle(32);
       Ias_when_PU[0]->SetMarkerStyle(32);
       Ias_when_PU_triple[0]->SetMarkerStyle(32);
       Ias_when_PU_ih_cut[0]->SetMarkerStyle(32);

       Ias_when_PU_ih_cut_triple[1]->SetMarkerStyle(2);
       Ias_when_PU[1]->SetMarkerStyle(2);
       Ias_when_PU_triple[1]->SetMarkerStyle(2);
       Ias_when_PU_ih_cut[1]->SetMarkerStyle(2);

       Ias_when_PU_ih_cut_triple[2]->SetMarkerStyle(5);
       Ias_when_PU[2]->SetMarkerStyle(5);
       Ias_when_PU_triple[2]->SetMarkerStyle(5);
       Ias_when_PU_ih_cut[2]->SetMarkerStyle(5);

       Ias_when_PU_ih_cut_triple[3]->SetMarkerStyle(49);
       Ias_when_PU[3]->SetMarkerStyle(49);
       Ias_when_PU_triple[3]->SetMarkerStyle(49);
       Ias_when_PU_ih_cut[3]->SetMarkerStyle(49);


       Ias_when_PU_ih_cut_triple[4]->SetMarkerStyle(27);
       Ias_when_PU[4]->SetMarkerStyle(27);
       Ias_when_PU_triple[4]->SetMarkerStyle(27);
       Ias_when_PU_ih_cut[4]->SetMarkerStyle(27);
       */

       IhStripVsP->Write();
       IhStripVsPtight->Write();
       IhBestStripVsP->Write();
       IhBestStripVsPtight->Write();

       IhStripVsP_presk->Write();
       IhStripVsPtight_presk->Write();

       IhBestStripVsP_presk->Write();
       IhBestStripVsPtight_presk->Write();

       PT_compute_ias->Write();
       P_compute_ias->Write();
       PT_hscp_compute_ias->Write();
       P_hscp_compute_ias->Write();



       Ias_1_VS_Delta_Ias_1_2->Write();


       Ias_vs_PU_above_40_triple->Write();

       Ias_all_triple_cutIH->Write();
       Ias_all_base_cutIH->Write();  
       Ias_all_triple_nocut->Write();
       Ias_all_base_nocut->Write();

       int prob_pt50 = 3;
       double p_pt50_single[3] = {0.8 ,0.9, 0.99}; double q_pt50_single[3];
       double p_pt50_triple[3] = {0.8 ,0.9, 0.99}; double q_pt50_triple[3];
       double pt_qt[3] = {80,90,99}; 
       IAS_single_pt50->GetQuantiles(prob_pt50,q_pt50_single,p_pt50_single);    
       IAS_triple_pt50->GetQuantiles(prob_pt50,q_pt50_triple,p_pt50_triple);    
         
       TCanvas *ias_qtl_pt50_single = new TCanvas("qtl_pt50_single","Ias 80-90-99% quantiles single pt > 50"); 
       auto ias_3_qtl_single = new TGraph(prob_pt50,pt_qt,q_pt50_single);

       auto ias_3_qtl_triple = new TGraph(prob_pt50,pt_qt,q_pt50_triple);


       ias_3_qtl_single->Draw("AC*");
       ias_3_qtl_single->SetLineColor(4);

       ias_3_qtl_triple->SetLineColor(2);
       ias_3_qtl_triple->Draw("SAME C*");
       ias_qtl_pt50_single->Write();

       for (int i = 0 ; i < ias_intervals ; i++){
           int prob_qtl = 3;
           double p[3] = {0.8 ,0.9, 0.99}; double q[3];
           Ias_when_PU[i]->GetQuantiles(prob_qtl,q,p);

           quantiles_PU_80[i] = q[0];
           quantiles_PU_90[i] = q[1];
           quantiles_PU_99[i] = q[2];

           double p_triple[3] = {0.8 , 0.9, 0.99 }; double q_triple[3];
           Ias_when_PU_triple[i]->GetQuantiles(prob_qtl,q_triple,p_triple);

           quantiles_PU_triple_80[i] = q_triple[0];
           quantiles_PU_triple_90[i] = q_triple[1];
           quantiles_PU_triple_99[i] = q_triple[2];
          
           double p_cut[3] = {0.8 ,0.9, 0.99}; double q_cut[3];
           Ias_when_PU_ih_cut[i]->GetQuantiles(prob_qtl,q_cut,p_cut);

           quantiles_PU_cut_80[i] = q_cut[0];
           quantiles_PU_cut_90[i] = q_cut[1];
           quantiles_PU_cut_99[i] = q_cut[2];
           
           double p_triple_cut[3] = {0.8 ,0.9,0.99}; double q_triple_cut[3];
           Ias_when_PU_ih_cut_triple[i]->GetQuantiles(prob_qtl,q_triple_cut,p_triple_cut);

           quantiles_PU_triple_cut_80[i] = q_triple_cut[0];
           quantiles_PU_triple_cut_90[i] = q_triple_cut[1];
           quantiles_PU_triple_cut_99[i] = q_triple_cut[2];


           Ias_when_PU[i]->Write();
           Ias_when_PU_ih_cut[i]->Write();
           Ias_when_PU_ih_cut_triple[i]->Write();
           Ias_when_PU_triple[i]->Write();
           Is_when_PU[i]->Write();
           Is_when_PU_ih_cut[i]->Write();
           Is_when_PU_ih_cut_triple[i]->Write();
           Is_when_PU_triple[i]->Write();

           PU_distrib[i]->Write();
       }



       TCanvas *ias_qtl_pu_80 = new TCanvas("80qtl","Ias 80% quantiles PU scenarios"); 
       auto ias_graph_qtl_80 = new TGraph(ias_intervals,ias_bins,quantiles_PU_80);
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_qtl_80->GetXaxis()->SetBinLabel(ias_graph_qtl_80->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }

       auto ias_graph_qtl_triple_80 = new TGraph(ias_intervals,ias_bins,quantiles_PU_triple_80);
       ias_graph_qtl_80->Draw("AC*");
       ias_graph_qtl_80->SetLineColor(4);

       ias_graph_qtl_triple_80->SetLineColor(2);
       ias_graph_qtl_triple_80->Draw("SAME C*");
       ias_qtl_pu_80->Write();
        
       TCanvas *ias_qtl_pu_cut_80 = new TCanvas("80qtl_cutih","Ias 80% quantiles PU scenarios + cut Ih > 3.47");
       auto ias_graph_qtl_cut_80 = new TGraph(ias_intervals,ias_bins,quantiles_PU_cut_80);
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_qtl_cut_80->GetXaxis()->SetBinLabel(ias_graph_qtl_80->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }

       auto ias_graph_qtl_cut_triple_80 = new TGraph(ias_intervals,ias_bins,quantiles_PU_triple_cut_80);
       ias_graph_qtl_cut_80->Draw("AC*");
       ias_graph_qtl_cut_80->SetLineColor(4);

       ias_graph_qtl_cut_triple_80->SetLineColor(2);
       ias_graph_qtl_cut_triple_80->Draw("SAME C*");
       ias_qtl_pu_cut_80->Write();  



       TCanvas *ias_qtl_pu = new TCanvas("90qtl","Ias 90% quantiles PU scenarios"); 
       auto ias_graph_qtl = new TGraph(ias_intervals,ias_bins,quantiles_PU_90);
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_qtl->GetXaxis()->SetBinLabel(ias_graph_qtl->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }

       auto ias_graph_qtl_triple = new TGraph(ias_intervals,ias_bins,quantiles_PU_triple_90);
       ias_graph_qtl->Draw("AC*");
       ias_graph_qtl->SetLineColor(4);

       ias_graph_qtl_triple->SetLineColor(2);
       ias_graph_qtl_triple->Draw("SAME C*");
       ias_qtl_pu->Write();

       TCanvas *ias_qtl_pu_cut = new TCanvas("90qtl_cutih","Ias 90% quantiles PU scenarios + cut Ih > 3.47");
       auto ias_graph_qtl_cut = new TGraph(ias_intervals,ias_bins,quantiles_PU_cut_90);
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_qtl_cut->GetXaxis()->SetBinLabel(ias_graph_qtl->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }

       auto ias_graph_qtl_cut_triple = new TGraph(ias_intervals,ias_bins,quantiles_PU_triple_cut_90);
       ias_graph_qtl_cut->Draw("AC*");
       ias_graph_qtl_cut->SetLineColor(4);

       ias_graph_qtl_cut_triple->SetLineColor(2);
       ias_graph_qtl_cut_triple->Draw("SAME C*");
       ias_qtl_pu_cut->Write();  


       TCanvas *ias_qtl_pu_99 = new TCanvas("99qtl","Ias 99% quantiles PU scenarios"); 
       auto ias_graph_qtl_99 = new TGraph(ias_intervals,ias_bins,quantiles_PU_99);
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_qtl_99->GetXaxis()->SetBinLabel(ias_graph_qtl->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }

       auto ias_graph_qtl_triple_99 = new TGraph(ias_intervals,ias_bins,quantiles_PU_triple_99);
       ias_graph_qtl_99->Draw("AC*");
       ias_graph_qtl_99->SetLineColor(4);

       ias_graph_qtl_triple_99->SetLineColor(2);
       ias_graph_qtl_triple_99->Draw("SAME C*");
       ias_qtl_pu_99->Write();

       TCanvas *ias_qtl_pu_cut_99 = new TCanvas("99qtl_cut","Ias 99% quantiles PU scenarios + cut Ih > 3.47");
       auto ias_graph_qtl_cut_99 = new TGraph(ias_intervals,ias_bins,quantiles_PU_cut_99);
       for (int u = 0 ; u < ias_intervals; u++){
          ias_graph_qtl_cut_99->GetXaxis()->SetBinLabel(ias_graph_qtl->GetXaxis()->FindBin(ias_bins[u]),ias_axis_names[u].c_str());
       }
       auto ias_graph_qtl_cut_triple_99 = new TGraph(ias_intervals,ias_bins,quantiles_PU_triple_cut_99);
       ias_graph_qtl_cut_99->Draw("AC*");
       ias_graph_qtl_cut_99->SetLineColor(4);

       ias_graph_qtl_cut_triple_99->SetLineColor(2);
       ias_graph_qtl_cut_triple_99->Draw("SAME C*");
       ias_qtl_pu_cut_99->Write();  

       for (int i = 0; i < ias_intervals ; i++){
           int binmin = Ias_when_PU[i]->GetXaxis()->FindBin(0.115);
           int binmax = Ias_when_PU[i]->GetXaxis()->FindBin(1);
           double itg = Ias_when_PU[i]->Integral(binmin,binmax);
           double itg2 = Ias_when_PU[i]->Integral();
           ratio_above_qtl[i] = (itg*1.0)/itg2;

           //

           Ias_when_PU[i]->Scale(1./Ias_when_PU[i]->Integral());   
           Ias_when_PU[i]->SetLineColor(i+2);
           Ias_when_PU[i]->SetLineWidth(2);

           int binmin_triple = Ias_when_PU_triple[i]->GetXaxis()->FindBin(0.115);
           int binmax_triple = Ias_when_PU_triple[i]->GetXaxis()->FindBin(1);
           double itg_triple = Ias_when_PU_triple[i]->Integral(binmin_triple,binmax_triple);
           double itg2_triple = Ias_when_PU_triple[i]->Integral();
           ratio_above_qtl_triple[i] = (itg_triple*1.0)/itg2_triple;

           Ias_when_PU_triple[i]->Scale(1./Ias_when_PU_triple[i]->Integral());
           Ias_when_PU_triple[i]->SetLineColor(i+2);
           Ias_when_PU_triple[i]->SetLineWidth(2);

           int binmin_ih_cut = Ias_when_PU_ih_cut[i]->GetXaxis()->FindBin(0.115);
           int binmax_ih_cut = Ias_when_PU_ih_cut[i]->GetXaxis()->FindBin(1);
           double itg_ih_cut = Ias_when_PU_ih_cut[i]->Integral(binmin_ih_cut,binmax_ih_cut);
           double itg2_ih_cut = Ias_when_PU_ih_cut[i]->Integral();
           ratio_above_qtl_cut[i] = (itg_ih_cut*1.0)/itg2_ih_cut;

           Ias_when_PU_ih_cut[i]->Scale(1./Ias_when_PU_ih_cut[i]->Integral());
           Ias_when_PU_ih_cut[i]->SetLineColor(i+2);
           Ias_when_PU_ih_cut[i]->SetLineWidth(2);

           int binmin_ih_cut_triple = Ias_when_PU_ih_cut_triple[i]->GetXaxis()->FindBin(0.115);
           int binmax_ih_cut_triple = Ias_when_PU_ih_cut_triple[i]->GetXaxis()->FindBin(1);
           double itg_ih_cut_triple = Ias_when_PU_ih_cut_triple[i]->Integral(binmin_ih_cut_triple,binmax_ih_cut_triple);
           double itg2_ih_cut_triple = Ias_when_PU_ih_cut_triple[i]->Integral();
           ratio_above_qtl_cut_triple[i] = (itg_ih_cut_triple*1.0)/itg2_ih_cut_triple;

           Ias_when_PU_ih_cut_triple[i]->Scale(1./Ias_when_PU_ih_cut_triple[i]->Integral());
           Ias_when_PU_ih_cut_triple[i]->SetLineColor(i+2);
           Ias_when_PU_ih_cut_triple[i]->SetLineWidth(2);
       }

       cout << " -------------- RATIOS OF INTEGRALS --------------" << endl;
       for (int i = 0 ; i < ias_intervals ; i++){
           cout << " NORMAL NO CUT : \n " << " PU bin " << i << " has ratio of integrals after 90% qtl = " << ratio_above_qtl[i] << endl;
           cout << " TRIPLE NO CUT : \n " << " PU bin " << i << " has ratio of integrals after 90% qtl = " << ratio_above_qtl_triple[i] << endl;
           cout << " NORMAL CUT IH: \n " << " PU bin " << i << " has ratio of integrals after 90% qtl = " << ratio_above_qtl_cut[i] << endl;
           cout << " TRIPLE CUT IH : \n " << " PU bin " << i << " has ratio of integrals after 90% qtl = " << ratio_above_qtl_cut_triple[i] << endl;
       }
       // ***** plot all normalized IAS distribution ******
       TCanvas *ias_pu_normal = new TCanvas("IAS_single","IAS from single template"); 
       Ias_when_PU[0]->Draw("C*");
       for (int k = 1; k < ias_intervals ; k++){
           Ias_when_PU[k]->Draw("SAME C*"); 
       }
       
       auto leg_norm = new TLegend(0.2, 0.6, 0.35, 0.8);
       for (int u = 0 ; u < ias_intervals ; u++){
           leg_norm->AddEntry(Ias_when_PU[u],names_pu[u].c_str(),"l"); 
       }
       leg_norm->Draw();
       ias_pu_normal->Write();

       TCanvas *ias_pu_normal_ratio = new TCanvas("IAS_single_RATIO","RATIO IAS from single template"); 
       auto rp = new TRatioPlot(Ias_when_PU[2], Ias_when_PU[0]);
       ias_pu_normal_ratio->SetTicks(0,1);

       rp->Draw();
       ias_pu_normal_ratio->Write(); 

       TCanvas *ias_pu_triple = new TCanvas("IAS_triple","IAS from multiple templates"); 
       Ias_when_PU_triple[0]->Draw("C*");
       for (int k = 1; k < ias_intervals ; k++){
           Ias_when_PU_triple[k]->Draw("SAME C*");
       }
       auto leg_triple = new TLegend(0.2, 0.6, 0.35, 0.8);
       for (int u = 0 ; u < ias_intervals ; u++){
           leg_triple->AddEntry(Ias_when_PU_triple[u],names_pu[u].c_str(),"l"); 
       }
       leg_triple->Draw();
       ias_pu_triple->Write();


       TCanvas *ias_pu_triple_ratio = new TCanvas("IAS_triple_RATIO","RATIO IAS from multiple templates"); 
       auto rp2 = new TRatioPlot(Ias_when_PU_triple[2], Ias_when_PU_triple[0]);
       ias_pu_triple_ratio->SetTicks(0,1);

       rp2->Draw();
       ias_pu_triple_ratio->Write();

       TCanvas *ias_pu_normal_cut = new TCanvas("IAS_single_cut_IH","IAS from single template, cut IH");  
       Ias_when_PU_ih_cut[0]->Draw("C*");
       for (int k = 1; k < ias_intervals ; k++){
           Ias_when_PU_ih_cut[k]->Draw("SAME C*");
       }
       auto leg_norm_cut = new TLegend(0.2, 0.6, 0.35, 0.8);
       for (int u = 0 ; u < ias_intervals ; u++){
           leg_norm_cut->AddEntry(Ias_when_PU_ih_cut[u],names_pu[u].c_str(),"l"); 
       }
       leg_norm_cut->Draw();
       ias_pu_normal_cut->Write();

       
       TCanvas *ias_pu_normal_cut_ratio = new TCanvas("IAS_single_cut_IH_RATIO","RATIO IAS from single template, cut IH");  
       auto rp3 = new TRatioPlot(Ias_when_PU_ih_cut[2], Ias_when_PU_ih_cut[0]);
       ias_pu_normal_cut_ratio->SetTicks(0,1);

       rp3->Draw();
       ias_pu_normal_cut_ratio->Write();

       TCanvas *ias_pu_triple_cut = new TCanvas("IAS_triple_cut_IH","IAS from multiple templates cut IH"); 
       Ias_when_PU_ih_cut_triple[0]->Draw("C*");
       for (int k = 1; k < ias_intervals ; k++){
           Ias_when_PU_ih_cut_triple[k]->Draw("SAME C*");
       }
       auto leg_triple_cut = new TLegend(0.2, 0.6, 0.35, 0.8);
       for (int u = 0 ; u < ias_intervals ; u++){
           leg_triple_cut->AddEntry(Ias_when_PU_ih_cut_triple[u],names_pu[u].c_str(),"l"); 
       }
       leg_triple_cut->Draw();
       ias_pu_triple_cut->Write();
       TCanvas *ias_pu_triple_cut_ratio = new TCanvas("IAS_triple_cut_IH_RATIO","RATIO IAS from multiple templates cut IH");
       auto rp4 = new TRatioPlot(Ias_when_PU_ih_cut_triple[2], Ias_when_PU_ih_cut_triple[0]);
       ias_pu_triple_cut_ratio->SetTicks(0,1);
       rp4->Draw();
       ias_pu_triple_cut_ratio->Write();

       Ias_toy->Write();
       
       
       PU_VS_NPV->Write();
       NPV_all->Write();
       NPV_presel->Write();

       Ias_below_5GeV_PU_below_20_base->Write();
       Ias_below_5GeV_PU_between_base->Write();
       Ias_below_5GeV_PU_above_40_base->Write(); 
    
       Ias_below_5GeV_PU_below_20_triple->Write();
       Ias_below_5GeV_PU_between_triple->Write();
       Ias_below_5GeV_PU_above_40_triple->Write();

       Ias_when_PU_in_16_18_base->Write();
       Ias_when_PU_in_16_18_triple->Write();

       Ias_when_PU_in_30_32_base->Write();
       Ias_when_PU_in_30_32_triple->Write();

       
        
       PATHLENGHT_BIN_ETA_0_01->Write();
       TptGivenPathlenght1->Write();


       if(compute_PE){ 
           for(int j = 0; j< ias_intervals; j++){
               for (int i = 0; i < nPE ; i++){
                   PE_Ias[j][i]->Write();
                   PE_Ias_cutih[j][i]->Write();
               }
           }
       }
       cout << "Number of events with npv > 99 : " << nb_above_99_npv << endl; 

       cout << " For cut PT > 50, ISO2 TK < 50 (and the frozen preselection) : " << endl; 
       cout << " IAS > 0.1 : " << endl;
       cout <<  "Ih > 3.47, # IAS single > 0.1 " << nb_ias_0p1_single << " , # IAS triple > 0.1 " << nb_ias_0p1_triple <<endl;

       cout << "Ih > 3.29 (val of C), # IAS single > 0.1 : " << nb_ias_0p1_single_cut3p29 << " , # IAS triple > 0.1 " << nb_ias_0p1_triple_cut3p29 <<endl;

       cout << "No cut on Ih, # IAS single > 0.1 : " << nb_ias_0p1_single_nocut << " , # IAS triple > 0.1 " << nb_ias_0p1_triple_nocut <<endl;

      
       cout << " IAS > 0.2 : " << endl;
       cout <<  "Ih > 3.47, # IAS single > 0.2 " << nb_ias_0p2_single << " , # IAS triple > 0.2 " << nb_ias_0p2_triple <<endl;

       cout << "Ih > 3.29 (val of C), # IAS single > 0.2 : " << nb_ias_0p2_single_cut3p29 << " , # IAS triple > 0.2 " << nb_ias_0p2_triple_cut3p29 <<endl;

       cout << "No cut on Ih, # IAS single > 0.2 : " << nb_ias_0p2_single_nocut << " , # IAS triple > 0.2 " << nb_ias_0p2_triple_nocut <<endl;


       cout << " IAS > 0.3 : " << endl;

       cout <<  "Ih > 3.47, # IAS single > 0.3 " << nb_ias_0p3_single << " , # IAS triple > 0.3 " << nb_ias_0p3_triple <<endl;

       cout << "Ih > 3.29 (val of C), # IAS single > 0.3 : " << nb_ias_0p3_single_cut3p29 << " , # IAS triple > 0.3 " << nb_ias_0p3_triple_cut3p29 <<endl;

       cout << "No cut on Ih, # IAS single > 0.3 : " << nb_ias_0p3_single_nocut << " , # IAS triple > 0.3 " << nb_ias_0p3_triple_nocut <<endl;


       cout << " IAS > 0.4 : " << endl;

       cout <<  "Ih > 3.47, # IAS single > 0.4 " << nb_ias_0p4_single << " , # IAS triple > 0.4 " << nb_ias_0p4_triple <<endl;

       cout << "Ih > 3.29 (val of C), # IAS single > 0.4 : " << nb_ias_0p4_single_cut3p29 << " , # IAS triple > 0.4 " << nb_ias_0p4_triple_cut3p29 <<endl;

       cout << "No cut on Ih, # IAS single > 0.4 : " << nb_ias_0p4_single_nocut << " , # IAS triple > 0.4 " << nb_ias_0p4_triple_nocut <<endl;


       OutIasStudy->Close();
  }
  

}


double run2analysis::getMassSpecial(float ih, float p, float K, float C, float N){

 float m1=sqrt((ih-C)*p*p/K);

 char fitfunc[2048];
 sprintf(fitfunc,"%6.4f*pow(x/%6.4f,2) + %6.4f + %6.4f*log(%6.4f/x) -%6.4f",K,p,C,N,p,ih);

 TF1 *fa1 = new TF1("fa1",fitfunc,0.1,m1*3.);

 float evalval = fa1->Eval(m1);
 float mbest=0;
 float ebest=1000;

 for (int i=0;i<2000;i++) {
   if (evalval<0) {
     float m2=m1+i*0.0001*m1;
     float evalval2 = fa1->Eval(m2);
     if (fabs(evalval2)<0.001) {
          if (fabs(evalval2)<ebest) {
                      mbest=m2;
                      ebest=fabs(evalval2);
          }
     }
     else if (evalval2>0.001) i+=2000; // this should allow to stop the search of solution
   }
   else {
     float m2=m1-i*0.0001*m1;
     if (m2>0) {
           float evalval2 = fa1->Eval(m2);
           if (fabs(evalval2)<0.001)  {
                   if (fabs(evalval2)<ebest) {
                         mbest=m2;
                         ebest=fabs(evalval2);
                    }
           }
           else if (evalval2<-0.001) i+=2000; // this should allow to stop the search of solution
      }
   }
  }

  return mbest;
}

double run2analysis::getdEdXIs(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, double dropHigherDeDxValue, int & nv, int & ns) {
  double result=-1;
//     double dropLowerDeDxValue=0.15;
     size_t MaxStripNOM=99;
     bool usePixel=true;
     bool useStrip=true;

     std::vector<double> vect;

     bool debugprint=false;
     unsigned int SiStripNOM = 0;
     ns=0;

     for(unsigned int h=0;h<charge.size();h++){
        if (debugprint) std::cout << "look on dedxHits in computedEdx " << h << std::endl;
        if(!usePixel && subdetId[h]<3)continue; // skip pixels
        if(!useStrip && subdetId[h]>=3)continue; // skip strips        
        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;

        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;
        if(useStrip && subdetId[h]>=3 && !mustBeInside[h])continue;
        if(useStrip && subdetId[h]>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

        int ClusterCharge = charge[h];
        if (subdetId[h]>=3 && charge[h]>=254) ns++;

        double scaleFactor = scaleFactors[0];
        if (subdetId[h]<3) scaleFactor *= scaleFactors[1]; // add pixel scaling
        if (debugprint) std::cout << " after SF " << std::endl;

        if(templateHisto){  //save discriminator probability
           double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?265:1));
           int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry[h]);
           int    BinY   = templateHisto->GetYaxis()->FindBin(pathlength[h]*10.0); //*10 because of cm-->mm
           int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
           double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
           vect.push_back(Prob); //save probability
           if (debugprint) std::cout << " after Prob vect.push_back " << std::endl;
        }else{
           double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
           double ChargeOverPathlength = scaleFactor*Norm*ClusterCharge/pathlength[h];
           vect.push_back(ChargeOverPathlength); //save charge
           if (debugprint) std::cout << " after ChargeOverPathlength vect.push_back " << std::endl;
        }
     }

     if(dropLowerDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
         int nTrunc = tmp.size()*dropLowerDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropLowerDeDxValue " << std::endl;

     if(dropHigherDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::less<double>() );
         int nTrunc = tmp.size()*dropHigherDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropHigherDeDxValue " << std::endl;






     int size = vect.size();
     nv = size;

     if(size>0){
        if(templateHisto){  //dEdx discriminator
          //Ias discriminator
          result = 1.0/(12*size);
           std::sort(vect.begin(), vect.end(), std::less<double>() );
           for(int i=1;i<=size;i++){
              result += pow(vect[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
           }
           result *= (3.0/size);
           if (debugprint) std::cout << " Ias discriminator " << result << std::endl;
        }else{  //dEdx estimator
           //harmonic2 estimator
           result=0;
//           double expo = -2;
           double expo = -1* n_estim;
           for(int i = 0; i< size; i ++){
              result+=pow(vect[i],expo);
           }
           result = pow(result/size,1./expo);
           if (debugprint) std::cout << " harmonic discriminator " << result << " with expo " << expo << std::endl;
        }
     }else{
        result = -1;
     }
     if (debugprint) std::cout << " ok finished computeDeDx " << std::endl;


  return result;
}

double run2analysis::getdEdXIs(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, int & nv, int & ns) {
  double result= getdEdXIs(charge, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, scaleFactors, templateHisto, n_estim, dropLowerDeDxValue, 0., nv, ns);
  return result;
}

double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, double dropHigherDeDxValue, int & nv, int & ns) {
  double result=-1;
//     double dropLowerDeDxValue=0.15;
     size_t MaxStripNOM=99;
     bool usePixel=true;
     bool useStrip=true;

     std::vector<double> vect;

     bool debugprint=false;
     unsigned int SiStripNOM = 0;
     ns=0;

     for(unsigned int h=0;h<charge.size();h++){
        if (debugprint) std::cout << "look on dedxHits in computedEdx " << h << std::endl;
        if(!usePixel && subdetId[h]<3)continue; // skip pixels
        if(!useStrip && subdetId[h]>=3)continue; // skip strips        
        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;

        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;
        if(useStrip && subdetId[h]>=3 && !mustBeInside[h])continue;
        if(useStrip && subdetId[h]>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

        int ClusterCharge = charge[h];
        if (subdetId[h]>=3 && charge[h]>=254) ns++;

        double scaleFactor = scaleFactors[0];
        if (subdetId[h]<3) scaleFactor *= scaleFactors[1]; // add pixel scaling
        if (debugprint) std::cout << " after SF " << std::endl;

        if(templateHisto){  //save discriminator probability
           double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?265:1));
           int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry[h]);
           int    BinY   = templateHisto->GetYaxis()->FindBin(pathlength[h]*10.0); //*10 because of cm-->mm
           int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
           double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
           vect.push_back(Prob); //save probability
           if (debugprint) std::cout << " after Prob vect.push_back " << std::endl;
        }else{
           double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
           double ChargeOverPathlength = scaleFactor*Norm*ClusterCharge/pathlength[h];
           vect.push_back(ChargeOverPathlength); //save charge
           if (debugprint) std::cout << " after ChargeOverPathlength vect.push_back " << std::endl;
        }
     }

     if(dropLowerDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
         int nTrunc = tmp.size()*dropLowerDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropLowerDeDxValue " << std::endl;

     if(dropHigherDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::less<double>() );
         int nTrunc = tmp.size()*dropHigherDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropHigherDeDxValue " << std::endl;






     int size = vect.size();
     nv = size;

     if(size>0){
        if(templateHisto){  //dEdx discriminator
          //Ias discriminator
          result = 1.0/(12*size);
           std::sort(vect.begin(), vect.end(), std::less<double>() );
           for(int i=1;i<=size;i++){
              result += vect[i-1] * pow(vect[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
           }
           result *= (3.0/size);
           if (debugprint) std::cout << " Ias discriminator " << result << std::endl;
        }else{  //dEdx estimator
           //harmonic2 estimator
           result=0;
//           double expo = -2;
           double expo = -1* n_estim;
           for(int i = 0; i< size; i ++){
              result+=pow(vect[i],expo);
           }
           result = pow(result/size,1./expo);
           if (debugprint) std::cout << " harmonic discriminator " << result << " with expo " << expo << std::endl;
        }
     }else{
        result = -1;
     }
     if (debugprint) std::cout << " ok finished computeDeDx " << std::endl;


  return result;
}

double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, int & nv, int & ns) {
  double result= getdEdX(charge, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, scaleFactors, templateHisto, n_estim, dropLowerDeDxValue, 0., nv, ns);
  return result;
}

//float run2analysis::FMIP(const vector<float>& charge, const<vector>& path, float thre = 4){
float run2analysis::FMIP(std::vector <float> charge, std::vector <float> path, float SFactor, float thre = 4){
   if(charge.size()!=path.size()) return -999; // error

   int nclusters = charge.size();
   if(nclusters==0) return -888;

   int nlow = 0; 
   for(unsigned int i=0;i<charge.size();i++){
       //hard-coded conversion factor to go for dEdx in MeV.cm2/g
       // ATTENTION : HERE THE CONVERSION FACTOR IS 247 WHILE ABOVE WE USE 265 !!!!!
       float dEdx = SFactor*charge[i]*(3.61*pow(10,-9)*247)*1000/path[i];
       if(dEdx<thre) nlow++;
   }
   return nlow*1./nclusters;
}


int run2analysis::GetLayerLabel(int subdetid_, UInt_t detid_, int year)
{
// from https://github.com/dapparu/HSCP/blob/c69cf1c71dd99289f72ab6d03077915776c85690/src/Cluster.cc
// and https://cmssdt.cern.ch/lxr/source/DataFormats/SiPixelDetId/interface/PXBDetId.h
// https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder
        if(subdetid_==1)
        {
             if (year==2016) {
                if(((detid_>>16)&0xF)==1) return 23;
                else if(((detid_>>16)&0xF)==2) return 24;
                else if(((detid_>>16)&0xF)==3) return 25;
                else if(((detid_>>16)&0xF)==4) return 26;  // do not exist in 2016
             }
             else {
                if(((detid_>>20)&0xF)==1) return 23;
                else if(((detid_>>20)&0xF)==2) return 24;
                else if(((detid_>>20)&0xF)==3) return 25;
                else if(((detid_>>20)&0xF)==4) return 26;
             }

        }
        else if(subdetid_==2)
        {
             if (year==2016) {
                if(((detid_>>16)&0xF)==1) return 27;
                else if(((detid_>>16)&0xF)==2) return 28;
                else if(((detid_>>16)&0xF)==3) return 29; // do not exist in 2016
             }
             else {
                if(((detid_>>18)&0xF)==1) return 27;
                else if(((detid_>>18)&0xF)==2) return 28;
                else if(((detid_>>18)&0xF)==3) return 29;
             }
/*  https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder for 2016
                cout << "  side " << int((detid_>>23)&0x3) <<
                        "  disk " << int((detid_>>16)&0xF) << 
                        "  blade "  << int((detid_>>10)&0x3F) <<
                        "  panel "  << int((detid_>>8)&0x3) <<
                        "  mod "  << int((detid_>>2)&0x3F) << endl;
*/
//                if (((detid_>>16)&0xF)==0) cout << " disk 0 ? " << std::endl;
        }
        else if(subdetid_==3)  // TIB
        {
                if(((detid_>>14)&0x7)==1) return 1;
                else if(((detid_>>14)&0x7)==2) return 2;
                else if(((detid_>>14)&0x7)==3) return 3;
                else if(((detid_>>14)&0x7)==4) return 4;
        }
        else if(subdetid_==5) // TOB
        {
                if(((detid_>>14)&0x7)==1) return 5;
                else if(((detid_>>14)&0x7)==2) return 6;
                else if(((detid_>>14)&0x7)==3) return 7;
                else if(((detid_>>14)&0x7)==4) return 8;
                else if(((detid_>>14)&0x7)==5) return 9;
                else if(((detid_>>14)&0x7)==6) return 10;
        }
        else if(subdetid_==4)  //TID
        {
                if(((detid_>>11)&0x3)==1) return 11;
                else if(((detid_>>11)&0x3)==2) return 12;
                else if(((detid_>>11)&0x3)==3) return 13;
        }
        else if(subdetid_==6) // TEC
        {
                if(((detid_>>14)&0xF)==1) return 14;
                else if(((detid_>>14)&0xF)==2) return 15;
                else if(((detid_>>14)&0xF)==3) return 16;
                else if(((detid_>>14)&0xF)==4) return 17;
                else if(((detid_>>14)&0xF)==5) return 18;
                else if(((detid_>>14)&0xF)==6) return 19;
                else if(((detid_>>14)&0xF)==7) return 20;
                else if(((detid_>>14)&0xF)==8) return 21;
                else if(((detid_>>14)&0xF)==9) return 22;
        }
        return -1;
}


void run2analysis::loadSFPixelCalib() {

   std::ifstream file_calib;
   file_calib.open("scale_for_cmssw2017.txt");
   icalib=0;
   if (!file_calib) { std::cerr << "cannot open file scale_for_cmssw2017.txt " << std::endl; }
   else
   {
      while (!file_calib.eof () && icalib<calmax) {
       // pix", "layerorside", "ladderorblade", "etaMin", "etaMax", "irunMin", "irunMax", "value
       file_calib >> pixVal[icalib] >> layerSideVal[icalib] >> ladderBladeVal[icalib] >> etaMinVal[icalib] >> etaMaxVal[icalib]
                  >> irunMinVal[icalib] >> irunMaxVal[icalib] >> scaleVal[icalib] ;
       if (icalib<10) { std::cout << pixVal[icalib] << " " << layerSideVal[icalib]  << " "  <<  ladderBladeVal[icalib] ; 
           std::cout << " " << etaMinVal[icalib] << " "  << etaMaxVal[icalib] << " " << irunMinVal[icalib] << " " << irunMaxVal[icalib] ;
           std::cout << " " << scaleVal[icalib] << std::endl ;
       }
       icalib++;
      }
   }
   std::cout << " file_calib : " << icalib << " entries in 2017" << std::endl;
   if (icalib>0) std::cout << " example runMin "<< irunMinVal[0] << " runMax "<< irunMaxVal[0] << " pix " << pixVal[0] << 
                " scale " << scaleVal[0] <<  std::endl;
   file_calib.close ();

   icalib2017=icalib;
   std::ifstream file_calib2;
   file_calib2.open("scale_for_cmssw2018.txt");
   if (!file_calib2) { std::cerr << "cannot open file scale_for_cmssw2018.txt " << std::endl; }
   else
   {  
      while (!file_calib2.eof () && icalib<calmax) {
       // pix", "layerorside", "ladderorblade", "etaMin", "etaMax", "irunMin", "irunMax", "value
       file_calib2 >> pixVal[icalib] >> layerSideVal[icalib] >> ladderBladeVal[icalib] >> etaMinVal[icalib] >> etaMaxVal[icalib]
                  >> irunMinVal[icalib] >> irunMaxVal[icalib] >> scaleVal[icalib] ;
       icalib++;
      }
   }
   icalib2018=icalib-icalib2017;
   std::cout << " file_calib : " << icalib2018 << " entries in 2018" << std::endl;
   if (icalib>icalib2017) std::cout << " example runMin "<< irunMinVal[icalib2017] << " runMax "<< irunMaxVal[icalib2017] 
                 << " pix " << pixVal[icalib2017] << " scale " << scaleVal[icalib2017] <<  std::endl;
   file_calib2.close ();
 



}

float run2analysis::GetSFPixel(int subdetid_, UInt_t detid_, int year, float eta, int run) {
  // https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder

   int pix = 0;
   int layerorside = 0;
   int ladderorblade = 0;

   if (subdetid_==1) {
        pix =1 ;
        if (year==2016) {
                layerorside = int((detid_>>16)&0xF);
                ladderorblade = int((detid_>>8)&0xFF);
        }
        else {
                layerorside = int((detid_>>20)&0xF);
                ladderorblade = int((detid_>>12)&0xFF);
        }
   }
   if (subdetid_==2) {
        pix =2 ;
        if (year==2016) {
                layerorside = int((detid_>>16)&0xF);
                ladderorblade = int((detid_>>10)&0x3F);
        }
        else {
                layerorside = int((detid_>>18)&0xF);
                ladderorblade = int((detid_>>12)&0x3F);
        }
   }

   float scale =1 ;

   if (year==2017) {
    for (int i=0; i<icalib2017; i++) {
     if (run >= irunMinVal[i] && run<irunMaxVal[i]) {
       if (eta >= etaMinVal[i] && eta < etaMaxVal[i]) { 
          if (pix == pixVal[i] ) {
             if (layerorside == layerSideVal[i]) {
               if (ladderorblade == ladderBladeVal[i]) {
                  scale = scaleVal[i];
                  return scale;
               }
             }
          }
       }
     }
    }
   }
   else if (year==2018) {
    for (int i=icalib2017; i<icalib; i++) {
     if (run >= irunMinVal[i] && run<irunMaxVal[i]) {
       if (eta >= etaMinVal[i] && eta < etaMaxVal[i]) { 
          if (pix == pixVal[i] ) {
             if (layerorside == layerSideVal[i]) {
               if (ladderorblade == ladderBladeVal[i]) {
                  scale = scaleVal[i];
                  return scale;
               }
             }
          }
       }
     }
    }
   }

   return scale;
   
}

void run2analysis::loadSFPixelTamas() {

   std::ifstream file_calib1;
   file_calib1.open("Tamas/CorrFactHistory2017L1.txt");
   icalibL1=0;
   if (!file_calib1) { std::cerr << "cannot open file CorrFactHistory2017L1.txt " << std::endl; }
   else
   {
      while (!file_calib1.eof () && icalibL1<calmax) {
       // from_run_number, corr_val, error_val
       file_calib1 >> irunMinValL1[icalibL1] >> scaleValL1[icalibL1] >> errorScaleValL1[icalibL1] ;
//       if (icalibL1<10) { std::cout << irunMinValL1[icalibL1] << " " <<  scaleValL1[icalibL1] << " " <<  errorScaleValL1[icalibL1] << std::endl ; }
       icalibL1++;
      }
   }
   std::cout << " file_calib1 : " << icalibL1 << " entries in L1 2017" << std::endl;
//   if (icalibL1>1) std::cout << " example runMin "<< irunMinValL1[0] << " runMax "<< irunMinValL1[1] << " scale " << scaleValL1[0] << " +/- " << errorScaleValL1[0] << std::endl;
   file_calib1.close ();

   icalibL1_2017=icalibL1;
   std::ifstream file_calib2;
   file_calib2.open("Tamas/CorrFactHistory2018L1.txt");
   if (!file_calib2) { std::cerr << "cannot open file CorrFactHistory2018L1.txt " << std::endl; }
   else
   {  
      while (!file_calib2.eof () && icalibL1<calmax) {
       file_calib2 >> irunMinValL1[icalibL1] >> scaleValL1[icalibL1] >> errorScaleValL1[icalibL1] ;
       icalibL1++;
      }
   }
   if (irunMinValL1[icalibL1-1]==0) icalibL1--;
   icalibL1_2018=icalibL1-icalibL1_2017;
   std::cout << " file_calib : " << icalibL1_2018 << " entries in L1 2018" << std::endl;
//   if (icalibL1>icalibL1_2017+1) std::cout << " example runMin "<< irunMinValL1[icalibL1_2017] << " runMax "<< irunMinValL1[icalibL1_2017+1] 
//                                           << " scale " << scaleValL1[icalibL1_2017] << " +/- " << errorScaleValL1[icalibL1_2017] << std::endl;
   file_calib2.close ();
 

   std::ifstream file_calib3;
   file_calib3.open("Tamas/CorrFactHistory2017L2.txt");
   icalibL2=0;
   if (!file_calib3) { std::cerr << "cannot open file CorrFactHistory2017L2.txt " << std::endl; }
   else
   {
      while (!file_calib3.eof () && icalibL2<calmax) {
       // from_run_number, corr_val, error_val
       file_calib3 >> irunMinValL2[icalibL2] >> scaleValL2[icalibL2] >> errorScaleValL2[icalibL2] ;
//       if (icalibL2<10) { std::cout << irunMinValL2[icalibL2] << " " <<  scaleValL2[icalibL2] << " " <<  errorScaleValL2[icalibL2] << std::endl ; }
       icalibL2++;
      }
   }
//   std::cout << " file_calib3 : " << icalibL2 << " entries in L2 2017" << std::endl;
//   if (icalibL2>1) std::cout << " example runMin "<< irunMinValL2[0] << " runMax "<< irunMinValL2[1] << " scale " << scaleValL2[0] << " +/- " << errorScaleValL2[0] << std::endl;
   file_calib3.close ();

   icalibL2_2017=icalibL2;
   std::ifstream file_calib4;
   file_calib4.open("Tamas/CorrFactHistory2018L2.txt");
   if (!file_calib4) { std::cerr << "cannot open file CorrFactHistory2018L2.txt " << std::endl; }
   else
   {  
      while (!file_calib4.eof () && icalibL2<calmax) {
       file_calib4 >> irunMinValL2[icalibL2] >> scaleValL2[icalibL2] >> errorScaleValL2[icalibL2] ;
       icalibL2++;
      }
   }
   if (irunMinValL2[icalibL2-1]==0) icalibL2--;
   icalibL2_2018=icalibL2-icalibL2_2017;
//   std::cout << " file_calib : " << icalibL2_2018 << " entries in L2 2018" << std::endl;
//   if (icalibL2>icalibL2_2017+1) std::cout << " example runMin "<< irunMinValL2[icalibL2_2017] << " runMax "<< irunMinValL2[icalibL2_2017+1] 
//                                           << " scale " << scaleValL2[icalibL2_2017] << " +/- " << errorScaleValL2[icalibL2_2017] << std::endl;
   file_calib4.close ();
 
   std::ifstream file_calib5;
   file_calib5.open("Tamas/CorrFactHistory2017L3.txt");
   icalibL3=0;
   if (!file_calib5) { std::cerr << "cannot open file CorrFactHistory2017L3.txt " << std::endl; }
   else
   {
      while (!file_calib5.eof () && icalibL3<calmax) {
       // from_run_number, corr_val, error_val
       file_calib5 >> irunMinValL3[icalibL3] >> scaleValL3[icalibL3] >> errorScaleValL3[icalibL3] ;
//       if (icalibL3<10) { std::cout << irunMinValL3[icalibL3] << " " <<  scaleValL3[icalibL3] << " " <<  errorScaleValL3[icalibL3] << std::endl ; }
       icalibL3++;
      }
   }
//   std::cout << " file_calib5 : " << icalibL3 << " entries in L3 2017" << std::endl;
//   if (icalibL3>1) std::cout << " example runMin "<< irunMinValL3[0] << " runMax "<< irunMinValL3[1] << " scale " << scaleValL3[0] << " +/- " << errorScaleValL3[0] << std::endl;
   file_calib5.close ();

   icalibL3_2017=icalibL3;
   std::ifstream file_calib6;
   file_calib6.open("Tamas/CorrFactHistory2018L3.txt");
   if (!file_calib6) { std::cerr << "cannot open file CorrFactHistory2018L3.txt " << std::endl; }
   else
   {  
      while (!file_calib6.eof () && icalibL3<calmax) {
       file_calib6 >> irunMinValL3[icalibL3] >> scaleValL3[icalibL3] >> errorScaleValL3[icalibL3] ;
       icalibL3++;
      }
   }
   if (irunMinValL3[icalibL3-1]==0) icalibL3--;
   icalibL3_2018=icalibL3-icalibL3_2017;
//   std::cout << " file_calib : " << icalibL3_2018 << " entries in L3 2018" << std::endl;
//   if (icalibL3>icalibL3_2017+1) std::cout << " example runMin "<< irunMinValL3[icalibL3_2017] << " runMax "<< irunMinValL3[icalibL3_2017+1] 
//                                           << " scale " << scaleValL3[icalibL3_2017] << " +/- " << errorScaleValL3[icalibL3_2017] << std::endl;
   file_calib6.close ();
 
   std::ifstream file_calib7;
   file_calib7.open("Tamas/CorrFactHistory2017L4.txt");
   icalibL4=0;
   if (!file_calib7) { std::cerr << "cannot open file CorrFactHistory2017L4.txt " << std::endl; }
   else
   {
      while (!file_calib7.eof () && icalibL4<calmax) {
       // from_run_number, corr_val, error_val
       file_calib7 >> irunMinValL4[icalibL4] >> scaleValL4[icalibL4] >> errorScaleValL4[icalibL4] ;
//       if (icalibL4<10) { std::cout << irunMinValL4[icalibL4] << " " <<  scaleValL4[icalibL4] << " " <<  errorScaleValL4[icalibL4] << std::endl ; }
       icalibL4++;
      }
   }
//   std::cout << " file_calib7 : " << icalibL4 << " entries in L4 2017" << std::endl;
//   if (icalibL4>1) std::cout << " example runMin "<< irunMinValL4[0] << " runMax "<< irunMinValL4[1] << " scale " << scaleValL4[0] << " +/- " << errorScaleValL4[0] << std::endl;
   file_calib7.close ();

   icalibL4_2017=icalibL4;
   std::ifstream file_calib8;
   file_calib8.open("Tamas/CorrFactHistory2018L4.txt");
   if (!file_calib8) { std::cerr << "cannot open file CorrFactHistory2018L4.txt " << std::endl; }
   else
   {  
      while (!file_calib8.eof () && icalibL4<calmax) {
       file_calib8 >> irunMinValL4[icalibL4] >> scaleValL4[icalibL4] >> errorScaleValL4[icalibL4] ;
       icalibL4++;
      }
   }
   if (irunMinValL4[icalibL4-1]==0) icalibL4--;
   icalibL4_2018=icalibL4-icalibL4_2017;
//   std::cout << " file_calib : " << icalibL4_2018 << " entries in L4 2018" << std::endl;
//   if (icalibL4>icalibL4_2017+1) std::cout << " example runMin "<< irunMinValL4[icalibL4_2017] << " runMax "<< irunMinValL4[icalibL4_2017+1] 
//                                           << " scale " << scaleValL4[icalibL4_2017] << " +/- " << errorScaleValL4[icalibL4_2017] << std::endl;
   file_calib8.close ();
 
   std::ifstream file_calib9;
   file_calib9.open("Tamas/CorrFactHistory2017R1.txt");
   icalibR1=0;
   if (!file_calib9) { std::cerr << "cannot open file CorrFactHistory2017R1.txt " << std::endl; }
   else
   {
      while (!file_calib9.eof () && icalibR1<calmax) {
       // from_run_number, corr_val, error_val
       file_calib9 >> irunMinValR1[icalibR1] >> scaleValR1[icalibR1] >> errorScaleValR1[icalibR1] ;
//       if (icalibR1<10) { std::cout << irunMinValR1[icalibR1] << " " <<  scaleValR1[icalibR1] << " " <<  errorScaleValR1[icalibR1] << std::endl ; }
       icalibR1++;
      }
   }
//   std::cout << " file_calib9 : " << icalibR1 << " entries in R1 2017" << std::endl;
//   if (icalibR1>1) std::cout << " example runMin "<< irunMinValR1[0] << " runMax "<< irunMinValR1[1] << " scale " << scaleValR1[0] << " +/- " << errorScaleValR1[0] << std::endl;
   file_calib9.close ();

   icalibR1_2017=icalibR1;
   std::ifstream file_calib10;
   file_calib10.open("Tamas/CorrFactHistory2018R1.txt");
   if (!file_calib10) { std::cerr << "cannot open file CorrFactHistory2018R1.txt " << std::endl; }
   else
   {  
      while (!file_calib10.eof () && icalibR1<calmax) {
       file_calib10 >> irunMinValR1[icalibR1] >> scaleValR1[icalibR1] >> errorScaleValR1[icalibR1] ;
       icalibR1++;
      }
   }
   if (irunMinValR1[icalibR1-1]==0) icalibR1--;
   icalibR1_2018=icalibR1-icalibR1_2017;
//   std::cout << " file_calib : " << icalibR1_2018 << " entries in R1 2018" << std::endl;
//   if (icalibR1>icalibR1_2017+1) std::cout << " example runMin "<< irunMinValR1[icalibR1_2017] << " runMax "<< irunMinValR1[icalibR1_2017+1] 
//                                           << " scale " << scaleValR1[icalibR1_2017] << " +/- " << errorScaleValR1[icalibR1_2017] << std::endl;
   file_calib10.close ();
 
   std::ifstream file_calib11;
   file_calib11.open("Tamas/CorrFactHistory2017R2.txt");
   icalibR2=0;
   if (!file_calib11) { std::cerr << "cannot open file CorrFactHistory2017R2.txt " << std::endl; }
   else
   {
      while (!file_calib11.eof () && icalibR2<calmax) {
       // from_run_number, corr_val, error_val
       file_calib11 >> irunMinValR2[icalibR2] >> scaleValR2[icalibR2] >> errorScaleValR2[icalibR2] ;
//       if (icalibR2<10) { std::cout << irunMinValR2[icalibR2] << " " <<  scaleValR2[icalibR2] << " " <<  errorScaleValR2[icalibR2] << std::endl ; }
       icalibR2++;
      }
   }
//   std::cout << " file_calib11 : " << icalibR2 << " entries in R2 2017" << std::endl;
//   if (icalibR2>1) std::cout << " example runMin "<< irunMinValR2[0] << " runMax "<< irunMinValR2[1] << " scale " << scaleValR2[0] << " +/- " << errorScaleValR2[0] << std::endl;
   file_calib11.close ();

   if (irunMinValR2[icalibR2-1]==0) icalibR2--;
   icalibR2_2017=icalibR2;
   std::ifstream file_calib12;
   file_calib12.open("Tamas/CorrFactHistory2018R2.txt");
   if (!file_calib12) { std::cerr << "cannot open file CorrFactHistory2018R2.txt " << std::endl; }
   else
   {  
      while (!file_calib12.eof () && icalibR2<calmax) {
       file_calib12 >> irunMinValR2[icalibR2] >> scaleValR2[icalibR2] >> errorScaleValR2[icalibR2] ;
       icalibR2++;
      }
   }
   icalibR2_2018=icalibR2-icalibR2_2017;
//   std::cout << " file_calib : " << icalibR2_2018 << " entries in R2 2018" << std::endl;
//   if (icalibR2>icalibR2_2017+1) std::cout << " example runMin "<< irunMinValR2[icalibR2_2017] << " runMax "<< irunMinValR2[icalibR2_2017+1] 
//                                           << " scale " << scaleValR2[icalibR2_2017] << " +/- " << errorScaleValR2[icalibR2_2017] << std::endl;
   file_calib12.close ();
 


}
float run2analysis::GetSFPixelTamas(int subdetid_, UInt_t detid_, int year, int run) {

   int pix = 0;
   int layerorside = 0;

   if (subdetid_==1) {
        pix =1 ;
        if (year==2016) {
                layerorside = int((detid_>>16)&0xF);
        }
        else {
                layerorside = int((detid_>>20)&0xF);
        }
   }
   if (subdetid_==2) {
        pix =2 ;
        layerorside = int((detid_>>23)&0x3); // 1=FPIX- 2=FPIX+
   }

   float scale =1 ;

   if (pix==1) {
     if (layerorside==1) {

       if (year==2017) {
        for (int i=0; i<icalibL1_2017; i++) {
         if (run >= irunMinValL1[i] && run<irunMinValL1[i+1]) {
                scale = scaleValL1[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL1_2017; i<icalibL1-1; i++) {
         if (run >= irunMinValL1[i] && run<irunMinValL1[i+1]) {
                scale = scaleValL1[i];
                return scale;
         }
        }
        if (run >= irunMinValL1[icalibL1-1]) {
                scale = scaleValL1[icalibL1-1];
                return scale;
        }
       }

     }
     else if (layerorside==2) {

       if (year==2017) {
        for (int i=0; i<icalibL2_2017; i++) {
         if (run >= irunMinValL2[i] && run<irunMinValL2[i+1]) {
                scale = scaleValL2[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL2_2017; i<icalibL2-1; i++) {
         if (run >= irunMinValL2[i] && run<irunMinValL2[i+1]) {
                scale = scaleValL2[i];
                return scale;
         }
        }
        if (run >= irunMinValL2[icalibL2-1]) {
                scale = scaleValL2[icalibL2-1];
                return scale;
        }
       }

     }
     else if (layerorside==3) {
     }

       if (year==2017) {
        for (int i=0; i<icalibL3_2017; i++) {
         if (run >= irunMinValL3[i] && run<irunMinValL3[i+1]) {
                scale = scaleValL3[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL3_2017; i<icalibL3-1; i++) {
         if (run >= irunMinValL3[i] && run<irunMinValL3[i+1]) {
                scale = scaleValL3[i];
                return scale;
         }
        }
        if (run >= irunMinValL3[icalibL3-1]) {
                scale = scaleValL3[icalibL3-1];
                return scale;
        }
       }

     else if (layerorside==4) {

       if (year==2017) {
        for (int i=0; i<icalibL4_2017; i++) {
         if (run >= irunMinValL4[i] && run<irunMinValL4[i+1]) {
                scale = scaleValL4[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL4_2017; i<icalibL4-1; i++) {
         if (run >= irunMinValL4[i] && run<irunMinValL4[i+1]) {
                scale = scaleValL4[i];
                return scale;
         }
        }
        if (run >= irunMinValL4[icalibL4-1]) {
                scale = scaleValL4[icalibL4-1];
                return scale;
        }
       }

     }
     else { 
        std::cout << "unknowm layer number in Barrel Pixel " << layerorside << std::endl; 
        return 0; 
     } 
     
   }
   else if (pix==2) {
     if (layerorside==1) {

       if (year==2017) {
        for (int i=0; i<icalibR1_2017; i++) {
         if (run >= irunMinValR1[i] && run<irunMinValR1[i+1]) {
                scale = scaleValR1[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibR1_2017; i<icalibR1-1; i++) {
         if (run >= irunMinValR1[i] && run<irunMinValR1[i+1]) {
                scale = scaleValR1[i];
                return scale;
         }
        }
        if (run >= irunMinValR1[icalibR1-1]) {
                scale = scaleValR1[icalibR1-1];
                return scale;
        }
       }

     }
     else if (layerorside==2) {

       if (year==2017) {
        for (int i=0; i<icalibR2_2017; i++) {
         if (run >= irunMinValR2[i] && run<irunMinValR2[i+1]) {
                scale = scaleValR2[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibR2_2017; i<icalibR2-1; i++) {
         if (run >= irunMinValR2[i] && run<irunMinValR2[i+1]) {
                scale = scaleValR2[i];
                return scale;
         }
        }
        if (run >= irunMinValR2[icalibR2-1]) {
                scale = scaleValR2[icalibR2-1];
                return scale;
        }
       }

     }
     else { 
        std::cout << "unknowm ring number in EndCap Pixel " << layerorside << std::endl; 
        return 0; 
     } 
   }

   return scale;
   
}




///
// Code from Dylan

std::vector<int> run2analysis::SaturationCorrection(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) {
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
//  TMatrix A(N,N);

//---  que pour 1 max bien net
 if(Q.size()<2 || Q.size()>8){
        for (unsigned int i=0;i<Q.size();i++){
                QII.push_back((int) Q[i]);
        }
        return QII;
  }
 if(way){
          vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
          if(*mQ>253){
                 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
                 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
                     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
          }
      else{
          return Q; // no saturation --> no x-talk inversion
      }
  }
//---
 // do nothing else
 return Q; 
}

///



bool run2analysis::clusterCleaning(const std::vector<int>&  Q, int crosstalkInv, uint8_t * exitCode)
{
     vector<int>  ampls = Q;
//   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true,20,25);


  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
  //
         Int_t NofMax=0; Int_t recur255=1; Int_t recur254=1;
         bool MaxOnStart=false;bool MaxInMiddle=false, MaxOnEnd =false;
         Int_t MaxPos=0;
        // D?but avec max
        if(ampls.size()!=1 && ((ampls[0]>ampls[1])
            || (ampls.size()>2 && ampls[0]==ampls[1] && ampls[1]>ampls[2] && ampls[0]!=254 && ampls[0]!=255)
            || (ampls.size()==2 && ampls[0]==ampls[1] && ampls[0]!=254 && ampls[0]!=255)) ){
          NofMax=NofMax+1;  MaxOnStart=true;  }

        // Maximum entour?
        if(ampls.size()>2){
          for (unsigned int i =1; i < ampls.size()-1; i++) {
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1])
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i;
                }
                if(ampls[i]==255 && ampls[i]==ampls[i-1]) {
                        recur255=recur255+1;
                        MaxPos=i-(recur255/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
                if(ampls[i]==254 && ampls[i]==ampls[i-1]) {
                        recur254=recur254+1;
                        MaxPos=i-(recur254/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
            }
        }
        // Fin avec un max
       if(ampls.size()>1){
          if(ampls[ampls.size()-1]>ampls[ampls.size()-2]
             || (ampls.size()>2 && ampls[ampls.size()-1]==ampls[ampls.size()-2] && ampls[ampls.size()-2]>ampls[ampls.size()-3] )
             ||  ampls[ampls.size()-1]==255){
           NofMax=NofMax+1;  MaxOnEnd=true;   }
         }
        // Si une seule strip touch?e
        if(ampls.size()==1){    NofMax=1;}

  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
  //
  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //

   bool shapecdtn=false;
   if (exitCode) *exitCode = 255;

      if(crosstalkInv==1){
        if(NofMax==1){shapecdtn=true; if (exitCode) *exitCode=0;}
        return shapecdtn;
      }

        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        Float_t coeff1=1.7;     Float_t coeff2=2.0;
        Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){ C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=3;
                                }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                }
                }
                if(MaxInMiddle==true){
                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*coeffn*C_M+coeff2*coeffnn*C_D+2*noise || C_M==255){
                                if(ampls.size()==3) shapecdtn=true ;
                                else if(ampls.size()>3){
                                        if(CDPos>MaxPos){
                                                if(ampls.size()-CDPos-1==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1==1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1>1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=(Float_t)ampls[CDPos+2];
                                                }
                                                if(MaxPos>=2){
                                                        C_Mnn=(Float_t)ampls[MaxPos-2];
                                                }
                                                else if(MaxPos<2) C_Mnn=0;
                                        }
                                        if(CDPos<MaxPos){
                                                if(CDPos==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(CDPos==1){
                                                        C_Dn=(Float_t)ampls[0];
                                                        C_Dnn=0;
                                                }
                                                if(CDPos>1){
                                                        C_Dn=(Float_t)ampls[CDPos-1];
                                                        C_Dnn=(Float_t)ampls[CDPos-2];
                                                }
                                                if(ampls.size()-LeftOfMaxPos>1 && MaxPos+2<(int)(ampls.size())-1){
                                                        C_Mnn=(Float_t)ampls[MaxPos+2];
                                                }else C_Mnn=0;
                                        }
                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                           && C_Mnn<=coeff1*coeffn*C_Mn+coeff2*coeffnn*C_M+2*noise
                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise) {
                                                shapecdtn=true;
                                        }

                                }
                        } else if (exitCode) *exitCode=6;
                }
        }
        else if (NofMax>1 && exitCode) *exitCode = 1; // more than one maximum
        if(ampls.size()==1){shapecdtn=true;}
        if(shapecdtn && exitCode) *exitCode=0;

   return shapecdtn;
}


TObject* run2analysis::GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
   if(pos < 256){
      std::string firstPart = Path.substr(0,pos);
      std::string endPart   = Path.substr(pos+1,Path.length());
      TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
      if(TMP!=NULL)return GetObjectFromPath(TMP,endPart,GetACopy);

      printf("ObjectNotFound: %s::%s\n",File->GetName(), Path.c_str());
      return NULL;
   }else{
      if(GetACopy){
         return (File->Get(Path.c_str()))->Clone();
      }else{
         return File->Get(Path.c_str());
      }
   }
}


TH3F* run2analysis::loadDeDxTemplate(std::string path, std::string nameHisto, bool splitByModuleType){
   cout << "  load DeDx template --> root file : " << path.c_str() << ", histo :  " << nameHisto.c_str() << endl;
   TFile* InputFile = new TFile(path.c_str());
//   TH3F* DeDxMap_ = (TH3F*)GetObjectFromPath(InputFile, "Charge_Vs_Path");
   TH3F* DeDxMap_ = (TH3F*)GetObjectFromPath(InputFile, nameHisto.c_str());
   cout << "  load DeDx template : " << nameHisto.c_str() << endl;
   if(!DeDxMap_){printf("dEdx templates in file %s can't be open\n", path.c_str()); exit(0);}

   TH3F* Prob_ChargePath  = (TH3F*)(DeDxMap_->Clone("Prob_ChargePath"));
   Prob_ChargePath->Reset();
   Prob_ChargePath->SetDirectory(0);

   if(!splitByModuleType){
      Prob_ChargePath->RebinX(Prob_ChargePath->GetNbinsX()-1); // <-- do not include pixel in the inclusive
   }

   for(int i=0;i<=Prob_ChargePath->GetXaxis()->GetNbins()+1;i++){
      for(int j=0;j<=Prob_ChargePath->GetYaxis()->GetNbins()+1;j++){
         double Ni = 0;
         for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){Ni+=DeDxMap_->GetBinContent(i,j,k);}

         for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){
            double tmp = 0;
            for(int l=0;l<=k;l++){ tmp+=DeDxMap_->GetBinContent(i,j,l);}

            if(Ni>0){
               Prob_ChargePath->SetBinContent (i, j, k, tmp/Ni);
            }else{
               Prob_ChargePath->SetBinContent (i, j, k, 0);
            }
         }
      }
   }
   InputFile->Close();
   return Prob_ChargePath;
}
