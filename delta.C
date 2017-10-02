#define delta_cxx
#include "delta.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <cmath>

void delta::Loop()
{
//   In a ROOT session, you can do:
//      root> .L delta.C
//      root> delta t
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

   Long64_t nentries = fChain->GetEntriesFast();

    TH1D* delta=new TH1D("delta","Histogram of delta BB",100,-0.5,2);
    TH1D* Hmass=new TH1D("Hmass","Histogram of Hcandidate mass",100,0.,600000);
    TH1D* Leadjetpt=new TH1D("Leadjetpt","Histogram of leadjet pt",100,0.,100000);
    TH1D* Leadjetm=new TH1D("Leadjetm","Histogram of leadjet mass",100,0.,100000);
    TH1D* subLeadjetpt=new TH1D("subLeadjetpt","Histogram of subleadjet pt",100,0.,600000);
    TH1D* subLeadjetm=new TH1D("subLeadjetm","Histogram of subleadjet mass",100,0.,25000);
    TH1D* subsubLeadjetpt=new TH1D("subsubLeadjetpt","Histogram of subsubleadjet pt",100,0.,100000);
    TH1D* subsubLeadjetm=new TH1D("subLeadjetm","Histogram of subsubleadjet mass",100,0.,10000);
    TH1D* numofH=new TH1D("numberofhiss","Histogram of pt of H",39,250000.,3000000.);
    TH1D* numofpass=new TH1D("numberofpass","R=0.2 Track Jet",39,250000.,3000000.);
    TH1D* numboostedjet=new TH1D("numberofboostedjet","Pt of bossted jet which match to truth higgs boson",39,250000.,3000000.);
    

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       double dR;
       double Hm;
       TLorentzVector Bquark;
       vector<TLorentzVector> vB;
       TLorentzVector subjet;
       vector<TLorentzVector> vsubjet;
       TLorentzVector Higgs;
       vector<TLorentzVector> vHiggs;
       TLorentzVector boostedjet;
       vector<TLorentzVector> vboostedjet;
       if(jentry%1000==0) cout<<jentry<<endl;
       TLorentzVector Bhadron;
       vector<TLorentzVector> vBhadron1;
       vector<TLorentzVector> vBhadron2;
       vector<vector <TLorentzVector> > vBhadron;
       
       //b hadron 1
       for(int i=0;i<truth_higgs_matchedbhad_pt->at(0).size();i++)
       {
           //cout<<"testok"<<"  ";
           Bhadron.SetPtEtaPhiM(truth_higgs_matchedbhad_pt->at(0)[i],truth_higgs_matchedbhad_eta->at(0)[i],truth_higgs_matchedbhad_phi->at(0)[i],truth_higgs_matchedbhad_m->at(0)[i]);
           vBhadron1.push_back(Bhadron);
       }
       vBhadron.push_back(vBhadron1);
       //b hadron 2
       for(int i=0;i<truth_higgs_matchedbhad_pt->at(1).size();i++)
       {
           //cout<<"testok"<<"  ";
           Bhadron.SetPtEtaPhiM(truth_higgs_matchedbhad_pt->at(1)[i],truth_higgs_matchedbhad_eta->at(1)[i],truth_higgs_matchedbhad_phi->at(1)[i],truth_higgs_matchedbhad_m->at(1)[i]);
           vBhadron2.push_back(Bhadron);
       }
       vBhadron.push_back(vBhadron2);

       
       //b quark
       for(int i=0;i<truth_part_pdgid->size();i++)
       {
           if(fabs(truth_part_pdgid->at(i))==5)
           {
               Bquark.SetPtEtaPhiM(truth_part_pt->at(i),truth_part_eta->at(i),truth_part_phi->at(i),truth_part_m->at(i));
               vB.push_back(Bquark);
               //cout<<"find one!"<<endl;
               //cout<<truth_part_eta->at(i)<<"  "<<truth_part_phi->at(i)<<endl;
           }
       }
       
       //truth higgs
       for(int i=0;i<truth_higgs_pt->size();i++)
       {
           Higgs.SetPtEtaPhiM(truth_higgs_pt->at(i),truth_higgs_eta->at(i),truth_higgs_phi->at(i),truth_higgs_m->at(i));
           vHiggs.push_back(Higgs);
           numofH->Fill(truth_higgs_pt->at(i));
       }
       
       int numberofb1=0;
       int numberofb2=0;
       vector<int> boostedjetmatchbhadron;
       vector<int> whichboostedjet;
       //boosted jets associated to a Higgs boson and two b hadron
       for(int i=0;i<hcand_boosted_pt->size();i++)
       {
           Hmass->Fill(hcand_boosted_m->at(i));
           boostedjet.SetPtEtaPhiM(hcand_boosted_pt->at(i),hcand_boosted_eta->at(i),hcand_boosted_phi->at(i),hcand_boosted_m->at(i));
           //cout<<hcand_boosted_eta->at(i)<<"  "<<hcand_boosted_phi->at(i)<<endl;

           //76Gev<mjet<146Gev
           if(hcand_boosted_m->at(i)<146000.&&hcand_boosted_m->at(i)>76000.)
           {
               if(boostedjet.DeltaR(vHiggs[0])<1.0||boostedjet.DeltaR(vHiggs[1])<1.0)
               {
                   //boostedjet.Eta();
                   for(int j=0;j<vBhadron1.size();j++)
                   {
                       if(boostedjet.DeltaR(vBhadron1[j])<1.)//I am not sure about this value(distance between b quark and large R jet)
                           numberofb1++;
                   }
                   for(int j=0;j<vBhadron2.size();j++)
                   {
                       if(boostedjet.DeltaR(vBhadron2[j])<1.)//I am not sure about this value(distance between b quark and large R jet)
                           numberofb2++;
                   }
               }
               if(numberofb1<2&&numberofb2<2)
               {
                   numberofb1=0;
                   numberofb2=0;
               }
               else if(numberofb1<2&&numberofb2>1)
               {
                   numberofb1=0;
                   numberofb2=-99;
                   numboostedjet->Fill(hcand_boosted_pt->at(i));
                   vboostedjet.push_back(boostedjet);
                   boostedjetmatchbhadron.push_back(1);
                   whichboostedjet.push_back(i);
               }
               else if(numberofb2<2&&numberofb1>1)
               {
                   numberofb2=0;
                   numberofb1=-99;
                   numboostedjet->Fill(hcand_boosted_pt->at(i));
                   vboostedjet.push_back(boostedjet);
                   boostedjetmatchbhadron.push_back(0);
                   whichboostedjet.push_back(i);
               }
               else if(numberofb2>1&&numberofb1>1)
               {
                   numberofb2=0;
                   numberofb1=0;
                   numboostedjet->Fill(hcand_boosted_pt->at(i));
                   vboostedjet.push_back(boostedjet);
                   //cout<<i<<"  "<<hcand_boosted_pt->size()<<endl;
                   boostedjetmatchbhadron.push_back(10);
                   whichboostedjet.push_back(i);
               }
           }
       }

       
       

       if(vB.size()!=4)
           cout<<vB.size()<<endl;
       //if(jet_ak2track_asso_lead_eta->size()!=jet_ak2track_asso_sublead_eta->size())
        //   cout<<jet_ak2track_asso_lead_eta->size()<<"  "<<jet_ak2track_asso_sublead_eta->size()<<endl;
       //cout<<jet_ak2track_asso_lead_phi->size();
  
       //subjet distance
       for(int i=0;i<jet_ak2track_asso_lead_eta->size();i++)
       {
           dR=sqrt(pow(jet_ak2track_asso_lead_eta->at(i)-jet_ak2track_asso_sublead_eta->at(i),2)+pow(jet_ak2track_asso_lead_phi->at(i)-jet_ak2track_asso_sublead_phi->at(i),2));
           delta->Fill(dR);
           Leadjetpt->Fill(jet_ak2track_asso_lead_pt->at(i));
           Leadjetm->Fill(jet_ak2track_asso_lead_m->at(i));
           subLeadjetpt->Fill(jet_ak2track_asso_sublead_pt->at(i));
           subLeadjetm->Fill(jet_ak2track_asso_sublead_m->at(i));
           subsubLeadjetpt->Fill(jet_ak2track_asso_subsublead_pt->at(i));
           subsubLeadjetm->Fill(jet_ak2track_asso_subsublead_m->at(i));
       }
       
       
       //find minimun
       if(vboostedjet.size()==1)
       {
          // if(jet_ak2track_asso_lead_m->size()!=1)
          //     cout<<"?"<<endl;
           if(boostedjetmatchbhadron[0]!=10)
           {
               TLorentzVector leadjet;
               TLorentzVector subleadjet;
               leadjet.SetPtEtaPhiM(jet_ak2track_asso_lead_pt->at(whichboostedjet[0]),jet_ak2track_asso_lead_eta->at(whichboostedjet[0]),jet_ak2track_asso_lead_phi->at(whichboostedjet[0]),jet_ak2track_asso_lead_m->at(whichboostedjet[0]));
               subleadjet.SetPtEtaPhiM(jet_ak2track_asso_sublead_pt->at(whichboostedjet[0]),jet_ak2track_asso_sublead_eta->at(whichboostedjet[0]),jet_ak2track_asso_sublead_phi->at(whichboostedjet[0]),jet_ak2track_asso_sublead_m->at(whichboostedjet[0]));
               int numberoflead=0;
               int numberofsublead=0;
               vector<int> leadwhichb1;
               vector<int> subleadwhichb1;
               for(int j=0;j<vBhadron[boostedjetmatchbhadron[0]].size();j++)
               {
                   //leadjet.DeltaR(vB[j]);
                   //leadjet.PseudoRapidity();
                   
                   //if (leadjet.CosTheta()*leadjet.CosTheta() < 1) continue;
                   //else cout<<leadjet.CosTheta()<<"  "<<leadjet.Pt()<<"  "<<jet_ak2track_asso_lead_pt->at(0)<<"  "<<jet_ak2track_asso_lead_eta->at(0)<<"  "<<jet_ak2track_asso_lead_phi->at(0)<<"  "<<jet_ak2track_asso_lead_m->at(0)<<endl;
                   if(leadjet.DeltaR(vBhadron[boostedjetmatchbhadron[0]][j])<0.3)
                   {
                       numberoflead++;
                       leadwhichb1.push_back(j);
                   }
               }
               for(int j=0;j<vBhadron[boostedjetmatchbhadron[0]].size();j++)
               {
                   if(subleadjet.DeltaR(vBhadron[boostedjetmatchbhadron[0]][j])<0.3)
                   {
                       numberofsublead++;
                       subleadwhichb1.push_back(j);
                   }
               }
               if(numberoflead>0&&numberofsublead>0)
               {
                   if(numberoflead==1&&numberofsublead==1&&leadwhichb1[0]==subleadwhichb1[0]) continue;
                   else numofpass->Fill(vboostedjet[0].Pt());
               }
           }
       }
       

       if(vboostedjet.size()==3)
           cout<<"???"<<endl;

       //vector<int> eachbquarkmin;
       vector<int> eachbhadronmin;
       int subjetmatch[4] ={0};
  
   
       //at least two large R jet associated to Higgs boson
       if(vboostedjet.size()==2)
       {
           for(int kk=0;kk<vboostedjet.size();kk++)
           {
               subjet.SetPtEtaPhiM(jet_ak2track_asso_lead_pt->at(whichboostedjet[kk]),jet_ak2track_asso_lead_eta->at(whichboostedjet[kk]),jet_ak2track_asso_lead_phi->at(whichboostedjet[kk]),jet_ak2track_asso_lead_m->at(whichboostedjet[kk]));
               vsubjet.push_back(subjet);
               subjet.SetPtEtaPhiM(jet_ak2track_asso_sublead_pt->at(whichboostedjet[kk]),jet_ak2track_asso_sublead_eta->at(whichboostedjet[kk]),jet_ak2track_asso_sublead_phi->at(whichboostedjet[kk]),jet_ak2track_asso_sublead_m->at(whichboostedjet[kk]));
               vsubjet.push_back(subjet);
           }
           //the first boostedjet
           for(int i=0;i<vBhadron[boostedjetmatchbhadron[0]].size();i++)
           {
               double DeltaRmin=9999;
               int test=0;
               int indexmin;
               for(int j=0;j<2;j++)
               {
                   if(DeltaRmin>vBhadron[boostedjetmatchbhadron[0]][i].DeltaR(vsubjet[j]))
                   {
                       DeltaRmin=vBhadron[boostedjetmatchbhadron[0]][i].DeltaR(vsubjet[j]);
                       indexmin=j;
                   }
               }
               if(DeltaRmin<0.3)
               {
                   for(int k=0;k<eachbhadronmin.size();k++)
                   {
                       if(indexmin==eachbhadronmin[k])
                       {
                           test=1;
                           break;
                       }
                   }
                   if(test==0)
                   {
                       subjetmatch[indexmin]=1;
                       eachbhadronmin.push_back(indexmin);
                   }
               }
           }
           //the second boostedjet
           for(int i=0;i<vBhadron[boostedjetmatchbhadron[1]].size();i++)
           {
               double DeltaRmin=9999;
               int test=0;
               int indexmin;
               for(int j=2;j<4;j++)
               {
                   if(DeltaRmin>vBhadron[boostedjetmatchbhadron[1]][i].DeltaR(vsubjet[j]))
                   {
                       DeltaRmin=vBhadron[boostedjetmatchbhadron[1]][i].DeltaR(vsubjet[j]);
                       indexmin=j;
                   }
               }
               if(DeltaRmin<0.3)
               {
                   for(int k=0;k<eachbhadronmin.size();k++)
                   {
                       if(indexmin==eachbhadronmin[k])
                       {
                           test=1;
                           break;
                       }
                   }
                   if(test==0)
                   {
                       subjetmatch[indexmin]=1;
                       eachbhadronmin.push_back(indexmin);
                   }
               }
           }
           //if B hadron match the subject
           if(subjetmatch[0]==1&&subjetmatch[1]==1&&subjetmatch[2]==1&&subjetmatch[3]==1)
           {
               for(int i=0;i<vboostedjet.size();i++)
               {
                   numofpass->Fill(vboostedjet[i].Pt());
               }
           }
           else if(subjetmatch[0]==1&&subjetmatch[1]==1)
           {
               numofpass->Fill(vboostedjet[0].Pt());
           }
           else if(subjetmatch[2]==1&&subjetmatch[3]==1)
           {
               numofpass->Fill(vboostedjet[1].Pt());
           }

           
           
           
           /* //using b quark but not b hadron
           for(int i=0;i<vB.size();i++)
           {
               double DeltaRmin=9999;
               int test=0;
               int indexmin;
               for(int j=0;j<vsubjet.size();j++)
               {
                   if(DeltaRmin>vB[i].DeltaR(vsubjet[j]))
                   {
                       DeltaRmin=vB[i].DeltaR(vsubjet[j]);
                       indexmin=j;
                   }
               }
               if(DeltaRmin<0.3)
               {
                   for(int k=0;k<eachbquarkmin.size();k++)
                   {
                       if(indexmin==eachbquarkmin[k])
                       {
                           test=1;
                           break;
                       }
                   }
                   if(test==0)
                   {
                       subjetmatch[indexmin]=1;
                       eachbquarkmin.push_back(indexmin);
                   }
               }
           }
           //if B quark match the subject
           if(subjetmatch[0]==1&&subjetmatch[1]==1&&subjetmatch[2]==1&&subjetmatch[3]==1)
           {
               for(int i=0;i<vboostedjet.size();i++)
               {
                   numofpass->Fill(vboostedjet[i].Pt());
               }
           }
           else if(subjetmatch[0]==1&&subjetmatch[1]==1)
           {
               double dRsubjetHiggs=100;
               int indexhiggs;
               for(int i=0;i<vboostedjet.size();i++)
               {
                   if(dRsubjetHiggs>(pow(vboostedjet[i].DeltaR(vsubjet[0]),2)+pow(vboostedjet[i].DeltaR(vsubjet[1]),2)))
                   {
                       dRsubjetHiggs=pow(vboostedjet[i].DeltaR(vsubjet[0]),2)+pow(vboostedjet[i].DeltaR(vsubjet[1]),2);
                       indexhiggs=i;
                   }
               }
               numofpass->Fill(vboostedjet[indexhiggs].Pt());
           }
           else if(subjetmatch[2]==1&&subjetmatch[3]==1)
           {
               double dRsubjetHiggs=100;
               int indexhiggs;
               for(int i=0;i<truth_higgs_pt->size();i++)
               {
                   if(dRsubjetHiggs>(pow(vboostedjet[i].DeltaR(vsubjet[2]),2)+pow(vboostedjet[i].DeltaR(vsubjet[3]),2)))
                   {
                       dRsubjetHiggs=pow(vboostedjet[i].DeltaR(vsubjet[2]),2)+pow(vboostedjet[i].DeltaR(vsubjet[3]),2);
                       indexhiggs=i;
                   }
               }
               numofpass->Fill(vboostedjet[indexhiggs].Pt());
           }*/
       }
       
       
       
   }
        TCanvas *c1=new TCanvas();
        delta->Draw();
        TCanvas *c2=new TCanvas();
        Hmass->Draw();
    Hmass->SetMarkerStyle(20);
    Hmass->GetXaxis()->SetTitle("Higgs Jet mass [Mev]");
    Hmass->GetYaxis()->SetTitle("N(Higgs jet)");
    Hmass->Draw("P0");
    
    gPad->BuildLegend();
        TCanvas *c3=new TCanvas();
        Leadjetpt->Draw();
        TCanvas *c4=new TCanvas();
        Leadjetm->Draw();
        TCanvas *c5=new TCanvas();
        subLeadjetpt->Draw();
        TCanvas *c6=new TCanvas();
        subLeadjetm->Draw();
        TCanvas *c7=new TCanvas();
        subsubLeadjetpt->Draw();
        TCanvas *c8=new TCanvas();
        subsubLeadjetm->Draw();
        TCanvas *c9=new TCanvas();
        numofH->Draw();
        TCanvas *c10=new TCanvas();
        gStyle->SetOptTitle(kFALSE);
        gStyle->SetOptStat(0);
        //numofpass->Divide(numboostedjet);
        numofpass->SetMarkerStyle(20);
    numofpass->GetXaxis()->SetTitle("Higgs Jet pt [Mev]");
    numofpass->GetYaxis()->SetTitle("N(Double subjet b-label|Higgs jet)");
        numofpass->Draw("P0");
    
        gPad->BuildLegend();
        //TCanvas *c11=new TCanvas();
        //numofpass->Divide(numofH);
        //numofpass->Draw();
        TCanvas *c12=new TCanvas();
    numboostedjet->SetMarkerStyle(20);
    numboostedjet->GetXaxis()->SetTitle("Higgs Jet pt [Mev]");
    numboostedjet->GetYaxis()->SetTitle("N(Higgs jet)");
    numboostedjet->Draw("P0");
    gPad->BuildLegend();
    
    
    //leadeta->Draw();
    //leadphi->Draw();
    //subleadeta->Draw();
    //subleadphi->Draw();

}
