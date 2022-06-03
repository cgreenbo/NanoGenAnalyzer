#define NanoGENanalyzer_cxx
#include "NanoGENanalyzer.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TTree.h>
#include <TObject.h>

#include <iostream>

using namespace std;

void NanoGENanalyzer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L NanoGENanalyzer.C
//      root> NanoGENanalyzer t
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


   TFile* fOutput = new TFile("output.root","RECREATE");
   
   TTree* tOutput = new TTree("Tree","Tree",0);
   TTree* tOutput_Real = new TTree("Tree_real","Events_real",0);
   TTree* tOutput_Complex = new TTree("Tree_complex","Events_complex",0);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;


   float weight;
   float sinTheta, cosTheta;
   float sinThetaStar, cosThetaStar;
   float sinPhiStar, cosPhiStar, PhiStar;
   float lepton_E_Wframe;
   float top_pt, W_pt, lepton_pt, neutrino_pt;
   float top_mass, W_mass, W_transverse_mass;
   float top_eta, neutrino_eta;
   float top_phi, neutrino_phi;
   float nature_lepton;

   /*________Variables to separate C and R solutions of W boson reconstruction________*/

   float sinTheta_Real, cosTheta_Real;
   float sinThetaStar_Real, cosThetaStar_Real;
   float sinPhiStar_Real, cosPhiStar_Real, PhiStar_Real;
   float lepton_E_Wframe_Real;
   float top_pt_Real, W_pt_Real, lepton_pt_Real, neutrino_pt_Real; 
   float top_eta_Real, neutrino_eta_Real;
   float top_phi_Real, neutrino_phi_Real;
   float top_mass_Real, W_mass_Real, W_transverse_mass_Real;

   float sinTheta_Complex, cosTheta_Complex;
   float sinThetaStar_Complex, cosThetaStar_Complex;
   float sinPhiStar_Complex, cosPhiStar_Complex, PhiStar_Complex;
   float lepton_E_Wframe_Complex;
   float top_pt_Complex, W_pt_Complex, lepton_pt_Complex, neutrino_pt_Complex;
   float top_eta_Complex, neutrino_eta_Complex;
   float top_phi_Complex, neutrino_phi_Complex;
   float top_mass_Complex, W_mass_Complex, W_transverse_mass_Complex;



   int nbOfRealSolutions = 0;
   int nbOfComplexSolutions = 0;

   int nbOfRealSolutions_Streco = 0;
   int nbOfComplexSolutions_Streco = 0;


   double met_pt;
   double met_phi;


   tOutput->Branch("nature_lepton",&nature_lepton,"nature_lepton/F");
   tOutput->Branch("weight",&weight,"weight/F");
   tOutput->Branch("cosTheta",&cosTheta,"cosTheta/F");
   tOutput->Branch("sinTheta",&sinTheta,"sinTheta/F");
   tOutput->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
   tOutput->Branch("sinThetaStar",&sinThetaStar,"sinThetaStar/F");
   tOutput->Branch("cosPhiStar",&cosPhiStar,"cosPhiStar/F");
   tOutput->Branch("sinPhiStar",&sinPhiStar,"sinPhiStar/F");
   tOutput->Branch("PhiStar",&PhiStar,"PhiStar/F");
   tOutput->Branch("lepton_E_Wframe",&lepton_E_Wframe,"lepton_E_Wframe/F");
   tOutput->Branch("top_pt",&top_pt,"top_pt/F");
   tOutput->Branch("top_eta",&top_eta,"top_eta/F");
   tOutput->Branch("top_phi",&top_phi,"top_phi/F");
   tOutput->Branch("top_mass", &top_mass, "top_mass/F");
   tOutput->Branch("W_pt",&W_pt,"W_pt/F");
   tOutput->Branch("lepton_pt",&lepton_pt,"lepton_pt/F");
   tOutput->Branch("neutrino_pt",&neutrino_pt,"neutrino_pt/F");
   tOutput->Branch("neutrino_eta",&neutrino_eta,"neutrino_eta/F");
   tOutput->Branch("neutrino_phi",&neutrino_phi,"neutrino_phi/F");
   tOutput->Branch("W_mass", &W_mass, "W_mass/F");
   tOutput->Branch("W_transverse_mass", &W_transverse_mass, "W_transverse_mass/F");

   tOutput_Real->Branch("cosTheta_Real",&cosTheta_Real,"cosTheta_Real/F");
   tOutput_Real->Branch("sinTheta_Real",&sinTheta_Real,"sinTheta_Real/F");
   tOutput_Real->Branch("cosThetaStar_Real",&cosThetaStar_Real,"cosThetaStar_Real/F");
   tOutput_Real->Branch("sinThetaStar_Real",&sinThetaStar_Real,"sinThetaStar_Real/F");
   tOutput_Real->Branch("cosPhiStar_Real",&cosPhiStar_Real,"cosPhiStar_Real/F");
   tOutput_Real->Branch("sinPhiStar_Real",&sinPhiStar_Real,"sinPhiStar_Real/F");
   tOutput_Real->Branch("PhiStar_Real",&PhiStar_Real,"PhiStar_Real/F");
   tOutput_Real->Branch("lepton_E_Wframe_Real",&lepton_E_Wframe_Real,"lepton_E_Wframe_Real/F");
   tOutput_Real->Branch("top_pt_Real",&top_pt_Real,"top_pt_Real/F");
   tOutput_Real->Branch("top_eta_Real",&top_eta_Real,"top_eta_Real/F");
   tOutput_Real->Branch("top_phi_Real",&top_phi_Real,"top_phi_Real/F");
   tOutput_Real->Branch("top_mass_Real", &top_mass_Real, "top_mass_Real/F");
   tOutput_Real->Branch("W_pt_Real",&W_pt_Real,"W_pt_Real/F");
   tOutput_Real->Branch("lepton_pt_Real",&lepton_pt_Real,"lepton_pt_Real/F");
   tOutput_Real->Branch("neutrino_pt_Real",&neutrino_pt_Real,"neutrino_pt_Real/F");
   tOutput_Real->Branch("neutrino_eta_Real",&neutrino_eta_Real,"neutrino_eta_Real/F");
   tOutput_Real->Branch("neutrino_phi_Real",&neutrino_phi_Real,"neutrino_phi_Real/F");
   tOutput_Real->Branch("W_mass_Real", &W_mass_Real, "W_mass_Real/F");
   tOutput_Real->Branch("W_transverse_mass_Real", &W_transverse_mass_Real, "W_transverse_mass_Real/F");

   tOutput_Complex->Branch("cosTheta_Complex",&cosTheta_Complex,"cosTheta_Complex/F");
   tOutput_Complex->Branch("sinTheta_Complex",&sinTheta_Complex,"sinTheta_Complex/F");
   tOutput_Complex->Branch("cosThetaStar_Complex",&cosThetaStar_Complex,"cosThetaStar_Complex/F");
   tOutput_Complex->Branch("sinThetaStar_Complex",&sinThetaStar_Complex,"sinThetaStar_Complex/F");
   tOutput_Complex->Branch("cosPhiStar_Complex",&cosPhiStar_Complex,"cosPhiStar_Complex/F");
   tOutput_Complex->Branch("sinPhiStar_Complex",&sinPhiStar_Complex,"sinPTreehiStar_Complex/F");
   tOutput_Complex->Branch("PhiStar_Complex",&PhiStar_Complex,"PhiStar_Complex/F");
   tOutput_Complex->Branch("lepton_E_Wframe_Complex",&lepton_E_Wframe_Complex,"lepton_E_Wframe_Complex/F");
   tOutput_Complex->Branch("top_pt_Complex",&top_pt_Complex,"top_pt_Complex/F");
   tOutput_Complex->Branch("top_eta_Complex",&top_eta_Complex,"top_eta_Complex/F");
   tOutput_Complex->Branch("top_phi_Complex",&top_phi_Complex,"top_phi_Complex/F");
   tOutput_Complex->Branch("top_mass_Complex", &top_mass_Complex, "top_mass_Complex/F");
   tOutput_Complex->Branch("W_pt_Complex",&W_pt_Complex,"W_pt_Complex/F");
   tOutput_Complex->Branch("lepton_pt_Complex",&lepton_pt_Complex,"lepton_pt_Complex/F");
   tOutput_Complex->Branch("neutrino_pt_Complex",&neutrino_pt_Complex,"neutrino_pt_Complex/F");
   tOutput_Complex->Branch("neutrino_eta_Complex",&neutrino_eta_Complex,"neutrino_eta_Complex/F");
   tOutput_Complex->Branch("neutrino_phi_Complex",&neutrino_phi_Complex,"neutrino_phi_Complex/F");
   tOutput_Complex->Branch("W_mass_Complex", &W_mass_Complex, "W_mass_Complex/F");
   tOutput_Complex->Branch("W_transverse_mass_Complex", &W_transverse_mass_Complex, "W_transverse_mass_Complex/F");



   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry%1000000==0)
      std::cout<<jentry<<std::endl;

      // if(jentry==1000000) break;

      met_pt = MET_fiducialGenPt;
      met_phi = MET_fiducialGenPhi;

      // std::cout<<"met_pt = "<<met_pt<<std::endl;
      // std::cout<<"met_phi = "<<met_phi<<std::endl;


      /*----- PARTICLE-LEVEL SELECTION -----*/

      unsigned int nSelectedGenDressedLepton = 0;
      TLorentzVector SelectedGenDressedLepton_sublead(0,0,0,0);
      unsigned int nSelectedGenJet = 0;
      std::vector<TLorentzVector> SelectedJet;
      std::vector<int> SelectedJet_hadronFlavour;
      int nbOfBjets = 0;
      int nbOfNonBjets = 0;
      double maxEtaNonBjet = 0.0;
      bool atLeastOneBjet = false;
      
      TLorentzVector SelectedGenDressedLepton_lead(0,0,0,0);
      TLorentzVector qSpec(0,0,0,0);
      TLorentzVector bJet(0,0,0,0);
      // Reco objects
      TLorentzVector Wboson(0,0,0,0);
      TLorentzVector neutrino(0,0,0,0);
      TLorentzVector top(0,0,0,0);


      TLorentzVector Wboson_Streco(0,0,0,0);
      TLorentzVector neutrino_Streco(0,0,0,0), neutrino_pz_plus_Streco(0,0,0,0), neutrino_pz_minus_Streco(0,0,0,0);
      TLorentzVector top_Streco(0,0,0,0);


      Bool_t isRealSolution;
      Bool_t isRealSolution_Streco;
      Bool_t isNotAValidReco;

      TLorentzVector p4nu_plus(0,0,0,0), p4nu_minus(0,0,0,0);


      for (unsigned int i=0; i<nGenDressedLepton; i++)
      {  
         // if(i==0) std::cout<<"Before selection on leptons"<<std::endl;
         // std::cout<<"lepton pt = "<<GenDressedLepton_pt[i]<<std::endl;
         
         if((GenDressedLepton_pt[i]<32 || fabs(GenDressedLepton_eta[i])>2.1 || (fabs(GenDressedLepton_eta[i])<1.5660 && fabs(GenDressedLepton_eta[i])>1.4442)) && fabs(GenDressedLepton_pdgId[i])==11) continue;
         if((GenDressedLepton_pt[i]<30 || fabs(GenDressedLepton_eta[i])>2.4) && fabs(GenDressedLepton_pdgId[i])==13) continue;

         if (GenDressedLepton_pt[i] > SelectedGenDressedLepton_lead.Pt()) 
         SelectedGenDressedLepton_lead.SetPtEtaPhiM(GenDressedLepton_pt[i], GenDressedLepton_eta[i], GenDressedLepton_phi[i], GenDressedLepton_mass[i]); 
         else if (GenDressedLepton_pt[i] > SelectedGenDressedLepton_sublead.Pt())
         SelectedGenDressedLepton_sublead.SetPtEtaPhiM(GenDressedLepton_pt[i], GenDressedLepton_eta[i], GenDressedLepton_phi[i], GenDressedLepton_mass[i]);

         nSelectedGenDressedLepton++;
	
      }

      // std::cout << "nSelectedGenDressedLepton="<<nSelectedGenDressedLepton<<endl;
      // std::cout << "SelectedGenDressedLepton_lead pt="<<SelectedGenDressedLepton_lead.Pt()<<" eta="<<SelectedGenDressedLepton_lead.Eta()<<endl;
      // std::cout << "SelectedGenDressedLepton_sublead pt="<<SelectedGenDressedLepton_sublead.Pt()<<" eta="<<SelectedGenDressedLepton_sublead.Eta()<<endl;


      for (unsigned int i=0; i<nGenJet; i++)
      {
	      // std::cout << "Jet "<<i<<" pt="<<GenJet_pt[i]<<" eta="<<GenJet_eta[i]<<" hadronFlavour="<< (int)GenJet_hadronFlavour[i] <<endl;
	      if(GenJet_pt[i]<40 || fabs(GenJet_eta[i])>4.7) continue;
         if(fabs(GenJet_eta[i])>=2.4 && GenJet_pt[i]<60) continue;

	      TLorentzVector jet;
	      jet.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);

	      bool keepJet = true;
	      for (unsigned int j=0; j<nGenDressedLepton; j++)
         {
            if((GenDressedLepton_pt[j]<32 || fabs(GenDressedLepton_eta[j])>2.1 || (fabs(GenDressedLepton_eta[j])<1.5660 && fabs(GenDressedLepton_eta[j])>1.4442)) && fabs(GenDressedLepton_pdgId[j])==11) continue;
            if((GenDressedLepton_pt[j]<30 || fabs(GenDressedLepton_eta[j])>2.4) && fabs(GenDressedLepton_pdgId[j])==13) continue;
	         TLorentzVector lepton;
	         lepton.SetPtEtaPhiM(GenDressedLepton_pt[j], GenDressedLepton_eta[j], GenDressedLepton_phi[j], GenDressedLepton_mass[j]);
	         if(jet.DeltaR(lepton)<0.4) keepJet = false;
	      }
         if(keepJet) 
         {
            SelectedJet.push_back(jet);
            SelectedJet_hadronFlavour.push_back(GenJet_hadronFlavour[i]);
            nSelectedGenJet++;
         }
      }

      // std::cout << "nGenJet="<< nGenJet<<" nSelectedGenJet="<<nSelectedGenJet<<endl;
      for (unsigned int i=0; i<nSelectedGenJet; i++) 
      {
         // std::cout << "SelectedJet pt="<<SelectedJet[i].Pt()<<" eta="<<SelectedJet[i].Eta()<<" hadronFlavour="<<SelectedJet_hadronFlavour[i]<<endl; 
	      if(SelectedJet_hadronFlavour[i]==5 && nbOfBjets<1)
         { 
            atLeastOneBjet = true;
            nbOfBjets += 1;
            bJet.SetPtEtaPhiM(SelectedJet[i].Pt(),SelectedJet[i].Eta(),SelectedJet[i].Phi(),SelectedJet[i].M());
         }
         else
         {
            nbOfNonBjets += 1;
            if(fabs(SelectedJet[i].Eta())>fabs(maxEtaNonBjet))
            {
               maxEtaNonBjet = SelectedJet[i].Eta();
               qSpec.SetPtEtaPhiM(SelectedJet[i].Pt(),SelectedJet[i].Eta(),SelectedJet[i].Phi(),SelectedJet[i].M());
            } 
         }
      }

      // std::cout<<"atLeastOneBjet="<<atLeastOneBjet<<endl;


      // std::cout<<"Passed 1st Selection"<<endl;
      // std::cout<<"Non b Jet to keep pt="<<qSpec.Pt()<<" eta="<<qSpec.Eta()<<endl;
      // std::cout<<"b Jet to keep pt="<<bJet.Pt()<<" eta="<<bJet.Eta()<<endl;
      // std::cout<<"Number of b jets = "<<nbOfBjets<<endl;
      // std::cout<<"Lepton to keep pT= "<<SelectedGenDressedLepton_lead.Pt()<<" eta = "<<SelectedGenDressedLepton_lead.Eta()<<endl;

      reconstructW(&Wboson, &neutrino, SelectedGenDressedLepton_lead, met_pt, met_phi, &isRealSolution, &isNotAValidReco);
         
      top = bJet + Wboson;

      if(isNotAValidReco == true) continue;
      if( nSelectedGenDressedLepton>2 || nSelectedGenDressedLepton==0 || nSelectedGenJet<2 || !atLeastOneBjet || nbOfBjets>1) continue;


      // ReconstructW_Streco(&Wboson_Streco, &neutrino_Streco, SelectedGenDressedLepton_lead, met_pt, met_phi, &isRealSolution_Streco, &p4nu_plus, &p4nu_minus);

      // top_Streco = bJet + Wboson_Streco;

      // if(isRealSolution_Streco) nbOfRealSolutions_Streco += 1;
      // else nbOfComplexSolutions_Streco +=1;


      // cout<<"W boson pt= "<<Wboson.Pt()<<" W boson Eta= "<<Wboson.Eta()<<" W boson Phi= "<<Wboson.Phi()<<" W boson mass= "<<Wboson.M()<<endl;
      // cout<<"neutrino pt= "<<neutrino.Pt()<<" neutrino Eta= "<<neutrino.Eta()<<" neutrino Phi= "<<neutrino.Phi()<<endl;
      // cout<<"top pt= "<<top.Pt()<<" top Eta= "<<top.Eta()<<" top Phi= "<<top.Phi()<<" top mass= "<<top.M()<<endl;

      // cout<<"W boson Streco pt= "<<Wboson_Streco.Pt()<<" W boson Streco Eta= "<<Wboson_Streco.Eta()<<" W boson Streco Phi= "<<Wboson_Streco.Phi()<<" W boson Streco mass= "<<Wboson_Streco.M()<<endl;
      // cout<<"neutrino Streco pt= "<<neutrino_Streco.Pt()<<" neutrino Streco Eta= "<<neutrino_Streco.Eta()<<" neutrino Streco Phi= "<<neutrino_Streco.Phi()<<endl;
      // cout<<"top Streco pt= "<<top_Streco.Pt()<<" top Streco Eta= "<<top_Streco.Eta()<<" top Streco Phi= "<<top_Streco.Phi()<<" top Streco mass= "<<top_Streco.M()<<endl;


      /*______________ANGLE RECONSTRUCTION______________*/

      
      
      sinTheta = calculate_sinTheta(Wboson, top, qSpec);
      cosTheta = calculate_cosTheta(Wboson, top, qSpec);

      cosThetaStar = calculate_cosThetaStar(SelectedGenDressedLepton_lead, Wboson);
      sinThetaStar = calculate_sinThetaStar(SelectedGenDressedLepton_lead, Wboson);

      PhiStar = calculate_PhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);

      sinPhiStar = calculate_sinPhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);
      cosPhiStar = calculate_cosPhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);

      if(isRealSolution)
      {
         nbOfRealSolutions += 1;
         // std::cout<<"real solution"<<std::endl;
         sinTheta_Real = calculate_sinTheta(Wboson, top, qSpec);
         cosTheta_Real = calculate_cosTheta(Wboson, top, qSpec);

         cosThetaStar_Real = calculate_cosThetaStar(SelectedGenDressedLepton_lead, Wboson);
         sinThetaStar_Real = calculate_sinThetaStar(SelectedGenDressedLepton_lead, Wboson);

         PhiStar_Real = calculate_PhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);

         sinPhiStar_Real = calculate_sinPhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);
         cosPhiStar_Real = calculate_cosPhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);

         top_pt_Real = top.Pt();
         top_eta_Real = top.Eta();
         top_phi_Real = top.Phi();
         top_mass_Real = top.M();

         neutrino_pt_Real = neutrino.Pt();
         neutrino_eta_Real = neutrino.Eta();
         neutrino_phi_Real = neutrino.Phi();

         tOutput_Real->Fill();
      }

      if(!isRealSolution)
      {
         nbOfComplexSolutions +=1;
         // std::cout<<"complex solution"<<std::endl;
         sinTheta_Complex = calculate_sinTheta(Wboson, top, qSpec);
         cosTheta_Complex = calculate_cosTheta(Wboson, top, qSpec);

         cosThetaStar_Complex = calculate_cosThetaStar(SelectedGenDressedLepton_lead, Wboson);
         sinThetaStar_Complex = calculate_sinThetaStar(SelectedGenDressedLepton_lead, Wboson);

         PhiStar_Complex = calculate_PhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);

         sinPhiStar_Complex = calculate_sinPhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);
         cosPhiStar_Complex = calculate_cosPhiStar(SelectedGenDressedLepton_lead, Wboson, qSpec);

         top_pt_Complex = top.Pt();
         top_eta_Complex = top.Eta();
         top_phi_Complex = top.Phi();
         top_mass_Complex = top.M();

         neutrino_pt_Complex = neutrino.Pt();
         neutrino_eta_Complex = neutrino.Eta();
         neutrino_phi_Complex = neutrino.Phi();

         tOutput_Complex->Fill();
      }
      
      top_pt = top.Pt();
      top_eta = top.Eta();
      top_phi = top.Phi();
      top_mass = top.M();

      neutrino_pt = neutrino.Pt();
      neutrino_eta = neutrino.Eta();
      neutrino_phi = neutrino.Phi();

      W_pt = Wboson.Pt();
      W_mass = Wboson.M();
      W_transverse_mass = Wboson.Mt();

      lepton_pt = SelectedGenDressedLepton_lead.Pt();


      tOutput->Fill();
   }

   // int fractionOfRealSolutions = nbOfRealSolutions / (nbOfComplexSolutions + nbOfRealSolutions);
   // std::cout<<"Fraction of Real soltuions: "<<fractionOfRealSolutions<<std::endl;
   std::cout<<"nbOfRealSolutions = "<<nbOfRealSolutions<<std::endl;
   std::cout<<"nbOfComplexSolutions = "<<nbOfComplexSolutions<<std::endl;
   // std::cout<<"\% of Real solutions = "<<nbOfRealSolutions/(nbOfComplexSolutions+nbOfRealSolutions)<<std::endl;
   // std::cout<<"\% of Real solutions Streco = "<<nbOfRealSolutions_Streco/(nbOfComplexSolutions_Streco+nbOfRealSolutions_Streco)<<std::endl;

   tOutput->Write(0,TObject::kWriteDelete);
   tOutput_Real->Write(0,TObject::kWriteDelete);
   tOutput_Complex->Write(0,TObject::kWriteDelete);
   fOutput->Close();
}
