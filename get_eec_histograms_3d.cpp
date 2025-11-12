#include "tTree.h"
#include "binning_histos_small.h"
#include "TH2.h"
#include "TH1.h"
#include <iostream>

//Skip MC events that have a too large weight to stabilize the distribution
bool skipMC(double jtpt, double pthat) {//double refpt
  //if (!(refpt>0)) return true;    
  if (pthat<0.35*jtpt) return true;
  return false;
}

//Print a vector (for debugging)
void print_vector(std::vector<Int_t> &vec){
  for(Int_t i = 0; i < vec.size(); i++){
    std::cout << vec.at(i) << " ";
  }
  std::cout << std::endl;
}

//Check that the same track is not twice in the list of track vectors (for debugging)
bool check_twice(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors){

  Int_t ntrk = trackVectors.size();

  for(Int_t i = 0; i < ntrk; i++){

    ROOT::Math::PtEtaPhiMVector v;
    v = trackVectors.at(i);

    for(Int_t j = 0; j < ntrk; j++){
      if((i != j) && (v == trackVectors.at(j))) return true;
    }
  }

  return false;
}


//Aggregate tracks using BDT score (reco level)
Float_t ReconstuctSingleB(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){

  //Define b hadron track vector
  ROOT::Math::PtEtaPhiMVector v;
  //Empty vector for a check
  ROOT::Math::PtEtaPhiMVector iszero;

  // Loop over all the tracks in the event
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
    //Check that the track is in the jet i 
    if (t.trkJetId[itrk] != ijet) continue; 

    //Define track vector
    ROOT::Math::PtEtaPhiMVector v1;

    //Get track information
    v1.SetEta(t.trkEta[itrk]);
    v1.SetPt(t.trkPt[itrk]);
    v1.SetPhi(t.trkPhi[itrk]);

    if(std::abs(t.trkPdgId[itrk])==211){
      v1.SetM(0.139570);
    }
    if(std::abs(t.trkPdgId[itrk])==13){
      v1.SetM(0.105658);
    }
    if(std::abs(t.trkPdgId[itrk])==11){
      v1.SetM(0.000510);
    }
    if(std::abs(t.trkPdgId[itrk])==2212){
      v1.SetM(0.938272);
    }
    if(std::abs(t.trkPdgId[itrk])==321){
      v1.SetM(0.493677);
    }

    //If the track is (likely) NOT a b hadron daughter 
    if (t.trkBdtScore[itrk]<=-0.95){
      //save as it is
      trackVectors.push_back(v1);
    }
    //If the track is a b daughter add it to the b hadron
    else{
      v += v1;
    }
	
  }

  //Check that a b hadron was reconstructed and save
  if(v != iszero) trackVectors.push_back(v);

  //Return the reconstrcuted b hadron mass
  return v.M();
}

//Aggregate with ideal aggregation (all true b daughters clustered, reco level)
Float_t ReconstuctSingleB_ideal(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){
  //Define b hadron track vector
  ROOT::Math::PtEtaPhiMVector v;
  //Empty vector for a check
  ROOT::Math::PtEtaPhiMVector iszero;

  // Loop over all the tracks in the event
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
    //Check that the tracks is in the jet i 
    if (t.trkJetId[itrk] != ijet) continue; 

    //Define track vector
    ROOT::Math::PtEtaPhiMVector v1;

    //Get track information
    v1.SetEta(t.trkEta[itrk]);
    v1.SetPt(t.trkPt[itrk]);
    v1.SetPhi(t.trkPhi[itrk]);

    if(std::abs(t.trkPdgId[itrk])==211){
      v1.SetM(0.139570);
    }
    if(std::abs(t.trkPdgId[itrk])==13){
      v1.SetM(0.105658);
    }
    if(std::abs(t.trkPdgId[itrk])==11){
      v1.SetM(0.000510);
    }
    if(std::abs(t.trkPdgId[itrk])==2212){
      v1.SetM(0.938272);
    }
    if(std::abs(t.trkPdgId[itrk])==321){
      v1.SetM(0.493677);
    }

    //If the track is NOT a true b hadron daughter
    if (t.trkMatchSta[itrk]<100){
      //save as it is
      trackVectors.push_back(v1);
    }
    //If the track is a true b daughter add to the b hadron
    else{
      v += v1;
    }
	
  }

  //Check that a b hadron was reconstructed and save
  if(v != iszero) trackVectors.push_back(v);

  //Return the reconstrcuted b hadron mass
  return v.M();
}

//Aggregation at gen level
Float_t ReconstuctSingleB_gen(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){
  //Define b hadron vector
  ROOT::Math::PtEtaPhiMVector v;
  //Empty vector for a check
  ROOT::Math::PtEtaPhiMVector iszero;

  // Loop over all the tracks in the event
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    //Check that the track is in jet i
    if (t.refTrkJetId[itrk] != ijet) continue; 

    //Define track vector
    ROOT::Math::PtEtaPhiMVector v1;

    //Get information
    v1.SetEta(t.refTrkEta[itrk]);
    v1.SetPt(t.refTrkPt[itrk]);
    v1.SetPhi(t.refTrkPhi[itrk]);

    if(std::abs(t.refTrkPdgId[itrk])==211){
      v1.SetM(0.139570);
    }
    if(std::abs(t.refTrkPdgId[itrk])==13){
      v1.SetM(0.105658);
    }
    if(std::abs(t.refTrkPdgId[itrk])==11){
      v1.SetM(0.000510);
    }
    if(std::abs(t.refTrkPdgId[itrk])==2212){
      v1.SetM(0.938272);
    }
    if(std::abs(t.refTrkPdgId[itrk])==321){
      v1.SetM(0.493677);
    }

    //If the track is NOT a b hadron daughter
    if (t.refTrkSta[itrk]<100){
      //save as it is
      trackVectors.push_back(v1);
    }
    //If it is a b daughter add to the b hadron
    else{
      v += v1;
    }
	
  }

  //Check that a b hadron was reconstructed and save
  if(v != iszero) trackVectors.push_back(v);

  //Return b hadron mass
  return v.M();
}

//Create a vector with all tracks without aggregation at reco
void create_trackVectors_reco(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){

  //Loop over the track in the event
  for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
    //Check that the track is in jet i
    if (t.trkJetId[itrk] != ijet) continue; 

    //Define the track vector
    ROOT::Math::PtEtaPhiMVector v1;

    //Get track information
    v1.SetEta(t.trkEta[itrk]);
    v1.SetPt(t.trkPt[itrk]);
    v1.SetPhi(t.trkPhi[itrk]);

    if(std::abs(t.trkPdgId[itrk])==211){
      v1.SetM(0.139570);
    }
    if(std::abs(t.trkPdgId[itrk])==13){
      v1.SetM(0.105658);
    }
    if(std::abs(t.trkPdgId[itrk])==11){
      v1.SetM(0.000510);
    }
    if(std::abs(t.trkPdgId[itrk])==2212){
      v1.SetM(0.938272);
    }
    if(std::abs(t.trkPdgId[itrk])==321){
      v1.SetM(0.493677);
    }

    //Save track
    trackVectors.push_back(v1);
	
  }

}

//Create a vector with all tracks without aggregation at reco
void create_trackVectors_gen(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){

  //Loop over all the tracks in the event
  for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
    //Check that the track is in jet i
    if (t.refTrkJetId[itrk] != ijet) continue; 

    //Define track vector
    ROOT::Math::PtEtaPhiMVector v1;
	
    //Save track information
    v1.SetEta(t.refTrkEta[itrk]);
    v1.SetPt(t.refTrkPt[itrk]);
    v1.SetPhi(t.refTrkPhi[itrk]);

    if(std::abs(t.refTrkPdgId[itrk])==211){
      v1.SetM(0.139570);
    }
    if(std::abs(t.refTrkPdgId[itrk])==13){
      v1.SetM(0.105658);
    }
    if(std::abs(t.refTrkPdgId[itrk])==11){
      v1.SetM(0.000510);
    }
    if(std::abs(t.refTrkPdgId[itrk])==2212){
      v1.SetM(0.938272);
    }
    if(std::abs(t.refTrkPdgId[itrk])==321){
      v1.SetM(0.493677);
    }

    //Save track
    trackVectors.push_back(v1);
	
  }
}


//_____________________________________________________________________________________________
//_______________________________Get EEC histograms at reco mod :)____________________________________
//_____________________________________________________________________________________________


//Create a 3D 2-point EEC(dr) histogram from tree
void do_hist(TString &filename,	TString folder, int &dataType, bool &isMC, Float_t &pT_low, Float_t &pT_high, Int_t &n, bool &btag, bool aggregated = true,  bool ideal_aggr = false){

  TString fin_name = filename;
  tTree t;
  t.Init(fin_name,isMC);
    
  // Turn off all branches and turn on only the interesting ones
  // Attention! If a branch is off, it will return bs without crashing 
  t.SetBranchStatus("*", 0);
  if(!isMC){ //so it's data
    std::vector<TString> active_branches = {
      "jtpt", "jteta", "jtphi", "jtm", "nref", "jtmB",
      "ntrk", "trkPt", "trkJetId","trkMatchSta",
      "trkEta", "trkPhi", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
      "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
    };
    t.SetBranchStatus(active_branches, 1);
  }
  else{
    std::vector<TString> active_branches = {"weight",
      "jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB",
      "ntrk", "trkPt", "trkJetId","trkMatchSta", "refTrkSta",
      "trkEta", "trkPhi", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
      "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
    };
    t.SetBranchStatus(active_branches, 1);
  }

  TH3D *h3D_all, *h3D_b1, *h3D_b2, *h3D_c, *h3D_light; 
    
  //Create a 3D histogram for the mB, dr and pT distributions
  h3D_all  = new TH3D("h3D_all", "h3D_all", bins_mb, mb_binsVector,
		      bins_dr, dr_binsVector,
		      jtpt_bins, jtpt_binsVector);
  h3D_b1 = new TH3D("h3D_b1", "h3D_b1", bins_mb, mb_binsVector,
		    bins_dr, dr_binsVector,
		    jtpt_bins, jtpt_binsVector);
  h3D_b2 = new TH3D("h3D_b2", "h3D_b2", bins_mb, mb_binsVector,
		    bins_dr, dr_binsVector,
		    jtpt_bins, jtpt_binsVector);
  h3D_c = new TH3D("h3D_c", "h3D_c", bins_mb, mb_binsVector,
		   bins_dr, dr_binsVector,
		   jtpt_bins, jtpt_binsVector);
  h3D_light = new TH3D("h3D_light", "h3D_light", bins_mb, mb_binsVector,
		       bins_dr, dr_binsVector,
		       jtpt_bins, jtpt_binsVector);

    
  //Save the prescale factor (only for 40 GeV trigger)
  double prescale_pf40 = 33.917210;

  std::cout << "Dataset = " << dataType << std::endl;
  std::cout << "Events = " << t.GetEntries() << std::endl;

  //looping over events
  std::cout << "Looping over events" << std::endl;

  int maxEvents = t.GetEntries();

  for (Long64_t ient = 0; ient < maxEvents && ient < t.GetEntries(); ient++) { //ATTENTION

    int mult = maxEvents/10;
      
    // Print progress
    if (ient % mult == 0) {
      std::cout << "entry nb = " << ient << std::endl;
    }
 
    //get tree entry ient
    t.GetEntry(ient); 

	
    //get the nr of tracks
    Int_t n_tracks = t.ntrk;

    //get the MC event weight
    Float_t weight_tree = t.weight;
    //if data no event weight
    if(!isMC) weight_tree = 1;

    //select HLT events in data depending if we are using the HighEG dataset or the LowEG dataset
    if(!isMC && dataType == 0){
      if(!(t.HLT_HIAK4PFJet100_v1 || t.HLT_HIAK4PFJet80_v1)) continue;
    }
	
	
    if(!isMC && dataType == -1){
      if(!((t.HLT_HIAK4PFJet60_v1 == 1 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) ||
	   (t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0))) continue;
    }

    //select HLT events with at least 40 GeV if MC
    if(isMC){
      if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
    }


    // Loop over jets 
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // Skip jets outside tracker 
      if (std::abs(t.jteta[ijet]) > 1.6) continue;

      // Skip large weight events in MC
      if ((isMC) && skipMC(t.jtpt[ijet], t.pthat)) continue;

      bool skip = false;


      //Select jets passing the b-jet tagging (if needed)
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;

      //Select jet phi
      /*if(t.jtphi[ijet] > -0.1) continue;
	if(t.jtphi[ijet] < 1.2) continue;*/


      //Select jet pt - low limit
      if (std::abs(t.jtpt[ijet]) < pT_low) continue;
      //Select jet pt - high limit
      if (std::abs(t.jtpt[ijet]) > pT_high) continue;

	
      //save jet pt, eta and phi
      Float_t jet_pt = t.jtpt[ijet];
      Float_t jet_eta = t.jteta[ijet];
      Float_t jet_phi = t.jtphi[ijet];

      //create a track vector for the i jet
      std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;


      //Define b hadron mass
      Float_t mB = 0;

      //Aggregate tracks (or not) with BDT or with gen information



      if (aggregated){
	mB = ReconstuctSingleB(trackVectors, t, ijet);
	if (ideal_aggr) mB = ReconstuctSingleB_ideal(trackVectors, t, ijet);
      }
      else{
	create_trackVectors_reco(trackVectors, t, ijet);
	//If not aggregated mB is the value saved in the tree
	mB = t.jtmB[ijet];
      }

      // Get the new nr of tracks
      n_tracks = trackVectors.size();
	    

      // Loop over the tracks 
      for (Int_t i = 0; i < n_tracks; i++) { 
		
	//Save track information
	Float_t etai = trackVectors[i].Eta();
	Float_t phii = trackVectors[i].Phi();
	Float_t ipt = trackVectors[i].Pt();


	//Select tracks with pT > 1 GeV
	if(ipt < 1) continue;
		

	// Loop over tracks (no double counting)
	for(Int_t j=0; j < i; j++){
		    
	  //Save track information
	  Float_t etaj = trackVectors[j].Eta();
	  Float_t phij = trackVectors[j].Phi();
	  Float_t jpt = trackVectors[j].Pt();

	  //Select tracks with pT > 1 GeV
	  if(jpt < 1) continue;
		    
		    
	  // Calculate dr
	  Float_t dr = t.calc_dr(etai, phii, etaj, phij);

	  // Calculate the eec weight (pTi*pTj)^n
	  Float_t eec = pow(ipt*jpt, n);

	  //Add prescale weight if necessary
	  if(t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) eec*=prescale_pf40;
		    
	  //Fix the under/overflow
	  if(dr < dr_min) dr = dr_min_fill;
	  if(dr >= dr_max) dr = dr_max_fill;
	  if(mB >= mb_max) mB = mb_max_fill;

	  //Fill

		   
	  h3D_all->Fill(mB, dr, jet_pt, eec*weight_tree);

	  if(isMC){
		    
	    if (t.jtHadFlav[ijet] == 5){
			
	      if (t.jtNbHad[ijet] == 1){
		h3D_b1->Fill(mB, dr, jet_pt, eec*weight_tree);   
	      }
	      else{
		h3D_b2->Fill(mB, dr, jet_pt, eec*weight_tree);
	      }
	    }		    
	    else if (t.jtHadFlav[ijet] == 4){
	      h3D_c->Fill(mB, dr, jet_pt, eec*weight_tree);
	    }
	    else{
	      h3D_light->Fill(mB, dr, jet_pt, eec*weight_tree);
	    }
	  }
	}
      }	   

    }
  }
    
  //Create the fout name depending on the selection
  TString fout_name = "hist_3d_";
 
  TString label;
  if (dataType <= 0) {
    label = "data";
    if(dataType==-1) label += "_lowEn";
    else label+= "_highEn";
  }
  else if (dataType == 1) {
    label = "MC_bjet";
  }
  else if (dataType == 2) {
    label = "MC_dijet";
  }
  else {
    label = "unknown"; // default case
  }

  if(aggregated){
    fout_name += "aggr_BDT_";
    if (ideal_aggr) fout_name += "ideal_aggr_";
  }
  else fout_name += "noaggr_";  


  if(!btag) label += "_notag"; 
 
  fout_name += TString(Form("n%i_",n)) + label + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + ".root";

  //Save histograms and close file
  std::cout << "Creating file: " << fout_name << std::endl;
  TFile *fhist = new TFile(folder+fout_name, "recreate");
  h3D_all->Write();
  h3D_b1->Write();
  h3D_b2->Write();
  h3D_c->Write();
  h3D_light->Write();
  fhist->Close();

}

//Create a 3D 3-point EEC(dr) histogram from tree

void do_hist_e3c(TString &filename, TString folder, int &dataType, bool &isMC, Float_t &pT_low, Float_t &pT_high, Int_t &n,  bool &btag, bool aggregated = true, bool ideal_aggr = false){
  TString fin_name = filename;
  tTree t;
  t.Init(fin_name,isMC);

  // !!!!!!!!
  // Turn off all branches and turn on only the interesting ones
  // Attention! If a branch is off, it will return bs without crashing 
  t.SetBranchStatus("*", 0);
  std::vector<TString> active_branches = {"weight",
    "jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB",
    "ntrk", "trkPt", "trkJetId","trkMatchSta", "refTrkSta",
    "trkEta", "trkPhi", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
    "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
  };
  t.SetBranchStatus(active_branches, 1);


  //create a 3D histogram for the mB, dr and pT distributions
  //  TH3D *h3D = new TH3D("h3D", "h3D", bins_mb, mb_binsVector,
  //	       bins_dr, dr_binsVector,
  //		       jtpt_bins, jtpt_binsVector)
  
  TH3D *h3D_all, *h3D_b1, *h3D_b2, *h3D_c, *h3D_light;                                                                                                                                                                                                                                                                 
  //Create a 3D histogram for the mB, dr and pT distributions
  h3D_all  = new TH3D("h3D_all", "h3D_all", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  h3D_b1 = new TH3D("h3D_b1", "h3D_b1", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  h3D_b2 = new TH3D("h3D_b2", "h3D_b2", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  h3D_c = new TH3D("h3D_c", "h3D_c", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);
  h3D_light = new TH3D("h3D_light", "h3D_light", bins_mb, mb_binsVector, bins_dr, dr_binsVector,jtpt_bins, jtpt_binsVector);

    
  //Save the prescale factor (only for 40 GeV trigger)
  double prescale_pf40 = 33.917210;

  std::cout << "Dataset = " << dataType << std::endl;
  std::cout << "Events = " << t.GetEntries() << std::endl;

  //looping over events
  std::cout << "Looping over events" << std::endl;

  int maxEvents = t.GetEntries();  
  for (Long64_t ient = 0; ient < maxEvents  && ient < t.GetEntries(); ient++) { //ATTENTION I ADDED IENT <1000


    int mult = maxEvents/10;                                                                                                                                                                                                                                                                                             // Print progress                                                                                                                                                                                                                                                                                                  
    if (ient % mult == 0) {                                                                                                                                                                                                                                                                                                  std::cout << "entry nb = " << ient << std::endl;                                                                                                                                                                                                                                                                    }  
    
    //get tree entry ient
    t.GetEntry(ient); 
	
    //get the nr of tracks
    Int_t n_tracks = t.ntrk;

    //get the MC event weight
    Float_t weight_tree = t.weight;
    //if data no event weight
    if(!isMC) weight_tree = 1;

    //select HLT events in data depending if we are using the HighEG dataset or the LowEG dataset
    if(!isMC && dataType == 0){
      if(!(t.HLT_HIAK4PFJet100_v1 || t.HLT_HIAK4PFJet80_v1)) continue;
    }
	
	
    if(!isMC && dataType == -1){
      if(!((t.HLT_HIAK4PFJet60_v1 == 1 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) ||
	   (t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0))) continue;
    }

    //select HLT events with at least 40 GeV if MC
    if(isMC){
      if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
    }


    // Loop over jets 
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // Skip jets outside tracker 
      if (std::abs(t.jteta[ijet]) > 1.6) continue;

      // Skip large weight events in MC
      if ((isMC) && skipMC(t.jtpt[ijet], t.pthat)) continue;

  
      //Select jets passing the b-jet tagging (if needed)
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;

      //Select jet phi
      /*if(t.jtphi[ijet] > -0.1) continue;
	if(t.jtphi[ijet] < 1.2) continue;*/


      //Select jet pt - low limit
      if (std::abs(t.jtpt[ijet]) < pT_low) continue;
      //Select jet pt - high limit
      if (std::abs(t.jtpt[ijet]) > pT_high) continue;

	
      //save jet pt, eta and phi
      Float_t jet_pt = t.jtpt[ijet];
      Float_t jet_eta = t.jteta[ijet];
      Float_t jet_phi = t.jtphi[ijet];

      //create a track vector for the i jet
      std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;


      //Define b hadron mass
      Float_t mB = 0;


      
      //Aggregate tracks (or not) with BDT or with gen information
      if (aggregated) {mB = ReconstuctSingleB(trackVectors, t, ijet);
	if (ideal_aggr) mB = ReconstuctSingleB_ideal(trackVectors, t, ijet);
      }

      else{
	create_trackVectors_reco(trackVectors, t, ijet);
	//If not aggregated mB is the value saved in the tree
	mB = t.jtmB[ijet];
      }

      // Get the new nr of tracks
      n_tracks = trackVectors.size();
	    

      // Loop over the tracks 
      for (Int_t i = 0; i < n_tracks; i++) { 
		
	Float_t etai = trackVectors[i].Eta();
	Float_t phii = trackVectors[i].Phi();
	Float_t ipt = trackVectors[i].Pt();


	//Select tracks with pT > 1 GeV
	if(ipt < 1) continue;
		

	// Loop over tracks (no double counting)
	for(Int_t j=0; j < i; j++){
		    
	  Float_t etaj = trackVectors[j].Eta();
	  Float_t phij = trackVectors[j].Phi();
	  Float_t jpt = trackVectors[j].Pt();

	  //Select tracks with pT > 1 GeV
	  if(jpt < 1) continue;
		    
	  // Loop over third track (no double counting)
	  for(Int_t k=0; k < j; k++){

	    Float_t etak = trackVectors[k].Eta();
	    Float_t phik = trackVectors[k].Phi();
	    Float_t kpt = trackVectors[k].Pt();

	    //Select tracks with pT > 1 GeV
	    if(kpt < 1) continue;
		    
	    //calculate dr and fill hisotgram
	    Float_t dr_ij = t.calc_dr(etai, phii, etaj, phij);
	    Float_t dr_ik = t.calc_dr(etai, phii, etak, phik);
	    Float_t dr_jk = t.calc_dr(etaj, phij, etak, phik);

	    //take the maximum as dr (see CMS paper)
	    Float_t dr = std::max({dr_ij, dr_ik, dr_jk});
		
	    // Calculate the eec
	    Float_t eec = pow(ipt*jpt*kpt, n);

		    
	    //Add prescale weight if necessary
	    if(t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) eec*=prescale_pf40;
		    
	    //Fix the under/overflow
	    if(dr < dr_min) dr = dr_min_fill;
	    if(dr >= dr_max) dr = dr_max_fill;
	    if(mB >= mb_max) mB = mb_max_fill;

	    //Fill
	    //	    h3D->Fill(mB, dr, jet_pt, eec*weight_tree);

	    h3D_all->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                                                                                                                                                                                        if(isMC){                                                                                                                                                                                                                                                                                                               if (t.jtHadFlav[ijet] == 5){                                                                                                                                                                                                                                                                                           if (t.jtNbHad[ijet] == 1){                                                                                                                                                                                                                                                                                             h3D_b1->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                                                                                                                                                                                    }                                                                                                                                                                                                                                                                                                        
		else{                                                                                                                                                                                                                                                                                                                   h3D_b2->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                                                                                                                                                                                      }                                                                                                                                                                                                                                                                                                                   }                                                                                                                                                                                                                                                                                                                     else if (t.jtHadFlav[ijet] == 4){                                                                                                                                                                                                                                                                                       h3D_c->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                                                                                                                                                                                       }                                                                                                                                                                                                                                                                                                                     else{                                                                                                                                                                                                                                                                                                                   h3D_light->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                                                                                                                                                                                   }                                                                                                                                                                                                                                                                                                       	        }                                                                                                                                                                                                                                                                                                                   }           
	}
      }
    }		   
  }
  
  

  //Create the fout name depending on the selection
  TString fout_name = "hist_3d_e3c_";
  
  TString label;                                                                                                                                                                                                                                                                                                       
  if (dataType == -1) {
    label = "data_lowEn"; 
  }
  else if (dataType == 0) {
    label = "data_highEn";
  }
  else if (dataType == 1) {
    label = "MC_bjet";
  }
  else if (dataType == 2) {
    label = "MC_dijet";
  }
  else {
    label = "unknown"; // default case
  }
  if(aggregated){
    fout_name += "aggr_BDT_";
    if (ideal_aggr) fout_name += "ideal_aggr_";
  }
  else fout_name += "noaggr_";
  
  if(!btag) label += "_notag"; 
  
  fout_name += TString(Form("n%i_",n)) + label + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + ".root";      

//Save histograms and close file                                                                                                                                                                                                                                                                                     
  std::cout << "Creating file: " << fout_name << std::endl;
  TFile *fhist = new TFile(folder+fout_name, "recreate");
  h3D_all->Write();
  h3D_b1->Write();
  h3D_b2->Write();
  h3D_c->Write();
  h3D_light->Write();
  fhist->Close();   
}

//_____________________________________________________________________________________________
//_______________________________Get EEC histograms at gen_____________________________________
//_____________________________________________________________________________________________

//Create a 3D 2-point EEC(dr) histogram from tree at gen level
void do_hist_gen(TString &filename, TString folder, int &dataType, bool &isMC, Float_t &pT_low, Float_t &pT_high, Int_t &n, bool &btag, bool aggregated = true,  bool ideal_aggr = false){
//void do_hist_gen(TString &filename,  TString &dataset, Float_t &pT_low, Float_t &pT_high, bool &aggregated, Int_t &cuts, TString &label, bool &btag, Int_t &n, TString folder){

  if (!isMC) {                                                                                                                                                                                                                                                                                                            std::cerr << "Error: This function should not be called with data. Process terminated." << std::endl;
    return;
  }
  
TString fin_name = filename;
tTree t;
t.Init(fin_name,isMC);
// Turn off all branches and turn on only the interesting ones
// Attention! If a branch is off, it will return bs without crashing 
t.SetBranchStatus("*", 0);
std::vector<TString> active_branches = {"weight",
  "nref", "refmB", 
  "refpt", "refeta", "refphi", "nref",
  "nrefTrk", "refTrkPt", "refTrkJetId",
  "refTrkEta", "refTrkPhi", "refTrkPdgId", "refTrkSta", "refTrkMass",
  "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore",
  "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
};
t.SetBranchStatus(active_branches, 1);

TH3D *h3D_all, *h3D_b1, *h3D_b2, *h3D_c, *h3D_light;                                                                                                               
                                                                                                                                                                       
//Create a 3D histogram for the mB, dr and pT distributions                                                                                                        
h3D_all  = new TH3D("h3D_all", "h3D_all", bins_mb, mb_binsVector,                                                                                                  
		    bins_dr, dr_binsVector,                                                                                                         
		    jtpt_bins, jtpt_binsVector);                                                                                                    
h3D_b1 = new TH3D("h3D_b1", "h3D_b1", bins_mb, mb_binsVector,                                                                                                      
		  bins_dr, dr_binsVector,                                                                                                         
		  jtpt_bins, jtpt_binsVector);                                                                                                    
h3D_b2 = new TH3D("h3D_b2", "h3D_b2", bins_mb, mb_binsVector,                                                                                                      
		  bins_dr, dr_binsVector,                                                                                                         
		  jtpt_bins, jtpt_binsVector);                                                                                                    
h3D_c = new TH3D("h3D_c", "h3D_c", bins_mb, mb_binsVector,                                                                                                         
		 bins_dr, dr_binsVector,                                                                                                         
		 jtpt_bins, jtpt_binsVector);                                                                                                    
h3D_light = new TH3D("h3D_light", "h3D_light", bins_mb, mb_binsVector,                                                                                             
		     bins_dr, dr_binsVector,                                                                                                         
		     jtpt_bins, jtpt_binsVector);                                                                                                    

//Save the prescale factor (only for 40 GeV trigger)                                                                                                               
double prescale_pf40 = 33.917210;                                                                                                                                  
std::cout << "Dataset = " << dataType << std::endl;                                                                                                                
std::cout << "Events = " << t.GetEntries() << std::endl;                                                                                                           
//looping over events                                                                                                                                              
std::cout << "Looping over events" << std::endl;                                                                                                                   
/*
    
//create a 3D histogram for the mB, dr and pT distributions
TH3D *h3D = new TH3D("h3D", "h3D", bins_mb, mb_binsVector,
bins_dr, dr_binsVector,
jtpt_bins, jtpt_binsVector);

    
//Save the prescale factor (only for 40 GeV trigger)
double prescale_pf40 = 33.917210;

std::cout << "Dataset = " << dataset << ", gen" << std::endl;
std::cout << "Selection = " << label << std::endl;
std::cout << "Events = " << t.GetEntries() << std::endl;

//looping over events
std::cout << "Looping over events" << std::endl;
*/    

int maxEvents = t.GetEntries();;                                                                                                                                             

for (Long64_t ient = 0; ient < maxEvents && ient < t.GetEntries(); ient++) { //ATTENTION                                                                           
  
  int mult = maxEvents/10;                                                                                                                                         
  
  // Print progress                                                                                                                                                
  if (ient % mult == 0) {                                                                                                                                          
    std::cout << "entry nb = " << ient << std::endl;                                                                                                           
  }
  t.GetEntry(ient); 
  
  //Get the MC event weight
  Float_t weight_tree = t.weight;
	
  //get the nr of tracks
  Int_t n_tracks = t.nrefTrk;
  
  //Select MC passing a trigger of at least 40 GeV
  if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
  
  // Loop over jets 
  for (Int_t ijet = 0; ijet < t.nref; ijet++) {
    
    // Skip jets outside tracker 
    if (std::abs(t.refeta[ijet]) > 1.6) continue;

    bool skip = false;

    //Select jets passing the b-jet tagging
    if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;

    //Select jet phi
    /*if(t.jtphi[ijet] > -0.1) continue;
      if(t.jtphi[ijet] < 1.2) continue;*/


    //Select jet pt - low limit
    if (std::abs(t.refpt[ijet]) < pT_low) continue;
    //Select jet pt - high limit
    if (std::abs(t.refpt[ijet]) > pT_high) continue;

    //save jet pt, eta and phi
    Float_t jet_pt = t.refpt[ijet];
    Float_t jet_eta = t.refeta[ijet];
    Float_t jet_phi = t.refphi[ijet];

    //create a track vector for the i jet
    std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;

    // Aggregate 1 b hadron, get mB
    Float_t mB = 0;
    if (aggregated) mB = ReconstuctSingleB_gen(trackVectors, t, ijet);
    // Otherwise fill all the tracks in the jet and use the mB saved in the tree
    else{
      create_trackVectors_gen(trackVectors, t, ijet);
      mB = t.jtmB[ijet];
    }

    // Set the new nr of tracks
    n_tracks = trackVectors.size();

    // Loop over the tracks 
    for (Int_t i = 0; i < n_tracks; i++) { 
		
      Float_t etai = trackVectors[i].Eta();
      Float_t phii = trackVectors[i].Phi();
      Float_t ipt = trackVectors[i].Pt();


      //Select tracks with pT > 1 GeV
      if(ipt < 1) continue;

      // Loop over tracks (no double counting)
      for(Int_t j=0; j < i; j++){
		    
	Float_t etaj = trackVectors[j].Eta();
	Float_t phij = trackVectors[j].Phi();
	Float_t jpt = trackVectors[j].Pt();

	//Select tracks with pT > 1 GeV
	if(jpt < 1) continue;
		    
	// Calculate dr
	Float_t dr = t.calc_dr(etai, phii, etaj, phij);

	// Calculate the eec
	Float_t eec = pow(ipt*jpt, n);

		    
	//Add prescale weight if necessary
	if(t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) eec*=prescale_pf40;
		    

	//Fix the under/overflow
	if(dr < dr_min) dr = dr_min_fill;
	if(dr >= dr_max) dr = dr_max_fill;
	if(mB >= mb_max) mB = mb_max_fill;

	//Fill
	h3D_all->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                    
	if(isMC){                                                                                                                                          
	  if (t.jtHadFlav[ijet] == 5){                                                                                                                     
	    if (t.jtNbHad[ijet] == 1){                                                                                                                     
	      h3D_b1->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                               
	    }                                                                                                                                              
	    else{                                                                                                                                          
	      h3D_b2->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                               
	    }                                                                                                                                              
	  }                                                                                                                                                
	  else if (t.jtHadFlav[ijet] == 4){                                                                                                                
	    h3D_c->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                                  
	  }                                                                                                                                                
	  else{                                                                                                                                            
	    h3D_light->Fill(mB, dr, jet_pt, eec*weight_tree);                                                                                              
	  }                                                                                                                                                
	}                                                                                                                                                  
      }

      
	//	h3D->Fill(mB, dr, jet_pt, eec*weight_tree);
    }
  }		   

 }



//Create the fout name depending on the selection
TString fout_name = "hist_3d_gen_";
TString label;                                                                                                                                                        
 if (dataType == 1) {
   label = "MC_bjet";
 }
 else if (dataType == 2) {
   label = "MC_dijet";
 }                                                                                                                                                                     
 else {
   label = "unknown";
 }     
 
 if(aggregated) fout_name += "aggr_";
 else fout_name += "noaggr_";
 
 if(!btag) label += "_notag"; 

 fout_name += TString(Form("n%i_",n)) + label + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + ".root";                                                      

//Save histograms and close file
 std::cout << "Creating file: " << fout_name << std::endl;
 TFile *fhist = new TFile(folder+fout_name, "recreate");
 h3D_all->Write();
 h3D_b1->Write();
 h3D_b2->Write();
 h3D_c->Write();
 h3D_light->Write();
 fhist->Close();
}   



//Save histograms and close file
//std::cout << "Creating file: " << fout_name << std::endl;
//TFile *fhist = new TFile(folder+fout_name, "recreate");
//h3D->Write();
//fhist->Close();
//}

//Create a 3D 3-point EEC(dr) histogram from tree at gen level
void do_hist_e3c_gen(TString &filename,	TString folder, int &dataType, bool isMC, Float_t &pT_low, Float_t &pT_high, Int_t &n, bool &btag, bool aggregated = true){
  if (!isMC) {
    std::cerr << "Error: This function should not be called with data. Process terminated." << std::endl;
    return;  
  }

  TString fin_name = filename;//
  tTree t;
  t.Init(fin_name,isMC);
      
  // Turn off all branches and turn on only the interesting ones
  // Attention! If a branch is off, it will return bs without crashing 
  t.SetBranchStatus("*", 0);
  std::vector<TString> active_branches = {"weight",
    "nref", "refmB", 
    "refpt", "refeta", "refphi", "nref",
    "nrefTrk", "refTrkPt", "refTrkJetId",
    "refTrkEta", "refTrkPhi", "refTrkPdgId", "refTrkSta", "refTrkMass",
    "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore",
    "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
  };
  t.SetBranchStatus(active_branches, 1);
    

  TH3D *h3D_all, *h3D_b1, *h3D_b2, *h3D_c, *h3D_light;                                                                                                                                                                                                                                                                   
  //Create a 3D histogram for the mB, dr and pT distributions
  h3D_all  = new TH3D("h3D_all", "h3D_all", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h3D_b1 = new TH3D("h3D_b1", "h3D_b1", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h3D_b2 = new TH3D("h3D_b2", "h3D_b2", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h3D_c = new TH3D("h3D_c", "h3D_c", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
  h3D_light = new TH3D("h3D_light", "h3D_light", bins_mb, mb_binsVector, bins_dr, dr_binsVector, jtpt_bins, jtpt_binsVector);
      
  //Save the prescale factor (only for 40 GeV trigger)
  double prescale_pf40 = 33.917210;

  std::cout << "Dataset = " << dataType << ", gen" << std::endl;
  std::cout << "Events = " << t.GetEntries() << std::endl;

  //looping over events
  std::cout << "Looping over events" << std::endl;


int maxEvents = t.GetEntries();;                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                       
for (Long64_t ient = 0; ient < maxEvents && ient < t.GetEntries(); ient++) { //ATTENTION
  int mult = maxEvents/10;
  if (ient % mult == 0) {
    std::cout << "entry nb = " << ient << std::endl;
  } 

    t.GetEntry(ient); 

    //Get the MC event weight
    Float_t weight_tree = t.weight;

    //get the nr of tracks
    Int_t n_tracks = t.nrefTrk;
	
    //Select MC passing a trigger of at least 40 GeV
    if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;

    // Loop over jets 
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // Skip jets outside tracker 
      if (std::abs(t.refeta[ijet]) > 1.6) continue;

       //Select jets passing the b-jet tagging
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;

      //Select jet phi
      /*if(t.jtphi[ijet] > -0.1) continue;
	if(t.jtphi[ijet] < 1.2) continue;*/


      //Select jet pt - low limit
      if (std::abs(t.refpt[ijet]) < pT_low) continue;
      //Select jet pt - high limit
      if (std::abs(t.refpt[ijet]) > pT_high) continue;

      //save jet pt, eta and phi
      Float_t jet_pt = t.refpt[ijet];
      Float_t jet_eta = t.refeta[ijet];
      Float_t jet_phi = t.refphi[ijet];

      //create a track vector for the i jet
      std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;
    
      // Aggregate 1 b hadron, get mB
      Float_t mB = 0;
      if (aggregated) mB = ReconstuctSingleB_gen(trackVectors, t, ijet);
      // Otherwise fill all the tracks in the jet and use the mB saved in the tree
      else{
	create_trackVectors_gen(trackVectors, t, ijet);
	mB = t.jtmB[ijet];
      }

      // Set the new nr of tracks
      n_tracks = trackVectors.size();

      // Loop over the tracks 
      for (Int_t i = 0; i < n_tracks; i++) { 
		
	Float_t etai = trackVectors[i].Eta();
	Float_t phii = trackVectors[i].Phi();
	Float_t ipt = trackVectors[i].Pt();


	//Select tracks with pT > 1 GeV
	if(ipt < 1) continue;
		

	// Loop over tracks (no double counting)
	for(Int_t j=0; j < i; j++){
		    
	  Float_t etaj = trackVectors[j].Eta();
	  Float_t phij = trackVectors[j].Phi();
	  Float_t jpt = trackVectors[j].Pt();

	  //Select tracks with pT > 1 GeV
	  if(jpt < 1) continue;
		    
	  // Loop over third track (no double counting)
	  for(Int_t k=0; k < j; k++){

	    Float_t etak = trackVectors[k].Eta();
	    Float_t phik = trackVectors[k].Phi();
	    Float_t kpt = trackVectors[k].Pt();

	    //Select tracks with pT > 1 GeV
	    if(kpt < 1) continue;
		    
	    //calculate dr and fill hisotgram
	    Float_t dr_ij = t.calc_dr(etai, phii, etaj, phij);
	    Float_t dr_ik = t.calc_dr(etai, phii, etak, phik);
	    Float_t dr_jk = t.calc_dr(etaj, phij, etak, phik);

	    //take the maximum as dr (see CMS paper)
	    Float_t dr = std::max({dr_ij, dr_ik, dr_jk});
		
	    // Calculate the eec
	    Float_t eec = pow(ipt*jpt*kpt, n);

		    
	    //Add prescale weight if necessary
	    if(t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) eec*=prescale_pf40;
		    
	    //Fix the under/overflow
	    if(dr < dr_min) dr = dr_min_fill;
	    if(dr >= dr_max) dr = dr_max_fill;
	    if(mB >= mb_max) mB = mb_max_fill;

	    //Fill
	    h3D_all->Fill(mB, dr, jet_pt, eec*weight_tree);
	    if(isMC){
	      if (t.jtHadFlav[ijet] == 5){
		if (t.jtNbHad[ijet] == 1){
		  h3D_b1->Fill(mB, dr, jet_pt, eec*weight_tree);
		}
		else{
		  h3D_b2->Fill(mB, dr, jet_pt, eec*weight_tree);
		}
	      }
	      else if (t.jtHadFlav[ijet] == 4){
		h3D_c->Fill(mB, dr, jet_pt, eec*weight_tree);
	      }
	      else{
		h3D_light->Fill(mB, dr, jet_pt, eec*weight_tree);
		
	      }                                                                                                                                                                                                                                                                                                            
	    }                                                                                                                                                                                                                                                                                                              
	  }  
	  
	}
      }
    }		      
 }
 
 
 
//Create the fout name depending on the selection
  TString fout_name = "hist_3d_e3c_gen_";
  TString label;
  if (dataType == 1) {
    label = "MC_bjet";
  }
  else if (dataType == 2) {
    label = "MC_dijet";
  }
  else {
    label = "unknown";
  }
  
  if(aggregated) fout_name += "aggr_";
  else fout_name += "noaggr_";

  if(!btag) label += "_notag"; 

 fout_name += TString(Form("n%i_",n)) + label + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + ".root";                                                                                                                                                                                                       
//Save histograms and close file
 std::cout << "Creating file: " << fout_name << std::endl;
 TFile *fhist = new TFile(folder+fout_name, "recreate");
 h3D_all->Write();
 h3D_b1->Write();
 h3D_b2->Write();
 h3D_c->Write();
 h3D_light->Write();
 fhist->Close();
}

//____________________________________________________________________________________________________________
//_______________________Functions to calculate efficiencies__________________________________________________
//____________________________________________________________________________________________________________


/*

//Get the signal/background rejection efficiency of the BDT score cut
void get_efficiency_histograms(TString &filename, TString &folder,  TString &dataType, bool &btag){
  TString fin_name = filename;
  tTree t;
  t.Init(fin_name, true);
  
  Double_t pT_high = 120;
  Double_t pT_low = 100;
  
  // Turn off all branches and turn on only the interesting ones
  // Attention! If a branch is off, it will return bs without crashing 
  t.SetBranchStatus("*", 0);
  std::vector<TString> active_branches = {"weight",
    "jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB",
    "ntrk", "trkPt", "trkJetId","trkMatchSta", "refTrkSta",
    "trkEta", "trkPhi", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
    "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
  };
  t.SetBranchStatus(active_branches, 1);
  

  //create the signal (b daughters) and background (non-b daughters) distributions
  TH1D *h_signal = new TH1D("h_signal", "h_signal", 40, -1,1);
  TH1D *h_background = new TH1D("h_background", "h_background", 40, -1,1);


//Save prescale (only for 40 GeV trigger)
  double prescale_pf40 = 33.917210;
  
  std::cout << "Calculating BDT efficiencies" << std::endl;
  std::cout << "Dataset = " << dataset << std::endl;
  std::cout << "Selection = " << label << std::endl;
  std::cout << "Events = " << t.GetEntries() << std::endl;
  
  //looping over events
  std::cout << "Looping over events" << std::endl;
  for (Long64_t ient = 0; ient < 10000 && ient < t.GetEntries(); ient++) {//ATTENTION
    
    // Print progress
    if (ient % 1000000 == 0) {
      std::cout << "entry nb = " << ient << std::endl;
    }
 
    t.GetEntry(ient); 


    //get the nr of tracks
    Int_t n_tracks = t.ntrk;


    //Select MC events passing the 40 GeV cut
    if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;


    // Loop over jets (apply the same cuts as when building the histograms)
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // Skip jets outside tracker 
      if (std::abs(t.jteta[ijet]) > 1.6) continue;

      // Skip large weight events in MC
      if (skipMC(t.jtpt[ijet], t.pthat)) continue;

      bool skip = false;

      // Select jet flavour and/or select on the number of b hadrons
      switch(cuts){
	//b-jet with one b hadron
      case 1: 
	if (t.jtHadFlav[ijet] < 5) skip = true;
	if (t.jtNbHad[ijet] != 1) skip = true;
	break;
	//b-jet with more than 1 b hadron
      case 2:
	if (t.jtHadFlav[ijet] < 5) skip = true;
	if (t.jtNbHad[ijet] < 2) skip = true;
	break;
	//non-b jets
      case 3:
	if (std::abs(t.jtHadFlav[ijet]) == 5) skip = true;
	break;
	//no flavour selection
      case 4:
	skip = false;
	break;
	//c-jets
      case 5:
	if(std::abs(t.jtHadFlav[ijet]) != 4) skip = true;
	break;
	//light (non-b non-c) jets
      case 6:
	if(std::abs(t.jtHadFlav[ijet]) >= 4) skip = true;
	break;
      }

      if (skip) continue;


      //Select jets passing the b-jet tagging (if needed)
      if (btag && std::abs(t.discr_particleNet_BvsAll[ijet]) <= 0.99) continue;

      //Select jet phi
      //if(t.jtphi[ijet] > -0.1) continue;
      //if(t.jtphi[ijet] < 1.2) continue;

      //Select jet pt - low limit
      if (std::abs(t.jtpt[ijet]) < pT_low) continue;
      //Select jet pt - high limit
      if (std::abs(t.jtpt[ijet]) > pT_high) continue;

      //fill the vector
      for (Int_t itrk = 0; itrk < n_tracks; itrk++) {
	//Track is in the jet
	if (t.trkJetId[itrk] != ijet) continue; 
	//Track has a match at gen
	if (t.trkMatchSta[itrk] < 1) continue;
	//Select only pT > 1 GeV
	if (t.trkPt[itrk] < 1) continue;

		
	Double_t bdt_score = t.trkBdtScore[itrk];

	//The track is or not a true b hadron decay daughter
	if(t.trkMatchSta[itrk]<100) h_background->Fill(bdt_score);
	else h_signal->Fill(bdt_score);

      }
    }
  }

  //Create fout name
  TString fout_name = "hist_efficiency_";

  if(!btag) label += "_notag"; 

  fout_name += label + "_" + dataset + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + + ".root";

  //Save histograms
  TFile *fhist = new TFile(folder+fout_name, "recreate");
  h_signal->Write();
  h_background->Write();

  //Print efficiencies (select a binning where one of the bin edges corresponds to your BDT cut value,
  //then integrate from that value to 1 and divide by the full integral)
  std::cout << "Signal efficiency for -0.9 = " << h_signal->Integral(3,40)/h_signal->Integral(1,40) << std::endl;
  std::cout << "Background rejection efficiency for -0.9 = " << 1-(h_background->Integral(3, 40)/h_background->Integral(1,40)) << std::endl;

  std::cout << "Signal efficiency for -0.95 = " << h_signal->Integral(2,40)/h_signal->Integral(1,40) << std::endl;
  std::cout << "Background rejection efficiency for -0.95 = " << 1-(h_background->Integral(2, 40)/h_background->Integral(1, 40)) << std::endl;


  fhist->Close();
}


*/

//_________________________________________Get the b-tagging efficiency
void get_tag_efficiency_histograms(TString &filename,  TString &folder, Int_t dataType){
  TString fin_name = filename;
  tTree t;
  t.Init(fin_name, true);
      
  Double_t pT_low = 100;
  Double_t pT_high = 120;

  // Turn off all branches and turn on only the interesting ones
  // Attention! If a branch is off, it will return bs without crashing 
  t.SetBranchStatus("*", 0);
  std::vector<TString> active_branches = {"weight",
    "jtpt", "pthat", "jteta", "jtphi", "jtm", "nref", "jtmB",
    "ntrk", "trkPt", "trkJetId","trkMatchSta", "refTrkSta",
    "trkEta", "trkPhi", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
    "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
  };
  t.SetBranchStatus(active_branches, 1);
    

  //Create signal (b-jets) and background (non-b jets) distributions
  TH1D *h_signal = new TH1D("h_signal", "h_signal", 100, 0, 1);
  TH1D *h_background = new TH1D("h_background", "h_background", 100, 0, 1);

  //Save prescale (only for 40 GeV trigger)
  double prescale_pf40 = 33.917210;

  std::cout << "Calculating tagging efficiencies" << std::endl;
  std::cout << "Dataset = " << dataType << std::endl;
  std::cout << "Events = " << t.GetEntries() << std::endl;

  //looping over events
  std::cout << "Looping over events" << std::endl;
  for (Long64_t ient = 0;  ient < t.GetEntries(); ient++) { //ATTENTION

    int mult = t.GetEntries()/10;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  // Print progress                                                                                                                                                                                                                                                                                                    
  if (ient % mult == 0) {
    std::cout << "entry nb = " << ient << std::endl;
  }        

    t.GetEntry(ient); 


    //get the nr of tracks
    Int_t n_tracks = t.ntrk;


    //Select MC events passing a trigger of 40 GeV
    if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;
	

    // Loop over jets 
    for (Int_t ijet = 0; ijet < t.nref; ijet++) {

      // Skip jets outside tracker 
      if (std::abs(t.jteta[ijet]) > 1.6) continue;

      //cut
      if (skipMC(t.jtpt[ijet], t.pthat)) continue;

	    
      //Select jet pt - low limit
      if (std::abs(t.jtpt[ijet]) < pT_low) continue;
      //Select jet pt - high limit
      if (std::abs(t.jtpt[ijet]) > pT_high) continue;


      Double_t tag_discr = t.discr_particleNet_BvsAll[ijet];

      //Fill signal and background histograms
      if (std::abs(t.jtHadFlav[ijet]) == 5){
	h_signal->Fill(tag_discr);
      }
      else h_background->Fill(tag_discr);

    }
  }

  TString label;
  if (dataType == 1) {                                                                                                                                                                                                                                                                                            
    label = "MC_bjet";                                                                                                                                                                                                                                                                                                
  }                                                                                                                                                                                                                                                                                                                   
  else if (dataType == 2) {                                                                                                                                                                                                                                                                                          
    label = "MC_dijet";                                                                                                                                                                                                                                                                                               
  }                                                                                                                                                                                                                                                                                                                   
  else {                                                                                                                                                                                                                                                                                                              
    label = "unknown"; // default case                                                                                                                                                                                                                                                                                
  }    

  
  //Create fout name
  TString fout_name = "hist_tag_efficiency_" +  label  + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high))) + + ".root";

  //Save histograms
  TFile *fhist = new TFile(folder+fout_name, "recreate");
  h_signal->Write();
  h_background->Write();


  //Print efficiencies (select a binning where one of the bin edges corresponds to your b-tagging cut value,
  //then integrate from that value to 1 and divide by the full integral)
  std::cout << "B jet tag efficiency for 0.99 = " << h_signal->Integral(99,100)/h_signal->Integral(1,100)*100 << "%" << std::endl;
  std::cout << "Background rejection efficiency for 0.99 = " << 1-(h_background->Integral(99,100)/h_background->Integral(1,100))*100 << "%" << std::endl;


  fhist->Close();
}












// ______ HIST _____
void get_eec_histograms_3d(int dataType = 2,
			   //-1 for data Low //0 for data High //1 for MC - bjet //2 for MC - dijet
			   
			   //Select pT range (later divided into 3 bins, central bin is the nominal one)
			   Float_t pT_low = 80, 
			   Float_t pT_high = 140, 

			   //select eec weight exponent
			   Int_t n=1, 

			   //apply b-tagging or aggregation	 
			   bool btag = false,
			   bool aggregated = true, 
			   bool ideal_aggr = false ){
  
    

  //Folder where to save files
  TString folder = "$mydata/test_for_code_mods/smaller_bins/";
  
  TString filename;
  TString dataset;
  
  bool isMC = true;
  

  if(dataType == -1){  //________________________________data High______________________________
    filename = "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_fixedMassBug/all_merged_HiForestMiniAOD.root";
    dataset = "LowEGJet";
    isMC = false;
    cout<<"you chose data Low" <<endl;
  }

  else if(dataType == 0) {
     filename = "/data_CMS/cms/kalipoliti/bJet2017G/HighEGJet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
     dataset = "HighEGJet";
     isMC = false;
     cout<<"you chose data High" <<endl;                                                                                                                                                                                                                                                                                  }         
  else if(dataType == 1){//________________________________bjet______________________________
    filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    dataset = "bjet";
    cout<<"you chose bjet MC" <<endl;
  }
  else if(dataType == 2){//________________________________dijet______________________________
    filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"; 
    dataset = "dijet";
    cout<<"you chose dijet MC" <<endl;
  }
  else{
    cout<<" undefined data type"<<endl;
    return; 
  }
  
 
    
  cout<<"let's start :)  " <<endl;
  do_hist(filename, folder, dataType, isMC, pT_low, pT_high, n, btag, aggregated, ideal_aggr);
  //do_hist_e3c(filenames, folder, dataType, isMC, pT_low, pT_high, n, btag, aggregated, ideal_aggr);

    if(dataType>0){
      cout<<" do_hist_gen "<<endl;
      do_hist_gen(filename, folder, dataType, isMC, pT_low, pT_high, n, btag, aggregated, ideal_aggr);
      //do_hist_e3c_gen(filenames, folder, dataType, isMC, pT_low, pT_high, n, btag, aggregated);
      //cout<<" get_efficiency_histograms "<<endl;
      //get_efficiency_histograms(filenames.at(i), datasets.at(i), cuts_vec.at(j), labels_vec.at(j), btag, folder);
      //cout<<" get_tag_efficiency_histograms "<<endl;
      //get_tag_efficiency_histograms(filenames.at(i), folder, dataType);
      cout<<" done "<<endl;
    }
}





