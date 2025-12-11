//Creates a tree where the information on jets and the dr and eec are store in separate arrays to
//build the response matrix (works for all dimensions of the unfolding)
#include "tTree.h"
#include "TH3.h"
#include "TSystem.h"
#include "TString.h"
#include "binning_histos_small.h"
#include "Math/GenVector/VectorUtil.h"
#pragma cling load("libGenVector.so")
#include "TSystem.h"

//Skip MC that have a too large event weight
bool skipMC(double jtpt, double pthat) {//double refpt
    //if (!(refpt>0)) return true;    
    if (pthat<0.35*jtpt) return true;
    return false;
}

// Prints a vector of Int_t (for debugging)
void print_vector(std::vector<Int_t> &vec){
    for(Int_t i = 0; i < (Int_t) vec.size(); i++){
        std::cout << vec.at(i) << " ";
    }
    std::cout << std::endl;
}

// Prints a vector of Float_t (for debugging)
void print_vector_float(std::vector<Float_t> &vec){
    for(Int_t i = 0; i < (Int_t) vec.size(); i++){
        std::cout << vec.at(i) << " ";
    }
    std::cout << std::endl;
}

// Prints a 2D vector (matrix) (for debugging)
void print_matrix(std::vector<std::vector<Float_t>> &v, Int_t &ncols, Int_t &nrows){
    for(Int_t r=0; r < nrows; r++){
        for(Int_t c=0; c < ncols; c++){
            std::cout << v.at(r).at(c) << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << "-------------" << std::endl;
}

// Checks if a Lorentz vector is included into the track vectors for aggregation
bool check_inclusion(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, ROOT::Math::PtEtaPhiMVector &v1){
    bool included = false;
    for (Int_t i = 0; i < (Int_t) trackVectors.size(); i++){
        if (trackVectors[i]==v1){
            included=true;
            break;
        }
    }

    return included;
}

//Aggregate tracks using BDT score (reco level)
Float_t ReconstuctSingleB_reco(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){
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

        //If the track is (likely) NOT a b hadron daughter 
        if (t.trkBdtScore[itrk]<=-0.9){
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
Float_t ReconstuctSingleB_reco_ideal(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, Int_t& ijet){
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

//Find the absolute minimum (i,j) in a matrix and set all other entries of row i and col j to 9999
void find_matrix_min(std::vector<std::vector<Float_t>> &v, Int_t &nrows, Int_t &ncols){
    Float_t min = 9999;
    Int_t min_col = 999;
    Int_t min_row = 999;

    //Find absolute minimum
    for(Int_t row_ind = 0; row_ind < nrows; row_ind++){
        for(Int_t col_ind = 0; col_ind < ncols; col_ind++){
            if(v[row_ind][col_ind] != -1 && v[row_ind][col_ind] < min){
                min = v[row_ind][col_ind];
                min_col = col_ind;
                min_row = row_ind;
            }
        }
    }

    //Set absolute minimum value to -1
    v[min_row][min_col] = -1;

    //Set all other values of row min_row and col min_col to 9999
    for(Int_t row_ind = 0; row_ind < nrows; row_ind++){
        if(row_ind != min_row) v[row_ind][min_col] = 9999;
    }

    for(Int_t col_ind = 0; col_ind < ncols; col_ind++){
        if(col_ind != min_col) v[min_row][col_ind] = 9999;
    }
}

//Check that the matrix is only filled with 9999 or -1 (i.e. that we are done with the matching)
bool check_all_matched(std::vector<std::vector<Float_t>> &v, Int_t &nrows,Int_t &ncols){
    
    for(Int_t row_ind = 0; row_ind < nrows; row_ind++){
        for(Int_t col_ind = 0; col_ind < ncols; col_ind++){
            if(v[row_ind][col_ind] != 9999 && v[row_ind][col_ind] != -1) return false;
        }
    }
    return true;
}

// Matching of reco and gen tracks
void match_tracks(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors_reco_return, std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors_reco_notmatched, std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors_gen_return, std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors_gen_notmatched, tTree& t, Int_t& ijet){
    
    //Save the vectors
    std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors_reco;
    std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors_gen;
    trackVectors_reco = trackVectors_reco_return;
    trackVectors_gen = trackVectors_gen_return;
    //Vectors to be returned are first emptied
    trackVectors_reco_return.clear();
    trackVectors_gen_return.clear();
    //Save the track corresponding to the b hadron
    ROOT::Math::PtEtaPhiMVector B_reco;
    ROOT::Math::PtEtaPhiMVector B_gen;

    //save sizes before matching
    Int_t reco_size = trackVectors_reco.size();
    Int_t gen_size = trackVectors_gen.size();

    //save sizes for debugging
    Int_t reco_size_save = reco_size;
    Int_t gen_size_save = gen_size;

    //store the matched indices (for debugging)
    std::vector<Int_t> indices_reco;
    std::vector<Int_t> indices_gen;

    //store the pT differences between all pairs in a matrix
    std::vector<std::vector<Float_t>> pT_diff;


    //Start the matching if there are some tracks
    if(reco_size != 0 && gen_size != 0){
        // Save b hadron tracks and lower the size to avoid matching them
        // (we don't match b hadrons since it is our assumption that we have only 
        //  one per jet and that it will match between reco and gen, b hadrons are
        //  always stored at the end of the vector, see RecostructSingleB)
        B_reco = trackVectors_reco.at(reco_size-1);
        B_gen = trackVectors_gen.at(gen_size-1);
        gen_size -= 1;
        reco_size -= 1;

        //Loop over the gen tracks - row index to build the pT difference matrix
        for(Int_t i_gen = 0; i_gen < gen_size; i_gen++){
            //temporary store
            std::vector<Float_t> pt_vec;

            //Loop over the reco tracks - colum index
            for(Int_t i_reco = 0; i_reco < reco_size; i_reco++){

                //save info
                Float_t eta_reco = trackVectors_reco[i_reco].Eta();
                Float_t eta_gen = trackVectors_gen[i_gen].Eta();
                Float_t phi_reco = trackVectors_reco[i_reco].Phi();
                Float_t phi_gen = trackVectors_gen[i_gen].Phi();
                Float_t pt_reco = trackVectors_reco[i_reco].Pt();
                Float_t pt_gen = trackVectors_gen[i_gen].Pt();
                Float_t pt_ratio = pt_reco/pt_gen;

                Float_t pt_diff = std::abs(pt_reco - pt_gen);

                //If the matching condition is satisfied, save the pT difference
                if(t.calc_dr(eta_reco, phi_reco, eta_gen, phi_gen) < 0.02 && pt_ratio > 0.8 && pt_ratio < 1.2){
                    pt_vec.push_back(pt_diff);
                }
                //else put it to a very high value that signifies no matching 
                //(no -999 because we later search for minima and this would give problems)
                else{
                    pt_vec.push_back(9999);
                }

            }

            //build the corresponding row
            pT_diff.push_back(pt_vec);
        }


        //Find the true minima for each pair by finding each time the absolute minimum until all are matched
        //The non-minimal entries of the matrix are set to 9999, for each pair of indices there will be only one 
        //entry left.
        while(!check_all_matched(pT_diff, gen_size, reco_size)){
            find_matrix_min(pT_diff, gen_size, reco_size);
        }
        
        //Save the reco and gen indices and the tracks matched for which a match was found
        for(Int_t ind_reco = 0; ind_reco < reco_size; ind_reco++){
            for(Int_t ind_gen = 0; ind_gen < gen_size; ind_gen++){
                if(pT_diff[ind_gen][ind_reco] < 9999){
                    indices_gen.push_back(ind_gen);
                    indices_reco.push_back(ind_reco);
                    trackVectors_reco_return.push_back(trackVectors_reco[ind_reco]);
                    trackVectors_gen_return.push_back(trackVectors_gen[ind_gen]);
                }
            }
        }

        //Put back the b hadron in the vectors
        trackVectors_reco_return.push_back(B_reco);
        trackVectors_gen_return.push_back(B_gen);


        //Save the non-matched reco vectors
        for(Int_t ind_reco = 0; ind_reco < reco_size; ind_reco++){

            ROOT::Math::PtEtaPhiMVector v = trackVectors_reco[ind_reco];

            if(!check_inclusion(trackVectors_reco_return, v)){
                trackVectors_reco_notmatched.push_back(v);
            }
        }

        //Save the non-matched gen vectors
        for(Int_t ind_gen = 0; ind_gen < gen_size; ind_gen++){

            ROOT::Math::PtEtaPhiMVector v = trackVectors_gen[ind_gen];

            if(!check_inclusion(trackVectors_gen_return, v)){
                trackVectors_gen_notmatched.push_back(v);
            }
        }

    }

    //check for difference in the matched number of tracks - debug
    reco_size = trackVectors_reco_return.size();
    gen_size = trackVectors_gen_return.size();
    if(reco_size != gen_size){
        std::cout << "reco_size = " << reco_size_save << ", gen_size = " << gen_size_save << std::endl;
        std::cout << "AFTER MATCHING = reco_size = " << reco_size << ", gen_size = " << gen_size << std::endl;
        /*std::cout << "indices_gen" << std::endl;
        print_vector(indices_gen);
        std::cout << "Indices gen erase" << std::endl;
        print_vector(indices_gen_erase);
        std::cout << "indices_reco" << std::endl;
        print_vector(indices_reco);
        std::cout << "Indices reco erase" << std::endl;
        print_vector(indices_reco_erase);*/
        
        std::cout << "-----------------------new jet--------------------------------------" << std::endl;

    }

}

//Create trees storing informations on 2-point EEC to build the response matrix
void do_trees(TString &filename,  Int_t &dataType, TString &label, TString &folder, Int_t &n, Float_t &pT_low, Float_t &pT_high, bool &aggregated,  bool &btag, bool &matching, Int_t &beg_event, Int_t &end_event, const char* output_suffix){


    bool isMC = true;
    if(dataType <= 0) {isMC = false;
    }
    TString fin_name = filename;//
    
    tTree t;
    t.Init(fin_name, isMC);

    //Create the fout name depending on the selection
    TString fout_name = "trees_nocuts_";

    if(matching) fout_name += "matched_";

    if(aggregated) fout_name += "aggr_BDT_";
    else fout_name += "noaggr_";
    if(!btag) label += "_notag"; 

    fout_name += TString(Form("n%i_",n))  + label + "_" + TString(Form("%i_%i",int(pT_low), int(pT_high)))+ output_suffix ;

    //Create output file and tree to store all the values
    TFile *fout = new TFile(folder+fout_name, "recreate");
    TTree *tree = new TTree("tree",   "tree_all_jets");
    
    //Define branches variables
    Int_t ndr_reco, ndr_gen, ndr_reco_tot, ndr_gen_tot;
    

    //for inclusive/large samples you might need to increase the array size if the code misteriously crashes
    Float_t eec_reco[4000], eec_gen[4000], dr_reco[4000], dr_gen[4000];
            
    Float_t jpt_reco, jpt_gen, weight,jt_eta_reco, jt_eta_gen, discr, pthat, mB_reco, mB_gen;

    Int_t jtHadFlav, jtNbHad;

    //Set branches
    //weight of the event
    tree->Branch("weight", &weight, "weight/F");
    //Number of matched dr and eec calculated (after all the cuts)
    tree->Branch("ndr_reco", &ndr_reco, "ndr_reco/I");
    tree->Branch("ndr_gen", &ndr_gen, "ndr_gen/I");
    //Total number of eec (matched and unmatched)
    tree->Branch("ndr_reco_tot", &ndr_reco_tot, "ndr_reco_tot/I");
    tree->Branch("ndr_gen_tot", &ndr_gen_tot, "ndr_gen_tot/I");
    //jet pt
    tree->Branch("jpt_reco", &jpt_reco, "jpt_reco/F");
    tree->Branch("jpt_gen", &jpt_gen, "jpt_gen/F");
    //mb
    tree->Branch("mB_reco", &mB_reco, "mB_reco/F");
    tree->Branch("mB_gen", &mB_gen, "mB_gen/F");
    //jet eta
    tree->Branch("jt_eta_reco", &jt_eta_reco, "jt_eta_reco/F");
    tree->Branch("jt_eta_gen", &jt_eta_gen, "jt_eta_gen/F");
    //btag
    tree->Branch("discr", &discr, "discr/F");
    tree->Branch("pthat", &pthat, "pthat/F");
    //dr array
    tree->Branch("dr_reco", dr_reco, "dr_reco[ndr_reco_tot]/F");
    tree->Branch("dr_gen", dr_gen, "dr_gen[ndr_gen_tot]/F");
    //eec array
    tree->Branch("eec_reco", eec_reco, "eec_reco[ndr_reco_tot]/F");
    tree->Branch("eec_gen", eec_gen, "eec_gen[ndr_gen_tot]/F");
    //jet flavour and number of b hadrons
    tree->Branch("jtHadFlav", &jtHadFlav, "jtHadFlav/I");
    tree->Branch("jtNbHad", &jtNbHad, "jtNbHad/I");


    // Turn off all branches and turn on only the interesting ones
    // Attention! If a branch is off, it will return bs without crashing 
    t.SetBranchStatus("*", 0);
    std::vector<TString> active_branches = {"weight", "pthat",
        "jtm", "nref", "refmB", 
        "refpt", "refeta", "refphi", "nref",
        "nrefTrk", "refTrkPt", "refTrkJetId",
        "refTrkEta", "refTrkPhi", "refTrkPdgId", "refTrkSta", "refTrkMass",
        "jtpt", "jteta", "jtphi", "jtm", "nref", "jtmB", "genpt",
        "ntrk", "trkPt", "trkJetId",
        "trkMatchSta", "refTrkSta",
        "trkEta", "trkPhi", "jtNbHad", "jtHadFlav", "discr_particleNet_BvsAll", "trkBdtScore", "trkPdgId",
        "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1", "HLT_HIAK4PFJet30_v1"
    };
    t.SetBranchStatus(active_branches, 1);
    
    //Prescale value (only for events passing a 40 GeV trigger)
    double prescale_pf40 = 33.917210;

    std::cout << "Dataset = " << dataType << std::endl;
    std::cout << "Selection = " << label << std::endl;
    std::cout << "Events = " << t.GetEntries() << std::endl;

    //For checks on the track efficiency of the matching
    Int_t tot_gen_tracks = 0;
    Int_t tot_reco_tracks = 0;
    
    Int_t tot_gen_matched_tracks = 0;
    Int_t tot_reco_matched_tracks = 0;

    Int_t tot_gen_matched_tracks_used = 0;
    Int_t tot_reco_matched_tracks_used = 0;

    //looping over events
    std::cout << "Looping over events" << std::endl;
    for (Long64_t ient = beg_event; ient <= end_event; ient++) { // 
        // Print progress

        int mult = end_event/10;
	    double percentage = round(ient * 100 / end_event);
	    // Print progress                                                                                                                                                        
	    if (ient % mult == 0) {                                                                                                                                                  
	    std::cout << "entry nb = " << ient << std::endl;                                                                                       
	    std::cout << percentage << "%" << std::endl;
	    }        
 
        t.GetEntry(ient); 

        //get the nr of tracks
        Int_t n_tracks = t.ntrk;
        Int_t n_tracks_gen = t.nrefTrk;

        //get the weight and pthat
        weight = t.weight;
        pthat = t.pthat;

        
        //see if MC pass the trigger of 40 GeV
        if(!(t.HLT_HIAK4PFJet40_v1 == 1)) continue;//Add prescale weight if necessary
        if(t.HLT_HIAK4PFJet40_v1 == 1 && t.HLT_HIAK4PFJet60_v1 == 0 && t.HLT_HIAK4PFJet80_v1 == 0 && t.HLT_HIAK4PFJet100_v1 == 0) weight*=prescale_pf40;
            


        // Loop over jets for reco
        for (Int_t ijet = 0; ijet < t.nref; ijet++) {
            jtHadFlav = t.jtHadFlav[ijet];
            jtNbHad = t.jtNbHad[ijet];
            
            //Save jet information
            jt_eta_gen = t.refeta[ijet];
            jt_eta_reco = t.jteta[ijet];
            jpt_gen = t.refpt[ijet];
            jpt_reco = t.jtpt[ijet];
            discr = t.discr_particleNet_BvsAll[ijet];


            //Select and build the track vectors

            //create a track vector for the i jet at GEN
            std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors_gen;
            std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors_gen_notmatched;


            //create a track vector for the i jet at RECO
            std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;
            std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors_notmatched;

            //Aggregate and match tracks
            
            // Aggregate 1 b hadron gen, get mB
            if (aggregated) mB_gen = ReconstuctSingleB_gen(trackVectors_gen, t, ijet);
            else{
                create_trackVectors_gen(trackVectors_gen, t, ijet);
                mB_gen = t.refmB[ijet];
            }

            // Aggregate 1 b hadron reco, get mB
            if (aggregated) mB_reco = ReconstuctSingleB_reco(trackVectors, t, ijet);
            else{
                create_trackVectors_reco(trackVectors, t, ijet);
                mB_reco = t.jtmB[ijet];
            } 

            //Get the total nr of tracks before matching
            tot_gen_tracks += trackVectors_gen.size();
            tot_reco_tracks += trackVectors.size();

            if(matching) match_tracks(trackVectors, trackVectors_notmatched, trackVectors_gen, trackVectors_gen_notmatched, t, ijet);
            
            //Save nr of tracks after the matching
            tot_gen_matched_tracks += trackVectors_gen.size();
            tot_reco_matched_tracks += trackVectors.size();

            //Create dr and eec weight

            //______________GEN_________________
            
            //save jet eta and phi
            Float_t jet_eta = t.refeta[ijet];
            Float_t jet_phi = t.refphi[ijet];

            // Set the new nr of tracks and keep record of the number of entries for dr
            n_tracks_gen = trackVectors_gen.size(); 
            Int_t count_dr_gen = 0;

            // Loop over the matched tracks gen
            for (Int_t i = 0; i < n_tracks_gen; i++) { 
                
                Float_t etai = trackVectors_gen[i].Eta();
                Float_t phii = trackVectors_gen[i].Phi();
                Float_t ipt = trackVectors_gen[i].Pt();

                //Select tracks with pT > 1 GeV
                if(ipt < 1) continue;
                //See how many tracks passed this cut
                tot_gen_matched_tracks_used += 1;
                

                // Loop over pairs
                for(Int_t j=0; j < i; j++){
                    
                    Float_t etaj = trackVectors_gen[j].Eta();
                    Float_t phij = trackVectors_gen[j].Phi();
                    Float_t jpt = trackVectors_gen[j].Pt();

                    if(jpt < 1) continue;
                    
                    // Calculate and save dr
                    dr_gen[count_dr_gen] = t.calc_dr(etai, phii, etaj, phij);

                    // Calculate and save the eec weight
                    eec_gen[count_dr_gen] = pow(ipt*jpt, n);

                    //One entry added to dr and eec array
                    count_dr_gen += 1;
                }
            } 

            //Save total matched pairs
            ndr_gen = count_dr_gen;

            //Pair the non-matched tracks withing themselves

            // Set the new nr of tracks and keep record of the number of entries for dr
            Int_t n_tracks_gen_tot = trackVectors_gen_notmatched.size();
            //Start filling the dr and eec array after the matched pairs
            Int_t count_dr_gen_tot = count_dr_gen;

            // Loop over all tracks gen
            for (Int_t i = 0; i < n_tracks_gen_tot; i++) { 
                
                Float_t etai = trackVectors_gen_notmatched[i].Eta();
                Float_t phii = trackVectors_gen_notmatched[i].Phi();
                Float_t ipt = trackVectors_gen_notmatched[i].Pt();


                if(ipt < 1) continue;

                // Loop over pairs
                for(Int_t j=0; j < i; j++){
                    
                    Float_t etaj = trackVectors_gen_notmatched[j].Eta();
                    Float_t phij = trackVectors_gen_notmatched[j].Phi();
                    Float_t jpt = trackVectors_gen_notmatched[j].Pt();

                    if(jpt < 1) continue;                    
                    
                    // Calculate and store the dr
                    dr_gen[count_dr_gen_tot] = t.calc_dr(etai, phii, etaj, phij);

                    // Calculate and store the eec weight
                    eec_gen[count_dr_gen_tot] = pow(ipt*jpt, n);

                    count_dr_gen_tot += 1;
                }
            } 

            // Pair non-matched tracks to the matched tracks
            for (Int_t i = 0; i < n_tracks_gen_tot; i++) { 
                
                Float_t etai = trackVectors_gen_notmatched[i].Eta();
                Float_t phii = trackVectors_gen_notmatched[i].Phi();
                Float_t ipt = trackVectors_gen_notmatched[i].Pt();

                if(ipt < 1) continue;

                // Loop over pairs
                for(Int_t j=0; j < n_tracks_gen; j++){
                    
                    Float_t etaj = trackVectors_gen[j].Eta();
                    Float_t phij = trackVectors_gen[j].Phi();
                    Float_t jpt = trackVectors_gen[j].Pt();

                    if(jpt < 1) continue;
                    
                    // Calculate and store the dr
                    dr_gen[count_dr_gen_tot] = t.calc_dr(etai, phii, etaj, phij);

                    // Calculate and store the eec
                    eec_gen[count_dr_gen_tot] = pow(ipt*jpt, n);

                    count_dr_gen_tot += 1;
                }
            }
            ndr_gen_tot = count_dr_gen_tot;

            
            

            //___________RECO_____________________

            //save jet eta and phi
            jet_eta = t.jteta[ijet];
            jet_phi = t.jtphi[ijet];

            
            // Set the new nr of tracks and keep record of the number of entries for dr
            n_tracks = trackVectors.size(); 
            Int_t count_dr_reco = 0;

            
            // Pair matched tracks within themselves
            // Loop over the tracks reco
            for (Int_t i = 0; i < n_tracks; i++) { 
                
                Float_t etai = trackVectors[i].Eta();
                Float_t phii = trackVectors[i].Phi();
                Float_t ipt = trackVectors[i].Pt();


                if(ipt < 1) continue;
                tot_reco_matched_tracks_used += 1;


                // Loop over pairs
                for(Int_t j=0; j < i; j++){
                    
                    Float_t etaj = trackVectors[j].Eta();
                    Float_t phij = trackVectors[j].Phi();
                    Float_t jpt = trackVectors[j].Pt();

                    if(jpt < 1) continue;

                    // Calculate and store the dr
                    dr_reco[count_dr_reco] = t.calc_dr(etai, phii, etaj, phij);

                    // Calculate and store the eec weight
                    eec_reco[count_dr_reco] = pow(ipt*jpt, n);
                    
                    //add one entry to the eec and dr arrays
                    count_dr_reco += 1;

                }
            }              

            //set the length of the eec and dr arrays
            ndr_reco = count_dr_reco;

            //Do the non-matched pairs

            // Set the new nr of tracks and keep record of the number of entries for dr
            Int_t n_tracks_reco_tot = trackVectors_notmatched.size(); 
            Int_t count_dr_reco_tot = count_dr_reco;

            // Match non-matched tracks within themselves
            for (Int_t i = 0; i < n_tracks_reco_tot; i++) { 
                
                Float_t etai = trackVectors_notmatched[i].Eta();
                Float_t phii = trackVectors_notmatched[i].Phi();
                Float_t ipt = trackVectors_notmatched[i].Pt();

                if(ipt < 1) continue;


                // Loop over pairs
                for(Int_t j=0; j < i; j++){
                    
                    Float_t etaj = trackVectors_notmatched[j].Eta();
                    Float_t phij = trackVectors_notmatched[j].Phi();
                    Float_t jpt = trackVectors_notmatched[j].Pt();

                    if(jpt < 1) continue;
                    
                    
                    // Calculate and store the dr
                    dr_reco[count_dr_reco_tot] = t.calc_dr(etai, phii, etaj, phij);

                    // Calculate and store the eec weight
                    eec_reco[count_dr_reco_tot] = pow(ipt*jpt, n);

                    count_dr_reco_tot += 1;
                }
            } 

            // Pair matched with non-matched tracks
            for (Int_t i = 0; i < n_tracks_reco_tot; i++) { 
                
                Float_t etai = trackVectors_notmatched[i].Eta();
                Float_t phii = trackVectors_notmatched[i].Phi();
                Float_t ipt = trackVectors_notmatched[i].Pt();

                if(ipt < 1) continue;

                // Loop over pairs
                for(Int_t j=0; j < n_tracks; j++){
                    
                    Float_t etaj = trackVectors[j].Eta();
                    Float_t phij = trackVectors[j].Phi();
                    Float_t jpt = trackVectors[j].Pt();

                    if(jpt < 1) continue;
                    
                    // Calculate and store the dr
                    dr_reco[count_dr_reco_tot] = t.calc_dr(etai, phii, etaj, phij);

                    // Calculate and store the eec weight
                    eec_reco[count_dr_reco_tot] = pow(ipt*jpt, n);

                    count_dr_reco_tot += 1;
                }
            } 

            ndr_reco_tot = count_dr_reco_tot;


            //Fill the tree entry
            tree->Fill();
        }
        
    }

tree->Write();

fout->Close();

//Print matching efficiency
std::cout << "Tot gen tracks = " << tot_gen_tracks << std::endl;
std::cout << "Tot reco tracks = " << tot_reco_tracks << std::endl;
std::cout << "Tot gen matched tracks = " << tot_gen_matched_tracks << std::endl;
std::cout << "Tot reco matched tracks = " << tot_reco_matched_tracks << std::endl;
std::cout << "Tot gen matched tracks used for pairs = " << tot_gen_matched_tracks_used << std::endl;
std::cout << "Tot reco matched tracks used for pairs = " << tot_reco_matched_tracks_used << std::endl;
std::cout << "---------------------" << std::endl;


std::cout << std::fixed << std::setprecision(4);
std::cout << "Gen matched tracks = " << 100.0 * tot_gen_matched_tracks / tot_gen_tracks << "%" << std::endl;
std::cout << "Reco matched tracks = " << 100.0 * tot_reco_matched_tracks / tot_reco_tracks << "%" << std::endl;

}

void create_trees_eec(int dataType = 1,                                                                                                                               //-1 for data Low //0 for data High //1 for MC - bjet //2 for MC - dijet //3 for bjet_herwig                                                                          
		      Float_t pT_low = 80,                                                                                                                              
		      Float_t pT_high = 140,                                                                                                                            
		      Int_t n=1,                                                                                                                                        
		      Int_t beg_event = 0,
		      Int_t end_event = 1000,

		      bool btag = false,
		      bool aggregated = false,
		      bool matching = true,
		      const char* output_suffix = ".root"){
   

    TString folder = "$mydata/eec_trees/fran_bins/unmerged_trees/";
    
    TString filename;
    TString label;

  if(dataType == 1){
    filename = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    label = "bjet";
}

  else if(dataType == 2){
    filename = "/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    label = "dijet";

}
  else if(dataType == 3){
    filename = "/data_CMS/cms/kalipoliti/herwigMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    label = "bjet_herwig";
  }

  else{std::cout << "error: data was chosen!!" << std::endl;
    return;
  } 

   label += "_inclusive";//,"moreb", "other","mc"};
    

    do_trees(filename, dataType, label, folder, n, pT_low, pT_high, aggregated, btag, matching, beg_event, end_event, output_suffix);
        }
  
