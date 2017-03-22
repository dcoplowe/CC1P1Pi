//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 20 06:36:15 2017 by ROOT version 5.34/05
// from TTree sel/Tuple created by an AnaTuple managed by AnaTupleManager
// found on file: CC1P1Pi_R13200_190317_5/grid/central_value/minerva/ana/v10r8p9/00/01/32/00/SIM_minerva_00013200_Subruns_0001-0002-0003-0004_CC1P1PiAnalysis_Ana_Tuple_v10r8p9-dcoplowe.root
//////////////////////////////////////////////////////////

#ifndef AnalysisReader_h
#define AnalysisReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "/grid/fermiapp/minerva/software_releases/lcgcmake/build/lcg_61/projects/ROOT-5.34.05/src/ROOT/5.34.05/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalysisReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        eventID;
   Int_t           sel_nuFlavor;
   Int_t           sel_nuHelicity;
   Int_t           sel_intCurrent;
   Int_t           sel_intType;
   Double_t        sel_E;
   Double_t        sel_Q2;
   Double_t        sel_x;
   Double_t        sel_y;
   Double_t        sel_W;
   Double_t        sel_score;
   Double_t        sel_leptonE[4];
   Double_t        sel_vtx[4];
   Int_t           sel_mu_PDG;
   Int_t           sel_mu_isKinked;
   Int_t           sel_mu_michel;
   Int_t           sel_pi_EX_michel;
   Int_t           sel_pi_EX_nNodes;
   Int_t           sel_pi_FSI;
   Int_t           sel_pi_LL_michel;
   Int_t           sel_pi_LL_nNodes;
   Int_t           sel_pi_PDG;
   Int_t           sel_pi_isKinked;
   Int_t           sel_pr_EX_michel;
   Int_t           sel_pr_EX_nNodes;
   Int_t           sel_pr_FSI;
   Int_t           sel_pr_LL_michel;
   Int_t           sel_pr_LL_nNodes;
   Int_t           sel_pr_PDG;
   Int_t           sel_pr_isKinked;
   Double_t        sel_Enu_EX;
   Double_t        sel_Enu_LL;
   Double_t        sel_Q2_EX;
   Double_t        sel_Q2_LL;
   Double_t        sel_dalphaT_EX;
   Double_t        sel_dalphaT_LL;
   Double_t        sel_dpTT_EX;
   Double_t        sel_dpTT_EX_tmumom;
   Double_t        sel_dpTT_EX_tnudir;
   Double_t        sel_dpTT_EX_tpimom;
   Double_t        sel_dpTT_EX_tprmom;
   Double_t        sel_dpTT_LL;
   Double_t        sel_dpTT_LL_tmumom;
   Double_t        sel_dpTT_LL_tnudir;
   Double_t        sel_dpTT_LL_tpimom;
   Double_t        sel_dpTT_LL_tprmom;
   Double_t        sel_dpTT_pi_EX;
   Double_t        sel_dpTT_pi_EX_tmumom;
   Double_t        sel_dpTT_pi_EX_tnudir;
   Double_t        sel_dpTT_pi_EX_tpimom;
   Double_t        sel_dpTT_pi_EX_tprmom;
   Double_t        sel_dpTT_pi_LL;
   Double_t        sel_dpTT_pi_LL_tmumom;
   Double_t        sel_dpTT_pi_LL_tnudir;
   Double_t        sel_dpTT_pi_LL_tpimom;
   Double_t        sel_dpTT_pi_LL_tprmom;
   Double_t        sel_dpTT_pi_dir_EX;
   Double_t        sel_dpTT_pi_dir_EX_tmumom;
   Double_t        sel_dpTT_pi_dir_EX_tnudir;
   Double_t        sel_dpTT_pi_dir_EX_tpidir;
   Double_t        sel_dpTT_pi_dir_EX_tprmom;
   Double_t        sel_dpTT_pi_dir_LL;
   Double_t        sel_dpTT_pi_dir_LL_tmumom;
   Double_t        sel_dpTT_pi_dir_LL_tnudir;
   Double_t        sel_dpTT_pi_dir_LL_tpidir;
   Double_t        sel_dpTT_pi_dir_LL_tprmom;
   Double_t        sel_dpTT_pr_EX;
   Double_t        sel_dpTT_pr_EX_tmumom;
   Double_t        sel_dpTT_pr_EX_tnudir;
   Double_t        sel_dpTT_pr_EX_tpimon;
   Double_t        sel_dpTT_pr_EX_tprmom;
   Double_t        sel_dpTT_pr_LL;
   Double_t        sel_dpTT_pr_LL_tmumom;
   Double_t        sel_dpTT_pr_LL_tnudir;
   Double_t        sel_dpTT_pr_LL_tpimon;
   Double_t        sel_dpTT_pr_LL_tprmom;
   Double_t        sel_dpTT_pr_dir_EX;
   Double_t        sel_dpTT_pr_dir_EX_tmumom;
   Double_t        sel_dpTT_pr_dir_EX_tnudir;
   Double_t        sel_dpTT_pr_dir_EX_tpimom;
   Double_t        sel_dpTT_pr_dir_EX_tprdir;
   Double_t        sel_dpTT_pr_dir_LL;
   Double_t        sel_dpTT_pr_dir_LL_tmumom;
   Double_t        sel_dpTT_pr_dir_LL_tnudir;
   Double_t        sel_dpTT_pr_dir_LL_tpimom;
   Double_t        sel_dpTT_pr_dir_LL_tprdir;
   Double_t        sel_dpT_EX;
   Double_t        sel_dpT_LL;
   Double_t        sel_dphiT_EX;
   Double_t        sel_dphiT_LL;
   Double_t        sel_mu_E;
   Double_t        sel_mu_KE;
   Double_t        sel_mu_Phi;
   Double_t        sel_mu_Theta;
   Double_t        sel_mu_XTheta;
   Double_t        sel_mu_YTheta;
   Double_t        sel_mu_chi2ndf;
   Double_t        sel_mu_det_frac;
   Double_t        sel_mu_det_otherE;
   Double_t        sel_mu_mom;
   Double_t        sel_mu_pTMag;
   Double_t        sel_mu_pTMag_tmumom;
   Double_t        sel_mu_pTMag_tnudir;
   Double_t        sel_mu_pTT;
   Double_t        sel_mu_pTT_tmumom;
   Double_t        sel_mu_pTT_tnudir;
   Double_t        sel_mu_score;
   Double_t        sel_mu_trueE;
   Double_t        sel_mu_trueKE;
   Double_t        sel_mu_truePhi;
   Double_t        sel_mu_trueTheta;
   Double_t        sel_mu_truemom;
   Double_t        sel_mu_truepTMag;
   Double_t        sel_mu_truepTT;
   Double_t        sel_pi_CaloE;
   Double_t        sel_pi_EX_E;
   Double_t        sel_pi_EX_E_altH;
   Double_t        sel_pi_EX_KE;
   Double_t        sel_pi_EX_Phi;
   Double_t        sel_pi_EX_Theta;
   Double_t        sel_pi_EX_XTheta;
   Double_t        sel_pi_EX_YTheta;
   Double_t        sel_pi_EX_ln_Q0;
   Double_t        sel_pi_EX_ln_Q1;
   Double_t        sel_pi_EX_ln_Q2;
   Double_t        sel_pi_EX_ln_Q3;
   Double_t        sel_pi_EX_ln_Q4;
   Double_t        sel_pi_EX_ln_Q5;
   Double_t        sel_pi_EX_mom;
   Double_t        sel_pi_EX_mom_altH;
   Double_t        sel_pi_EX_pTMag;
   Double_t        sel_pi_EX_pTMag_tnudir;
   Double_t        sel_pi_EX_pTMag_tpimom;
   Double_t        sel_pi_EX_pTT;
   Double_t        sel_pi_EX_pTT_tnudir;
   Double_t        sel_pi_EX_pTT_tpimom;
   Double_t        sel_pi_EX_score;
   Double_t        sel_pi_EX_score_altH;
   Double_t        sel_pi_LL_E;
   Double_t        sel_pi_LL_E_altH;
   Double_t        sel_pi_LL_KE;
   Double_t        sel_pi_LL_Phi;
   Double_t        sel_pi_LL_Theta;
   Double_t        sel_pi_LL_XTheta;
   Double_t        sel_pi_LL_YTheta;
   Double_t        sel_pi_LL_ln_Q0;
   Double_t        sel_pi_LL_ln_Q1;
   Double_t        sel_pi_LL_ln_Q2;
   Double_t        sel_pi_LL_ln_Q3;
   Double_t        sel_pi_LL_ln_Q4;
   Double_t        sel_pi_LL_ln_Q5;
   Double_t        sel_pi_LL_mom;
   Double_t        sel_pi_LL_mom_altH;
   Double_t        sel_pi_LL_pTMag;
   Double_t        sel_pi_LL_pTMag_tnudir;
   Double_t        sel_pi_LL_pTMag_tpimom;
   Double_t        sel_pi_LL_pTT;
   Double_t        sel_pi_LL_pTT_tnudir;
   Double_t        sel_pi_LL_pTT_tpimom;
   Double_t        sel_pi_LL_score;
   Double_t        sel_pi_LL_score_altH;
   Double_t        sel_pi_chi2ndf;
   Double_t        sel_pi_det_frac;
   Double_t        sel_pi_det_otherE;
   Double_t        sel_pi_trueE;
   Double_t        sel_pi_trueKE;
   Double_t        sel_pi_truePhi;
   Double_t        sel_pi_trueTheta;
   Double_t        sel_pi_truemom;
   Double_t        sel_pi_truepTMag;
   Double_t        sel_pi_truepTT;
   Double_t        sel_pr_CaloE;
   Double_t        sel_pr_EX_E;
   Double_t        sel_pr_EX_E_altH;
   Double_t        sel_pr_EX_KE;
   Double_t        sel_pr_EX_Phi;
   Double_t        sel_pr_EX_Theta;
   Double_t        sel_pr_EX_XTheta;
   Double_t        sel_pr_EX_YTheta;
   Double_t        sel_pr_EX_ln_Q0;
   Double_t        sel_pr_EX_ln_Q1;
   Double_t        sel_pr_EX_ln_Q2;
   Double_t        sel_pr_EX_ln_Q3;
   Double_t        sel_pr_EX_ln_Q4;
   Double_t        sel_pr_EX_ln_Q5;
   Double_t        sel_pr_EX_mom;
   Double_t        sel_pr_EX_mom_altH;
   Double_t        sel_pr_EX_pTMag;
   Double_t        sel_pr_EX_pTMag_tnudir;
   Double_t        sel_pr_EX_pTMag_tprmom;
   Double_t        sel_pr_EX_pTT;
   Double_t        sel_pr_EX_pTT_tnudir;
   Double_t        sel_pr_EX_pTT_tprmom;
   Double_t        sel_pr_EX_score;
   Double_t        sel_pr_EX_score_altH;
   Double_t        sel_pr_LL_E;
   Double_t        sel_pr_LL_E_altH;
   Double_t        sel_pr_LL_KE;
   Double_t        sel_pr_LL_Phi;
   Double_t        sel_pr_LL_Theta;
   Double_t        sel_pr_LL_XTheta;
   Double_t        sel_pr_LL_YTheta;
   Double_t        sel_pr_LL_ln_Q0;
   Double_t        sel_pr_LL_ln_Q1;
   Double_t        sel_pr_LL_ln_Q2;
   Double_t        sel_pr_LL_ln_Q3;
   Double_t        sel_pr_LL_ln_Q4;
   Double_t        sel_pr_LL_ln_Q5;
   Double_t        sel_pr_LL_mom;
   Double_t        sel_pr_LL_mom_altH;
   Double_t        sel_pr_LL_pTMag;
   Double_t        sel_pr_LL_pTMag_tnudir;
   Double_t        sel_pr_LL_pTMag_tprmom;
   Double_t        sel_pr_LL_pTT;
   Double_t        sel_pr_LL_pTT_tnudir;
   Double_t        sel_pr_LL_pTT_tprmom;
   Double_t        sel_pr_LL_score;
   Double_t        sel_pr_LL_score_altH;
   Double_t        sel_pr_chi2ndf;
   Double_t        sel_pr_det_frac;
   Double_t        sel_pr_det_otherE;
   Double_t        sel_pr_trueE;
   Double_t        sel_pr_trueKE;
   Double_t        sel_pr_truePhi;
   Double_t        sel_pr_trueTheta;
   Double_t        sel_pr_truemom;
   Double_t        sel_pr_truepTMag;
   Double_t        sel_pr_truepTT;
   Double_t        sel_trueEnu;
   Double_t        sel_trueQ2;
   Double_t        sel_truedalphaT;
   Double_t        sel_truedpT;
   Double_t        sel_truedpTT;
   Double_t        sel_truedpTT_pi;
   Double_t        sel_truedpTT_pi_dir;
   Double_t        sel_truedpTT_pr;
   Double_t        sel_truedpTT_pr_dir;
   Double_t        sel_truedphiT;
   Int_t           n_iso_blobs;
   Int_t           sel_iso_blob_nclusters[1];   //[n_iso_blobs]
   Double_t        sel_VTX[3];
   Double_t        sel_dpT_vec_EX[3];
   Double_t        sel_dpT_vec_LL[3];
   Double_t        sel_iso_blob_energy[1];   //[n_iso_blobs]
   Double_t        sel_meanPDP[3];
   Double_t        sel_mu_4mom[4];
   Double_t        sel_mu_endpos[3];
   Double_t        sel_mu_pT[3];
   Double_t        sel_mu_pT_tmumom[3];
   Double_t        sel_mu_pT_tnudir[3];
   Double_t        sel_mu_startdir[3];
   Double_t        sel_mu_startpos[3];
   Double_t        sel_mu_true4mom[4];
   Double_t        sel_mu_trueendpos[3];
   Double_t        sel_mu_truepT[3];
   Double_t        sel_mu_truestartdir[3];
   Double_t        sel_mu_truestartpos[3];
   Double_t        sel_nu_dir_001[3];
   Double_t        sel_nu_dir_PDP[3];
   Double_t        sel_pi_EX_4mom[4];
   Double_t        sel_pi_EX_pT[3];
   Double_t        sel_pi_EX_pT_tnudir[3];
   Double_t        sel_pi_EX_pT_tpimom[3];
   Double_t        sel_pi_LL_4mom[4];
   Double_t        sel_pi_LL_pT[3];
   Double_t        sel_pi_LL_pT_tnudir[3];
   Double_t        sel_pi_LL_pT_tpimom[3];
   Double_t        sel_pi_endpos[3];
   Double_t        sel_pi_startdir[3];
   Double_t        sel_pi_startpos[3];
   Double_t        sel_pi_true4mom[4];
   Double_t        sel_pi_trueendpos[3];
   Double_t        sel_pi_truepT[3];
   Double_t        sel_pi_truestartdir[3];
   Double_t        sel_pi_truestartpos[3];
   Double_t        sel_pr_EX_4mom[4];
   Double_t        sel_pr_EX_pT[3];
   Double_t        sel_pr_EX_pT_tnudir[3];
   Double_t        sel_pr_EX_pT_tprmom[3];
   Double_t        sel_pr_LL_4mom[4];
   Double_t        sel_pr_LL_pT[3];
   Double_t        sel_pr_LL_pT_tnudir[3];
   Double_t        sel_pr_LL_pT_tprmom[3];
   Double_t        sel_pr_endpos[3];
   Double_t        sel_pr_startdir[3];
   Double_t        sel_pr_startpos[3];
   Double_t        sel_pr_true4mom[4];
   Double_t        sel_pr_trueendpos[3];
   Double_t        sel_pr_truepT[3];
   Double_t        sel_pr_truestartdir[3];
   Double_t        sel_pr_truestartpos[3];
   Double_t        sel_truePDP[3];
   Double_t        sel_trueVTX[3];
   Double_t        sel_true_nu_dir_PDP[3];
   Double_t        sel_truedpT_vec[3];
   Int_t           physEvtNum;
   Int_t           n_hyps;
   Int_t           processType;
   Int_t           primaryPart;
   Int_t           n_slices;
   Int_t           slice_numbers[1];   //[n_slices]
   Int_t           shared_slice;
   Double_t        vtx[4];
   Double_t        vtxErr[4];
   Double_t        E[4];
   Bool_t          found_truth;
   Bool_t          isMinosMatchTrack;
   Bool_t          isMinosMatchStub;
   Int_t           contained_evt;
   Int_t           muon_charge;
   Int_t           n_anchored_long_trk_prongs;
   Int_t           n_anchored_short_trk_prongs;
   Int_t           n_iso_trk_prongs;
   Int_t           n_prongs;
   Int_t           n_tracks3;
   Int_t           ncuts;
   Int_t           new_tracks;
   Int_t           nsplits;
   Int_t           target_region;
   Int_t           true_target_region;
   Int_t           vert_exists;
   Int_t           accum_level[2];
   Bool_t          truth_reco_isMinosMatch;
   Int_t           truth_muon_charge;
   Int_t           truth_n_ele;
   Int_t           truth_n_kPM;
   Int_t           truth_n_kaO;
   Int_t           truth_n_muo;
   Int_t           truth_n_ntn;
   Int_t           truth_n_pho;
   Int_t           truth_n_pi0;
   Int_t           truth_n_piM;
   Int_t           truth_n_piP;
   Int_t           truth_n_pro;
   Int_t           truth_n_tau;
   Int_t           truth_ncuts;
   Int_t           truth_nsplits;
   Int_t           truth_pi_EX_michel;
   Int_t           truth_pi_LL_michel;
   Int_t           truth_pr_EX_michel;
   Int_t           truth_pr_LL_michel;
   Int_t           truth_reco_target;
   Int_t           truth_should_be_accepted;
   Int_t           truth_true_target_region;
   Double_t        truth_mu_E;
   Double_t        truth_mu_KE;
   Double_t        truth_mu_Phi;
   Double_t        truth_mu_Theta;
   Double_t        truth_mu_mom;
   Double_t        truth_mu_pTMag;
   Double_t        truth_mu_pTT;
   Double_t        truth_pi_E;
   Double_t        truth_pi_EX_score;
   Double_t        truth_pi_EX_score_altH;
   Double_t        truth_pi_KE;
   Double_t        truth_pi_LL_score;
   Double_t        truth_pi_LL_score_altH;
   Double_t        truth_pi_Phi;
   Double_t        truth_pi_Theta;
   Double_t        truth_pi_mom;
   Double_t        truth_pi_pTMag;
   Double_t        truth_pi_pTT;
   Double_t        truth_pr_E;
   Double_t        truth_pr_EX_score;
   Double_t        truth_pr_EX_score_altH;
   Double_t        truth_pr_KE;
   Double_t        truth_pr_LL_score;
   Double_t        truth_pr_LL_score_altH;
   Double_t        truth_pr_Phi;
   Double_t        truth_pr_Theta;
   Double_t        truth_pr_mom;
   Double_t        truth_pr_pTMag;
   Double_t        truth_pr_pTT;
   Double_t        truth_trueEnu;
   Double_t        truth_trueQ2;
   Double_t        truth_truedalphaT;
   Double_t        truth_truedpT;
   Double_t        truth_truedpTT;
   Double_t        truth_truedpTT_pi;
   Double_t        truth_truedpTT_pi_dir;
   Double_t        truth_truedpTT_pr;
   Double_t        truth_truedpTT_pr_dir;
   Double_t        truth_truedphiT;
   Int_t           truth_accum_level[2];
   Double_t        truth_mu_4mom[4];
   Double_t        truth_mu_pT[3];
   Double_t        truth_pi_4mom[4];
   Double_t        truth_pi_pT[3];
   Double_t        truth_pr_4mom[4];
   Double_t        truth_pr_pT[3];
   Double_t        truth_truedpT_vec[3];
   Int_t           ev_run;
   Int_t           ev_subrun;
   Int_t           ev_detector;
   Int_t           ev_triggerType;
   Int_t           ev_gate;
   Int_t           ev_global_gate;
   Int_t           ev_gps_time_sec;
   Int_t           ev_gps_time_usec;
   Int_t           mc_run;
   Int_t           mc_subrun;
   Int_t           mc_nInteractions;
   Int_t           mc_MIState;
   Double_t        mc_pot;
   Int_t           mc_beamConfig;
   Int_t           mc_processType;
   Int_t           mc_nthEvtInSpill;
   Int_t           mc_nthEvtInFile;
   Int_t           mc_intType;
   Int_t           mc_current;
   Int_t           mc_charm;
   Double_t        mc_weight;
   Double_t        mc_XSec;
   Double_t        mc_diffXSec;
   Int_t           mc_incoming;
   Double_t        mc_fluxDriverProb;
   Int_t           mc_targetNucleus;
   Int_t           mc_targetZ;
   Int_t           mc_targetA;
   Int_t           mc_targetNucleon;
   Int_t           mc_struckQuark;
   Int_t           mc_seaQuark;
   Int_t           mc_resID;
   Int_t           mc_primaryLepton;
   Double_t        mc_incomingE;
   Double_t        mc_Bjorkenx;
   Double_t        mc_Bjorkeny;
   Double_t        mc_Q2;
   Double_t        mc_nuT;
   Double_t        mc_w;
   Double_t        mc_vtx[4];
   Double_t        mc_incomingPartVec[4];
   Double_t        mc_initNucVec[4];
   Double_t        mc_primFSLepton[4];
   Int_t           mc_nFSPart;
   Double_t        mc_FSPartPx[159];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[159];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[159];   //[mc_nFSPart]
   Double_t        mc_FSPartE[159];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[159];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[187];   //[mc_er_nPart]
   Int_t           mc_er_status[187];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[187];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[187];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[187];   //[mc_er_nPart]
   Double_t        mc_er_Px[187];   //[mc_er_nPart]
   Double_t        mc_er_Py[187];   //[mc_er_nPart]
   Double_t        mc_er_Pz[187];   //[mc_er_nPart]
   Double_t        mc_er_E[187];   //[mc_er_nPart]
   Int_t           mc_er_FD[187];   //[mc_er_nPart]
   Int_t           mc_er_LD[187];   //[mc_er_nPart]
   Int_t           mc_er_mother[187];   //[mc_er_nPart]
   Int_t           mc_fr_nNuAncestorIDs;
   Int_t           mc_fr_nuAncestorIDs[8];   //[mc_fr_nNuAncestorIDs]
   Int_t           mc_fr_nuParentID;
   Int_t           mc_fr_decMode;
   Double_t        mc_fr_primProtonVtx[3];
   Double_t        mc_fr_primProtonP[4];
   Double_t        mc_fr_nuParentDecVtx[3];
   Double_t        mc_fr_nuParentProdVtx[3];
   Double_t        mc_fr_nuParentProdP[4];
   Double_t        mc_cvweight_total;
   Double_t        wgt;
   Double_t        mc_cvweight_totalFlux;
   Double_t        mc_cvweight_totalXsec;
   Double_t        mc_ppfx1_cvweight;
   Double_t        mc_hornCurrent_cvweight;
   Double_t        mc_gen1_cvweight_total;
   Double_t        gen1_wgt;
   Double_t        mc_gen1_cvweight_totalFlux;
   Double_t        mc_gen1_cvweight_NA49;
   Int_t           mc_wgt_Flux_BeamFocus_sz;
   Double_t        mc_wgt_Flux_BeamFocus[1];   //[mc_wgt_Flux_BeamFocus_sz]
   Int_t           mc_wgt_gen1_Flux_Tertiary_sz;
   Double_t        mc_wgt_gen1_Flux_Tertiary[1];   //[mc_wgt_gen1_Flux_Tertiary_sz]
   Int_t           mc_wgt_gen1_Flux_NA49_sz;
   Double_t        mc_wgt_gen1_Flux_NA49[1];   //[mc_wgt_gen1_Flux_NA49_sz]
   Int_t           mc_wgt_Norm_sz;
   Double_t        mc_wgt_Norm[1];   //[mc_wgt_Norm_sz]
   Int_t           mc_wgt_ppfx1_Total_sz;
   Double_t        mc_wgt_ppfx1_Total[1];   //[mc_wgt_ppfx1_Total_sz]
   Int_t           prong_nParticles[5];   //[n_prongs]
   Double_t        prong_part_score[5];   //[n_prongs]
   Double_t        prong_part_mass[5];   //[n_prongs]
   Int_t           prong_part_charge[5];   //[n_prongs]
   Int_t           prong_part_pid[5];   //[n_prongs]
   // vector<vector<double> > *prong_part_E;
   // vector<vector<double> > *prong_part_pos;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_sel_nuFlavor;   //!
   TBranch        *b_sel_nuHelicity;   //!
   TBranch        *b_sel_intCurrent;   //!
   TBranch        *b_sel_intType;   //!
   TBranch        *b_sel_E;   //!
   TBranch        *b_sel_Q2;   //!
   TBranch        *b_sel_x;   //!
   TBranch        *b_sel_y;   //!
   TBranch        *b_sel_W;   //!
   TBranch        *b_sel_score;   //!
   TBranch        *b_sel_leptonE;   //!
   TBranch        *b_sel_vtx;   //!
   TBranch        *b_sel_mu_PDG;   //!
   TBranch        *b_sel_mu_isKinked;   //!
   TBranch        *b_sel_mu_michel;   //!
   TBranch        *b_sel_pi_EX_michel;   //!
   TBranch        *b_sel_pi_EX_nNodes;   //!
   TBranch        *b_sel_pi_FSI;   //!
   TBranch        *b_sel_pi_LL_michel;   //!
   TBranch        *b_sel_pi_LL_nNodes;   //!
   TBranch        *b_sel_pi_PDG;   //!
   TBranch        *b_sel_pi_isKinked;   //!
   TBranch        *b_sel_pr_EX_michel;   //!
   TBranch        *b_sel_pr_EX_nNodes;   //!
   TBranch        *b_sel_pr_FSI;   //!
   TBranch        *b_sel_pr_LL_michel;   //!
   TBranch        *b_sel_pr_LL_nNodes;   //!
   TBranch        *b_sel_pr_PDG;   //!
   TBranch        *b_sel_pr_isKinked;   //!
   TBranch        *b_sel_Enu_EX;   //!
   TBranch        *b_sel_Enu_LL;   //!
   TBranch        *b_sel_Q2_EX;   //!
   TBranch        *b_sel_Q2_LL;   //!
   TBranch        *b_sel_dalphaT_EX;   //!
   TBranch        *b_sel_dalphaT_LL;   //!
   TBranch        *b_sel_dpTT_EX;   //!
   TBranch        *b_sel_dpTT_EX_tmumom;   //!
   TBranch        *b_sel_dpTT_EX_tnudir;   //!
   TBranch        *b_sel_dpTT_EX_tpimom;   //!
   TBranch        *b_sel_dpTT_EX_tprmom;   //!
   TBranch        *b_sel_dpTT_LL;   //!
   TBranch        *b_sel_dpTT_LL_tmumom;   //!
   TBranch        *b_sel_dpTT_LL_tnudir;   //!
   TBranch        *b_sel_dpTT_LL_tpimom;   //!
   TBranch        *b_sel_dpTT_LL_tprmom;   //!
   TBranch        *b_sel_dpTT_pi_EX;   //!
   TBranch        *b_sel_dpTT_pi_EX_tmumom;   //!
   TBranch        *b_sel_dpTT_pi_EX_tnudir;   //!
   TBranch        *b_sel_dpTT_pi_EX_tpimom;   //!
   TBranch        *b_sel_dpTT_pi_EX_tprmom;   //!
   TBranch        *b_sel_dpTT_pi_LL;   //!
   TBranch        *b_sel_dpTT_pi_LL_tmumom;   //!
   TBranch        *b_sel_dpTT_pi_LL_tnudir;   //!
   TBranch        *b_sel_dpTT_pi_LL_tpimom;   //!
   TBranch        *b_sel_dpTT_pi_LL_tprmom;   //!
   TBranch        *b_sel_dpTT_pi_dir_EX;   //!
   TBranch        *b_sel_dpTT_pi_dir_EX_tmumom;   //!
   TBranch        *b_sel_dpTT_pi_dir_EX_tnudir;   //!
   TBranch        *b_sel_dpTT_pi_dir_EX_tpidir;   //!
   TBranch        *b_sel_dpTT_pi_dir_EX_tprmom;   //!
   TBranch        *b_sel_dpTT_pi_dir_LL;   //!
   TBranch        *b_sel_dpTT_pi_dir_LL_tmumom;   //!
   TBranch        *b_sel_dpTT_pi_dir_LL_tnudir;   //!
   TBranch        *b_sel_dpTT_pi_dir_LL_tpidir;   //!
   TBranch        *b_sel_dpTT_pi_dir_LL_tprmom;   //!
   TBranch        *b_sel_dpTT_pr_EX;   //!
   TBranch        *b_sel_dpTT_pr_EX_tmumom;   //!
   TBranch        *b_sel_dpTT_pr_EX_tnudir;   //!
   TBranch        *b_sel_dpTT_pr_EX_tpimon;   //!
   TBranch        *b_sel_dpTT_pr_EX_tprmom;   //!
   TBranch        *b_sel_dpTT_pr_LL;   //!
   TBranch        *b_sel_dpTT_pr_LL_tmumom;   //!
   TBranch        *b_sel_dpTT_pr_LL_tnudir;   //!
   TBranch        *b_sel_dpTT_pr_LL_tpimon;   //!
   TBranch        *b_sel_dpTT_pr_LL_tprmom;   //!
   TBranch        *b_sel_dpTT_pr_dir_EX;   //!
   TBranch        *b_sel_dpTT_pr_dir_EX_tmumom;   //!
   TBranch        *b_sel_dpTT_pr_dir_EX_tnudir;   //!
   TBranch        *b_sel_dpTT_pr_dir_EX_tpimom;   //!
   TBranch        *b_sel_dpTT_pr_dir_EX_tprdir;   //!
   TBranch        *b_sel_dpTT_pr_dir_LL;   //!
   TBranch        *b_sel_dpTT_pr_dir_LL_tmumom;   //!
   TBranch        *b_sel_dpTT_pr_dir_LL_tnudir;   //!
   TBranch        *b_sel_dpTT_pr_dir_LL_tpimom;   //!
   TBranch        *b_sel_dpTT_pr_dir_LL_tprdir;   //!
   TBranch        *b_sel_dpT_EX;   //!
   TBranch        *b_sel_dpT_LL;   //!
   TBranch        *b_sel_dphiT_EX;   //!
   TBranch        *b_sel_dphiT_LL;   //!
   TBranch        *b_sel_mu_E;   //!
   TBranch        *b_sel_mu_KE;   //!
   TBranch        *b_sel_mu_Phi;   //!
   TBranch        *b_sel_mu_Theta;   //!
   TBranch        *b_sel_mu_XTheta;   //!
   TBranch        *b_sel_mu_YTheta;   //!
   TBranch        *b_sel_mu_chi2ndf;   //!
   TBranch        *b_sel_mu_det_frac;   //!
   TBranch        *b_sel_mu_det_otherE;   //!
   TBranch        *b_sel_mu_mom;   //!
   TBranch        *b_sel_mu_pTMag;   //!
   TBranch        *b_sel_mu_pTMag_tmumom;   //!
   TBranch        *b_sel_mu_pTMag_tnudir;   //!
   TBranch        *b_sel_mu_pTT;   //!
   TBranch        *b_sel_mu_pTT_tmumom;   //!
   TBranch        *b_sel_mu_pTT_tnudir;   //!
   TBranch        *b_sel_mu_score;   //!
   TBranch        *b_sel_mu_trueE;   //!
   TBranch        *b_sel_mu_trueKE;   //!
   TBranch        *b_sel_mu_truePhi;   //!
   TBranch        *b_sel_mu_trueTheta;   //!
   TBranch        *b_sel_mu_truemom;   //!
   TBranch        *b_sel_mu_truepTMag;   //!
   TBranch        *b_sel_mu_truepTT;   //!
   TBranch        *b_sel_pi_CaloE;   //!
   TBranch        *b_sel_pi_EX_E;   //!
   TBranch        *b_sel_pi_EX_E_altH;   //!
   TBranch        *b_sel_pi_EX_KE;   //!
   TBranch        *b_sel_pi_EX_Phi;   //!
   TBranch        *b_sel_pi_EX_Theta;   //!
   TBranch        *b_sel_pi_EX_XTheta;   //!
   TBranch        *b_sel_pi_EX_YTheta;   //!
   TBranch        *b_sel_pi_EX_ln_Q0;   //!
   TBranch        *b_sel_pi_EX_ln_Q1;   //!
   TBranch        *b_sel_pi_EX_ln_Q2;   //!
   TBranch        *b_sel_pi_EX_ln_Q3;   //!
   TBranch        *b_sel_pi_EX_ln_Q4;   //!
   TBranch        *b_sel_pi_EX_ln_Q5;   //!
   TBranch        *b_sel_pi_EX_mom;   //!
   TBranch        *b_sel_pi_EX_mom_altH;   //!
   TBranch        *b_sel_pi_EX_pTMag;   //!
   TBranch        *b_sel_pi_EX_pTMag_tnudir;   //!
   TBranch        *b_sel_pi_EX_pTMag_tpimom;   //!
   TBranch        *b_sel_pi_EX_pTT;   //!
   TBranch        *b_sel_pi_EX_pTT_tnudir;   //!
   TBranch        *b_sel_pi_EX_pTT_tpimom;   //!
   TBranch        *b_sel_pi_EX_score;   //!
   TBranch        *b_sel_pi_EX_score_altH;   //!
   TBranch        *b_sel_pi_LL_E;   //!
   TBranch        *b_sel_pi_LL_E_altH;   //!
   TBranch        *b_sel_pi_LL_KE;   //!
   TBranch        *b_sel_pi_LL_Phi;   //!
   TBranch        *b_sel_pi_LL_Theta;   //!
   TBranch        *b_sel_pi_LL_XTheta;   //!
   TBranch        *b_sel_pi_LL_YTheta;   //!
   TBranch        *b_sel_pi_LL_ln_Q0;   //!
   TBranch        *b_sel_pi_LL_ln_Q1;   //!
   TBranch        *b_sel_pi_LL_ln_Q2;   //!
   TBranch        *b_sel_pi_LL_ln_Q3;   //!
   TBranch        *b_sel_pi_LL_ln_Q4;   //!
   TBranch        *b_sel_pi_LL_ln_Q5;   //!
   TBranch        *b_sel_pi_LL_mom;   //!
   TBranch        *b_sel_pi_LL_mom_altH;   //!
   TBranch        *b_sel_pi_LL_pTMag;   //!
   TBranch        *b_sel_pi_LL_pTMag_tnudir;   //!
   TBranch        *b_sel_pi_LL_pTMag_tpimom;   //!
   TBranch        *b_sel_pi_LL_pTT;   //!
   TBranch        *b_sel_pi_LL_pTT_tnudir;   //!
   TBranch        *b_sel_pi_LL_pTT_tpimom;   //!
   TBranch        *b_sel_pi_LL_score;   //!
   TBranch        *b_sel_pi_LL_score_altH;   //!
   TBranch        *b_sel_pi_chi2ndf;   //!
   TBranch        *b_sel_pi_det_frac;   //!
   TBranch        *b_sel_pi_det_otherE;   //!
   TBranch        *b_sel_pi_trueE;   //!
   TBranch        *b_sel_pi_trueKE;   //!
   TBranch        *b_sel_pi_truePhi;   //!
   TBranch        *b_sel_pi_trueTheta;   //!
   TBranch        *b_sel_pi_truemom;   //!
   TBranch        *b_sel_pi_truepTMag;   //!
   TBranch        *b_sel_pi_truepTT;   //!
   TBranch        *b_sel_pr_CaloE;   //!
   TBranch        *b_sel_pr_EX_E;   //!
   TBranch        *b_sel_pr_EX_E_altH;   //!
   TBranch        *b_sel_pr_EX_KE;   //!
   TBranch        *b_sel_pr_EX_Phi;   //!
   TBranch        *b_sel_pr_EX_Theta;   //!
   TBranch        *b_sel_pr_EX_XTheta;   //!
   TBranch        *b_sel_pr_EX_YTheta;   //!
   TBranch        *b_sel_pr_EX_ln_Q0;   //!
   TBranch        *b_sel_pr_EX_ln_Q1;   //!
   TBranch        *b_sel_pr_EX_ln_Q2;   //!
   TBranch        *b_sel_pr_EX_ln_Q3;   //!
   TBranch        *b_sel_pr_EX_ln_Q4;   //!
   TBranch        *b_sel_pr_EX_ln_Q5;   //!
   TBranch        *b_sel_pr_EX_mom;   //!
   TBranch        *b_sel_pr_EX_mom_altH;   //!
   TBranch        *b_sel_pr_EX_pTMag;   //!
   TBranch        *b_sel_pr_EX_pTMag_tnudir;   //!
   TBranch        *b_sel_pr_EX_pTMag_tprmom;   //!
   TBranch        *b_sel_pr_EX_pTT;   //!
   TBranch        *b_sel_pr_EX_pTT_tnudir;   //!
   TBranch        *b_sel_pr_EX_pTT_tprmom;   //!
   TBranch        *b_sel_pr_EX_score;   //!
   TBranch        *b_sel_pr_EX_score_altH;   //!
   TBranch        *b_sel_pr_LL_E;   //!
   TBranch        *b_sel_pr_LL_E_altH;   //!
   TBranch        *b_sel_pr_LL_KE;   //!
   TBranch        *b_sel_pr_LL_Phi;   //!
   TBranch        *b_sel_pr_LL_Theta;   //!
   TBranch        *b_sel_pr_LL_XTheta;   //!
   TBranch        *b_sel_pr_LL_YTheta;   //!
   TBranch        *b_sel_pr_LL_ln_Q0;   //!
   TBranch        *b_sel_pr_LL_ln_Q1;   //!
   TBranch        *b_sel_pr_LL_ln_Q2;   //!
   TBranch        *b_sel_pr_LL_ln_Q3;   //!
   TBranch        *b_sel_pr_LL_ln_Q4;   //!
   TBranch        *b_sel_pr_LL_ln_Q5;   //!
   TBranch        *b_sel_pr_LL_mom;   //!
   TBranch        *b_sel_pr_LL_mom_altH;   //!
   TBranch        *b_sel_pr_LL_pTMag;   //!
   TBranch        *b_sel_pr_LL_pTMag_tnudir;   //!
   TBranch        *b_sel_pr_LL_pTMag_tprmom;   //!
   TBranch        *b_sel_pr_LL_pTT;   //!
   TBranch        *b_sel_pr_LL_pTT_tnudir;   //!
   TBranch        *b_sel_pr_LL_pTT_tprmom;   //!
   TBranch        *b_sel_pr_LL_score;   //!
   TBranch        *b_sel_pr_LL_score_altH;   //!
   TBranch        *b_sel_pr_chi2ndf;   //!
   TBranch        *b_sel_pr_det_frac;   //!
   TBranch        *b_sel_pr_det_otherE;   //!
   TBranch        *b_sel_pr_trueE;   //!
   TBranch        *b_sel_pr_trueKE;   //!
   TBranch        *b_sel_pr_truePhi;   //!
   TBranch        *b_sel_pr_trueTheta;   //!
   TBranch        *b_sel_pr_truemom;   //!
   TBranch        *b_sel_pr_truepTMag;   //!
   TBranch        *b_sel_pr_truepTT;   //!
   TBranch        *b_sel_trueEnu;   //!
   TBranch        *b_sel_trueQ2;   //!
   TBranch        *b_sel_truedalphaT;   //!
   TBranch        *b_sel_truedpT;   //!
   TBranch        *b_sel_truedpTT;   //!
   TBranch        *b_sel_truedpTT_pi;   //!
   TBranch        *b_sel_truedpTT_pi_dir;   //!
   TBranch        *b_sel_truedpTT_pr;   //!
   TBranch        *b_sel_truedpTT_pr_dir;   //!
   TBranch        *b_sel_truedphiT;   //!
   TBranch        *b_n_iso_blobs;   //!
   TBranch        *b_sel_iso_blob_nclusters;   //!
   TBranch        *b_sel_VTX;   //!
   TBranch        *b_sel_dpT_vec_EX;   //!
   TBranch        *b_sel_dpT_vec_LL;   //!
   TBranch        *b_sel_iso_blob_energy;   //!
   TBranch        *b_sel_meanPDP;   //!
   TBranch        *b_sel_mu_4mom;   //!
   TBranch        *b_sel_mu_endpos;   //!
   TBranch        *b_sel_mu_pT;   //!
   TBranch        *b_sel_mu_pT_tmumom;   //!
   TBranch        *b_sel_mu_pT_tnudir;   //!
   TBranch        *b_sel_mu_startdir;   //!
   TBranch        *b_sel_mu_startpos;   //!
   TBranch        *b_sel_mu_true4mom;   //!
   TBranch        *b_sel_mu_trueendpos;   //!
   TBranch        *b_sel_mu_truepT;   //!
   TBranch        *b_sel_mu_truestartdir;   //!
   TBranch        *b_sel_mu_truestartpos;   //!
   TBranch        *b_sel_nu_dir_001;   //!
   TBranch        *b_sel_nu_dir_PDP;   //!
   TBranch        *b_sel_pi_EX_4mom;   //!
   TBranch        *b_sel_pi_EX_pT;   //!
   TBranch        *b_sel_pi_EX_pT_tnudir;   //!
   TBranch        *b_sel_pi_EX_pT_tpimom;   //!
   TBranch        *b_sel_pi_LL_4mom;   //!
   TBranch        *b_sel_pi_LL_pT;   //!
   TBranch        *b_sel_pi_LL_pT_tnudir;   //!
   TBranch        *b_sel_pi_LL_pT_tpimom;   //!
   TBranch        *b_sel_pi_endpos;   //!
   TBranch        *b_sel_pi_startdir;   //!
   TBranch        *b_sel_pi_startpos;   //!
   TBranch        *b_sel_pi_true4mom;   //!
   TBranch        *b_sel_pi_trueendpos;   //!
   TBranch        *b_sel_pi_truepT;   //!
   TBranch        *b_sel_pi_truestartdir;   //!
   TBranch        *b_sel_pi_truestartpos;   //!
   TBranch        *b_sel_pr_EX_4mom;   //!
   TBranch        *b_sel_pr_EX_pT;   //!
   TBranch        *b_sel_pr_EX_pT_tnudir;   //!
   TBranch        *b_sel_pr_EX_pT_tprmom;   //!
   TBranch        *b_sel_pr_LL_4mom;   //!
   TBranch        *b_sel_pr_LL_pT;   //!
   TBranch        *b_sel_pr_LL_pT_tnudir;   //!
   TBranch        *b_sel_pr_LL_pT_tprmom;   //!
   TBranch        *b_sel_pr_endpos;   //!
   TBranch        *b_sel_pr_startdir;   //!
   TBranch        *b_sel_pr_startpos;   //!
   TBranch        *b_sel_pr_true4mom;   //!
   TBranch        *b_sel_pr_trueendpos;   //!
   TBranch        *b_sel_pr_truepT;   //!
   TBranch        *b_sel_pr_truestartdir;   //!
   TBranch        *b_sel_pr_truestartpos;   //!
   TBranch        *b_sel_truePDP;   //!
   TBranch        *b_sel_trueVTX;   //!
   TBranch        *b_sel_true_nu_dir_PDP;   //!
   TBranch        *b_sel_truedpT_vec;   //!
   TBranch        *b_physEvtNum;   //!
   TBranch        *b_n_hyps;   //!
   TBranch        *b_processType;   //!
   TBranch        *b_primaryPart;   //!
   TBranch        *b_n_slices;   //!
   TBranch        *b_slice_numbers;   //!
   TBranch        *b_shared_slice;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vtxErr;   //!
   TBranch        *b_E;   //!
   TBranch        *b_found_truth;   //!
   TBranch        *b_isMinosMatchTrack;   //!
   TBranch        *b_isMinosMatchStub;   //!
   TBranch        *b_contained_evt;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_n_anchored_long_trk_prongs;   //!
   TBranch        *b_n_anchored_short_trk_prongs;   //!
   TBranch        *b_n_iso_trk_prongs;   //!
   TBranch        *b_n_prongs;   //!
   TBranch        *b_n_tracks3;   //!
   TBranch        *b_ncuts;   //!
   TBranch        *b_new_tracks;   //!
   TBranch        *b_nsplits;   //!
   TBranch        *b_target_region;   //!
   TBranch        *b_true_target_region;   //!
   TBranch        *b_vert_exists;   //!
   TBranch        *b_accum_level;   //!
   TBranch        *b_truth_reco_isMinosMatch;   //!
   TBranch        *b_truth_muon_charge;   //!
   TBranch        *b_truth_n_ele;   //!
   TBranch        *b_truth_n_kPM;   //!
   TBranch        *b_truth_n_kaO;   //!
   TBranch        *b_truth_n_muo;   //!
   TBranch        *b_truth_n_ntn;   //!
   TBranch        *b_truth_n_pho;   //!
   TBranch        *b_truth_n_pi0;   //!
   TBranch        *b_truth_n_piM;   //!
   TBranch        *b_truth_n_piP;   //!
   TBranch        *b_truth_n_pro;   //!
   TBranch        *b_truth_n_tau;   //!
   TBranch        *b_truth_ncuts;   //!
   TBranch        *b_truth_nsplits;   //!
   TBranch        *b_truth_pi_EX_michel;   //!
   TBranch        *b_truth_pi_LL_michel;   //!
   TBranch        *b_truth_pr_EX_michel;   //!
   TBranch        *b_truth_pr_LL_michel;   //!
   TBranch        *b_truth_reco_target;   //!
   TBranch        *b_truth_should_be_accepted;   //!
   TBranch        *b_truth_true_target_region;   //!
   TBranch        *b_truth_mu_E;   //!
   TBranch        *b_truth_mu_KE;   //!
   TBranch        *b_truth_mu_Phi;   //!
   TBranch        *b_truth_mu_Theta;   //!
   TBranch        *b_truth_mu_mom;   //!
   TBranch        *b_truth_mu_pTMag;   //!
   TBranch        *b_truth_mu_pTT;   //!
   TBranch        *b_truth_pi_E;   //!
   TBranch        *b_truth_pi_EX_score;   //!
   TBranch        *b_truth_pi_EX_score_altH;   //!
   TBranch        *b_truth_pi_KE;   //!
   TBranch        *b_truth_pi_LL_score;   //!
   TBranch        *b_truth_pi_LL_score_altH;   //!
   TBranch        *b_truth_pi_Phi;   //!
   TBranch        *b_truth_pi_Theta;   //!
   TBranch        *b_truth_pi_mom;   //!
   TBranch        *b_truth_pi_pTMag;   //!
   TBranch        *b_truth_pi_pTT;   //!
   TBranch        *b_truth_pr_E;   //!
   TBranch        *b_truth_pr_EX_score;   //!
   TBranch        *b_truth_pr_EX_score_altH;   //!
   TBranch        *b_truth_pr_KE;   //!
   TBranch        *b_truth_pr_LL_score;   //!
   TBranch        *b_truth_pr_LL_score_altH;   //!
   TBranch        *b_truth_pr_Phi;   //!
   TBranch        *b_truth_pr_Theta;   //!
   TBranch        *b_truth_pr_mom;   //!
   TBranch        *b_truth_pr_pTMag;   //!
   TBranch        *b_truth_pr_pTT;   //!
   TBranch        *b_truth_trueEnu;   //!
   TBranch        *b_truth_trueQ2;   //!
   TBranch        *b_truth_truedalphaT;   //!
   TBranch        *b_truth_truedpT;   //!
   TBranch        *b_truth_truedpTT;   //!
   TBranch        *b_truth_truedpTT_pi;   //!
   TBranch        *b_truth_truedpTT_pi_dir;   //!
   TBranch        *b_truth_truedpTT_pr;   //!
   TBranch        *b_truth_truedpTT_pr_dir;   //!
   TBranch        *b_truth_truedphiT;   //!
   TBranch        *b_truth_accum_level;   //!
   TBranch        *b_truth_mu_4mom;   //!
   TBranch        *b_truth_mu_pT;   //!
   TBranch        *b_truth_pi_4mom;   //!
   TBranch        *b_truth_pi_pT;   //!
   TBranch        *b_truth_pr_4mom;   //!
   TBranch        *b_truth_pr_pT;   //!
   TBranch        *b_truth_truedpT_vec;   //!
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_subrun;   //!
   TBranch        *b_ev_detector;   //!
   TBranch        *b_ev_triggerType;   //!
   TBranch        *b_ev_gate;   //!
   TBranch        *b_ev_global_gate;   //!
   TBranch        *b_ev_gps_time_sec;   //!
   TBranch        *b_ev_gps_time_usec;   //!
   TBranch        *b_mc_run;   //!
   TBranch        *b_mc_subrun;   //!
   TBranch        *b_mc_nInteractions;   //!
   TBranch        *b_mc_MIState;   //!
   TBranch        *b_mc_pot;   //!
   TBranch        *b_mc_beamConfig;   //!
   TBranch        *b_mc_processType;   //!
   TBranch        *b_mc_nthEvtInSpill;   //!
   TBranch        *b_mc_nthEvtInFile;   //!
   TBranch        *b_mc_intType;   //!
   TBranch        *b_mc_current;   //!
   TBranch        *b_mc_charm;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_XSec;   //!
   TBranch        *b_mc_diffXSec;   //!
   TBranch        *b_mc_incoming;   //!
   TBranch        *b_mc_fluxDriverProb;   //!
   TBranch        *b_mc_targetNucleus;   //!
   TBranch        *b_mc_targetZ;   //!
   TBranch        *b_mc_targetA;   //!
   TBranch        *b_mc_targetNucleon;   //!
   TBranch        *b_mc_struckQuark;   //!
   TBranch        *b_mc_seaQuark;   //!
   TBranch        *b_mc_resID;   //!
   TBranch        *b_mc_primaryLepton;   //!
   TBranch        *b_mc_incomingE;   //!
   TBranch        *b_mc_Bjorkenx;   //!
   TBranch        *b_mc_Bjorkeny;   //!
   TBranch        *b_mc_Q2;   //!
   TBranch        *b_mc_nuT;   //!
   TBranch        *b_mc_w;   //!
   TBranch        *b_mc_vtx;   //!
   TBranch        *b_mc_incomingPartVec;   //!
   TBranch        *b_mc_initNucVec;   //!
   TBranch        *b_mc_primFSLepton;   //!
   TBranch        *b_mc_nFSPart;   //!
   TBranch        *b_mc_FSPartPx;   //!
   TBranch        *b_mc_FSPartPy;   //!
   TBranch        *b_mc_FSPartPz;   //!
   TBranch        *b_mc_FSPartE;   //!
   TBranch        *b_mc_FSPartPDG;   //!
   TBranch        *b_mc_er_nPart;   //!
   TBranch        *b_mc_er_ID;   //!
   TBranch        *b_mc_er_status;   //!
   TBranch        *b_mc_er_posInNucX;   //!
   TBranch        *b_mc_er_posInNucY;   //!
   TBranch        *b_mc_er_posInNucZ;   //!
   TBranch        *b_mc_er_Px;   //!
   TBranch        *b_mc_er_Py;   //!
   TBranch        *b_mc_er_Pz;   //!
   TBranch        *b_mc_er_E;   //!
   TBranch        *b_mc_er_FD;   //!
   TBranch        *b_mc_er_LD;   //!
   TBranch        *b_mc_er_mother;   //!
   TBranch        *b_mc_fr_nNuAncestorIDs;   //!
   TBranch        *b_mc_fr_nuAncestorIDs;   //!
   TBranch        *b_mc_fr_nuParentID;   //!
   TBranch        *b_mc_fr_decMode;   //!
   TBranch        *b_mc_fr_primProtonVtx;   //!
   TBranch        *b_mc_fr_primProtonP;   //!
   TBranch        *b_mc_fr_nuParentDecVtx;   //!
   TBranch        *b_mc_fr_nuParentProdVtx;   //!
   TBranch        *b_mc_fr_nuParentProdP;   //!
   TBranch        *b_mc_cvweight_total;   //!
   TBranch        *b_wgt;   //!
   TBranch        *b_mc_cvweight_totalFlux;   //!
   TBranch        *b_mc_cvweight_totalXsec;   //!
   TBranch        *b_mc_ppfx1_cvweight;   //!
   TBranch        *b_mc_hornCurrent_cvweight;   //!
   TBranch        *b_mc_gen1_cvweight_total;   //!
   TBranch        *b_gen1_wgt;   //!
   TBranch        *b_mc_gen1_cvweight_totalFlux;   //!
   TBranch        *b_mc_gen1_cvweight_NA49;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus_sz;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus;   //!
   TBranch        *b_mc_wgt_gen1_Flux_Tertiary_sz;   //!
   TBranch        *b_mc_wgt_gen1_Flux_Tertiary;   //!
   TBranch        *b_mc_wgt_gen1_Flux_NA49_sz;   //!
   TBranch        *b_mc_wgt_gen1_Flux_NA49;   //!
   TBranch        *b_mc_wgt_Norm_sz;   //!
   TBranch        *b_mc_wgt_Norm;   //!
   TBranch        *b_mc_wgt_ppfx1_Total_sz;   //!
   TBranch        *b_mc_wgt_ppfx1_Total;   //!
   TBranch        *b_prong_nParticles;   //!
   TBranch        *b_prong_part_score;   //!
   TBranch        *b_prong_part_mass;   //!
   TBranch        *b_prong_part_charge;   //!
   TBranch        *b_prong_part_pid;   //!
   // TBranch        *b_prong_part_E;   //!
   // TBranch        *b_prong_part_pos;   //!

 AnalysisReader(TTree *tree, TTree *savetree);
    Int_t GetEntry(Long64_t entry);
    Int_t GetEntries();
    virtual ~AnalysisReader();
    virtual void     Init(TTree *tree);

    void SetOutTree();
    void FillOutTree(){ m_savetree->Fill(); }

 private:

TTree * m_savetree;
};

#endif

#ifdef AnalysisReader_cxx
AnalysisReader::AnalysisReader(TTree *tree, TTree *savetree) : fChain(0)
{
   Init(tree);
   m_savetree = savetree;
}

AnalysisReader::~AnalysisReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Int_t AnalysisReader::GetEntries(){
    if(!fChain) return 0;
    return fChain->GetEntries();
}

void AnalysisReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   // prong_part_E = 0;
   // prong_part_pos = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("sel_nuFlavor", &sel_nuFlavor, &b_sel_nuFlavor);
   fChain->SetBranchAddress("sel_nuHelicity", &sel_nuHelicity, &b_sel_nuHelicity);
   fChain->SetBranchAddress("sel_intCurrent", &sel_intCurrent, &b_sel_intCurrent);
   fChain->SetBranchAddress("sel_intType", &sel_intType, &b_sel_intType);
   fChain->SetBranchAddress("sel_E", &sel_E, &b_sel_E);
   fChain->SetBranchAddress("sel_Q2", &sel_Q2, &b_sel_Q2);
   fChain->SetBranchAddress("sel_x", &sel_x, &b_sel_x);
   fChain->SetBranchAddress("sel_y", &sel_y, &b_sel_y);
   fChain->SetBranchAddress("sel_W", &sel_W, &b_sel_W);
   fChain->SetBranchAddress("sel_score", &sel_score, &b_sel_score);
   fChain->SetBranchAddress("sel_leptonE", sel_leptonE, &b_sel_leptonE);
   fChain->SetBranchAddress("sel_vtx", sel_vtx, &b_sel_vtx);
   fChain->SetBranchAddress("sel_mu_PDG", &sel_mu_PDG, &b_sel_mu_PDG);
   fChain->SetBranchAddress("sel_mu_isKinked", &sel_mu_isKinked, &b_sel_mu_isKinked);
   fChain->SetBranchAddress("sel_mu_michel", &sel_mu_michel, &b_sel_mu_michel);
   fChain->SetBranchAddress("sel_pi_EX_michel", &sel_pi_EX_michel, &b_sel_pi_EX_michel);
   fChain->SetBranchAddress("sel_pi_EX_nNodes", &sel_pi_EX_nNodes, &b_sel_pi_EX_nNodes);
   fChain->SetBranchAddress("sel_pi_FSI", &sel_pi_FSI, &b_sel_pi_FSI);
   fChain->SetBranchAddress("sel_pi_LL_michel", &sel_pi_LL_michel, &b_sel_pi_LL_michel);
   fChain->SetBranchAddress("sel_pi_LL_nNodes", &sel_pi_LL_nNodes, &b_sel_pi_LL_nNodes);
   fChain->SetBranchAddress("sel_pi_PDG", &sel_pi_PDG, &b_sel_pi_PDG);
   fChain->SetBranchAddress("sel_pi_isKinked", &sel_pi_isKinked, &b_sel_pi_isKinked);
   fChain->SetBranchAddress("sel_pr_EX_michel", &sel_pr_EX_michel, &b_sel_pr_EX_michel);
   fChain->SetBranchAddress("sel_pr_EX_nNodes", &sel_pr_EX_nNodes, &b_sel_pr_EX_nNodes);
   fChain->SetBranchAddress("sel_pr_FSI", &sel_pr_FSI, &b_sel_pr_FSI);
   fChain->SetBranchAddress("sel_pr_LL_michel", &sel_pr_LL_michel, &b_sel_pr_LL_michel);
   fChain->SetBranchAddress("sel_pr_LL_nNodes", &sel_pr_LL_nNodes, &b_sel_pr_LL_nNodes);
   fChain->SetBranchAddress("sel_pr_PDG", &sel_pr_PDG, &b_sel_pr_PDG);
   fChain->SetBranchAddress("sel_pr_isKinked", &sel_pr_isKinked, &b_sel_pr_isKinked);
   fChain->SetBranchAddress("sel_Enu_EX", &sel_Enu_EX, &b_sel_Enu_EX);
   fChain->SetBranchAddress("sel_Enu_LL", &sel_Enu_LL, &b_sel_Enu_LL);
   fChain->SetBranchAddress("sel_Q2_EX", &sel_Q2_EX, &b_sel_Q2_EX);
   fChain->SetBranchAddress("sel_Q2_LL", &sel_Q2_LL, &b_sel_Q2_LL);
   fChain->SetBranchAddress("sel_dalphaT_EX", &sel_dalphaT_EX, &b_sel_dalphaT_EX);
   fChain->SetBranchAddress("sel_dalphaT_LL", &sel_dalphaT_LL, &b_sel_dalphaT_LL);
   fChain->SetBranchAddress("sel_dpTT_EX", &sel_dpTT_EX, &b_sel_dpTT_EX);
   fChain->SetBranchAddress("sel_dpTT_EX_tmumom", &sel_dpTT_EX_tmumom, &b_sel_dpTT_EX_tmumom);
   fChain->SetBranchAddress("sel_dpTT_EX_tnudir", &sel_dpTT_EX_tnudir, &b_sel_dpTT_EX_tnudir);
   fChain->SetBranchAddress("sel_dpTT_EX_tpimom", &sel_dpTT_EX_tpimom, &b_sel_dpTT_EX_tpimom);
   fChain->SetBranchAddress("sel_dpTT_EX_tprmom", &sel_dpTT_EX_tprmom, &b_sel_dpTT_EX_tprmom);
   fChain->SetBranchAddress("sel_dpTT_LL", &sel_dpTT_LL, &b_sel_dpTT_LL);
   fChain->SetBranchAddress("sel_dpTT_LL_tmumom", &sel_dpTT_LL_tmumom, &b_sel_dpTT_LL_tmumom);
   fChain->SetBranchAddress("sel_dpTT_LL_tnudir", &sel_dpTT_LL_tnudir, &b_sel_dpTT_LL_tnudir);
   fChain->SetBranchAddress("sel_dpTT_LL_tpimom", &sel_dpTT_LL_tpimom, &b_sel_dpTT_LL_tpimom);
   fChain->SetBranchAddress("sel_dpTT_LL_tprmom", &sel_dpTT_LL_tprmom, &b_sel_dpTT_LL_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pi_EX", &sel_dpTT_pi_EX, &b_sel_dpTT_pi_EX);
   fChain->SetBranchAddress("sel_dpTT_pi_EX_tmumom", &sel_dpTT_pi_EX_tmumom, &b_sel_dpTT_pi_EX_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pi_EX_tnudir", &sel_dpTT_pi_EX_tnudir, &b_sel_dpTT_pi_EX_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pi_EX_tpimom", &sel_dpTT_pi_EX_tpimom, &b_sel_dpTT_pi_EX_tpimom);
   fChain->SetBranchAddress("sel_dpTT_pi_EX_tprmom", &sel_dpTT_pi_EX_tprmom, &b_sel_dpTT_pi_EX_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pi_LL", &sel_dpTT_pi_LL, &b_sel_dpTT_pi_LL);
   fChain->SetBranchAddress("sel_dpTT_pi_LL_tmumom", &sel_dpTT_pi_LL_tmumom, &b_sel_dpTT_pi_LL_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pi_LL_tnudir", &sel_dpTT_pi_LL_tnudir, &b_sel_dpTT_pi_LL_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pi_LL_tpimom", &sel_dpTT_pi_LL_tpimom, &b_sel_dpTT_pi_LL_tpimom);
   fChain->SetBranchAddress("sel_dpTT_pi_LL_tprmom", &sel_dpTT_pi_LL_tprmom, &b_sel_dpTT_pi_LL_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_EX", &sel_dpTT_pi_dir_EX, &b_sel_dpTT_pi_dir_EX);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_EX_tmumom", &sel_dpTT_pi_dir_EX_tmumom, &b_sel_dpTT_pi_dir_EX_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_EX_tnudir", &sel_dpTT_pi_dir_EX_tnudir, &b_sel_dpTT_pi_dir_EX_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_EX_tpidir", &sel_dpTT_pi_dir_EX_tpidir, &b_sel_dpTT_pi_dir_EX_tpidir);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_EX_tprmom", &sel_dpTT_pi_dir_EX_tprmom, &b_sel_dpTT_pi_dir_EX_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_LL", &sel_dpTT_pi_dir_LL, &b_sel_dpTT_pi_dir_LL);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_LL_tmumom", &sel_dpTT_pi_dir_LL_tmumom, &b_sel_dpTT_pi_dir_LL_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_LL_tnudir", &sel_dpTT_pi_dir_LL_tnudir, &b_sel_dpTT_pi_dir_LL_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_LL_tpidir", &sel_dpTT_pi_dir_LL_tpidir, &b_sel_dpTT_pi_dir_LL_tpidir);
   fChain->SetBranchAddress("sel_dpTT_pi_dir_LL_tprmom", &sel_dpTT_pi_dir_LL_tprmom, &b_sel_dpTT_pi_dir_LL_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pr_EX", &sel_dpTT_pr_EX, &b_sel_dpTT_pr_EX);
   fChain->SetBranchAddress("sel_dpTT_pr_EX_tmumom", &sel_dpTT_pr_EX_tmumom, &b_sel_dpTT_pr_EX_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pr_EX_tnudir", &sel_dpTT_pr_EX_tnudir, &b_sel_dpTT_pr_EX_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pr_EX_tpimon", &sel_dpTT_pr_EX_tpimon, &b_sel_dpTT_pr_EX_tpimon);
   fChain->SetBranchAddress("sel_dpTT_pr_EX_tprmom", &sel_dpTT_pr_EX_tprmom, &b_sel_dpTT_pr_EX_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pr_LL", &sel_dpTT_pr_LL, &b_sel_dpTT_pr_LL);
   fChain->SetBranchAddress("sel_dpTT_pr_LL_tmumom", &sel_dpTT_pr_LL_tmumom, &b_sel_dpTT_pr_LL_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pr_LL_tnudir", &sel_dpTT_pr_LL_tnudir, &b_sel_dpTT_pr_LL_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pr_LL_tpimon", &sel_dpTT_pr_LL_tpimon, &b_sel_dpTT_pr_LL_tpimon);
   fChain->SetBranchAddress("sel_dpTT_pr_LL_tprmom", &sel_dpTT_pr_LL_tprmom, &b_sel_dpTT_pr_LL_tprmom);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_EX", &sel_dpTT_pr_dir_EX, &b_sel_dpTT_pr_dir_EX);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_EX_tmumom", &sel_dpTT_pr_dir_EX_tmumom, &b_sel_dpTT_pr_dir_EX_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_EX_tnudir", &sel_dpTT_pr_dir_EX_tnudir, &b_sel_dpTT_pr_dir_EX_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_EX_tpimom", &sel_dpTT_pr_dir_EX_tpimom, &b_sel_dpTT_pr_dir_EX_tpimom);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_EX_tprdir", &sel_dpTT_pr_dir_EX_tprdir, &b_sel_dpTT_pr_dir_EX_tprdir);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_LL", &sel_dpTT_pr_dir_LL, &b_sel_dpTT_pr_dir_LL);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_LL_tmumom", &sel_dpTT_pr_dir_LL_tmumom, &b_sel_dpTT_pr_dir_LL_tmumom);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_LL_tnudir", &sel_dpTT_pr_dir_LL_tnudir, &b_sel_dpTT_pr_dir_LL_tnudir);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_LL_tpimom", &sel_dpTT_pr_dir_LL_tpimom, &b_sel_dpTT_pr_dir_LL_tpimom);
   fChain->SetBranchAddress("sel_dpTT_pr_dir_LL_tprdir", &sel_dpTT_pr_dir_LL_tprdir, &b_sel_dpTT_pr_dir_LL_tprdir);
   fChain->SetBranchAddress("sel_dpT_EX", &sel_dpT_EX, &b_sel_dpT_EX);
   fChain->SetBranchAddress("sel_dpT_LL", &sel_dpT_LL, &b_sel_dpT_LL);
   fChain->SetBranchAddress("sel_dphiT_EX", &sel_dphiT_EX, &b_sel_dphiT_EX);
   fChain->SetBranchAddress("sel_dphiT_LL", &sel_dphiT_LL, &b_sel_dphiT_LL);
   fChain->SetBranchAddress("sel_mu_E", &sel_mu_E, &b_sel_mu_E);
   fChain->SetBranchAddress("sel_mu_KE", &sel_mu_KE, &b_sel_mu_KE);
   fChain->SetBranchAddress("sel_mu_Phi", &sel_mu_Phi, &b_sel_mu_Phi);
   fChain->SetBranchAddress("sel_mu_Theta", &sel_mu_Theta, &b_sel_mu_Theta);
   fChain->SetBranchAddress("sel_mu_XTheta", &sel_mu_XTheta, &b_sel_mu_XTheta);
   fChain->SetBranchAddress("sel_mu_YTheta", &sel_mu_YTheta, &b_sel_mu_YTheta);
   fChain->SetBranchAddress("sel_mu_chi2ndf", &sel_mu_chi2ndf, &b_sel_mu_chi2ndf);
   fChain->SetBranchAddress("sel_mu_det_frac", &sel_mu_det_frac, &b_sel_mu_det_frac);
   fChain->SetBranchAddress("sel_mu_det_otherE", &sel_mu_det_otherE, &b_sel_mu_det_otherE);
   fChain->SetBranchAddress("sel_mu_mom", &sel_mu_mom, &b_sel_mu_mom);
   fChain->SetBranchAddress("sel_mu_pTMag", &sel_mu_pTMag, &b_sel_mu_pTMag);
   fChain->SetBranchAddress("sel_mu_pTMag_tmumom", &sel_mu_pTMag_tmumom, &b_sel_mu_pTMag_tmumom);
   fChain->SetBranchAddress("sel_mu_pTMag_tnudir", &sel_mu_pTMag_tnudir, &b_sel_mu_pTMag_tnudir);
   fChain->SetBranchAddress("sel_mu_pTT", &sel_mu_pTT, &b_sel_mu_pTT);
   fChain->SetBranchAddress("sel_mu_pTT_tmumom", &sel_mu_pTT_tmumom, &b_sel_mu_pTT_tmumom);
   fChain->SetBranchAddress("sel_mu_pTT_tnudir", &sel_mu_pTT_tnudir, &b_sel_mu_pTT_tnudir);
   fChain->SetBranchAddress("sel_mu_score", &sel_mu_score, &b_sel_mu_score);
   fChain->SetBranchAddress("sel_mu_trueE", &sel_mu_trueE, &b_sel_mu_trueE);
   fChain->SetBranchAddress("sel_mu_trueKE", &sel_mu_trueKE, &b_sel_mu_trueKE);
   fChain->SetBranchAddress("sel_mu_truePhi", &sel_mu_truePhi, &b_sel_mu_truePhi);
   fChain->SetBranchAddress("sel_mu_trueTheta", &sel_mu_trueTheta, &b_sel_mu_trueTheta);
   fChain->SetBranchAddress("sel_mu_truemom", &sel_mu_truemom, &b_sel_mu_truemom);
   fChain->SetBranchAddress("sel_mu_truepTMag", &sel_mu_truepTMag, &b_sel_mu_truepTMag);
   fChain->SetBranchAddress("sel_mu_truepTT", &sel_mu_truepTT, &b_sel_mu_truepTT);
   fChain->SetBranchAddress("sel_pi_CaloE", &sel_pi_CaloE, &b_sel_pi_CaloE);
   fChain->SetBranchAddress("sel_pi_EX_E", &sel_pi_EX_E, &b_sel_pi_EX_E);
   fChain->SetBranchAddress("sel_pi_EX_E_altH", &sel_pi_EX_E_altH, &b_sel_pi_EX_E_altH);
   fChain->SetBranchAddress("sel_pi_EX_KE", &sel_pi_EX_KE, &b_sel_pi_EX_KE);
   fChain->SetBranchAddress("sel_pi_EX_Phi", &sel_pi_EX_Phi, &b_sel_pi_EX_Phi);
   fChain->SetBranchAddress("sel_pi_EX_Theta", &sel_pi_EX_Theta, &b_sel_pi_EX_Theta);
   fChain->SetBranchAddress("sel_pi_EX_XTheta", &sel_pi_EX_XTheta, &b_sel_pi_EX_XTheta);
   fChain->SetBranchAddress("sel_pi_EX_YTheta", &sel_pi_EX_YTheta, &b_sel_pi_EX_YTheta);
   fChain->SetBranchAddress("sel_pi_EX_ln_Q0", &sel_pi_EX_ln_Q0, &b_sel_pi_EX_ln_Q0);
   fChain->SetBranchAddress("sel_pi_EX_ln_Q1", &sel_pi_EX_ln_Q1, &b_sel_pi_EX_ln_Q1);
   fChain->SetBranchAddress("sel_pi_EX_ln_Q2", &sel_pi_EX_ln_Q2, &b_sel_pi_EX_ln_Q2);
   fChain->SetBranchAddress("sel_pi_EX_ln_Q3", &sel_pi_EX_ln_Q3, &b_sel_pi_EX_ln_Q3);
   fChain->SetBranchAddress("sel_pi_EX_ln_Q4", &sel_pi_EX_ln_Q4, &b_sel_pi_EX_ln_Q4);
   fChain->SetBranchAddress("sel_pi_EX_ln_Q5", &sel_pi_EX_ln_Q5, &b_sel_pi_EX_ln_Q5);
   fChain->SetBranchAddress("sel_pi_EX_mom", &sel_pi_EX_mom, &b_sel_pi_EX_mom);
   fChain->SetBranchAddress("sel_pi_EX_mom_altH", &sel_pi_EX_mom_altH, &b_sel_pi_EX_mom_altH);
   fChain->SetBranchAddress("sel_pi_EX_pTMag", &sel_pi_EX_pTMag, &b_sel_pi_EX_pTMag);
   fChain->SetBranchAddress("sel_pi_EX_pTMag_tnudir", &sel_pi_EX_pTMag_tnudir, &b_sel_pi_EX_pTMag_tnudir);
   fChain->SetBranchAddress("sel_pi_EX_pTMag_tpimom", &sel_pi_EX_pTMag_tpimom, &b_sel_pi_EX_pTMag_tpimom);
   fChain->SetBranchAddress("sel_pi_EX_pTT", &sel_pi_EX_pTT, &b_sel_pi_EX_pTT);
   fChain->SetBranchAddress("sel_pi_EX_pTT_tnudir", &sel_pi_EX_pTT_tnudir, &b_sel_pi_EX_pTT_tnudir);
   fChain->SetBranchAddress("sel_pi_EX_pTT_tpimom", &sel_pi_EX_pTT_tpimom, &b_sel_pi_EX_pTT_tpimom);
   fChain->SetBranchAddress("sel_pi_EX_score", &sel_pi_EX_score, &b_sel_pi_EX_score);
   fChain->SetBranchAddress("sel_pi_EX_score_altH", &sel_pi_EX_score_altH, &b_sel_pi_EX_score_altH);
   fChain->SetBranchAddress("sel_pi_LL_E", &sel_pi_LL_E, &b_sel_pi_LL_E);
   fChain->SetBranchAddress("sel_pi_LL_E_altH", &sel_pi_LL_E_altH, &b_sel_pi_LL_E_altH);
   fChain->SetBranchAddress("sel_pi_LL_KE", &sel_pi_LL_KE, &b_sel_pi_LL_KE);
   fChain->SetBranchAddress("sel_pi_LL_Phi", &sel_pi_LL_Phi, &b_sel_pi_LL_Phi);
   fChain->SetBranchAddress("sel_pi_LL_Theta", &sel_pi_LL_Theta, &b_sel_pi_LL_Theta);
   fChain->SetBranchAddress("sel_pi_LL_XTheta", &sel_pi_LL_XTheta, &b_sel_pi_LL_XTheta);
   fChain->SetBranchAddress("sel_pi_LL_YTheta", &sel_pi_LL_YTheta, &b_sel_pi_LL_YTheta);
   fChain->SetBranchAddress("sel_pi_LL_ln_Q0", &sel_pi_LL_ln_Q0, &b_sel_pi_LL_ln_Q0);
   fChain->SetBranchAddress("sel_pi_LL_ln_Q1", &sel_pi_LL_ln_Q1, &b_sel_pi_LL_ln_Q1);
   fChain->SetBranchAddress("sel_pi_LL_ln_Q2", &sel_pi_LL_ln_Q2, &b_sel_pi_LL_ln_Q2);
   fChain->SetBranchAddress("sel_pi_LL_ln_Q3", &sel_pi_LL_ln_Q3, &b_sel_pi_LL_ln_Q3);
   fChain->SetBranchAddress("sel_pi_LL_ln_Q4", &sel_pi_LL_ln_Q4, &b_sel_pi_LL_ln_Q4);
   fChain->SetBranchAddress("sel_pi_LL_ln_Q5", &sel_pi_LL_ln_Q5, &b_sel_pi_LL_ln_Q5);
   fChain->SetBranchAddress("sel_pi_LL_mom", &sel_pi_LL_mom, &b_sel_pi_LL_mom);
   fChain->SetBranchAddress("sel_pi_LL_mom_altH", &sel_pi_LL_mom_altH, &b_sel_pi_LL_mom_altH);
   fChain->SetBranchAddress("sel_pi_LL_pTMag", &sel_pi_LL_pTMag, &b_sel_pi_LL_pTMag);
   fChain->SetBranchAddress("sel_pi_LL_pTMag_tnudir", &sel_pi_LL_pTMag_tnudir, &b_sel_pi_LL_pTMag_tnudir);
   fChain->SetBranchAddress("sel_pi_LL_pTMag_tpimom", &sel_pi_LL_pTMag_tpimom, &b_sel_pi_LL_pTMag_tpimom);
   fChain->SetBranchAddress("sel_pi_LL_pTT", &sel_pi_LL_pTT, &b_sel_pi_LL_pTT);
   fChain->SetBranchAddress("sel_pi_LL_pTT_tnudir", &sel_pi_LL_pTT_tnudir, &b_sel_pi_LL_pTT_tnudir);
   fChain->SetBranchAddress("sel_pi_LL_pTT_tpimom", &sel_pi_LL_pTT_tpimom, &b_sel_pi_LL_pTT_tpimom);
   fChain->SetBranchAddress("sel_pi_LL_score", &sel_pi_LL_score, &b_sel_pi_LL_score);
   fChain->SetBranchAddress("sel_pi_LL_score_altH", &sel_pi_LL_score_altH, &b_sel_pi_LL_score_altH);
   fChain->SetBranchAddress("sel_pi_chi2ndf", &sel_pi_chi2ndf, &b_sel_pi_chi2ndf);
   fChain->SetBranchAddress("sel_pi_det_frac", &sel_pi_det_frac, &b_sel_pi_det_frac);
   fChain->SetBranchAddress("sel_pi_det_otherE", &sel_pi_det_otherE, &b_sel_pi_det_otherE);
   fChain->SetBranchAddress("sel_pi_trueE", &sel_pi_trueE, &b_sel_pi_trueE);
   fChain->SetBranchAddress("sel_pi_trueKE", &sel_pi_trueKE, &b_sel_pi_trueKE);
   fChain->SetBranchAddress("sel_pi_truePhi", &sel_pi_truePhi, &b_sel_pi_truePhi);
   fChain->SetBranchAddress("sel_pi_trueTheta", &sel_pi_trueTheta, &b_sel_pi_trueTheta);
   fChain->SetBranchAddress("sel_pi_truemom", &sel_pi_truemom, &b_sel_pi_truemom);
   fChain->SetBranchAddress("sel_pi_truepTMag", &sel_pi_truepTMag, &b_sel_pi_truepTMag);
   fChain->SetBranchAddress("sel_pi_truepTT", &sel_pi_truepTT, &b_sel_pi_truepTT);
   fChain->SetBranchAddress("sel_pr_CaloE", &sel_pr_CaloE, &b_sel_pr_CaloE);
   fChain->SetBranchAddress("sel_pr_EX_E", &sel_pr_EX_E, &b_sel_pr_EX_E);
   fChain->SetBranchAddress("sel_pr_EX_E_altH", &sel_pr_EX_E_altH, &b_sel_pr_EX_E_altH);
   fChain->SetBranchAddress("sel_pr_EX_KE", &sel_pr_EX_KE, &b_sel_pr_EX_KE);
   fChain->SetBranchAddress("sel_pr_EX_Phi", &sel_pr_EX_Phi, &b_sel_pr_EX_Phi);
   fChain->SetBranchAddress("sel_pr_EX_Theta", &sel_pr_EX_Theta, &b_sel_pr_EX_Theta);
   fChain->SetBranchAddress("sel_pr_EX_XTheta", &sel_pr_EX_XTheta, &b_sel_pr_EX_XTheta);
   fChain->SetBranchAddress("sel_pr_EX_YTheta", &sel_pr_EX_YTheta, &b_sel_pr_EX_YTheta);
   fChain->SetBranchAddress("sel_pr_EX_ln_Q0", &sel_pr_EX_ln_Q0, &b_sel_pr_EX_ln_Q0);
   fChain->SetBranchAddress("sel_pr_EX_ln_Q1", &sel_pr_EX_ln_Q1, &b_sel_pr_EX_ln_Q1);
   fChain->SetBranchAddress("sel_pr_EX_ln_Q2", &sel_pr_EX_ln_Q2, &b_sel_pr_EX_ln_Q2);
   fChain->SetBranchAddress("sel_pr_EX_ln_Q3", &sel_pr_EX_ln_Q3, &b_sel_pr_EX_ln_Q3);
   fChain->SetBranchAddress("sel_pr_EX_ln_Q4", &sel_pr_EX_ln_Q4, &b_sel_pr_EX_ln_Q4);
   fChain->SetBranchAddress("sel_pr_EX_ln_Q5", &sel_pr_EX_ln_Q5, &b_sel_pr_EX_ln_Q5);
   fChain->SetBranchAddress("sel_pr_EX_mom", &sel_pr_EX_mom, &b_sel_pr_EX_mom);
   fChain->SetBranchAddress("sel_pr_EX_mom_altH", &sel_pr_EX_mom_altH, &b_sel_pr_EX_mom_altH);
   fChain->SetBranchAddress("sel_pr_EX_pTMag", &sel_pr_EX_pTMag, &b_sel_pr_EX_pTMag);
   fChain->SetBranchAddress("sel_pr_EX_pTMag_tnudir", &sel_pr_EX_pTMag_tnudir, &b_sel_pr_EX_pTMag_tnudir);
   fChain->SetBranchAddress("sel_pr_EX_pTMag_tprmom", &sel_pr_EX_pTMag_tprmom, &b_sel_pr_EX_pTMag_tprmom);
   fChain->SetBranchAddress("sel_pr_EX_pTT", &sel_pr_EX_pTT, &b_sel_pr_EX_pTT);
   fChain->SetBranchAddress("sel_pr_EX_pTT_tnudir", &sel_pr_EX_pTT_tnudir, &b_sel_pr_EX_pTT_tnudir);
   fChain->SetBranchAddress("sel_pr_EX_pTT_tprmom", &sel_pr_EX_pTT_tprmom, &b_sel_pr_EX_pTT_tprmom);
   fChain->SetBranchAddress("sel_pr_EX_score", &sel_pr_EX_score, &b_sel_pr_EX_score);
   fChain->SetBranchAddress("sel_pr_EX_score_altH", &sel_pr_EX_score_altH, &b_sel_pr_EX_score_altH);
   fChain->SetBranchAddress("sel_pr_LL_E", &sel_pr_LL_E, &b_sel_pr_LL_E);
   fChain->SetBranchAddress("sel_pr_LL_E_altH", &sel_pr_LL_E_altH, &b_sel_pr_LL_E_altH);
   fChain->SetBranchAddress("sel_pr_LL_KE", &sel_pr_LL_KE, &b_sel_pr_LL_KE);
   fChain->SetBranchAddress("sel_pr_LL_Phi", &sel_pr_LL_Phi, &b_sel_pr_LL_Phi);
   fChain->SetBranchAddress("sel_pr_LL_Theta", &sel_pr_LL_Theta, &b_sel_pr_LL_Theta);
   fChain->SetBranchAddress("sel_pr_LL_XTheta", &sel_pr_LL_XTheta, &b_sel_pr_LL_XTheta);
   fChain->SetBranchAddress("sel_pr_LL_YTheta", &sel_pr_LL_YTheta, &b_sel_pr_LL_YTheta);
   fChain->SetBranchAddress("sel_pr_LL_ln_Q0", &sel_pr_LL_ln_Q0, &b_sel_pr_LL_ln_Q0);
   fChain->SetBranchAddress("sel_pr_LL_ln_Q1", &sel_pr_LL_ln_Q1, &b_sel_pr_LL_ln_Q1);
   fChain->SetBranchAddress("sel_pr_LL_ln_Q2", &sel_pr_LL_ln_Q2, &b_sel_pr_LL_ln_Q2);
   fChain->SetBranchAddress("sel_pr_LL_ln_Q3", &sel_pr_LL_ln_Q3, &b_sel_pr_LL_ln_Q3);
   fChain->SetBranchAddress("sel_pr_LL_ln_Q4", &sel_pr_LL_ln_Q4, &b_sel_pr_LL_ln_Q4);
   fChain->SetBranchAddress("sel_pr_LL_ln_Q5", &sel_pr_LL_ln_Q5, &b_sel_pr_LL_ln_Q5);
   fChain->SetBranchAddress("sel_pr_LL_mom", &sel_pr_LL_mom, &b_sel_pr_LL_mom);
   fChain->SetBranchAddress("sel_pr_LL_mom_altH", &sel_pr_LL_mom_altH, &b_sel_pr_LL_mom_altH);
   fChain->SetBranchAddress("sel_pr_LL_pTMag", &sel_pr_LL_pTMag, &b_sel_pr_LL_pTMag);
   fChain->SetBranchAddress("sel_pr_LL_pTMag_tnudir", &sel_pr_LL_pTMag_tnudir, &b_sel_pr_LL_pTMag_tnudir);
   fChain->SetBranchAddress("sel_pr_LL_pTMag_tprmom", &sel_pr_LL_pTMag_tprmom, &b_sel_pr_LL_pTMag_tprmom);
   fChain->SetBranchAddress("sel_pr_LL_pTT", &sel_pr_LL_pTT, &b_sel_pr_LL_pTT);
   fChain->SetBranchAddress("sel_pr_LL_pTT_tnudir", &sel_pr_LL_pTT_tnudir, &b_sel_pr_LL_pTT_tnudir);
   fChain->SetBranchAddress("sel_pr_LL_pTT_tprmom", &sel_pr_LL_pTT_tprmom, &b_sel_pr_LL_pTT_tprmom);
   fChain->SetBranchAddress("sel_pr_LL_score", &sel_pr_LL_score, &b_sel_pr_LL_score);
   fChain->SetBranchAddress("sel_pr_LL_score_altH", &sel_pr_LL_score_altH, &b_sel_pr_LL_score_altH);
   fChain->SetBranchAddress("sel_pr_chi2ndf", &sel_pr_chi2ndf, &b_sel_pr_chi2ndf);
   fChain->SetBranchAddress("sel_pr_det_frac", &sel_pr_det_frac, &b_sel_pr_det_frac);
   fChain->SetBranchAddress("sel_pr_det_otherE", &sel_pr_det_otherE, &b_sel_pr_det_otherE);
   fChain->SetBranchAddress("sel_pr_trueE", &sel_pr_trueE, &b_sel_pr_trueE);
   fChain->SetBranchAddress("sel_pr_trueKE", &sel_pr_trueKE, &b_sel_pr_trueKE);
   fChain->SetBranchAddress("sel_pr_truePhi", &sel_pr_truePhi, &b_sel_pr_truePhi);
   fChain->SetBranchAddress("sel_pr_trueTheta", &sel_pr_trueTheta, &b_sel_pr_trueTheta);
   fChain->SetBranchAddress("sel_pr_truemom", &sel_pr_truemom, &b_sel_pr_truemom);
   fChain->SetBranchAddress("sel_pr_truepTMag", &sel_pr_truepTMag, &b_sel_pr_truepTMag);
   fChain->SetBranchAddress("sel_pr_truepTT", &sel_pr_truepTT, &b_sel_pr_truepTT);
   fChain->SetBranchAddress("sel_trueEnu", &sel_trueEnu, &b_sel_trueEnu);
   fChain->SetBranchAddress("sel_trueQ2", &sel_trueQ2, &b_sel_trueQ2);
   fChain->SetBranchAddress("sel_truedalphaT", &sel_truedalphaT, &b_sel_truedalphaT);
   fChain->SetBranchAddress("sel_truedpT", &sel_truedpT, &b_sel_truedpT);
   fChain->SetBranchAddress("sel_truedpTT", &sel_truedpTT, &b_sel_truedpTT);
   fChain->SetBranchAddress("sel_truedpTT_pi", &sel_truedpTT_pi, &b_sel_truedpTT_pi);
   fChain->SetBranchAddress("sel_truedpTT_pi_dir", &sel_truedpTT_pi_dir, &b_sel_truedpTT_pi_dir);
   fChain->SetBranchAddress("sel_truedpTT_pr", &sel_truedpTT_pr, &b_sel_truedpTT_pr);
   fChain->SetBranchAddress("sel_truedpTT_pr_dir", &sel_truedpTT_pr_dir, &b_sel_truedpTT_pr_dir);
   fChain->SetBranchAddress("sel_truedphiT", &sel_truedphiT, &b_sel_truedphiT);
   fChain->SetBranchAddress("n_iso_blobs", &n_iso_blobs, &b_n_iso_blobs);
   fChain->SetBranchAddress("sel_iso_blob_nclusters", &sel_iso_blob_nclusters, &b_sel_iso_blob_nclusters);
   fChain->SetBranchAddress("sel_VTX", sel_VTX, &b_sel_VTX);
   fChain->SetBranchAddress("sel_dpT_vec_EX", sel_dpT_vec_EX, &b_sel_dpT_vec_EX);
   fChain->SetBranchAddress("sel_dpT_vec_LL", sel_dpT_vec_LL, &b_sel_dpT_vec_LL);
   fChain->SetBranchAddress("sel_iso_blob_energy", &sel_iso_blob_energy, &b_sel_iso_blob_energy);
   fChain->SetBranchAddress("sel_meanPDP", sel_meanPDP, &b_sel_meanPDP);
   fChain->SetBranchAddress("sel_mu_4mom", sel_mu_4mom, &b_sel_mu_4mom);
   fChain->SetBranchAddress("sel_mu_endpos", sel_mu_endpos, &b_sel_mu_endpos);
   fChain->SetBranchAddress("sel_mu_pT", sel_mu_pT, &b_sel_mu_pT);
   fChain->SetBranchAddress("sel_mu_pT_tmumom", sel_mu_pT_tmumom, &b_sel_mu_pT_tmumom);
   fChain->SetBranchAddress("sel_mu_pT_tnudir", sel_mu_pT_tnudir, &b_sel_mu_pT_tnudir);
   fChain->SetBranchAddress("sel_mu_startdir", sel_mu_startdir, &b_sel_mu_startdir);
   fChain->SetBranchAddress("sel_mu_startpos", sel_mu_startpos, &b_sel_mu_startpos);
   fChain->SetBranchAddress("sel_mu_true4mom", sel_mu_true4mom, &b_sel_mu_true4mom);
   fChain->SetBranchAddress("sel_mu_trueendpos", sel_mu_trueendpos, &b_sel_mu_trueendpos);
   fChain->SetBranchAddress("sel_mu_truepT", sel_mu_truepT, &b_sel_mu_truepT);
   fChain->SetBranchAddress("sel_mu_truestartdir", sel_mu_truestartdir, &b_sel_mu_truestartdir);
   fChain->SetBranchAddress("sel_mu_truestartpos", sel_mu_truestartpos, &b_sel_mu_truestartpos);
   fChain->SetBranchAddress("sel_nu_dir_001", sel_nu_dir_001, &b_sel_nu_dir_001);
   fChain->SetBranchAddress("sel_nu_dir_PDP", sel_nu_dir_PDP, &b_sel_nu_dir_PDP);
   fChain->SetBranchAddress("sel_pi_EX_4mom", sel_pi_EX_4mom, &b_sel_pi_EX_4mom);
   fChain->SetBranchAddress("sel_pi_EX_pT", sel_pi_EX_pT, &b_sel_pi_EX_pT);
   fChain->SetBranchAddress("sel_pi_EX_pT_tnudir", sel_pi_EX_pT_tnudir, &b_sel_pi_EX_pT_tnudir);
   fChain->SetBranchAddress("sel_pi_EX_pT_tpimom", sel_pi_EX_pT_tpimom, &b_sel_pi_EX_pT_tpimom);
   fChain->SetBranchAddress("sel_pi_LL_4mom", sel_pi_LL_4mom, &b_sel_pi_LL_4mom);
   fChain->SetBranchAddress("sel_pi_LL_pT", sel_pi_LL_pT, &b_sel_pi_LL_pT);
   fChain->SetBranchAddress("sel_pi_LL_pT_tnudir", sel_pi_LL_pT_tnudir, &b_sel_pi_LL_pT_tnudir);
   fChain->SetBranchAddress("sel_pi_LL_pT_tpimom", sel_pi_LL_pT_tpimom, &b_sel_pi_LL_pT_tpimom);
   fChain->SetBranchAddress("sel_pi_endpos", sel_pi_endpos, &b_sel_pi_endpos);
   fChain->SetBranchAddress("sel_pi_startdir", sel_pi_startdir, &b_sel_pi_startdir);
   fChain->SetBranchAddress("sel_pi_startpos", sel_pi_startpos, &b_sel_pi_startpos);
   fChain->SetBranchAddress("sel_pi_true4mom", sel_pi_true4mom, &b_sel_pi_true4mom);
   fChain->SetBranchAddress("sel_pi_trueendpos", sel_pi_trueendpos, &b_sel_pi_trueendpos);
   fChain->SetBranchAddress("sel_pi_truepT", sel_pi_truepT, &b_sel_pi_truepT);
   fChain->SetBranchAddress("sel_pi_truestartdir", sel_pi_truestartdir, &b_sel_pi_truestartdir);
   fChain->SetBranchAddress("sel_pi_truestartpos", sel_pi_truestartpos, &b_sel_pi_truestartpos);
   fChain->SetBranchAddress("sel_pr_EX_4mom", sel_pr_EX_4mom, &b_sel_pr_EX_4mom);
   fChain->SetBranchAddress("sel_pr_EX_pT", sel_pr_EX_pT, &b_sel_pr_EX_pT);
   fChain->SetBranchAddress("sel_pr_EX_pT_tnudir", sel_pr_EX_pT_tnudir, &b_sel_pr_EX_pT_tnudir);
   fChain->SetBranchAddress("sel_pr_EX_pT_tprmom", sel_pr_EX_pT_tprmom, &b_sel_pr_EX_pT_tprmom);
   fChain->SetBranchAddress("sel_pr_LL_4mom", sel_pr_LL_4mom, &b_sel_pr_LL_4mom);
   fChain->SetBranchAddress("sel_pr_LL_pT", sel_pr_LL_pT, &b_sel_pr_LL_pT);
   fChain->SetBranchAddress("sel_pr_LL_pT_tnudir", sel_pr_LL_pT_tnudir, &b_sel_pr_LL_pT_tnudir);
   fChain->SetBranchAddress("sel_pr_LL_pT_tprmom", sel_pr_LL_pT_tprmom, &b_sel_pr_LL_pT_tprmom);
   fChain->SetBranchAddress("sel_pr_endpos", sel_pr_endpos, &b_sel_pr_endpos);
   fChain->SetBranchAddress("sel_pr_startdir", sel_pr_startdir, &b_sel_pr_startdir);
   fChain->SetBranchAddress("sel_pr_startpos", sel_pr_startpos, &b_sel_pr_startpos);
   fChain->SetBranchAddress("sel_pr_true4mom", sel_pr_true4mom, &b_sel_pr_true4mom);
   fChain->SetBranchAddress("sel_pr_trueendpos", sel_pr_trueendpos, &b_sel_pr_trueendpos);
   fChain->SetBranchAddress("sel_pr_truepT", sel_pr_truepT, &b_sel_pr_truepT);
   fChain->SetBranchAddress("sel_pr_truestartdir", sel_pr_truestartdir, &b_sel_pr_truestartdir);
   fChain->SetBranchAddress("sel_pr_truestartpos", sel_pr_truestartpos, &b_sel_pr_truestartpos);
   fChain->SetBranchAddress("sel_truePDP", sel_truePDP, &b_sel_truePDP);
   fChain->SetBranchAddress("sel_trueVTX", sel_trueVTX, &b_sel_trueVTX);
   fChain->SetBranchAddress("sel_true_nu_dir_PDP", sel_true_nu_dir_PDP, &b_sel_true_nu_dir_PDP);
   fChain->SetBranchAddress("sel_truedpT_vec", sel_truedpT_vec, &b_sel_truedpT_vec);
   fChain->SetBranchAddress("physEvtNum", &physEvtNum, &b_physEvtNum);
   fChain->SetBranchAddress("n_hyps", &n_hyps, &b_n_hyps);
   fChain->SetBranchAddress("processType", &processType, &b_processType);
   fChain->SetBranchAddress("primaryPart", &primaryPart, &b_primaryPart);
   fChain->SetBranchAddress("n_slices", &n_slices, &b_n_slices);
   fChain->SetBranchAddress("slice_numbers", slice_numbers, &b_slice_numbers);
   fChain->SetBranchAddress("shared_slice", &shared_slice, &b_shared_slice);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("vtxErr", vtxErr, &b_vtxErr);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("found_truth", &found_truth, &b_found_truth);
   fChain->SetBranchAddress("isMinosMatchTrack", &isMinosMatchTrack, &b_isMinosMatchTrack);
   fChain->SetBranchAddress("isMinosMatchStub", &isMinosMatchStub, &b_isMinosMatchStub);
   fChain->SetBranchAddress("contained_evt", &contained_evt, &b_contained_evt);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("n_anchored_long_trk_prongs", &n_anchored_long_trk_prongs, &b_n_anchored_long_trk_prongs);
   fChain->SetBranchAddress("n_anchored_short_trk_prongs", &n_anchored_short_trk_prongs, &b_n_anchored_short_trk_prongs);
   fChain->SetBranchAddress("n_iso_trk_prongs", &n_iso_trk_prongs, &b_n_iso_trk_prongs);
   fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
   fChain->SetBranchAddress("n_tracks3", &n_tracks3, &b_n_tracks3);
   fChain->SetBranchAddress("ncuts", &ncuts, &b_ncuts);
   fChain->SetBranchAddress("new_tracks", &new_tracks, &b_new_tracks);
   fChain->SetBranchAddress("nsplits", &nsplits, &b_nsplits);
   fChain->SetBranchAddress("target_region", &target_region, &b_target_region);
   fChain->SetBranchAddress("true_target_region", &true_target_region, &b_true_target_region);
   fChain->SetBranchAddress("vert_exists", &vert_exists, &b_vert_exists);
   fChain->SetBranchAddress("accum_level", accum_level, &b_accum_level);
   fChain->SetBranchAddress("truth_reco_isMinosMatch", &truth_reco_isMinosMatch, &b_truth_reco_isMinosMatch);
   fChain->SetBranchAddress("truth_muon_charge", &truth_muon_charge, &b_truth_muon_charge);
   fChain->SetBranchAddress("truth_n_ele", &truth_n_ele, &b_truth_n_ele);
   fChain->SetBranchAddress("truth_n_kPM", &truth_n_kPM, &b_truth_n_kPM);
   fChain->SetBranchAddress("truth_n_kaO", &truth_n_kaO, &b_truth_n_kaO);
   fChain->SetBranchAddress("truth_n_muo", &truth_n_muo, &b_truth_n_muo);
   fChain->SetBranchAddress("truth_n_ntn", &truth_n_ntn, &b_truth_n_ntn);
   fChain->SetBranchAddress("truth_n_pho", &truth_n_pho, &b_truth_n_pho);
   fChain->SetBranchAddress("truth_n_pi0", &truth_n_pi0, &b_truth_n_pi0);
   fChain->SetBranchAddress("truth_n_piM", &truth_n_piM, &b_truth_n_piM);
   fChain->SetBranchAddress("truth_n_piP", &truth_n_piP, &b_truth_n_piP);
   fChain->SetBranchAddress("truth_n_pro", &truth_n_pro, &b_truth_n_pro);
   fChain->SetBranchAddress("truth_n_tau", &truth_n_tau, &b_truth_n_tau);
   fChain->SetBranchAddress("truth_ncuts", &truth_ncuts, &b_truth_ncuts);
   fChain->SetBranchAddress("truth_nsplits", &truth_nsplits, &b_truth_nsplits);
   fChain->SetBranchAddress("truth_pi_EX_michel", &truth_pi_EX_michel, &b_truth_pi_EX_michel);
   fChain->SetBranchAddress("truth_pi_LL_michel", &truth_pi_LL_michel, &b_truth_pi_LL_michel);
   fChain->SetBranchAddress("truth_pr_EX_michel", &truth_pr_EX_michel, &b_truth_pr_EX_michel);
   fChain->SetBranchAddress("truth_pr_LL_michel", &truth_pr_LL_michel, &b_truth_pr_LL_michel);
   fChain->SetBranchAddress("truth_reco_target", &truth_reco_target, &b_truth_reco_target);
   fChain->SetBranchAddress("truth_should_be_accepted", &truth_should_be_accepted, &b_truth_should_be_accepted);
   fChain->SetBranchAddress("truth_true_target_region", &truth_true_target_region, &b_truth_true_target_region);
   fChain->SetBranchAddress("truth_mu_E", &truth_mu_E, &b_truth_mu_E);
   fChain->SetBranchAddress("truth_mu_KE", &truth_mu_KE, &b_truth_mu_KE);
   fChain->SetBranchAddress("truth_mu_Phi", &truth_mu_Phi, &b_truth_mu_Phi);
   fChain->SetBranchAddress("truth_mu_Theta", &truth_mu_Theta, &b_truth_mu_Theta);
   fChain->SetBranchAddress("truth_mu_mom", &truth_mu_mom, &b_truth_mu_mom);
   fChain->SetBranchAddress("truth_mu_pTMag", &truth_mu_pTMag, &b_truth_mu_pTMag);
   fChain->SetBranchAddress("truth_mu_pTT", &truth_mu_pTT, &b_truth_mu_pTT);
   fChain->SetBranchAddress("truth_pi_E", &truth_pi_E, &b_truth_pi_E);
   fChain->SetBranchAddress("truth_pi_EX_score", &truth_pi_EX_score, &b_truth_pi_EX_score);
   fChain->SetBranchAddress("truth_pi_EX_score_altH", &truth_pi_EX_score_altH, &b_truth_pi_EX_score_altH);
   fChain->SetBranchAddress("truth_pi_KE", &truth_pi_KE, &b_truth_pi_KE);
   fChain->SetBranchAddress("truth_pi_LL_score", &truth_pi_LL_score, &b_truth_pi_LL_score);
   fChain->SetBranchAddress("truth_pi_LL_score_altH", &truth_pi_LL_score_altH, &b_truth_pi_LL_score_altH);
   fChain->SetBranchAddress("truth_pi_Phi", &truth_pi_Phi, &b_truth_pi_Phi);
   fChain->SetBranchAddress("truth_pi_Theta", &truth_pi_Theta, &b_truth_pi_Theta);
   fChain->SetBranchAddress("truth_pi_mom", &truth_pi_mom, &b_truth_pi_mom);
   fChain->SetBranchAddress("truth_pi_pTMag", &truth_pi_pTMag, &b_truth_pi_pTMag);
   fChain->SetBranchAddress("truth_pi_pTT", &truth_pi_pTT, &b_truth_pi_pTT);
   fChain->SetBranchAddress("truth_pr_E", &truth_pr_E, &b_truth_pr_E);
   fChain->SetBranchAddress("truth_pr_EX_score", &truth_pr_EX_score, &b_truth_pr_EX_score);
   fChain->SetBranchAddress("truth_pr_EX_score_altH", &truth_pr_EX_score_altH, &b_truth_pr_EX_score_altH);
   fChain->SetBranchAddress("truth_pr_KE", &truth_pr_KE, &b_truth_pr_KE);
   fChain->SetBranchAddress("truth_pr_LL_score", &truth_pr_LL_score, &b_truth_pr_LL_score);
   fChain->SetBranchAddress("truth_pr_LL_score_altH", &truth_pr_LL_score_altH, &b_truth_pr_LL_score_altH);
   fChain->SetBranchAddress("truth_pr_Phi", &truth_pr_Phi, &b_truth_pr_Phi);
   fChain->SetBranchAddress("truth_pr_Theta", &truth_pr_Theta, &b_truth_pr_Theta);
   fChain->SetBranchAddress("truth_pr_mom", &truth_pr_mom, &b_truth_pr_mom);
   fChain->SetBranchAddress("truth_pr_pTMag", &truth_pr_pTMag, &b_truth_pr_pTMag);
   fChain->SetBranchAddress("truth_pr_pTT", &truth_pr_pTT, &b_truth_pr_pTT);
   fChain->SetBranchAddress("truth_trueEnu", &truth_trueEnu, &b_truth_trueEnu);
   fChain->SetBranchAddress("truth_trueQ2", &truth_trueQ2, &b_truth_trueQ2);
   fChain->SetBranchAddress("truth_truedalphaT", &truth_truedalphaT, &b_truth_truedalphaT);
   fChain->SetBranchAddress("truth_truedpT", &truth_truedpT, &b_truth_truedpT);
   fChain->SetBranchAddress("truth_truedpTT", &truth_truedpTT, &b_truth_truedpTT);
   fChain->SetBranchAddress("truth_truedpTT_pi", &truth_truedpTT_pi, &b_truth_truedpTT_pi);
   fChain->SetBranchAddress("truth_truedpTT_pi_dir", &truth_truedpTT_pi_dir, &b_truth_truedpTT_pi_dir);
   fChain->SetBranchAddress("truth_truedpTT_pr", &truth_truedpTT_pr, &b_truth_truedpTT_pr);
   fChain->SetBranchAddress("truth_truedpTT_pr_dir", &truth_truedpTT_pr_dir, &b_truth_truedpTT_pr_dir);
   fChain->SetBranchAddress("truth_truedphiT", &truth_truedphiT, &b_truth_truedphiT);
   fChain->SetBranchAddress("truth_accum_level", truth_accum_level, &b_truth_accum_level);
   fChain->SetBranchAddress("truth_mu_4mom", truth_mu_4mom, &b_truth_mu_4mom);
   fChain->SetBranchAddress("truth_mu_pT", truth_mu_pT, &b_truth_mu_pT);
   fChain->SetBranchAddress("truth_pi_4mom", truth_pi_4mom, &b_truth_pi_4mom);
   fChain->SetBranchAddress("truth_pi_pT", truth_pi_pT, &b_truth_pi_pT);
   fChain->SetBranchAddress("truth_pr_4mom", truth_pr_4mom, &b_truth_pr_4mom);
   fChain->SetBranchAddress("truth_pr_pT", truth_pr_pT, &b_truth_pr_pT);
   fChain->SetBranchAddress("truth_truedpT_vec", truth_truedpT_vec, &b_truth_truedpT_vec);
   fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   fChain->SetBranchAddress("ev_subrun", &ev_subrun, &b_ev_subrun);
   fChain->SetBranchAddress("ev_detector", &ev_detector, &b_ev_detector);
   fChain->SetBranchAddress("ev_triggerType", &ev_triggerType, &b_ev_triggerType);
   fChain->SetBranchAddress("ev_gate", &ev_gate, &b_ev_gate);
   fChain->SetBranchAddress("ev_global_gate", &ev_global_gate, &b_ev_global_gate);
   fChain->SetBranchAddress("ev_gps_time_sec", &ev_gps_time_sec, &b_ev_gps_time_sec);
   fChain->SetBranchAddress("ev_gps_time_usec", &ev_gps_time_usec, &b_ev_gps_time_usec);
   fChain->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
   fChain->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
   fChain->SetBranchAddress("mc_nInteractions", &mc_nInteractions, &b_mc_nInteractions);
   fChain->SetBranchAddress("mc_MIState", &mc_MIState, &b_mc_MIState);
   fChain->SetBranchAddress("mc_pot", &mc_pot, &b_mc_pot);
   fChain->SetBranchAddress("mc_beamConfig", &mc_beamConfig, &b_mc_beamConfig);
   fChain->SetBranchAddress("mc_processType", &mc_processType, &b_mc_processType);
   fChain->SetBranchAddress("mc_nthEvtInSpill", &mc_nthEvtInSpill, &b_mc_nthEvtInSpill);
   fChain->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile, &b_mc_nthEvtInFile);
   fChain->SetBranchAddress("mc_intType", &mc_intType, &b_mc_intType);
   fChain->SetBranchAddress("mc_current", &mc_current, &b_mc_current);
   fChain->SetBranchAddress("mc_charm", &mc_charm, &b_mc_charm);
   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   fChain->SetBranchAddress("mc_XSec", &mc_XSec, &b_mc_XSec);
   fChain->SetBranchAddress("mc_diffXSec", &mc_diffXSec, &b_mc_diffXSec);
   fChain->SetBranchAddress("mc_incoming", &mc_incoming, &b_mc_incoming);
   fChain->SetBranchAddress("mc_fluxDriverProb", &mc_fluxDriverProb, &b_mc_fluxDriverProb);
   fChain->SetBranchAddress("mc_targetNucleus", &mc_targetNucleus, &b_mc_targetNucleus);
   fChain->SetBranchAddress("mc_targetZ", &mc_targetZ, &b_mc_targetZ);
   fChain->SetBranchAddress("mc_targetA", &mc_targetA, &b_mc_targetA);
   fChain->SetBranchAddress("mc_targetNucleon", &mc_targetNucleon, &b_mc_targetNucleon);
   fChain->SetBranchAddress("mc_struckQuark", &mc_struckQuark, &b_mc_struckQuark);
   fChain->SetBranchAddress("mc_seaQuark", &mc_seaQuark, &b_mc_seaQuark);
   fChain->SetBranchAddress("mc_resID", &mc_resID, &b_mc_resID);
   fChain->SetBranchAddress("mc_primaryLepton", &mc_primaryLepton, &b_mc_primaryLepton);
   fChain->SetBranchAddress("mc_incomingE", &mc_incomingE, &b_mc_incomingE);
   fChain->SetBranchAddress("mc_Bjorkenx", &mc_Bjorkenx, &b_mc_Bjorkenx);
   fChain->SetBranchAddress("mc_Bjorkeny", &mc_Bjorkeny, &b_mc_Bjorkeny);
   fChain->SetBranchAddress("mc_Q2", &mc_Q2, &b_mc_Q2);
   fChain->SetBranchAddress("mc_nuT", &mc_nuT, &b_mc_nuT);
   fChain->SetBranchAddress("mc_w", &mc_w, &b_mc_w);
   fChain->SetBranchAddress("mc_vtx", mc_vtx, &b_mc_vtx);
   fChain->SetBranchAddress("mc_incomingPartVec", mc_incomingPartVec, &b_mc_incomingPartVec);
   fChain->SetBranchAddress("mc_initNucVec", mc_initNucVec, &b_mc_initNucVec);
   fChain->SetBranchAddress("mc_primFSLepton", mc_primFSLepton, &b_mc_primFSLepton);
   fChain->SetBranchAddress("mc_nFSPart", &mc_nFSPart, &b_mc_nFSPart);
   fChain->SetBranchAddress("mc_FSPartPx", mc_FSPartPx, &b_mc_FSPartPx);
   fChain->SetBranchAddress("mc_FSPartPy", mc_FSPartPy, &b_mc_FSPartPy);
   fChain->SetBranchAddress("mc_FSPartPz", mc_FSPartPz, &b_mc_FSPartPz);
   fChain->SetBranchAddress("mc_FSPartE", mc_FSPartE, &b_mc_FSPartE);
   fChain->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG, &b_mc_FSPartPDG);
   fChain->SetBranchAddress("mc_er_nPart", &mc_er_nPart, &b_mc_er_nPart);
   fChain->SetBranchAddress("mc_er_ID", mc_er_ID, &b_mc_er_ID);
   fChain->SetBranchAddress("mc_er_status", mc_er_status, &b_mc_er_status);
   fChain->SetBranchAddress("mc_er_posInNucX", mc_er_posInNucX, &b_mc_er_posInNucX);
   fChain->SetBranchAddress("mc_er_posInNucY", mc_er_posInNucY, &b_mc_er_posInNucY);
   fChain->SetBranchAddress("mc_er_posInNucZ", mc_er_posInNucZ, &b_mc_er_posInNucZ);
   fChain->SetBranchAddress("mc_er_Px", mc_er_Px, &b_mc_er_Px);
   fChain->SetBranchAddress("mc_er_Py", mc_er_Py, &b_mc_er_Py);
   fChain->SetBranchAddress("mc_er_Pz", mc_er_Pz, &b_mc_er_Pz);
   fChain->SetBranchAddress("mc_er_E", mc_er_E, &b_mc_er_E);
   fChain->SetBranchAddress("mc_er_FD", mc_er_FD, &b_mc_er_FD);
   fChain->SetBranchAddress("mc_er_LD", mc_er_LD, &b_mc_er_LD);
   fChain->SetBranchAddress("mc_er_mother", mc_er_mother, &b_mc_er_mother);
   fChain->SetBranchAddress("mc_fr_nNuAncestorIDs", &mc_fr_nNuAncestorIDs, &b_mc_fr_nNuAncestorIDs);
   fChain->SetBranchAddress("mc_fr_nuAncestorIDs", mc_fr_nuAncestorIDs, &b_mc_fr_nuAncestorIDs);
   fChain->SetBranchAddress("mc_fr_nuParentID", &mc_fr_nuParentID, &b_mc_fr_nuParentID);
   fChain->SetBranchAddress("mc_fr_decMode", &mc_fr_decMode, &b_mc_fr_decMode);
   fChain->SetBranchAddress("mc_fr_primProtonVtx", mc_fr_primProtonVtx, &b_mc_fr_primProtonVtx);
   fChain->SetBranchAddress("mc_fr_primProtonP", mc_fr_primProtonP, &b_mc_fr_primProtonP);
   fChain->SetBranchAddress("mc_fr_nuParentDecVtx", mc_fr_nuParentDecVtx, &b_mc_fr_nuParentDecVtx);
   fChain->SetBranchAddress("mc_fr_nuParentProdVtx", mc_fr_nuParentProdVtx, &b_mc_fr_nuParentProdVtx);
   fChain->SetBranchAddress("mc_fr_nuParentProdP", mc_fr_nuParentProdP, &b_mc_fr_nuParentProdP);
   fChain->SetBranchAddress("mc_cvweight_total", &mc_cvweight_total, &b_mc_cvweight_total);
   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("mc_cvweight_totalFlux", &mc_cvweight_totalFlux, &b_mc_cvweight_totalFlux);
   fChain->SetBranchAddress("mc_cvweight_totalXsec", &mc_cvweight_totalXsec, &b_mc_cvweight_totalXsec);
   fChain->SetBranchAddress("mc_ppfx1_cvweight", &mc_ppfx1_cvweight, &b_mc_ppfx1_cvweight);
   fChain->SetBranchAddress("mc_hornCurrent_cvweight", &mc_hornCurrent_cvweight, &b_mc_hornCurrent_cvweight);
   fChain->SetBranchAddress("mc_gen1_cvweight_total", &mc_gen1_cvweight_total, &b_mc_gen1_cvweight_total);
   fChain->SetBranchAddress("gen1_wgt", &gen1_wgt, &b_gen1_wgt);
   fChain->SetBranchAddress("mc_gen1_cvweight_totalFlux", &mc_gen1_cvweight_totalFlux, &b_mc_gen1_cvweight_totalFlux);
   fChain->SetBranchAddress("mc_gen1_cvweight_NA49", &mc_gen1_cvweight_NA49, &b_mc_gen1_cvweight_NA49);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", &mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary_sz", &mc_wgt_gen1_Flux_Tertiary_sz, &b_mc_wgt_gen1_Flux_Tertiary_sz);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary", &mc_wgt_gen1_Flux_Tertiary, &b_mc_wgt_gen1_Flux_Tertiary);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49_sz", &mc_wgt_gen1_Flux_NA49_sz, &b_mc_wgt_gen1_Flux_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49", &mc_wgt_gen1_Flux_NA49, &b_mc_wgt_gen1_Flux_NA49);
   fChain->SetBranchAddress("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, &b_mc_wgt_Norm_sz);
   fChain->SetBranchAddress("mc_wgt_Norm", &mc_wgt_Norm, &b_mc_wgt_Norm);
   fChain->SetBranchAddress("mc_wgt_ppfx1_Total_sz", &mc_wgt_ppfx1_Total_sz, &b_mc_wgt_ppfx1_Total_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx1_Total", &mc_wgt_ppfx1_Total, &b_mc_wgt_ppfx1_Total);
   fChain->SetBranchAddress("prong_nParticles", prong_nParticles, &b_prong_nParticles);
   fChain->SetBranchAddress("prong_part_score", prong_part_score, &b_prong_part_score);
   fChain->SetBranchAddress("prong_part_mass", prong_part_mass, &b_prong_part_mass);
   fChain->SetBranchAddress("prong_part_charge", prong_part_charge, &b_prong_part_charge);
   fChain->SetBranchAddress("prong_part_pid", prong_part_pid, &b_prong_part_pid);
   // fChain->SetBranchAddress("prong_part_E", &prong_part_E, &b_prong_part_E);
   // fChain->SetBranchAddress("prong_part_pos", &prong_part_pos, &b_prong_part_pos);
}

void AnalysisReader::SetOutTree(){

   m_savetree->Branch("eventID", &eventID, "eventID/D");
   m_savetree->Branch("sel_nuFlavor", &sel_nuFlavor, "sel_nuFlavor/I");
   m_savetree->Branch("sel_nuHelicity", &sel_nuHelicity, "sel_nuHelicity/I");
   m_savetree->Branch("sel_intCurrent", &sel_intCurrent, "sel_intCurrent/I");
   m_savetree->Branch("sel_intType", &sel_intType, "sel_intType/I");
   m_savetree->Branch("sel_E", &sel_E, "sel_E/D");
   m_savetree->Branch("sel_Q2", &sel_Q2, "sel_Q2/D");
   m_savetree->Branch("sel_x", &sel_x, "sel_x/D");
   m_savetree->Branch("sel_y", &sel_y, "sel_y/D");
   m_savetree->Branch("sel_W", &sel_W, "sel_W/D");
   m_savetree->Branch("sel_score", &sel_score, "sel_score/D");
   m_savetree->Branch("sel_leptonE", sel_leptonE, "sel_leptonE[4]/D");
   m_savetree->Branch("sel_vtx", sel_vtx, "sel_vtx[3]/D");
   m_savetree->Branch("sel_mu_PDG", &sel_mu_PDG, "sel_mu_PDG/I");
   m_savetree->Branch("sel_mu_isKinked", &sel_mu_isKinked, "sel_mu_isKinked/I");
   m_savetree->Branch("sel_mu_michel", &sel_mu_michel, "sel_mu_michel/I");
   m_savetree->Branch("sel_pi_EX_michel", &sel_pi_EX_michel, "sel_pi_EX_michel/I");
   m_savetree->Branch("sel_pi_EX_nNodes", &sel_pi_EX_nNodes, "sel_pi_EX_nNodes/I");
   m_savetree->Branch("sel_pi_FSI", &sel_pi_FSI, "sel_pi_FSI/I");
   m_savetree->Branch("sel_pi_LL_michel", &sel_pi_LL_michel, "sel_pi_LL_michel/I");
   m_savetree->Branch("sel_pi_LL_nNodes", &sel_pi_LL_nNodes, "sel_pi_LL_nNodes/I");
   m_savetree->Branch("sel_pi_PDG", &sel_pi_PDG, "sel_pi_PDG/I");
   m_savetree->Branch("sel_pi_isKinked", &sel_pi_isKinked, "sel_pi_isKinked/I");
   m_savetree->Branch("sel_pr_EX_michel", &sel_pr_EX_michel, "sel_pr_EX_michel/I");
   m_savetree->Branch("sel_pr_EX_nNodes", &sel_pr_EX_nNodes, "sel_pr_EX_nNodes/I");
   m_savetree->Branch("sel_pr_FSI", &sel_pr_FSI, "sel_pr_FSI/I");
   m_savetree->Branch("sel_pr_LL_michel", &sel_pr_LL_michel, "sel_pr_LL_michel/I");
   m_savetree->Branch("sel_pr_LL_nNodes", &sel_pr_LL_nNodes, "sel_pr_LL_nNodes/I");
   m_savetree->Branch("sel_pr_PDG", &sel_pr_PDG, "sel_pr_PDG/I");
   m_savetree->Branch("sel_pr_isKinked", &sel_pr_isKinked, "sel_pr_isKinked/I");
   m_savetree->Branch("sel_Enu_EX", &sel_Enu_EX, "sel_Enu_EX/D");
   m_savetree->Branch("sel_Enu_LL", &sel_Enu_LL, "sel_Enu_LL/D");
   m_savetree->Branch("sel_Q2_EX", &sel_Q2_EX, "sel_Q2_EX/D");
   m_savetree->Branch("sel_Q2_LL", &sel_Q2_LL, "sel_Q2_LL/D");
   m_savetree->Branch("sel_dalphaT_EX", &sel_dalphaT_EX, "sel_dalphaT_EX/D");
   m_savetree->Branch("sel_dalphaT_LL", &sel_dalphaT_LL, "sel_dalphaT_LL/D");
   m_savetree->Branch("sel_dpTT_EX", &sel_dpTT_EX, "sel_dpTT_EX/D");
   m_savetree->Branch("sel_dpTT_EX_tmumom", &sel_dpTT_EX_tmumom, "sel_dpTT_EX_tmumom/D");
   m_savetree->Branch("sel_dpTT_EX_tnudir", &sel_dpTT_EX_tnudir, "sel_dpTT_EX_tnudir/D");
   m_savetree->Branch("sel_dpTT_EX_tpimom", &sel_dpTT_EX_tpimom, "sel_dpTT_EX_tpimom/D");
   m_savetree->Branch("sel_dpTT_EX_tprmom", &sel_dpTT_EX_tprmom, "sel_dpTT_EX_tprmom/D");
   m_savetree->Branch("sel_dpTT_LL", &sel_dpTT_LL, "sel_dpTT_LL/D");
   m_savetree->Branch("sel_dpTT_LL_tmumom", &sel_dpTT_LL_tmumom, "sel_dpTT_LL_tmumom/D");
   m_savetree->Branch("sel_dpTT_LL_tnudir", &sel_dpTT_LL_tnudir, "sel_dpTT_LL_tnudir/D");
   m_savetree->Branch("sel_dpTT_LL_tpimom", &sel_dpTT_LL_tpimom, "sel_dpTT_LL_tpimom/D");
   m_savetree->Branch("sel_dpTT_LL_tprmom", &sel_dpTT_LL_tprmom, "sel_dpTT_LL_tprmom/D");
   m_savetree->Branch("sel_dpTT_pi_EX", &sel_dpTT_pi_EX, "sel_dpTT_pi_EX/D");
   m_savetree->Branch("sel_dpTT_pi_EX_tmumom", &sel_dpTT_pi_EX_tmumom, "sel_dpTT_pi_EX_tmumom/D");
   m_savetree->Branch("sel_dpTT_pi_EX_tnudir", &sel_dpTT_pi_EX_tnudir, "sel_dpTT_pi_EX_tnudir/D");
   m_savetree->Branch("sel_dpTT_pi_EX_tpimom", &sel_dpTT_pi_EX_tpimom, "sel_dpTT_pi_EX_tpimom/D");
   m_savetree->Branch("sel_dpTT_pi_EX_tprmom", &sel_dpTT_pi_EX_tprmom, "sel_dpTT_pi_EX_tprmom/D");
   m_savetree->Branch("sel_dpTT_pi_LL", &sel_dpTT_pi_LL, "sel_dpTT_pi_LL/D");
   m_savetree->Branch("sel_dpTT_pi_LL_tmumom", &sel_dpTT_pi_LL_tmumom, "sel_dpTT_pi_LL_tmumom/D");
   m_savetree->Branch("sel_dpTT_pi_LL_tnudir", &sel_dpTT_pi_LL_tnudir, "sel_dpTT_pi_LL_tnudir/D");
   m_savetree->Branch("sel_dpTT_pi_LL_tpimom", &sel_dpTT_pi_LL_tpimom, "sel_dpTT_pi_LL_tpimom/D");
   m_savetree->Branch("sel_dpTT_pi_LL_tprmom", &sel_dpTT_pi_LL_tprmom, "sel_dpTT_pi_LL_tprmom/D");
   m_savetree->Branch("sel_dpTT_pi_dir_EX", &sel_dpTT_pi_dir_EX, "sel_dpTT_pi_dir_EX/D");
   m_savetree->Branch("sel_dpTT_pi_dir_EX_tmumom", &sel_dpTT_pi_dir_EX_tmumom, "sel_dpTT_pi_dir_EX_tmumom/D");
   m_savetree->Branch("sel_dpTT_pi_dir_EX_tnudir", &sel_dpTT_pi_dir_EX_tnudir, "sel_dpTT_pi_dir_EX_tnudir/D");
   m_savetree->Branch("sel_dpTT_pi_dir_EX_tpidir", &sel_dpTT_pi_dir_EX_tpidir, "sel_dpTT_pi_dir_EX_tpidir/D");
   m_savetree->Branch("sel_dpTT_pi_dir_EX_tprmom", &sel_dpTT_pi_dir_EX_tprmom, "sel_dpTT_pi_dir_EX_tprmom/D");
   m_savetree->Branch("sel_dpTT_pi_dir_LL", &sel_dpTT_pi_dir_LL, "sel_dpTT_pi_dir_LL/D");
   m_savetree->Branch("sel_dpTT_pi_dir_LL_tmumom", &sel_dpTT_pi_dir_LL_tmumom, "sel_dpTT_pi_dir_LL_tmumom/D");
   m_savetree->Branch("sel_dpTT_pi_dir_LL_tnudir", &sel_dpTT_pi_dir_LL_tnudir, "sel_dpTT_pi_dir_LL_tnudir/D");
   m_savetree->Branch("sel_dpTT_pi_dir_LL_tpidir", &sel_dpTT_pi_dir_LL_tpidir, "sel_dpTT_pi_dir_LL_tpidir/D");
   m_savetree->Branch("sel_dpTT_pi_dir_LL_tprmom", &sel_dpTT_pi_dir_LL_tprmom, "sel_dpTT_pi_dir_LL_tprmom/D");
   m_savetree->Branch("sel_dpTT_pr_EX", &sel_dpTT_pr_EX, "sel_dpTT_pr_EX/D");
   m_savetree->Branch("sel_dpTT_pr_EX_tmumom", &sel_dpTT_pr_EX_tmumom, "sel_dpTT_pr_EX_tmumom/D");
   m_savetree->Branch("sel_dpTT_pr_EX_tnudir", &sel_dpTT_pr_EX_tnudir, "sel_dpTT_pr_EX_tnudir/D");
   m_savetree->Branch("sel_dpTT_pr_EX_tpimon", &sel_dpTT_pr_EX_tpimon, "sel_dpTT_pr_EX_tpimon/D");
   m_savetree->Branch("sel_dpTT_pr_EX_tprmom", &sel_dpTT_pr_EX_tprmom, "sel_dpTT_pr_EX_tprmom/D");
   m_savetree->Branch("sel_dpTT_pr_LL", &sel_dpTT_pr_LL, "sel_dpTT_pr_LL/D");
   m_savetree->Branch("sel_dpTT_pr_LL_tmumom", &sel_dpTT_pr_LL_tmumom, "sel_dpTT_pr_LL_tmumom/D");
   m_savetree->Branch("sel_dpTT_pr_LL_tnudir", &sel_dpTT_pr_LL_tnudir, "sel_dpTT_pr_LL_tnudir/D");
   m_savetree->Branch("sel_dpTT_pr_LL_tpimon", &sel_dpTT_pr_LL_tpimon, "sel_dpTT_pr_LL_tpimon/D");
   m_savetree->Branch("sel_dpTT_pr_LL_tprmom", &sel_dpTT_pr_LL_tprmom, "sel_dpTT_pr_LL_tprmom/D");
   m_savetree->Branch("sel_dpTT_pr_dir_EX", &sel_dpTT_pr_dir_EX, "sel_dpTT_pr_dir_EX/D");
   m_savetree->Branch("sel_dpTT_pr_dir_EX_tmumom", &sel_dpTT_pr_dir_EX_tmumom, "sel_dpTT_pr_dir_EX_tmumom/D");
   m_savetree->Branch("sel_dpTT_pr_dir_EX_tnudir", &sel_dpTT_pr_dir_EX_tnudir, "sel_dpTT_pr_dir_EX_tnudir/D");
   m_savetree->Branch("sel_dpTT_pr_dir_EX_tpimom", &sel_dpTT_pr_dir_EX_tpimom, "sel_dpTT_pr_dir_EX_tpimom/D");
   m_savetree->Branch("sel_dpTT_pr_dir_EX_tprdir", &sel_dpTT_pr_dir_EX_tprdir, "sel_dpTT_pr_dir_EX_tprdir/D");
   m_savetree->Branch("sel_dpTT_pr_dir_LL", &sel_dpTT_pr_dir_LL, "sel_dpTT_pr_dir_LL/D");
   m_savetree->Branch("sel_dpTT_pr_dir_LL_tmumom", &sel_dpTT_pr_dir_LL_tmumom, "sel_dpTT_pr_dir_LL_tmumom/D");
   m_savetree->Branch("sel_dpTT_pr_dir_LL_tnudir", &sel_dpTT_pr_dir_LL_tnudir, "sel_dpTT_pr_dir_LL_tnudir/D");
   m_savetree->Branch("sel_dpTT_pr_dir_LL_tpimom", &sel_dpTT_pr_dir_LL_tpimom, "sel_dpTT_pr_dir_LL_tpimom/D");
   m_savetree->Branch("sel_dpTT_pr_dir_LL_tprdir", &sel_dpTT_pr_dir_LL_tprdir, "sel_dpTT_pr_dir_LL_tprdir/D");
   m_savetree->Branch("sel_dpT_EX", &sel_dpT_EX, "sel_dpT_EX/D");
   m_savetree->Branch("sel_dpT_LL", &sel_dpT_LL, "sel_dpT_LL/D");
   m_savetree->Branch("sel_dphiT_EX", &sel_dphiT_EX, "sel_dphiT_EX/D");
   m_savetree->Branch("sel_dphiT_LL", &sel_dphiT_LL, "sel_dphiT_LL/D");
   m_savetree->Branch("sel_mu_E", &sel_mu_E, "sel_mu_E/D");
   m_savetree->Branch("sel_mu_KE", &sel_mu_KE, "sel_mu_KE/D");
   m_savetree->Branch("sel_mu_Phi", &sel_mu_Phi, "sel_mu_Phi/D");
   m_savetree->Branch("sel_mu_Theta", &sel_mu_Theta, "sel_mu_Theta/D");
   m_savetree->Branch("sel_mu_XTheta", &sel_mu_XTheta, "sel_mu_XTheta/D");
   m_savetree->Branch("sel_mu_YTheta", &sel_mu_YTheta, "sel_mu_YTheta/D");
   m_savetree->Branch("sel_mu_chi2ndf", &sel_mu_chi2ndf, "sel_mu_chi2ndf/D");
   m_savetree->Branch("sel_mu_det_frac", &sel_mu_det_frac, "sel_mu_det_frac/D");
   m_savetree->Branch("sel_mu_det_otherE", &sel_mu_det_otherE, "sel_mu_det_otherE/D");
   m_savetree->Branch("sel_mu_mom", &sel_mu_mom, "sel_mu_mom/D");
   m_savetree->Branch("sel_mu_pTMag", &sel_mu_pTMag, "sel_mu_pTMag/D");
   m_savetree->Branch("sel_mu_pTMag_tmumom", &sel_mu_pTMag_tmumom, "sel_mu_pTMag_tmumom/D");
   m_savetree->Branch("sel_mu_pTMag_tnudir", &sel_mu_pTMag_tnudir, "sel_mu_pTMag_tnudir/D");
   m_savetree->Branch("sel_mu_pTT", &sel_mu_pTT, "sel_mu_pTT/D");
   m_savetree->Branch("sel_mu_pTT_tmumom", &sel_mu_pTT_tmumom, "sel_mu_pTT_tmumom/D");
   m_savetree->Branch("sel_mu_pTT_tnudir", &sel_mu_pTT_tnudir, "sel_mu_pTT_tnudir/D");
   m_savetree->Branch("sel_mu_score", &sel_mu_score, "sel_mu_score/D");
   m_savetree->Branch("sel_mu_trueE", &sel_mu_trueE, "sel_mu_trueE/D");
   m_savetree->Branch("sel_mu_trueKE", &sel_mu_trueKE, "sel_mu_trueKE/D");
   m_savetree->Branch("sel_mu_truePhi", &sel_mu_truePhi, "sel_mu_truePhi/D");
   m_savetree->Branch("sel_mu_trueTheta", &sel_mu_trueTheta, "sel_mu_trueTheta/D");
   m_savetree->Branch("sel_mu_truemom", &sel_mu_truemom, "sel_mu_truemom/D");
   m_savetree->Branch("sel_mu_truepTMag", &sel_mu_truepTMag, "sel_mu_truepTMag/D");
   m_savetree->Branch("sel_mu_truepTT", &sel_mu_truepTT, "sel_mu_truepTT/D");
   m_savetree->Branch("sel_pi_CaloE", &sel_pi_CaloE, "sel_pi_CaloE/D");
   m_savetree->Branch("sel_pi_EX_E", &sel_pi_EX_E, "sel_pi_EX_E/D");
   m_savetree->Branch("sel_pi_EX_E_altH", &sel_pi_EX_E_altH, "sel_pi_EX_E_altH/D");
   m_savetree->Branch("sel_pi_EX_KE", &sel_pi_EX_KE, "sel_pi_EX_KE/D");
   m_savetree->Branch("sel_pi_EX_Phi", &sel_pi_EX_Phi, "sel_pi_EX_Phi/D");
   m_savetree->Branch("sel_pi_EX_Theta", &sel_pi_EX_Theta, "sel_pi_EX_Theta/D");
   m_savetree->Branch("sel_pi_EX_XTheta", &sel_pi_EX_XTheta, "sel_pi_EX_XTheta/D");
   m_savetree->Branch("sel_pi_EX_YTheta", &sel_pi_EX_YTheta, "sel_pi_EX_YTheta/D");
   m_savetree->Branch("sel_pi_EX_ln_Q0", &sel_pi_EX_ln_Q0, "sel_pi_EX_ln_Q0/D");
   m_savetree->Branch("sel_pi_EX_ln_Q1", &sel_pi_EX_ln_Q1, "sel_pi_EX_ln_Q1/D");
   m_savetree->Branch("sel_pi_EX_ln_Q2", &sel_pi_EX_ln_Q2, "sel_pi_EX_ln_Q2/D");
   m_savetree->Branch("sel_pi_EX_ln_Q3", &sel_pi_EX_ln_Q3, "sel_pi_EX_ln_Q3/D");
   m_savetree->Branch("sel_pi_EX_ln_Q4", &sel_pi_EX_ln_Q4, "sel_pi_EX_ln_Q4/D");
   m_savetree->Branch("sel_pi_EX_ln_Q5", &sel_pi_EX_ln_Q5, "sel_pi_EX_ln_Q5/D");
   m_savetree->Branch("sel_pi_EX_mom", &sel_pi_EX_mom, "sel_pi_EX_mom/D");
   m_savetree->Branch("sel_pi_EX_mom_altH", &sel_pi_EX_mom_altH, "sel_pi_EX_mom_altH/D");
   m_savetree->Branch("sel_pi_EX_pTMag", &sel_pi_EX_pTMag, "sel_pi_EX_pTMag/D");
   m_savetree->Branch("sel_pi_EX_pTMag_tnudir", &sel_pi_EX_pTMag_tnudir, "sel_pi_EX_pTMag_tnudir/D");
   m_savetree->Branch("sel_pi_EX_pTMag_tpimom", &sel_pi_EX_pTMag_tpimom, "sel_pi_EX_pTMag_tpimom/D");
   m_savetree->Branch("sel_pi_EX_pTT", &sel_pi_EX_pTT, "sel_pi_EX_pTT/D");
   m_savetree->Branch("sel_pi_EX_pTT_tnudir", &sel_pi_EX_pTT_tnudir, "sel_pi_EX_pTT_tnudir/D");
   m_savetree->Branch("sel_pi_EX_pTT_tpimom", &sel_pi_EX_pTT_tpimom, "sel_pi_EX_pTT_tpimom/D");
   m_savetree->Branch("sel_pi_EX_score", &sel_pi_EX_score, "sel_pi_EX_score/D");
   m_savetree->Branch("sel_pi_EX_score_altH", &sel_pi_EX_score_altH, "sel_pi_EX_score_altH/D");
   m_savetree->Branch("sel_pi_LL_E", &sel_pi_LL_E, "sel_pi_LL_E/D");
   m_savetree->Branch("sel_pi_LL_E_altH", &sel_pi_LL_E_altH, "sel_pi_LL_E_altH/D");
   m_savetree->Branch("sel_pi_LL_KE", &sel_pi_LL_KE, "sel_pi_LL_KE/D");
   m_savetree->Branch("sel_pi_LL_Phi", &sel_pi_LL_Phi, "sel_pi_LL_Phi/D");
   m_savetree->Branch("sel_pi_LL_Theta", &sel_pi_LL_Theta, "sel_pi_LL_Theta/D");
   m_savetree->Branch("sel_pi_LL_XTheta", &sel_pi_LL_XTheta, "sel_pi_LL_XTheta/D");
   m_savetree->Branch("sel_pi_LL_YTheta", &sel_pi_LL_YTheta, "sel_pi_LL_YTheta/D");
   m_savetree->Branch("sel_pi_LL_ln_Q0", &sel_pi_LL_ln_Q0, "sel_pi_LL_ln_Q0/D");
   m_savetree->Branch("sel_pi_LL_ln_Q1", &sel_pi_LL_ln_Q1, "sel_pi_LL_ln_Q1/D");
   m_savetree->Branch("sel_pi_LL_ln_Q2", &sel_pi_LL_ln_Q2, "sel_pi_LL_ln_Q2/D");
   m_savetree->Branch("sel_pi_LL_ln_Q3", &sel_pi_LL_ln_Q3, "sel_pi_LL_ln_Q3/D");
   m_savetree->Branch("sel_pi_LL_ln_Q4", &sel_pi_LL_ln_Q4, "sel_pi_LL_ln_Q4/D");
   m_savetree->Branch("sel_pi_LL_ln_Q5", &sel_pi_LL_ln_Q5, "sel_pi_LL_ln_Q5/D");
   m_savetree->Branch("sel_pi_LL_mom", &sel_pi_LL_mom, "sel_pi_LL_mom/D");
   m_savetree->Branch("sel_pi_LL_mom_altH", &sel_pi_LL_mom_altH, "sel_pi_LL_mom_altH/D");
   m_savetree->Branch("sel_pi_LL_pTMag", &sel_pi_LL_pTMag, "sel_pi_LL_pTMag/D");
   m_savetree->Branch("sel_pi_LL_pTMag_tnudir", &sel_pi_LL_pTMag_tnudir, "sel_pi_LL_pTMag_tnudir/D");
   m_savetree->Branch("sel_pi_LL_pTMag_tpimom", &sel_pi_LL_pTMag_tpimom, "sel_pi_LL_pTMag_tpimom/D");
   m_savetree->Branch("sel_pi_LL_pTT", &sel_pi_LL_pTT, "sel_pi_LL_pTT/D");
   m_savetree->Branch("sel_pi_LL_pTT_tnudir", &sel_pi_LL_pTT_tnudir, "sel_pi_LL_pTT_tnudir/D");
   m_savetree->Branch("sel_pi_LL_pTT_tpimom", &sel_pi_LL_pTT_tpimom, "sel_pi_LL_pTT_tpimom/D");
   m_savetree->Branch("sel_pi_LL_score", &sel_pi_LL_score, "sel_pi_LL_score/D");
   m_savetree->Branch("sel_pi_LL_score_altH", &sel_pi_LL_score_altH, "sel_pi_LL_score_altH/D");
   m_savetree->Branch("sel_pi_chi2ndf", &sel_pi_chi2ndf, "sel_pi_chi2ndf/D");
   m_savetree->Branch("sel_pi_det_frac", &sel_pi_det_frac, "sel_pi_det_frac/D");
   m_savetree->Branch("sel_pi_det_otherE", &sel_pi_det_otherE, "sel_pi_det_otherE/D");
   m_savetree->Branch("sel_pi_trueE", &sel_pi_trueE, "sel_pi_trueE/D");
   m_savetree->Branch("sel_pi_trueKE", &sel_pi_trueKE, "sel_pi_trueKE/D");
   m_savetree->Branch("sel_pi_truePhi", &sel_pi_truePhi, "sel_pi_truePhi/D");
   m_savetree->Branch("sel_pi_trueTheta", &sel_pi_trueTheta, "sel_pi_trueTheta/D");
   m_savetree->Branch("sel_pi_truemom", &sel_pi_truemom, "sel_pi_truemom/D");
   m_savetree->Branch("sel_pi_truepTMag", &sel_pi_truepTMag, "sel_pi_truepTMag/D");
   m_savetree->Branch("sel_pi_truepTT", &sel_pi_truepTT, "sel_pi_truepTT/D");
   m_savetree->Branch("sel_pr_CaloE", &sel_pr_CaloE, "sel_pr_CaloE/D");
   m_savetree->Branch("sel_pr_EX_E", &sel_pr_EX_E, "sel_pr_EX_E/D");
   m_savetree->Branch("sel_pr_EX_E_altH", &sel_pr_EX_E_altH, "sel_pr_EX_E_altH/D");
   m_savetree->Branch("sel_pr_EX_KE", &sel_pr_EX_KE, "sel_pr_EX_KE/D");
   m_savetree->Branch("sel_pr_EX_Phi", &sel_pr_EX_Phi, "sel_pr_EX_Phi/D");
   m_savetree->Branch("sel_pr_EX_Theta", &sel_pr_EX_Theta, "sel_pr_EX_Theta/D");
   m_savetree->Branch("sel_pr_EX_XTheta", &sel_pr_EX_XTheta, "sel_pr_EX_XTheta/D");
   m_savetree->Branch("sel_pr_EX_YTheta", &sel_pr_EX_YTheta, "sel_pr_EX_YTheta/D");
   m_savetree->Branch("sel_pr_EX_ln_Q0", &sel_pr_EX_ln_Q0, "sel_pr_EX_ln_Q0/D");
   m_savetree->Branch("sel_pr_EX_ln_Q1", &sel_pr_EX_ln_Q1, "sel_pr_EX_ln_Q1/D");
   m_savetree->Branch("sel_pr_EX_ln_Q2", &sel_pr_EX_ln_Q2, "sel_pr_EX_ln_Q2/D");
   m_savetree->Branch("sel_pr_EX_ln_Q3", &sel_pr_EX_ln_Q3, "sel_pr_EX_ln_Q3/D");
   m_savetree->Branch("sel_pr_EX_ln_Q4", &sel_pr_EX_ln_Q4, "sel_pr_EX_ln_Q4/D");
   m_savetree->Branch("sel_pr_EX_ln_Q5", &sel_pr_EX_ln_Q5, "sel_pr_EX_ln_Q5/D");
   m_savetree->Branch("sel_pr_EX_mom", &sel_pr_EX_mom, "sel_pr_EX_mom/D");
   m_savetree->Branch("sel_pr_EX_mom_altH", &sel_pr_EX_mom_altH, "sel_pr_EX_mom_altH/D");
   m_savetree->Branch("sel_pr_EX_pTMag", &sel_pr_EX_pTMag, "sel_pr_EX_pTMag/D");
   m_savetree->Branch("sel_pr_EX_pTMag_tnudir", &sel_pr_EX_pTMag_tnudir, "sel_pr_EX_pTMag_tnudir/D");
   m_savetree->Branch("sel_pr_EX_pTMag_tprmom", &sel_pr_EX_pTMag_tprmom, "sel_pr_EX_pTMag_tprmom/D");
   m_savetree->Branch("sel_pr_EX_pTT", &sel_pr_EX_pTT, "sel_pr_EX_pTT/D");
   m_savetree->Branch("sel_pr_EX_pTT_tnudir", &sel_pr_EX_pTT_tnudir, "sel_pr_EX_pTT_tnudir/D");
   m_savetree->Branch("sel_pr_EX_pTT_tprmom", &sel_pr_EX_pTT_tprmom, "sel_pr_EX_pTT_tprmom/D");
   m_savetree->Branch("sel_pr_EX_score", &sel_pr_EX_score, "sel_pr_EX_score/D");
   m_savetree->Branch("sel_pr_EX_score_altH", &sel_pr_EX_score_altH, "sel_pr_EX_score_altH/D");
   m_savetree->Branch("sel_pr_LL_E", &sel_pr_LL_E, "sel_pr_LL_E/D");
   m_savetree->Branch("sel_pr_LL_E_altH", &sel_pr_LL_E_altH, "sel_pr_LL_E_altH/D");
   m_savetree->Branch("sel_pr_LL_KE", &sel_pr_LL_KE, "sel_pr_LL_KE/D");
   m_savetree->Branch("sel_pr_LL_Phi", &sel_pr_LL_Phi, "sel_pr_LL_Phi/D");
   m_savetree->Branch("sel_pr_LL_Theta", &sel_pr_LL_Theta, "sel_pr_LL_Theta/D");
   m_savetree->Branch("sel_pr_LL_XTheta", &sel_pr_LL_XTheta, "sel_pr_LL_XTheta/D");
   m_savetree->Branch("sel_pr_LL_YTheta", &sel_pr_LL_YTheta, "sel_pr_LL_YTheta/D");
   m_savetree->Branch("sel_pr_LL_ln_Q0", &sel_pr_LL_ln_Q0, "sel_pr_LL_ln_Q0/D");
   m_savetree->Branch("sel_pr_LL_ln_Q1", &sel_pr_LL_ln_Q1, "sel_pr_LL_ln_Q1/D");
   m_savetree->Branch("sel_pr_LL_ln_Q2", &sel_pr_LL_ln_Q2, "sel_pr_LL_ln_Q2/D");
   m_savetree->Branch("sel_pr_LL_ln_Q3", &sel_pr_LL_ln_Q3, "sel_pr_LL_ln_Q3/D");
   m_savetree->Branch("sel_pr_LL_ln_Q4", &sel_pr_LL_ln_Q4, "sel_pr_LL_ln_Q4/D");
   m_savetree->Branch("sel_pr_LL_ln_Q5", &sel_pr_LL_ln_Q5, "sel_pr_LL_ln_Q5/D");
   m_savetree->Branch("sel_pr_LL_mom", &sel_pr_LL_mom, "sel_pr_LL_mom/D");
   m_savetree->Branch("sel_pr_LL_mom_altH", &sel_pr_LL_mom_altH, "sel_pr_LL_mom_altH/D");
   m_savetree->Branch("sel_pr_LL_pTMag", &sel_pr_LL_pTMag, "sel_pr_LL_pTMag/D");
   m_savetree->Branch("sel_pr_LL_pTMag_tnudir", &sel_pr_LL_pTMag_tnudir, "sel_pr_LL_pTMag_tnudir/D");
   m_savetree->Branch("sel_pr_LL_pTMag_tprmom", &sel_pr_LL_pTMag_tprmom, "sel_pr_LL_pTMag_tprmom/D");
   m_savetree->Branch("sel_pr_LL_pTT", &sel_pr_LL_pTT, "sel_pr_LL_pTT/D");
   m_savetree->Branch("sel_pr_LL_pTT_tnudir", &sel_pr_LL_pTT_tnudir, "sel_pr_LL_pTT_tnudir/D");
   m_savetree->Branch("sel_pr_LL_pTT_tprmom", &sel_pr_LL_pTT_tprmom, "sel_pr_LL_pTT_tprmom/D");
   m_savetree->Branch("sel_pr_LL_score", &sel_pr_LL_score, "sel_pr_LL_score/D");
   m_savetree->Branch("sel_pr_LL_score_altH", &sel_pr_LL_score_altH, "sel_pr_LL_score_altH/D");
   m_savetree->Branch("sel_pr_chi2ndf", &sel_pr_chi2ndf, "sel_pr_chi2ndf/D");
   m_savetree->Branch("sel_pr_det_frac", &sel_pr_det_frac, "sel_pr_det_frac/D");
   m_savetree->Branch("sel_pr_det_otherE", &sel_pr_det_otherE, "sel_pr_det_otherE/D");
   m_savetree->Branch("sel_pr_trueE", &sel_pr_trueE, "sel_pr_trueE/D");
   m_savetree->Branch("sel_pr_trueKE", &sel_pr_trueKE, "sel_pr_trueKE/D");
   m_savetree->Branch("sel_pr_truePhi", &sel_pr_truePhi, "sel_pr_truePhi/D");
   m_savetree->Branch("sel_pr_trueTheta", &sel_pr_trueTheta, "sel_pr_trueTheta/D");
   m_savetree->Branch("sel_pr_truemom", &sel_pr_truemom, "sel_pr_truemom/D");
   m_savetree->Branch("sel_pr_truepTMag", &sel_pr_truepTMag, "sel_pr_truepTMag/D");
   m_savetree->Branch("sel_pr_truepTT", &sel_pr_truepTT, "sel_pr_truepTT/D");
   m_savetree->Branch("sel_trueEnu", &sel_trueEnu, "sel_trueEnu/D");
   m_savetree->Branch("sel_trueQ2", &sel_trueQ2, "sel_trueQ2/D");
   m_savetree->Branch("sel_truedalphaT", &sel_truedalphaT, "sel_truedalphaT/D");
   m_savetree->Branch("sel_truedpT", &sel_truedpT, "sel_truedpT/D");
   m_savetree->Branch("sel_truedpTT", &sel_truedpTT, "sel_truedpTT/D");
   m_savetree->Branch("sel_truedpTT_pi", &sel_truedpTT_pi, "sel_truedpTT_pi/D");
   m_savetree->Branch("sel_truedpTT_pi_dir", &sel_truedpTT_pi_dir, "sel_truedpTT_pi_dir/D");
   m_savetree->Branch("sel_truedpTT_pr", &sel_truedpTT_pr, "sel_truedpTT_pr/D");
   m_savetree->Branch("sel_truedpTT_pr_dir", &sel_truedpTT_pr_dir, "sel_truedpTT_pr_dir/D");
   m_savetree->Branch("sel_truedphiT", &sel_truedphiT, "sel_truedphiT/D");
   m_savetree->Branch("n_iso_blobs", &n_iso_blobs, "n_iso_blobs/I");
   m_savetree->Branch("sel_iso_blob_nclusters", &sel_iso_blob_nclusters, "sel_iso_blob_nclusters[n_iso_blobs]/I");
   m_savetree->Branch("sel_VTX", sel_VTX, "sel_VTX[3]/D");
   m_savetree->Branch("sel_dpT_vec_EX", sel_dpT_vec_EX, "sel_dpT_vec_EX[3]/D");
   m_savetree->Branch("sel_dpT_vec_LL", sel_dpT_vec_LL, "sel_dpT_vec_LL[3]/D");
   m_savetree->Branch("sel_iso_blob_energy", &sel_iso_blob_energy, "sel_iso_blob_energy[n_iso_blobs]/D");
   m_savetree->Branch("sel_meanPDP", sel_meanPDP, "sel_meanPDP[3]/D");
   m_savetree->Branch("sel_mu_4mom", sel_mu_4mom, "sel_mu_4mom[4]/D");
   m_savetree->Branch("sel_mu_endpos", sel_mu_endpos, "sel_mu_endpos[3]/D");
   m_savetree->Branch("sel_mu_pT", sel_mu_pT, "sel_mu_pT[3]/D");
   m_savetree->Branch("sel_mu_pT_tmumom", sel_mu_pT_tmumom, "sel_mu_pT_tmumom[3]/D");
   m_savetree->Branch("sel_mu_pT_tnudir", sel_mu_pT_tnudir, "sel_mu_pT_tnudir[3]/D");
   m_savetree->Branch("sel_mu_startdir", sel_mu_startdir, "sel_mu_startdir[3]/D");
   m_savetree->Branch("sel_mu_startpos", sel_mu_startpos, "sel_mu_startpos[3]/D");
   m_savetree->Branch("sel_mu_true4mom", sel_mu_true4mom, "sel_mu_true4mom[4]/D");
   m_savetree->Branch("sel_mu_trueendpos", sel_mu_trueendpos, "sel_mu_trueendpos[3]/D");
   m_savetree->Branch("sel_mu_truepT", sel_mu_truepT, "sel_mu_truepT[3]/D");
   m_savetree->Branch("sel_mu_truestartdir", sel_mu_truestartdir, "sel_mu_truestartdir[3]/D");
   m_savetree->Branch("sel_mu_truestartpos", sel_mu_truestartpos, "sel_mu_truestartpos[3]/D");
   m_savetree->Branch("sel_nu_dir_001", sel_nu_dir_001, "sel_nu_dir_001[3]/D");
   m_savetree->Branch("sel_nu_dir_PDP", sel_nu_dir_PDP, "sel_nu_dir_PDP[3]/D");
   m_savetree->Branch("sel_pi_EX_4mom", sel_pi_EX_4mom, "sel_pi_EX_4mom[4]/D");
   m_savetree->Branch("sel_pi_EX_pT", sel_pi_EX_pT, "sel_pi_EX_pT[3]/D");
   m_savetree->Branch("sel_pi_EX_pT_tnudir", sel_pi_EX_pT_tnudir, "sel_pi_EX_pT_tnudir[3]/D");
   m_savetree->Branch("sel_pi_EX_pT_tpimom", sel_pi_EX_pT_tpimom, "sel_pi_EX_pT_tpimom[3]/D");
   m_savetree->Branch("sel_pi_LL_4mom", sel_pi_LL_4mom, "sel_pi_LL_4mom[4]/D");
   m_savetree->Branch("sel_pi_LL_pT", sel_pi_LL_pT, "sel_pi_LL_pT[3]/D");
   m_savetree->Branch("sel_pi_LL_pT_tnudir", sel_pi_LL_pT_tnudir, "sel_pi_LL_pT_tnudir[3]/D");
   m_savetree->Branch("sel_pi_LL_pT_tpimom", sel_pi_LL_pT_tpimom, "sel_pi_LL_pT_tpimom[3]/D");
   m_savetree->Branch("sel_pi_endpos", sel_pi_endpos, "sel_pi_endpos[3]/D");
   m_savetree->Branch("sel_pi_startdir", sel_pi_startdir, "sel_pi_startdir[3]/D");
   m_savetree->Branch("sel_pi_startpos", sel_pi_startpos, "sel_pi_startpos[3]/D");
   m_savetree->Branch("sel_pi_true4mom", sel_pi_true4mom, "sel_pi_true4mom[4]/D");
   m_savetree->Branch("sel_pi_trueendpos", sel_pi_trueendpos, "sel_pi_trueendpos[3]/D");
   m_savetree->Branch("sel_pi_truepT", sel_pi_truepT, "sel_pi_truepT[3]/D");
   m_savetree->Branch("sel_pi_truestartdir", sel_pi_truestartdir, "sel_pi_truestartdir[3]/D");
   m_savetree->Branch("sel_pi_truestartpos", sel_pi_truestartpos, "sel_pi_truestartpos[3]/D");
   m_savetree->Branch("sel_pr_EX_4mom", sel_pr_EX_4mom, "sel_pr_EX_4mom[4]/D");
   m_savetree->Branch("sel_pr_EX_pT", sel_pr_EX_pT, "sel_pr_EX_pT[3]/D");
   m_savetree->Branch("sel_pr_EX_pT_tnudir", sel_pr_EX_pT_tnudir, "sel_pr_EX_pT_tnudir[3]/D");
   m_savetree->Branch("sel_pr_EX_pT_tprmom", sel_pr_EX_pT_tprmom, "sel_pr_EX_pT_tprmom[3]/D");
   m_savetree->Branch("sel_pr_LL_4mom", sel_pr_LL_4mom, "sel_pr_LL_4mom[4]/D");
   m_savetree->Branch("sel_pr_LL_pT", sel_pr_LL_pT, "sel_pr_LL_pT[3]/D");
   m_savetree->Branch("sel_pr_LL_pT_tnudir", sel_pr_LL_pT_tnudir, "sel_pr_LL_pT_tnudir[3]/D");
   m_savetree->Branch("sel_pr_LL_pT_tprmom", sel_pr_LL_pT_tprmom, "sel_pr_LL_pT_tprmom[3]/D");
   m_savetree->Branch("sel_pr_endpos", sel_pr_endpos, "sel_pr_endpos[3]/D");
   m_savetree->Branch("sel_pr_startdir", sel_pr_startdir, "sel_pr_startdir[3]/D");
   m_savetree->Branch("sel_pr_startpos", sel_pr_startpos, "sel_pr_startpos[3]/D");
   m_savetree->Branch("sel_pr_true4mom", sel_pr_true4mom, "sel_pr_true4mom[4]/D");
   m_savetree->Branch("sel_pr_trueendpos", sel_pr_trueendpos, "sel_pr_trueendpos[3]/D");
   m_savetree->Branch("sel_pr_truepT", sel_pr_truepT, "sel_pr_truepT[3]/D");
   m_savetree->Branch("sel_pr_truestartdir", sel_pr_truestartdir, "sel_pr_truestartdir[3]/D");
   m_savetree->Branch("sel_pr_truestartpos", sel_pr_truestartpos, "sel_pr_truestartpos[3]/D");
   m_savetree->Branch("sel_truePDP", sel_truePDP, "sel_truePDP[3]/D");
   m_savetree->Branch("sel_trueVTX", sel_trueVTX, "sel_trueVTX[3]/D");
   m_savetree->Branch("sel_true_nu_dir_PDP", sel_true_nu_dir_PDP, "sel_true_nu_dir_PDP[3]/D");
   m_savetree->Branch("sel_truedpT_vec", sel_truedpT_vec, "sel_truedpT_vec[3]/D");
   m_savetree->Branch("physEvtNum", &physEvtNum, "physEvtNum/I");
   m_savetree->Branch("n_hyps", &n_hyps, "n_hyps/I");
   m_savetree->Branch("processType", &processType, "processType/I");
   m_savetree->Branch("primaryPart", &primaryPart, "primaryPart/I");
   m_savetree->Branch("n_slices", &n_slices, "n_slices/I");
   m_savetree->Branch("slice_numbers", slice_numbers, "slice_numbers[n_slices]/I");
   m_savetree->Branch("shared_slice", &shared_slice, "shared_slice/I");
   m_savetree->Branch("vtx", vtx, "vtx[4]/D");
   m_savetree->Branch("vtxErr", vtxErr, "vtxErr[4]/D");
   m_savetree->Branch("E", E, "E[4]/D");
   m_savetree->Branch("found_truth", &found_truth, "found_truth/O");
   m_savetree->Branch("isMinosMatchTrack", &isMinosMatchTrack, "isMinosMatchTrack/O");
   m_savetree->Branch("isMinosMatchStub", &isMinosMatchStub, "isMinosMatchStub/O");
   m_savetree->Branch("contained_evt", &contained_evt, "contained_evt/I");
   m_savetree->Branch("muon_charge", &muon_charge, "muon_charge/I");
   m_savetree->Branch("n_anchored_long_trk_prongs", &n_anchored_long_trk_prongs, "n_anchored_long_trk_prongs/I");
   m_savetree->Branch("n_anchored_short_trk_prongs", &n_anchored_short_trk_prongs, "n_anchored_short_trk_prongs/I");
   m_savetree->Branch("n_iso_trk_prongs", &n_iso_trk_prongs, "n_iso_trk_prongs/I");
   m_savetree->Branch("n_prongs", &n_prongs, "n_prongs/I");
   m_savetree->Branch("n_tracks3", &n_tracks3, "n_tracks3/I");
   m_savetree->Branch("ncuts", &ncuts, "ncuts/I");
   m_savetree->Branch("new_tracks", &new_tracks, "new_tracks/I");
   m_savetree->Branch("nsplits", &nsplits, "nsplits/I");
   m_savetree->Branch("target_region", &target_region, "target_region/I");
   m_savetree->Branch("true_target_region", &true_target_region, "true_target_region/I");
   m_savetree->Branch("vert_exists", &vert_exists, "vert_exists/I");
   m_savetree->Branch("accum_level", accum_level, "accum_level[2]/I");
   m_savetree->Branch("truth_reco_isMinosMatch", &truth_reco_isMinosMatch, "truth_reco_isMinosMatch/O");
   m_savetree->Branch("truth_muon_charge", &truth_muon_charge, "truth_muon_charge/I");
   m_savetree->Branch("truth_n_ele", &truth_n_ele, "truth_n_ele/I");
   m_savetree->Branch("truth_n_kPM", &truth_n_kPM, "truth_n_kPM/I");
   m_savetree->Branch("truth_n_kaO", &truth_n_kaO, "truth_n_kaO/I");
   m_savetree->Branch("truth_n_muo", &truth_n_muo, "truth_n_muo/I");
   m_savetree->Branch("truth_n_ntn", &truth_n_ntn, "truth_n_ntn/I");
   m_savetree->Branch("truth_n_pho", &truth_n_pho, "truth_n_pho/I");
   m_savetree->Branch("truth_n_pi0", &truth_n_pi0, "truth_n_pi0/I");
   m_savetree->Branch("truth_n_piM", &truth_n_piM, "truth_n_piM/I");
   m_savetree->Branch("truth_n_piP", &truth_n_piP, "truth_n_piP/I");
   m_savetree->Branch("truth_n_pro", &truth_n_pro, "truth_n_pro/I");
   m_savetree->Branch("truth_n_tau", &truth_n_tau, "truth_n_tau/I");
   m_savetree->Branch("truth_ncuts", &truth_ncuts, "truth_ncuts/I");
   m_savetree->Branch("truth_nsplits", &truth_nsplits, "truth_nsplits/I");
   m_savetree->Branch("truth_pi_EX_michel", &truth_pi_EX_michel, "truth_pi_EX_michel/I");
   m_savetree->Branch("truth_pi_LL_michel", &truth_pi_LL_michel, "truth_pi_LL_michel/I");
   m_savetree->Branch("truth_pr_EX_michel", &truth_pr_EX_michel, "truth_pr_EX_michel/I");
   m_savetree->Branch("truth_pr_LL_michel", &truth_pr_LL_michel, "truth_pr_LL_michel/I");
   m_savetree->Branch("truth_reco_target", &truth_reco_target, "truth_reco_target/I");
   m_savetree->Branch("truth_should_be_accepted", &truth_should_be_accepted, "truth_should_be_accepted/I");
   m_savetree->Branch("truth_true_target_region", &truth_true_target_region, "truth_true_target_region/I");
   m_savetree->Branch("truth_mu_E", &truth_mu_E, "truth_mu_E/D");
   m_savetree->Branch("truth_mu_KE", &truth_mu_KE, "truth_mu_KE/D");
   m_savetree->Branch("truth_mu_Phi", &truth_mu_Phi, "truth_mu_Phi/D");
   m_savetree->Branch("truth_mu_Theta", &truth_mu_Theta, "truth_mu_Theta/D");
   m_savetree->Branch("truth_mu_mom", &truth_mu_mom, "truth_mu_mom/D");
   m_savetree->Branch("truth_mu_pTMag", &truth_mu_pTMag, "truth_mu_pTMag/D");
   m_savetree->Branch("truth_mu_pTT", &truth_mu_pTT, "truth_mu_pTT/D");
   m_savetree->Branch("truth_pi_E", &truth_pi_E, "truth_pi_E/D");
   m_savetree->Branch("truth_pi_EX_score", &truth_pi_EX_score, "truth_pi_EX_score/D");
   m_savetree->Branch("truth_pi_EX_score_altH", &truth_pi_EX_score_altH, "truth_pi_EX_score_altH/D");
   m_savetree->Branch("truth_pi_KE", &truth_pi_KE, "truth_pi_KE/D");
   m_savetree->Branch("truth_pi_LL_score", &truth_pi_LL_score, "truth_pi_LL_score/D");
   m_savetree->Branch("truth_pi_LL_score_altH", &truth_pi_LL_score_altH, "truth_pi_LL_score_altH/D");
   m_savetree->Branch("truth_pi_Phi", &truth_pi_Phi, "truth_pi_Phi/D");
   m_savetree->Branch("truth_pi_Theta", &truth_pi_Theta, "truth_pi_Theta/D");
   m_savetree->Branch("truth_pi_mom", &truth_pi_mom, "truth_pi_mom/D");
   m_savetree->Branch("truth_pi_pTMag", &truth_pi_pTMag, "truth_pi_pTMag/D");
   m_savetree->Branch("truth_pi_pTT", &truth_pi_pTT, "truth_pi_pTT/D");
   m_savetree->Branch("truth_pr_E", &truth_pr_E, "truth_pr_E/D");
   m_savetree->Branch("truth_pr_EX_score", &truth_pr_EX_score, "truth_pr_EX_score/D");
   m_savetree->Branch("truth_pr_EX_score_altH", &truth_pr_EX_score_altH, "truth_pr_EX_score_altH/D");
   m_savetree->Branch("truth_pr_KE", &truth_pr_KE, "truth_pr_KE/D");
   m_savetree->Branch("truth_pr_LL_score", &truth_pr_LL_score, "truth_pr_LL_score/D");
   m_savetree->Branch("truth_pr_LL_score_altH", &truth_pr_LL_score_altH, "truth_pr_LL_score_altH/D");
   m_savetree->Branch("truth_pr_Phi", &truth_pr_Phi, "truth_pr_Phi/D");
   m_savetree->Branch("truth_pr_Theta", &truth_pr_Theta, "truth_pr_Theta/D");
   m_savetree->Branch("truth_pr_mom", &truth_pr_mom, "truth_pr_mom/D");
   m_savetree->Branch("truth_pr_pTMag", &truth_pr_pTMag, "truth_pr_pTMag/D");
   m_savetree->Branch("truth_pr_pTT", &truth_pr_pTT, "truth_pr_pTT/D");
   m_savetree->Branch("truth_trueEnu", &truth_trueEnu, "truth_trueEnu/D");
   m_savetree->Branch("truth_trueQ2", &truth_trueQ2, "truth_trueQ2/D");
   m_savetree->Branch("truth_truedalphaT", &truth_truedalphaT, "truth_truedalphaT/D");
   m_savetree->Branch("truth_truedpT", &truth_truedpT, "truth_truedpT/D");
   m_savetree->Branch("truth_truedpTT", &truth_truedpTT, "truth_truedpTT/D");
   m_savetree->Branch("truth_truedpTT_pi", &truth_truedpTT_pi, "truth_truedpTT_pi/D");
   m_savetree->Branch("truth_truedpTT_pi_dir", &truth_truedpTT_pi_dir, "truth_truedpTT_pi_dir/D");
   m_savetree->Branch("truth_truedpTT_pr", &truth_truedpTT_pr, "truth_truedpTT_pr/D");
   m_savetree->Branch("truth_truedpTT_pr_dir", &truth_truedpTT_pr_dir, "truth_truedpTT_pr_dir/D");
   m_savetree->Branch("truth_truedphiT", &truth_truedphiT, "truth_truedphiT/D");
   m_savetree->Branch("truth_accum_level", truth_accum_level, "truth_accum_level[2]/I");
   m_savetree->Branch("truth_mu_4mom", truth_mu_4mom, "truth_mu_4mom[4]/D");
   m_savetree->Branch("truth_mu_pT", truth_mu_pT, "truth_mu_pT[3]/D");
   m_savetree->Branch("truth_pi_4mom", truth_pi_4mom, "truth_pi_4mom[4]/D");
   m_savetree->Branch("truth_pi_pT", truth_pi_pT, "truth_pi_pT[3]/D");
   m_savetree->Branch("truth_pr_4mom", truth_pr_4mom, "truth_pr_4mom[4]/D");
   m_savetree->Branch("truth_pr_pT", truth_pr_pT, "truth_pr_pT[3]/D");
   m_savetree->Branch("truth_truedpT_vec", truth_truedpT_vec, "truth_truedpT_vec[3]/D");
   m_savetree->Branch("ev_run", &ev_run, "ev_run/I");
   m_savetree->Branch("ev_subrun", &ev_subrun, "ev_subrun/I");
   m_savetree->Branch("ev_detector", &ev_detector, "ev_detector/I");
   m_savetree->Branch("ev_triggerType", &ev_triggerType, "ev_triggerType/I");
   m_savetree->Branch("ev_gate", &ev_gate, "ev_gate/I");
   m_savetree->Branch("ev_global_gate", &ev_global_gate, "ev_global_gate/I");
   m_savetree->Branch("ev_gps_time_sec", &ev_gps_time_sec, "ev_gps_time_sec/I");
   m_savetree->Branch("ev_gps_time_usec", &ev_gps_time_usec, "ev_gps_time_usec/I");
   m_savetree->Branch("mc_run", &mc_run, "mc_run/I");
   m_savetree->Branch("mc_subrun", &mc_subrun, "mc_subrun/I");
   m_savetree->Branch("mc_nInteractions", &mc_nInteractions, "mc_nInteractions/I");
   m_savetree->Branch("mc_MIState", &mc_MIState, "mc_MIState/I");
   m_savetree->Branch("mc_pot", &mc_pot, "mc_pot/D");
   m_savetree->Branch("mc_beamConfig", &mc_beamConfig, "mc_beamConfig/I");
   m_savetree->Branch("mc_processType", &mc_processType, "mc_processType/I");
   m_savetree->Branch("mc_nthEvtInSpill", &mc_nthEvtInSpill, "mc_nthEvtInSpill/I");
   m_savetree->Branch("mc_nthEvtInFile", &mc_nthEvtInFile, "mc_nthEvtInFile/I");
   m_savetree->Branch("mc_intType", &mc_intType, "mc_intType/I");
   m_savetree->Branch("mc_current", &mc_current, "mc_current/I");
   m_savetree->Branch("mc_charm", &mc_charm, "mc_charm/I");
   m_savetree->Branch("mc_weight", &mc_weight, "mc_weight/D");
   m_savetree->Branch("mc_XSec", &mc_XSec, "mc_XSec/D");
   m_savetree->Branch("mc_diffXSec", &mc_diffXSec, "mc_diffXSec/D");
   m_savetree->Branch("mc_incoming", &mc_incoming, "mc_incoming/I");
   m_savetree->Branch("mc_fluxDriverProb", &mc_fluxDriverProb, "mc_fluxDriverProb/D");
   m_savetree->Branch("mc_targetNucleus", &mc_targetNucleus, "mc_targetNucleus/I");
   m_savetree->Branch("mc_targetZ", &mc_targetZ, "mc_targetZ/I");
   m_savetree->Branch("mc_targetA", &mc_targetA, "mc_targetA/I");
   m_savetree->Branch("mc_targetNucleon", &mc_targetNucleon, "mc_targetNucleon/I");
   m_savetree->Branch("mc_struckQuark", &mc_struckQuark, "mc_struckQuark/I");
   m_savetree->Branch("mc_seaQuark", &mc_seaQuark, "mc_seaQuark/I");
   m_savetree->Branch("mc_resID", &mc_resID, "mc_resID/I");
   m_savetree->Branch("mc_primaryLepton", &mc_primaryLepton, "mc_primaryLepton/I");
   m_savetree->Branch("mc_incomingE", &mc_incomingE, "mc_incomingE/D");
   m_savetree->Branch("mc_Bjorkenx", &mc_Bjorkenx, "mc_Bjorkenx/D");
   m_savetree->Branch("mc_Bjorkeny", &mc_Bjorkeny, "mc_Bjorkeny/D");
   m_savetree->Branch("mc_Q2", &mc_Q2, "mc_Q2/D");
   m_savetree->Branch("mc_nuT", &mc_nuT, "mc_nuT/D");
   m_savetree->Branch("mc_w", &mc_w, "mc_w/D");
   m_savetree->Branch("mc_vtx", mc_vtx, "mc_vtx[4]/D");
   m_savetree->Branch("mc_incomingPartVec", mc_incomingPartVec, "mc_incomingPartVec[4]/D");
   m_savetree->Branch("mc_initNucVec", mc_initNucVec, "mc_initNucVec[4]/D");
   m_savetree->Branch("mc_primFSLepton", mc_primFSLepton, "mc_primFSLepton[4]/D");
   m_savetree->Branch("mc_nFSPart", &mc_nFSPart, "mc_nFSPart/I");
   m_savetree->Branch("mc_FSPartPx", mc_FSPartPx, "mc_FSPartPx[mc_nFSPart]/D");
   m_savetree->Branch("mc_FSPartPy", mc_FSPartPy, "mc_FSPartPy[mc_nFSPart]/D");
   m_savetree->Branch("mc_FSPartPz", mc_FSPartPz, "mc_FSPartPz[mc_nFSPart]/D");
   m_savetree->Branch("mc_FSPartE", mc_FSPartE, "mc_FSPartE[mc_nFSPart]/D");
   m_savetree->Branch("mc_FSPartPDG", mc_FSPartPDG, "mc_FSPartPDG[mc_nFSPart]/I");
   m_savetree->Branch("mc_er_nPart", &mc_er_nPart, "mc_er_nPart/I");
   m_savetree->Branch("mc_er_ID", mc_er_ID, "mc_er_ID[mc_er_nPart]/I");
   m_savetree->Branch("mc_er_status", mc_er_status, "mc_er_status[mc_er_nPart]/I");
   m_savetree->Branch("mc_er_posInNucX", mc_er_posInNucX, "mc_er_posInNucX[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_posInNucY", mc_er_posInNucY, "mc_er_posInNucY[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_posInNucZ", mc_er_posInNucZ, "mc_er_posInNucZ[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_Px", mc_er_Px, "mc_er_Px[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_Py", mc_er_Py, "mc_er_Py[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_Pz", mc_er_Pz, "mc_er_Pz[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_E", mc_er_E, "mc_er_E[mc_er_nPart]/D");
   m_savetree->Branch("mc_er_FD", mc_er_FD, "mc_er_FD[mc_er_nPart]/I");
   m_savetree->Branch("mc_er_LD", mc_er_LD, "mc_er_LD[mc_er_nPart]/I");
   m_savetree->Branch("mc_er_mother", mc_er_mother, "mc_er_mother[mc_er_nPart]/I");
   m_savetree->Branch("mc_fr_nNuAncestorIDs", &mc_fr_nNuAncestorIDs, "mc_fr_nNuAncestorIDs/I");
   m_savetree->Branch("mc_fr_nuAncestorIDs", mc_fr_nuAncestorIDs, "mc_fr_nuAncestorIDs[mc_fr_nNuAncestorIDs]/I");
   m_savetree->Branch("mc_fr_nuParentID", &mc_fr_nuParentID, "mc_fr_nuParentID/I");
   m_savetree->Branch("mc_fr_decMode", &mc_fr_decMode, "mc_fr_decMode/I");
   m_savetree->Branch("mc_fr_primProtonVtx", mc_fr_primProtonVtx, "mc_fr_primProtonVtx[3]/D");
   m_savetree->Branch("mc_fr_primProtonP", mc_fr_primProtonP, "mc_fr_primProtonP[4]/D");
   m_savetree->Branch("mc_fr_nuParentDecVtx", mc_fr_nuParentDecVtx, "mc_fr_nuParentDecVtx[3]/D");
   m_savetree->Branch("mc_fr_nuParentProdVtx", mc_fr_nuParentProdVtx, "mc_fr_nuParentProdVtx[3]/D");
   m_savetree->Branch("mc_fr_nuParentProdP", mc_fr_nuParentProdP, "mc_fr_nuParentProdP[4]/D");
   m_savetree->Branch("mc_cvweight_total", &mc_cvweight_total, "mc_cvweight_total/D");
   m_savetree->Branch("wgt", &wgt, "wgt/D");
   m_savetree->Branch("mc_cvweight_totalFlux", &mc_cvweight_totalFlux, "mc_cvweight_totalFlux/D");
   m_savetree->Branch("mc_cvweight_totalXsec", &mc_cvweight_totalXsec, "mc_cvweight_totalXsec/D");
   m_savetree->Branch("mc_ppfx1_cvweight", &mc_ppfx1_cvweight, "mc_ppfx1_cvweight/D");
   m_savetree->Branch("mc_hornCurrent_cvweight", &mc_hornCurrent_cvweight, "mc_hornCurrent_cvweight/D");
   m_savetree->Branch("mc_gen1_cvweight_total", &mc_gen1_cvweight_total, "mc_gen1_cvweight_total/D");
   m_savetree->Branch("gen1_wgt", &gen1_wgt, "gen1_wgt/D");
   m_savetree->Branch("mc_gen1_cvweight_totalFlux", &mc_gen1_cvweight_totalFlux, "mc_gen1_cvweight_totalFlux/D");
   m_savetree->Branch("mc_gen1_cvweight_NA49", &mc_gen1_cvweight_NA49, "mc_gen1_cvweight_NA49/D");
   m_savetree->Branch("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, "mc_wgt_Flux_BeamFocus_sz/I");
   m_savetree->Branch("mc_wgt_Flux_BeamFocus", &mc_wgt_Flux_BeamFocus, "mc_wgt_Flux_BeamFocus[mc_wgt_Flux_BeamFocus_sz]/D");
   m_savetree->Branch("mc_wgt_gen1_Flux_Tertiary_sz", &mc_wgt_gen1_Flux_Tertiary_sz, "mc_wgt_gen1_Flux_Tertiary_sz/I");
   m_savetree->Branch("mc_wgt_gen1_Flux_Tertiary", &mc_wgt_gen1_Flux_Tertiary, "mc_wgt_gen1_Flux_Tertiary[mc_wgt_gen1_Flux_Tertiary_sz]/D");
   m_savetree->Branch("mc_wgt_gen1_Flux_NA49_sz", &mc_wgt_gen1_Flux_NA49_sz, "mc_wgt_gen1_Flux_NA49_sz/I");
   m_savetree->Branch("mc_wgt_gen1_Flux_NA49", &mc_wgt_gen1_Flux_NA49, "mc_wgt_gen1_Flux_NA49[mc_wgt_gen1_Flux_NA49_sz]/D");
   m_savetree->Branch("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, "mc_wgt_Norm_sz/I");
   m_savetree->Branch("mc_wgt_Norm", &mc_wgt_Norm, "mc_wgt_Norm[mc_wgt_Norm_sz]/D");
   m_savetree->Branch("mc_wgt_ppfx1_Total_sz", &mc_wgt_ppfx1_Total_sz, "mc_wgt_ppfx1_Total_sz/I");
   m_savetree->Branch("mc_wgt_ppfx1_Total", &mc_wgt_ppfx1_Total, "mc_wgt_ppfx1_Total[mc_wgt_ppfx1_Total_sz]/D");
   m_savetree->Branch("prong_nParticles", prong_nParticles, "prong_nParticles[n_prongs]/I");
   m_savetree->Branch("prong_part_score", prong_part_score, "prong_part_score[n_prongs]/D");
   m_savetree->Branch("prong_part_mass", prong_part_mass, "prong_part_mass[n_prongs]/D");
   m_savetree->Branch("prong_part_charge", prong_part_charge, "prong_part_charge[n_prongs]/I");
   m_savetree->Branch("prong_part_pid", prong_part_pid, "prong_part_pid[n_prongs]/I");

   // m_savetree->Branch("prong_part_E", &prong_part_E, "prong_part_E");
   // m_savetree->Branch("prong_part_pos", &prong_part_pos, "prong_part_pos/");


}
#endif // #ifdef AnalysisReader_cxx
