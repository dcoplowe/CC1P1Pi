//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 20 06:36:15 2017 by ROOT version 5.34/05
// from TTree sel/Tuple created by an AnaTuple managed by AnaTupleManager
// found on file: CC1P1Pi_R13200_190317_5/grid/central_value/minerva/ana/v10r8p9/00/01/32/00/SIM_minerva_00013200_Subruns_0001-0002-0003-0004_CC1P1PiAnalysis_Ana_Tuple_v10r8p9-dcoplowe.root
//////////////////////////////////////////////////////////

#ifndef AnalysisReader_h
#define AnalysisReader_h

#include <TROOT.h>

// Forward declarations:
class TTree;
class TChain;
class TBranch;

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
   void FillOutTree();

 private:

TTree * m_savetree;
};

#endif

