#ifndef CC1P1PIANALYSIS_H
#define CC1P1PIANALYSIS_H 1

#include "AnaUtils/MinervaAnalysisTool.h"

//Forward declarations:
class ICCPionIncUtils;
class IMuonUtils;
class IMinervaCoordSysTool;
class INuclearTargetTool;
class IProtonUtils;
//class IParticleMakerTool;
class IParticleTool;//This is for LL PID
class ITruthMatcher;
class TString;
class TVector3;

class CC1P1PiAnalysis : public MinervaAnalysisTool
{
public:
    
    CC1P1PiAnalysis( const std::string& type, const std::string& name, const IInterface* parent );
    
    ~CC1P1PiAnalysis(){};
    
    StatusCode initialize();

    StatusCode finalize();
    
    StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const;
    
    StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInts ) const;
    
    StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;
    
    bool truthIsPlausible( const Minerva::PhysicsEvent * event ) const;
    
protected:
    
    
private:
    
    //---------------Particles being passed around for analysis-------------------//
    //Particles being passed around for analysis:
    //Making them mutable objects enables them to be changed in const functions -- this coding practice is taken from the CCProtonPi0 analysis
    mutable SmartRef<Minerva::Prong> m_MuonProng;
    mutable SmartRef<Minerva::Particle> m_MuonParticle;
    double * m_Muon4Mom;
    double * m_Muontrue4Mom;
    
    mutable SmartRef<Minerva::Prong> m_EX_ProtonProng;
    mutable SmartRef<Minerva::Particle> m_EX_ProtonParticle;
    mutable SmartRef<Minerva::Particle> m_EX_ProtonParticle_AltH;
    double * m_EX_Proton4Mom;

    mutable SmartRef<Minerva::Prong> m_EX_PionProng;
    mutable SmartRef<Minerva::Particle> m_EX_PionParticle;
    mutable SmartRef<Minerva::Particle> m_EX_PionParticle_AltH;
    double * m_EX_Pion4Mom;
    
    mutable SmartRef<Minerva::Prong> m_LL_ProtonProng;
    mutable SmartRef<Minerva::Particle> m_LL_ProtonParticle;
    mutable SmartRef<Minerva::Particle> m_LL_ProtonParticle_AltH;
    double * m_LL_Proton4Mom;
    
    mutable SmartRef<Minerva::Prong> m_LL_PionProng;
    mutable SmartRef<Minerva::Particle> m_LL_PionParticle;
    mutable SmartRef<Minerva::Particle> m_LL_PionParticle_AltH;
    double * m_LL_Pion4Mom;

    double * m_Protontrue4Mom;
    double * m_Piontrue4Mom;
    
    double * m_Muon_dir;
    double * m_Proton_dir;
    double * m_Pion_dir;
    
    double * m_Muontrue_dir;
    double * m_Protontrue_dir;
    double * m_Piontrue_dir;
    
    void ResetParticles() const;
    //----------------------------------------------------------------------------//

    //------ Find short tracks originating from the vertex ------//
    ICCPionIncUtils * m_ccPionIncUtils;
    void FindShortTracks(Minerva::PhysicsEvent * event) const;

    //------ Find Muon ------//
    bool FindMuon(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth, SmartRef<Minerva::Prong>& muonProng, SmartRef<Minerva::Particle>& muonPart) const;
    IMuonUtils * m_muonUtils;
    std::string  m_muonUtilsAlias; //What are these for?
    double       m_minMuonScore;
    //-----------------------//
    
    //------ Determine if vertex is in FV of either scintillator or carbon ------//
    bool VertIsIn(TString targetRegion, Minerva::PhysicsEvent* event) const;
    IMinervaCoordSysTool * m_coordSysTool;
    INuclearTargetTool *   m_nuclearTargetTool;
    std::string m_nuclearTargetToolAlias;
    
    //Default is the full detector fiducial volume
    double m_default_apothem;
    double m_default_upZ;
    double m_default_downZ;
    
    //Scintillator fiducial volume taken from CCQE Analysis
    double m_scint_apothem;
    double m_scint_upZ;
    double m_scint_downZ;
    
    //Carbon fiducial volume taken from XXX 3rd target info
    double m_carbon_apothem;
    double m_carbon_upZ;
    double m_carbon_downZ;
    //-----------------------//
    
    //Find Proton/Pion:
    bool FindParticles(Minerva::PhysicsEvent * event, std::string method = "LL") const;
    
    bool IsEventContained(Minerva::PhysicsEvent * event) const;
    
    double m_det_apothem;
    double m_det_upZ;
    double m_det_downZ;
    
    //dEdX method:
    bool EXMethod(Minerva::PhysicsEvent * event) const;
    IProtonUtils * m_protonUtils;
    //IParticleMakerTool * m_particleMaker;
    //std::string          m_particleMakerAlias;
    
    //Likelihood method:
    bool LLMethod(Minerva::PhysicsEvent * event) const;
    IParticleTool * m_LikelihoodPIDTool;
    
    double m_minProtonScore;
    double m_maxProtonChi2;
    
    double m_minPionScore;
    double m_maxPionChi2;
    
    double * m_ProtonScore; //[0] - Proton score, [1] - Pion score
    double * m_ProtonChi2ndf;
    double * m_PionScore; //[0] - Proton score, [1] - Pion score
    double * m_PionChi2ndf; //[0] - Proton Chi2, [1] - Pion Chi2
    
    int m_Proton_PDG;// = 2212;//proton
    int m_Pion_PDG;// = 211;//pi+
    
    //Generic Particle information builder:
    void SetPartInfo(std::string name);
    
    void FillPartInfo(std::string name, const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, Minerva::NeutrinoInt* cc1p1piHyp) const;
    
    void SetCommonBranches();// const;
    
    void FillCommonBranches(const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, Minerva::NeutrinoInt* cc1p1piHyp) const;

    void FillMomDepVars(std::string name, SmartRef<Minerva::Particle> particle, double mass, const Minerva::PhysicsEvent *event, Minerva::NeutrinoInt* cc1p1piHyp) const;
    
    void DefineTruthTree();
    void FillTruthTree(Minerva::GenMinInteraction* truth) const;
    
    void Rotate2BeamCoords(std::vector<double> val) const;

    //Accumulation level counter:
    int m_ncuts;// = 5;
    int m_nsplits;
    mutable int * m_accum_level;
    int m_accum_level_to_save;
    void SetAccumLevel(int split = -999) const;
    void ResetAccumLevel() const;
    void SaveAccumLevel(Minerva::PhysicsEvent * event, Minerva::GenMinInteraction* truth) const;
    
    int m_PID_method;
    std::string m_PID_tool;
    
    //Determine the truth information from the track:
    ITruthMatcher * m_truthMatcher;
    
    //Set four vecs for the final state 'name' particle:
    void SetGlobal4Vec(std::string name, Gaudi::LorentzVector vec, bool truth = false) const;
    void SetGlobal4Vec(std::string name, std::vector<double> vec, bool truth = false) const;
    
    
    //Transverse variables:
    //Parent Decay Point Var:
    double m_PDP_X;
    double m_PDP_Y;
    double m_PDP_Z;
    
    TVector3 * m_PDP;
    TVector3 * GetNuDirRec(double vtx[]) const;
    
    TVector3 * GetTransverseVars(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT, bool is_truth = false) const;
    TVector3 * GetPT(double vtx[], const TVector3 *& mom, bool is_truth = false) const;
    void SetDPT(TVector3 * deltapt, const TVector3 *& ptmuon, const TVector3 *& ptproton, const TVector3 *& ptpion) const;
    TVector3 * GetVecT(const TVector3 *& refdir, const TVector3 *& mom) const;
    double GetDPTT(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, bool is_truth = false) const;
    
    void PrintInfo(std::string var, bool print = true) const;
    bool m_print_acc_level;
    bool m_print_cuts;
    bool m_print_other;
    bool m_print_cut_verbose;
    
    void EventFinished() const;
};

#endif

