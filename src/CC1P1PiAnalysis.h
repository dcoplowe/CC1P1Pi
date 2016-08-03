#ifndef CC1P1PIANALYSIS_H
#define CC1P1PIANALYSIS_H 1

#include "AnaUtils/MinervaAnalysisTool.h"

//Forward declarations:
class IMuonUtils;
class IMinervaCoordSysTool;
class INuclearTargetTool;
class IParticleMakerTool;
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
    
    mutable SmartRef<Minerva::Prong> m_ProtonProng;
    mutable SmartRef<Minerva::Particle> m_ProtonParticle;
    double * m_Proton4Mom;
    double * m_Protontrue4Mom;
    
    mutable SmartRef<Minerva::Prong> m_PionProng;
    mutable SmartRef<Minerva::Particle> m_PionParticle;
    double * m_Pion4Mom;
    double * m_Piontrue4Mom;
    
    void ResetParticles() const;
    //----------------------------------------------------------------------------//


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
    bool FindParticles(Minerva::PhysicsEvent* event) const;
    double m_det_apothem;
    double m_det_upZ;
    double m_det_downZ;
    
    IParticleMakerTool * m_particleMaker;
    std::string          m_particleMakerAlias;
    
    double m_minProtonScore;
    double m_maxProtonChi2;
    
    double m_minPionScore;
    double m_maxPionChi2;
    
    double * m_ProtonScore; //[0] - Proton score, [1] - Pion score
    double * m_PionScore; //[0] - Proton score, [1] - Pion score
    double * m_Chi2NDF; //[0] - Proton Chi2, [1] - Pion Chi2
    
    int m_Proton_PDG;// = 2212;//proton
    int m_Pion_PDG;// = 211;//pi+
    
    //Generic Particle information builder:
    void SetPartInfo(std::string name);
    
    void FillPartInfo(std::string name, const Minerva::GenMinInteraction *truth, Minerva::NeutrinoInt* cc1p1piHyp) const;
    
    void SetCommonBranches();// const;
    
    void FillCommonBranches(const Minerva::GenMinInteraction *truth, Minerva::NeutrinoInt* cc1p1piHyp) const;
    
    void Rotate2BeamCoords(std::vector<double> val) const;

    //Accumulation level counter:
    int m_ncuts;// = 5;
    int * m_accum_level;
    void SetAccumLevel(int cut) const;
    void ResetAccumLevel() const;
    void SaveAccumLevel(Minerva::PhysicsEvent * event) const;
    
    //Determine the truth information from the track:
    ITruthMatcher * m_truthMatcher;
    
    //Set four vecs for the final state 'name' particle:
    void SetGlobal4Vec(std::string name, Gaudi::LorentzVector vec, bool truth = false) const;
    
    //Transverse variables:
    TVector3 * GetTransverseVars(double &vtx[3], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT, bool is_truth = false) const;
    //TVector3 * GetPT(double vtx[], TVector3 mom, bool is_truth = false) const;
    //void SetDPT(TVector3 * deltapt, TVector3 * ptmuon, TVector3 * ptproton, TVector3 * ptpion) const;
    //static const TVector3 * GetVecT(const TVector3 * refdir, const TVector3 * mom);
    //double GetDPTT(double vtx[], TVector3 * mumom, TVector3 * prmom, TVector3 * pimom, bool is_truth = false) const;
    
};

#endif

