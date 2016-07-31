#ifndef CC1P1PIANALYSIS_H
#define CC1P1PIANALYSIS_H 1

#include "AnaUtils/MinervaAnalysisTool.h"

//Forward declarations:
class IMuonUtils;
class IMinervaCoordSysTool;
class INuclearTargetTool;
class IParticleMakerTool;

//Test to see if I can forward declare TString if it is being used in a function -- It worked!!
class TString;

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
    
    mutable SmartRef<Minerva::Prong> m_ProtonProng;
    mutable SmartRef<Minerva::Particle> m_ProtonParticle;
    
    mutable SmartRef<Minerva::Prong> m_PionProng;
    mutable SmartRef<Minerva::Particle> m_PionParticle;
    
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
    
    int m_Proton_PDG;// = 2212;//proton
    int m_Pion_PDG;// = 211;//pi+
    
    //Generic Particle information builder:
    void SetPartInfo(std::string name);
    
    void FillPartInfo(SmartRef<Minerva::Prong> prong, SmartRef<Minerva::Particle> particle);
    
};

#endif

