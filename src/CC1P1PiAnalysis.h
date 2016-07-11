#ifndef CC1P1PIANALYSIS_H
#define CC1P1PIANALYSIS_H 1


#include "AnaUtils/MinervaAnalysisTool.h"

// ************** Forward Declarations ************** //

class CC1P1PiAnalysis : public MinervaAnalysisTool{

public:
    
    CC1P1PiAnalysis( const std::string& type, const std::string& name, const IInterface* parent );
    ~CC1P1PiAnalysis(){};
    
    StatusCode initialize();
    
    StatusCode finalize();
    
    StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const;
    
    StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInts ) const;
    
    StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;
    
private:
    bool truthIsPlausible( const Minerva::PhysicsEvent * event ) const;

    
    
};

#endif