#include "CC1P1PiAnalysis.h"

// ************** Forward Declared Headers ************** //



// ************** Other Headers ************** //



ECLARE_TOOL_FACTORY( CC1P1PiAnalysis );

using namespace Minerva;

CC1P1PiAnalysis::CC1P1PiAnalysis( const std::string& type, const std::string& name, const IInterface* parent ) : MinervaAnalysisTool( type, name, parent ){
    
    // ************** Mandatory Code: ************** //
    declareInterface<IInteractionHypothesis>(this);
    m_anaSignature = "numuCC1P1PiAnalysis";
    m_hypMeths.push_back( "numuCC1P1Pi" );
    declareProperty("HypothesisMethods", m_hypMeths);
    // ********************************************* //

    
    
}

StatusCode CC1P1PiAnalysis::initialize(){
    
}

StatusCode CC1P1PiAnalysis::finalize(){
    
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const{
    
}

StatusCode CC1P1PiAnalysis::interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInts ) const{
    
}

StatusCode CC1P1PiAnalysis::tagTruth( Minerva::GenMinInteraction* truth ) const{
    
}