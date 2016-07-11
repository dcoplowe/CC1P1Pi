#include "CC1P1PiAnalysis.h"

// ************** Forward Declared Headers ************** //



// ************** Other Headers ************** //



DECLARE_TOOL_FACTORY( CC1P1PiAnalysis );

using namespace Minerva;

CC1P1PiAnalysis::CC1P1PiAnalysis( const std::string& type, const std::string& name, const IInterface* parent ) : MinervaAnalysisTool( type, name, parent ) {
    
    // ************** Mandatory Code: ************** //
    declareInterface<IInteractionHypothesis>(this);
    m_anaSignature = "numuCC1P1PiAnalysis";
    m_hypMeths.push_back( "numuCC1P1Pi" );
    declareProperty("HypothesisMethods", m_hypMeths);
    // ********************************************* //

    
    
}

StatusCode CC1P1PiAnalysis::initialize(){
    
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    
    if( sc.isFailure() ){
        return Error( "Failed to finalize!", sc );
    }
    
    
    return sc;
}

StatusCode CC1P1PiAnalysis::finalize(){
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    
    if( sc.isFailure() ){
        return Error( "Failed to finalize!", sc );
    }
    
    return sc;
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth ) const{
    
    debug() << "CC1P1PiAnalysis::reconstructEvent : Called." << endmsg;
    
    
    
    markEvent( event );
    // Now interpret the event and add NeutrinoInts
    std::vector<Minerva::NeutrinoInt*> nuInts;
    interpretEvent( event, truth, nuInts );
    
    // Add the newly create NeutrinoInts to this PhysicsEvent
    StatusCode sc = addInteractionHyp( event, nuInts );
    
    return sc;
    
}

StatusCode CC1P1PiAnalysis::interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInts ) const{
    
    (void)event;
    (void)truth;
    (void)nuInts;
    
    debug() << "Exit interpretEvent()" << endmsg;
    
    return StatusCode::SUCCESS;
}

StatusCode CC1P1PiAnalysis::tagTruth( Minerva::GenMinInteraction* truth ) const{
    truth->setIntData("should_be_accepted", 1);
    
    return StatusCode::SUCCESS;
}