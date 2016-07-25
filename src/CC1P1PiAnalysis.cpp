#include "CC1P1PiAnalysis.h"

// ************** Forward Declared Headers ************** //

// ************** Other Headers ************** //

DECLARE_TOOL_FACTORY( CC1P1PiAnalysis );

//using namespace Minerva;

CC1P1PiAnalysis::CC1P1PiAnalysis( const std::string& type, const std::string& name, const IInterface* parent ) : MinervaAnalysisTool( type, name, parent ) {
    
    // ************** Mandatory Code: ************** //
    declareInterface<IInteractionHypothesis>(this);
    m_anaSignature = "CC1P1PiAnalysis";
    m_hypMeths.push_back( "CC1P1Pi" );
    declareProperty("HypothesisMethods", m_hypMeths);
    // ********************************************* //
}

StatusCode CC1P1PiAnalysis::initialize(){
    info()<<"CC1P1PiAnalysis::initialize()"<<endmsg;

    StatusCode sc = this->MinervaAnalysisTool::initialize();
    
    if( sc.isFailure() ){
        return Error( "Failed to finalize!", sc );
    }
    
    
    return sc;
}

StatusCode CC1P1PiAnalysis::finalize(){
    info()<<"CC1P1PiAnalysis::finalize()"<<endmsg;
    
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    
    if( sc.isFailure() ){
        return Error( "Failed to finalize!", sc );
    }
    
    return sc;
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth ) const{
    
    debug() << "***** CHIPS AND CHEESE ********* CC1P1PiAnalysis::reconstructEvent : Called. ************* " << endmsg;
    debug() << gateData() << "Event Number: " << event->physicsEventNumber() << ", no. of slices " << event->sliceNumbers() << endmsg;
    debug() << "Event time: " << event->time() / 1000.0 << " microseconds" << endmsg;

    //*********** 1 : Find vertex              ***********//
    debug() << "1) Find vertex" << endmsg;
    
    if(!event->hasInteractionVertex()){
        debug() << "No event vertex. Quitting..." << endmsg;
        return StatusCode::SUCCESS;
    }

    //*********** 2 : Vertex has only 3 tracks ***********//
    //Only want a total of three outgoing tracks therefore total number of
    //tracks is equal to no. of outgoing tracks.
    debug() << "2) Three tracks" << endmsg;
    SmartRef<Minerva::Vertex> reco_vertex = event->interactionVertex();

    unsigned int ntot_tracks = reco_vertex->getNTracks();
    unsigned int nout_tracks = reco_vertex->getNOutgoingTracks();
    
    if(!(ntot_tracks == nout_tracks && ntot_tracks == 3)){
        debug() << "Event doesn't contain extactly three tracks." << endmsg;
        return StatusCode::SUCCESS;
    }

    //*********** 3 : Vertex in active tracker or carbon target ***********//

    
    //*********** 4 : Muon track coming from common vertex ***********//
    
    //*********** 5 : PID on p/pi+ ***********//

    
    
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

bool CC1P1PiAnalysis::truthIsPlausible( const Minerva::PhysicsEvent * ) const
{
    // truthIsPlausible is now a pure virtual function of MinervaAnalysisTool, so you must implement it.
    // It is called automatically by PhysicsEventAnalysisAlg AFTER reconstructEvent() and interpretEvent() are run.
    // If it returns false, the event is not included in the analysis ntuple
    throw MinervaException( "You need to implement truthIsPlausible" );
    return false;
}