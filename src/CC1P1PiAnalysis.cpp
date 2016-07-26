//include the header file...
#include "CC1P1PiAnalysis.h"

#include "Event/TG4PrimaryTrajectory.h"
#include <AnaUtils/IMuonUtils.h>

//this command allows other parts of Gaudi to use the tool
DECLARE_TOOL_FACTORY( CC1P1PiAnalysis );

CC1P1PiAnalysis::CC1P1PiAnalysis(const std::string& type, const std::string& name, const IInterface* parent ) : MinervaAnalysisTool( type, name, parent )
{
    
    declareInterface<IInteractionHypothesis>(this);
    
    m_anaSignature = "CC1P1Pi";
    
    m_hypMeths.push_back( "CC1P1Pi" );
    //m_hypMeths.push_back( "InterpretationB" );
    declareProperty("HypothesisMethods", m_hypMeths);
    
    // Declare other properties you can set from an options file.
    declareProperty( "SomeProperty", m_someProperty = 1 );
    
}

//! Initialize
StatusCode CC1P1PiAnalysis::initialize()
{
    debug() << "DAVID : : CC1P1PiAnalysis::initialize()" << endmsg;
    
    // Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() )
        return Error( "Failed to initialize!", sc );
    
    //---------------------------------------------------------------------
    // Declare recon vars
    //---------------------------------------------------------------------
    //declareDoubleEventBranch( "time_width", -999.0 ); // Inherited from Template
    //declareIntEventBranch( "n_orig_prongs", -1 ); // Inherited from Template
    
    declareIntEventBranch( "n_tracks3", -999);
    declareIntEventBranch( "vert_exists", -999);
    //  declareContainerDoubleEventBranch( "shower_momentum", 4, -999. );
    
    //---------------------------------------------------------------------
    // Declare the Interpretations block branches
    //---------------------------------------------------------------------
    
    // The size of the iso_blob_energy branch is controlled by the n_iso_blobs branch, which you should not manually create.
    declareContainerDoubleBranch( m_hypMeths, "iso_blob_energy", "n_iso_blobs" ); // Inherited from Template
    // variable-sized containerDouble branches for each hypMeth.
    
    // Another example of variable length container branch using the same sizer: n_iso_blobs
    declareContainerIntBranch( m_hypMeths, "iso_blob_nclusters", "n_iso_blobs" ); // Inherited from Template
    
    //---------------------------------------------------------------------
    // Declare the Truth block branches.
    // Truth branches contain information matched to a GenMinInteraction
    //---------------------------------------------------------------------
    declareIntTruthBranch( "should_be_accepted", 0 ); // Inherited from Template
    
    return sc;
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{
    debug() << "CC1P1PiAnalysis::reconstructEvent" << endmsg;
    
    //*********** 1 : Find vertex              ***********//
    debug()<< "1) Find vertex" << endmsg;
    
    if( !event->hasInteractionVertex() ){
        debug() << "No event vertex. Quitting..." << endmsg;
        //event->setIntData("vert_exists", 0);
        return StatusCode::SUCCESS;
    }
    //else{
        debug() << "Found vertex!" << endmsg;
        event->setIntData("vert_exists", 1);
    //}
    
    //*********** 2 : Vertex has only 3 tracks ***********//
    //Only want a total of three outgoing tracks therefore total number of
    //tracks is equal to no. of outgoing tracks.
    debug()<< "2) Three tracks" << endmsg;
    SmartRef<Minerva::Vertex> reco_vertex = event->interactionVertex();
    
    unsigned int ntot_tracks = reco_vertex->getNTracks();
    unsigned int nout_tracks = reco_vertex->getNOutgoingTracks();
    
    if(!(ntot_tracks == nout_tracks && ntot_tracks == 3)){
        debug() << "Event doesn't contain extactly three tracks." << endmsg;
        return StatusCode::SUCCESS;
    }
    
    debug()<< "Has 3 tracks!" << endmsg;
    
    event->setIntData( "n_tracks3", 3);
    
    //*********** 3 : Vertex in active tracker or carbon target ***********//
    
    
    //*********** 4 : Muon track coming from common vertex ***********//
    
    //*********** 5 : PID on p/pi+ ***********//

    
    
    
    
    // Set the PhysicsEvent reconstructionSignature to m_anaSignature, so I know that this tool reconstructed this event.
    // If you mark the event it will go to your analysis DST.  If you don't want it to go there, don't mark it!
    markEvent( event );
    
    
    // Now interpret the event and add NeutrinoInts
    std::vector<Minerva::NeutrinoInt*> nuInts;
    interpretEvent( event, truth, nuInts );
    
    // You can also use other analysis tools to interpret the event.
    //m_ccInclusive->interpretEvent( event, interaction, nuInts );
    
    // Add the newly create NeutrinoInts to this PhysicsEvent
    StatusCode sc = addInteractionHyp( event, nuInts );
    
    return sc;
}

StatusCode CC1P1PiAnalysis::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *interaction, std::vector<Minerva::NeutrinoInt*>& nuInts ) const
{
    debug() << "CC1P1PiAnalysis::interpretEvent" << endmsg;
    
    counter("N Primary Prongs") += event->primaryProngs().size();
    // If you decide you want to interpret the event, create a new NeutrinoInt.
    Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( "InterpretationA" );
    
    // Add the NeutrinoInt to the vector in return value
    nuInts.push_back( nuInt );
    
    // Set the NeutrinoInt data
    nuInt->setEnergy( 5.5 );
    
    // Add your own extra data to the NeutrinoInt.  Remember that you declared branches so that extra data could be written to an analysis DST...
    std::vector<double> isoEnergyVec;
    isoEnergyVec.push_back( 1.2 );
    isoEnergyVec.push_back( 5.5 );
    isoEnergyVec.push_back( 3.6 );
    nuInt->setContainerDoubleData( "iso_blob_energy", isoEnergyVec );
    
    
    if( interaction )
    {
        debug() << " There's an interaction so make some truth plots." << endmsg;
    }
    
    return StatusCode::SUCCESS;
}

//! Attach information to the GenMinInteraction
StatusCode CC1P1PiAnalysis::tagTruth( Minerva::GenMinInteraction *truth ) const
{
    debug() << "CC1P1PiAnalysis::tagTruth" << endmsg;
    
    truth->setIntData("should_be_accepted", 1);
    
    return StatusCode::SUCCESS;
}


//! Finalize
StatusCode CC1P1PiAnalysis::finalize()
{
    debug() << "CC1P1PiAnalysis::finalize()" << endmsg;
    
    // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() )
        return Error( "Failed to finalize!", sc );
    
    return sc;
    
}

//! Implement the pure virtual function truthIsPlausible (inherited from IInteractionHypothesis):
//! did this PhysicsEvent come (primarily) from the MC?  (See DocDB 10471.)
bool CC1P1PiAnalysis::truthIsPlausible( const Minerva::PhysicsEvent * event ) const
{
    
    debug() << "CC1P1PiAnalysis::truthIsPlausible" << endmsg;

    // Here you need to decide whether the things that you require for the event
    // to pass your signal selection were made up primarily of MC.
    // SEE DOCDB 10471 IF YOU ARE UNSURE HOW TO IMPLEMENT THIS METHOD.
    
    // in a MINOS-matched-muon analysis, for example, you usually just want the muon to be plausible
    SmartRef<Minerva::Prong> muonProng;
    SmartRef<Minerva::Particle> muonPart;
    if ( ! MuonUtils->findMuonProng(event, muonProng, muonPart) )  // returns false if it can't find a muon
        return false;
    return muonIsPlausible(muonProng);
}

