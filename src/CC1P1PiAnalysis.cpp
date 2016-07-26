//include the header file...
#include "CC1P1PiAnalysis.h"

#include "Event/TG4PrimaryTrajectory.h"
#include <AnaUtils/IMuonUtils.h>

//Forward declared headers:
#include "AnaUtils/IMuonUtils.h"

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
    declareProperty("MuonUtilsAlias", m_muonUtilsAlias = "CC1P1PiMuonUtils");
    declareProperty("MinMuonScore", m_minMuonScore = 0.9);
}

//! Initialize
StatusCode CC1P1PiAnalysis::initialize()
{
    debug() << "CC1P1PiAnalysis::initialize()" << endmsg;
    
    // Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() )
        return Error( "Failed to initialize!", sc );
    
    
    try{ m_muonUtils = tool<IMuonUtils>("MuonUtils", m_muonUtilsAlias); }
    catch( GaudiException& e){
        error() << "Could not find MuonUtils with alias: " << m_muonUtilsAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    //---------------------------------------------------------------------
    // Declare the Interpretations block branches
    //---------------------------------------------------------------------
    
    // The size of the iso_blob_energy branch is controlled by the n_iso_blobs branch, which you should not manually create.
    declareContainerDoubleBranch( m_hypMeths, "iso_blob_energy", "n_iso_blobs" ); // Inherited from Template
    // variable-sized containerDouble branches for each hypMeth.
    
    // Another example of variable length container branch using the same sizer: n_iso_blobs
    declareContainerIntBranch( m_hypMeths, "iso_blob_nclusters", "n_iso_blobs" ); // Inherited from Template
    
    //---------------------------------------------------------------------
    // Declare recon vars
    //---------------------------------------------------------------------

    declareIntEventBranch( "n_tracks3", -999);
    declareIntEventBranch( "vert_exists", -999);
    //  declareContainerDoubleEventBranch( "shower_momentum", 4, -999. );
    declareBoolEventBranch( "isMinosMatchTrack");
    declareBoolEventBranch( "isMinosMatchStub");
    
    //---------------------------------------------------------------------
    // Declare the Truth block branches.
    // Truth branches contain information matched to a GenMinInteraction
    //---------------------------------------------------------------------
    
    declareBoolTruthBranch("reco_isMinosMatch");
    declareIntTruthBranch( "should_be_accepted", 0 ); // Inherited from Template
    
    return sc;
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{
    debug() << "CC1P1PiAnalysis::reconstructEvent" << endmsg;
    
    //--------------------------------------------------------------
    // Initialize truth reco booleans
    //--------------------------------------------------------------
    if (truth) {
        truth->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatch", false );
    }
    
    //--------------------------------------------------------------
    // Initialize reco booleans
    //--------------------------------------------------------------
    event->filtertaglist()->setOrAddFilterTag( "isMinosMatchTrack", false );
    event->filtertaglist()->setOrAddFilterTag( "isMinosMatchStub", false );
    
    //----------- 1 : Find vertex              -----------//
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
    
    //----------- 2 : Vertex has only 3 tracks -----------//
    //Only want a total of three outgoing tracks therefore total number of
    //tracks is equal to no. of outgoing tracks.
    debug()<< "2) Three tracks" << endmsg;
    SmartRef<Minerva::Vertex> reco_vertex = event->interactionVertex();
    
    if( !reco_vertex ) {
        bool pass = true; std::string tag = "BadObject";
        event->filtertaglist()->addFilterTag(tag,pass);
        error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    unsigned int ntot_tracks = reco_vertex->getNTracks();
    unsigned int nout_tracks = reco_vertex->getNOutgoingTracks();
    
    if(!(ntot_tracks == nout_tracks && ntot_tracks == 3)){
        debug() << "Event doesn't contain extactly three tracks." << endmsg;
        return StatusCode::SUCCESS;
    }
    
    debug()<< "Has 3 tracks!" << endmsg;
    
    event->setIntData( "n_tracks3", 3);
    
    //----------- 3 : Muon track coming from common vertex -----------//
    debug()<< "3) Muon Track" << endmsg;
    SmartRef<Minerva::Prong>    muonProng = (Minerva::Prong*)NULL;
    SmartRef<Minerva::Particle> muonPart = (Minerva::Particle*)NULL;
    
    if(!FindMuon(event, truth, muonProng, muonPart )){
        debug() << "Muon not found..." << endmsg;
        return StatusCode::SUCCESS;
    }
    debug()<< "Muon track found!" << endmsg;
    
    //----------- 4 : Vertex in active tracker or carbon target -----------//
    
    //----------- 5 : PID on p/pi+ -----------//

    
    
    
    
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


//Selection Functions:

bool CC1P1PiAnalysis::FindMuon(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth, SmartRef<Minerva::Prong>& muonProng, SmartRef<Minerva::Particle>& muonPart ) const
{
    
    bool is_minos_track = false;
    bool is_minos_stub = false;
    
    if(m_muonUtils->findMuonProng(event, muonProng, muonPart)){
        
        if ( !muonProng ) {
            warning() << "Identified a muon Prong, but it is NULL!" << endmsg;
            return false; // We sort of did crash...
        }
        
       /* double mc_frac = -1.0;
        if ( m_doPlausibilityCuts && !muonIsPlausible( muonProng, mc_frac) ) {
            debug()<<"Muon is not plausible"<<endmsg;
            return false;
        }*/
        
        debug() << " Muon Particle Score: " << muonPart->score() << endmsg;
        if (muonPart->score() >= m_minMuonScore) {
            
            muonProng->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
            muonPart->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
            
            if (muonProng->MinosTrack()) is_minos_track = true;
            if (muonProng->MinosStub()) is_minos_stub = true;
            
            if (is_minos_stub && is_minos_track) counter("MuonHasMinosStubAndTrack")++;
            else counter("MuonHasMinosStubAndTrack")+=0;
            if (!is_minos_stub && !is_minos_track) counter("MuonIsNotMinosMatched")++;
            else counter("MuonIsNotMinosMatched")+=0;
            
            event->filtertaglist()->setOrAddFilterTag("isMinosMatchTrack", is_minos_track );
            event->filtertaglist()->setOrAddFilterTag("isMinosMatchStub", is_minos_stub );
            if (truth) truth->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatch", true );
        }
        else {
            debug()<<"Muon prong does not pass score cut"<<endmsg;
            return false;
        }
        
    } 
    else {
        debug() << "Did not find a muon prong!" << endmsg;
        return false;
    }
    
    return true;
}
