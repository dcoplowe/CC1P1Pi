//include the header file...
#include "CC1P1PiAnalysis.h"

#include "Event/TG4PrimaryTrajectory.h"

//For transverse calcs:
//#include "MINERVAUtils.h"

//Forward declared headers:
#include "AnaUtils/IMuonUtils.h"
#include "GeoUtils/IMinervaCoordSysTool.h"
#include "DetDesc/Material.h"
#include "GeoUtils/INuclearTargetTool.h"
#include "ParticleMaker/IParticleMakerTool.h"
#include "TruthMatcher/ITruthMatcher.h"

//Root headers:
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>

//#include <stdio.h>
//#include <stdlib.h>

#ifndef EPSILON
#define EPSILON  1e-10
#endif

//this command allows other parts of Gaudi to use the tool
DECLARE_TOOL_FACTORY( CC1P1PiAnalysis );

CC1P1PiAnalysis::CC1P1PiAnalysis(const std::string& type, const std::string& name, const IInterface* parent ) : MinervaAnalysisTool( type, name, parent )
{
    
    declareInterface<IInteractionHypothesis>(this);
    
    m_anaSignature = "CC1P1Pi";
    
    m_hypMeths.push_back( "CC1P1Pi" );
    //m_hypMeths.push_back( "InterpretationB" );
    declareProperty("HypothesisMethods", m_hypMeths);
    
    //For muon PID:
    // Declare other properties you can set from an options file.
    declareProperty("MuonUtilsAlias", m_muonUtilsAlias = "CC1P1PiMuonUtils");
    declareProperty("MinMuonScore",   m_minMuonScore = 0.9);
    
    //For Nuclear Target/Scintillator ID:
    //Taken from CCQE Two Track code
    declareProperty("default_apothem", m_default_apothem = 850.0*CLHEP::mm); //Taken from CCQETwoTrack
    declareProperty("default_upZ",    m_default_upZ = 4284.46*CLHEP::mm);//Taken from NukeCCQETwoTrack
    declareProperty("default_downZ",  m_default_downZ = 8300.00*CLHEP::mm);//Taken from CCQETwoTrack
    
    declareProperty("scint_apothem", m_scint_apothem = 850.0*CLHEP::mm);//Taken from CCQETwoTrack
    declareProperty("scint_upZ" ,    m_scint_upZ = 5990.00*CLHEP::mm);//Taken from CCQETwoTrack
    declareProperty("scint_downZ" ,  m_scint_downZ = 8300.00*CLHEP::mm);//Taken from CCQETwoTrack
    
    //Taken from NukeCQETwoTrack - taken as front/back face nuclear target
    //May want to change this to Target3 vertex cut in NukeCCQETwoTrack
    declareProperty("carbon_apothem", m_carbon_apothem = 900.0*CLHEP::mm);
    declareProperty("carbon_upZ",     m_carbon_upZ = 4815.04*CLHEP::mm);
    declareProperty("carbon_downZ",   m_carbon_downZ = 5073.59*CLHEP::mm);

    declareProperty("NuclearTargetToolAlias", m_nuclearTargetToolAlias  = "CC1P1PiTargetTool");
    
    //For hadron PID:
    declareProperty("det_apothem", m_det_apothem = 1200.0*CLHEP::mm);//Same as Proton utils:
    declareProperty("det_upZ", m_det_upZ = 4000.0*CLHEP::mm);//Same as Proton utils:
    declareProperty("det_downZ", m_det_downZ = 10000.0*CLHEP::mm);//Same as Proton utils:
    
    declareProperty("ParticleMakerAlias", m_particleMakerAlias = "CC1P1PiParticleMaker");
    
    //These values are taken from ProtonUtils:
    declareProperty("minProtonScore", m_minProtonScore = 0.05);
    declareProperty("maxProtonChi2",  m_maxProtonChi2  = 50.);
    
    //Need to determine better scores (currently same as those in ProtonUtils for Proton Hyp:
    declareProperty("minPionScore", m_minPionScore = 0.05);
    declareProperty("maxPionChi2",  m_maxPionChi2  = 50.);

    declareProperty("Proton_PDG", m_Proton_PDG = 2212);//proton
    declareProperty("Pion_PDG", m_Pion_PDG = 211);//pi+
    
    declareProperty("n_cuts", m_ncuts = 5);
    m_accum_level = 0;// = new int [ m_ncuts ]; -- no need for N dim. as we only have one branch.

    //Mean Parent Decay Point
    declareProperty("PDP_X", m_PDP_X = 0.231135);
    declareProperty("PDP_Y", m_PDP_Y = 45.368069);
    declareProperty("PDP_Z", m_PDP_Z = 766.384058);
    
    declareProperty("accum_level_to_save", m_accum_level_to_save = 5);//Defualt to no of cuts so that we only save interesting events.
    
    declareProperty("print_cuts", m_print_cuts = false);
    declareProperty("print_cut_verbose", m_print_cut_verbose = false);
    declareProperty("print_acc_level", m_print_acc_level = true);
    declareProperty("print_other", m_print_other = false);
    
    m_PDP = new TVector3(m_PDP_X, m_PDP_Y, m_PDP_X);
    
    m_ProtonScore = new double [2];
    m_PionScore = new double [2];
    m_Chi2NDF = new double [2];
    
    //Want to pass 4-mom vectors in order to produce Q2, Enu and transverse variables:
    m_Pion4Mom = new double [4];//TXYZ
    m_Piontrue4Mom = new double [4];
    
    m_Proton4Mom = new double [4];
    m_Protontrue4Mom = new double [4];
    
    m_Muon4Mom = new double [4];
    m_Muontrue4Mom = new double [4];
    
}

//! Initialize
StatusCode CC1P1PiAnalysis::initialize()
{
    PrintInfo("CC1P1PiAnalysis::initialize()", m_print_other);
    
    
    // Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() )
        return Error( "Failed to initialize!", sc );
    
    
    try{ m_muonUtils = tool<IMuonUtils>("MuonUtils", m_muonUtilsAlias); }
    catch( GaudiException& e){
        error() << "Could not find MuonUtils with alias: " << m_muonUtilsAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ m_coordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool"); }
    catch( GaudiException& e ) {
        error() << "Could not obtain MinervaCoordSysTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ m_nuclearTargetTool = tool<INuclearTargetTool>("NuclearTargetTool", m_nuclearTargetToolAlias);
        m_nuclearTargetTool->m_locked = false;
        m_nuclearTargetTool->addAllPassiveNuclearTargets();
        m_nuclearTargetTool->lock();
    } catch( GaudiException& e ) {
        error() << "Could not obtain NuclearTargetTool: " << m_nuclearTargetToolAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try { m_particleMaker = tool<IParticleMakerTool>("ParticleMakerTool", m_particleMakerAlias); }
    catch( GaudiException& e){
        error() << "Could not obtain ParticleMakerTool: " << m_particleMakerAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try { m_truthMatcher = tool<ITruthMatcher>("TruthMatcher"); }
    catch( GaudiException& e){
        error() << "Could not obtain TruthMather! " << endmsg;
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

    declareIntEventBranch("n_tracks3", -999);
    declareIntEventBranch("vert_exists", -999);
    declareIntEventBranch("target_region", -999);//1 - Scint, 2 - carbon, 3 - other - There shouldn't be any of these as these events will be cut.
    //  declareContainerDoubleEventBranch( "shower_momentum", 4, -999. );
    declareBoolEventBranch("isMinosMatchTrack");
    declareBoolEventBranch("isMinosMatchStub");
    
    SetCommonBranches();
    
    SetPartInfo("mu");
    SetPartInfo("pi");
    SetPartInfo("pr");

    //---------------------------------------------------------------------
    // Declare the Truth block branches.
    // Truth branches contain information matched to a GenMinInteraction
    //---------------------------------------------------------------------
    declareIntEventBranch("accum_level", 0);
    declareIntTruthBranch("accum_level", 0);
    declareIntEventBranch("ncuts", m_ncuts);
    declareIntTruthBranch("ncuts", m_ncuts);
    
    declareBoolTruthBranch("reco_isMinosMatch");
    declareIntTruthBranch("should_be_accepted", 0); // Inherited from Template
    
    return sc;
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{
    PrintInfo("CC1P1PiAnalysis::reconstructEvent", m_print_other);
    
    PrintInfo("-------------------------------------------------------------------------------------------------", m_print_acc_level);
    PrintInfo("----------------------------           New Event             ------------------------------------", m_print_acc_level);
    PrintInfo("-------------------------------------------------------------------------------------------------", m_print_acc_level);
    
    //Clear Particle Prongs and Particle objects:
    ResetParticles();
    ResetAccumLevel();

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
    PrintInfo("1) Find vertex", m_print_cuts);
    PrintInfo("AL should be 0", m_print_acc_level);
    //PrintInfo(Form("***** Accum. Level %d *****", m_accum_level), m_print_acc_level);
    PrintInfo( ("***** Accum. Level " + std::to_string(m_accum_level) + " *****").c_str(), m_print_acc_level);
    
    if( !event->hasInteractionVertex() ){
        PrintInfo("No event vertex. Quitting...", m_print_cuts);
        //event->setIntData("vert_exists", 0);
        SaveAccumLevel(event, truth);
        EventFinished();
        return StatusCode::SUCCESS;
    }
    //else{
    PrintInfo("Found vertex!", m_print_cuts);
    event->setIntData("vert_exists", 1);
    counter("c_vertex")++;
    PrintInfo("AL should be 1", m_print_acc_level);
    SetAccumLevel();
    
    //}
    
    //----------- 2 : Vertex has only 3 tracks -----------//
    //Only want a total of three outgoing tracks therefore total number of
    //tracks is equal to no. of outgoing tracks.
    PrintInfo("2) Three tracks", m_print_cuts);
    SmartRef<Minerva::Vertex> reco_vertex = event->interactionVertex();
    
    if( !reco_vertex ) {
        bool pass = true; std::string tag = "BadObject";
        event->filtertaglist()->addFilterTag(tag,pass);
        error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
        PrintInfo("AL save 1 ?", m_print_acc_level);
        SaveAccumLevel(event, truth);
        EventFinished();
        return StatusCode::SUCCESS;
    }
    
    unsigned int ntot_tracks = reco_vertex->getNTracks();
    unsigned int nout_tracks = reco_vertex->getNOutgoingTracks();
    
    //Check if there are the same number of prongs as tracks:
    Minerva::ProngVect n_prong_check = event->primaryProngs();
    unsigned int n_prongs = n_prong_check.size();
    
    PrintInfo(Form("n_tracks = %d n_prongs = %d", ntot_tracks, n_prongs), m_print_cuts);
    if(ntot_tracks == n_prongs){
        PrintInfo(" !! EQUAL !!", m_print_cuts);
    }
    //debug() << " " << endmsg;
    
    if(!(ntot_tracks == nout_tracks && ntot_tracks == 3)){
        PrintInfo("Event doesn't contain extactly three tracks.", m_print_cuts);
        PrintInfo("AL save 1 ?", m_print_acc_level);
        SaveAccumLevel(event, truth);
        EventFinished();
        return StatusCode::SUCCESS;
    }
    
    counter("c_3tracks")++;
    PrintInfo("AL should be 2", m_print_acc_level);
    SetAccumLevel();
    
    PrintInfo("Has 3 tracks!", m_print_cuts);
    
    event->setIntData("n_tracks3", 3);
    event->setIntData("n_prongs", n_prongs);
    
    //----------- 3 : Muon track coming from common vertex -----------//
    PrintInfo("3) Muon Track", m_print_cuts);
   // SmartRef<Minerva::Prong>    muonProng = (Minerva::Prong*)NULL;
   // SmartRef<Minerva::Particle> muonPart = (Minerva::Particle*)NULL;
    
    if(!FindMuon(event, truth, m_MuonProng, m_MuonParticle)){
        PrintInfo("Muon not found...", m_print_cuts);
        PrintInfo("AL save 2 ?", m_print_acc_level);
        SaveAccumLevel(event, truth);
        EventFinished();
        return StatusCode::SUCCESS;
    }
    PrintInfo("Muon track found!", m_print_cuts);
    
    counter("c_muon_trk")++;
    PrintInfo("AL should be 3", m_print_acc_level);
    SetAccumLevel();
    
    //----------- 4 : Vertex in active tracker or carbon target -----------//
    //Not a cut but an action to determine the location of the vertex.
    PrintInfo("4) Vertex in Carbon or Scintillator", m_print_cuts);
    if(VertIsIn("Scint", event)){
        PrintInfo("Yes in SCINTILLATOR", m_print_cuts);
        event->setIntData("target_region", 1);
        counter("c_tar_scint")++;
    }
    else if (VertIsIn("Carbon", event)){
        PrintInfo("Yes in CARBON TARGET", m_print_cuts);
        event->setIntData("target_region", 2);
        counter("c_tar_carbon")++;
    }
    else{
        PrintInfo("Event not in either...", m_print_cuts);
        event->setIntData("target_region", 3);//Probably don't need this...
        counter("c_tar_other")++;
        PrintInfo("AL save 3 ?", m_print_acc_level);
        SaveAccumLevel(event, truth);
        EventFinished();
        return StatusCode::SUCCESS;
    }
    
    PrintInfo("AL should be 4", m_print_acc_level);
    SetAccumLevel();
    
    //----------- 5 : PID on p/pi+ -----------//
    PrintInfo("5) PID: p/pi+", m_print_cuts);
    
    //HadronSystem hadrons;
    
    bool tFinPar = FindParticles(event);
    
    if(!tFinPar){
        PrintInfo("Failed to identify particles...", m_print_cuts);
        PrintInfo("AL save 4 ?", m_print_acc_level);
        SaveAccumLevel(event, truth);
        EventFinished();
        return StatusCode::SUCCESS;
    }
    else{
        PrintInfo("Finished Selection Successfully. Pheeewwww ;)", m_print_cuts);
    }
    
    PrintInfo("AL should be 5", m_print_acc_level);
    SetAccumLevel();
    
    SaveAccumLevel(event, truth);//markEvent is called in SaveAccumLevel.
    
    // Set the PhysicsEvent reconstructionSignature to m_anaSignature, so I know that this tool reconstructed this event.
    // If you mark the event it will go to your analysis DST.  If you don't want it to go there, don't mark it!
    //markEvent( event );
    
    
    // Now interpret the event and add NeutrinoInts
    std::vector<Minerva::NeutrinoInt*> nuInts;
    interpretEvent( event, truth, nuInts );
    
    // You can also use other analysis tools to interpret the event.
    //m_ccInclusive->interpretEvent( event, interaction, nuInts );
    
    // Add the newly create NeutrinoInts to this PhysicsEvent
    StatusCode sc = addInteractionHyp( event, nuInts );
    
    EventFinished();
    
    return sc;
}

StatusCode CC1P1PiAnalysis::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *interaction, std::vector<Minerva::NeutrinoInt*>& nuInts ) const
{
    //debug() << "CC1P1PiAnalysis::interpretEvent" << endmsg;
    
    Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( m_anaSignature );
    
    FillPartInfo("mu", event, interaction, nuInt);
    FillPartInfo("pr", event, interaction, nuInt);
    FillPartInfo("pi", event, interaction, nuInt);
    
    FillCommonBranches(event, interaction, nuInt);
    
    // Add the NeutrinoInt to the vector in return value
    nuInts.push_back( nuInt );
    
    return StatusCode::SUCCESS;
}

//! Attach information to the GenMinInteraction
StatusCode CC1P1PiAnalysis::tagTruth( Minerva::GenMinInteraction *truth ) const
{
    //debug() << "CC1P1PiAnalysis::tagTruth" << endmsg;
    
    truth->setIntData("should_be_accepted", 1);
    
    return StatusCode::SUCCESS;
}


//! Finalize
StatusCode CC1P1PiAnalysis::finalize()
{
    PrintInfo("CC1P1PiAnalysis::finalize()", m_print_other);
    
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
    
    //debug() << "CC1P1PiAnalysis::truthIsPlausible" << endmsg;

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
        
        PrintInfo(Form("Muon Particle Score: %f", muonPart->score()), m_print_cut_verbose);
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
            PrintInfo("Muon prong does not pass score cut", m_print_cut_verbose);
            return false;
        }
        
    } 
    else {
        PrintInfo("Did not find a muon prong!", m_print_cut_verbose);
        return false;
    }
    
    return true;
}

bool CC1P1PiAnalysis::VertIsIn(TString targetRegion, Minerva::PhysicsEvent* event) const
{
    //This function checks if the vertex is in the target region specified by the string and in the fiducial volume. Currently this works for only
    //carbon and scintillator but can be fixed to work with any target.
    
    double apothem = m_default_apothem;
    double upZ = m_default_upZ;
    double downZ = m_default_downZ;
    
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    
    if(targetRegion.Contains("Scint", TString::kIgnoreCase)){
        apothem = m_scint_apothem;
        upZ = m_scint_upZ;
        downZ = m_scint_downZ;
    }
    else if(targetRegion.Contains("Carbon", TString::kIgnoreCase)){
        
        const Gaudi::XYZPoint p_vert = vertex->position();
        const Material * material = m_nuclearTargetTool->getSectionMaterial(p_vert);
        int materialZ = -999;
        
        if(material){
            materialZ = (int)(material->Z()+0.5);
            debug() << "  retrieve the target's section material name = " << material->name() << ", and Z = " << material->Z() << endmsg;
        }
        //Carbon is in target region 3 and has Z == 6.
        
        if(materialZ != 6) return false;
        
        apothem = m_carbon_apothem;
        upZ = m_carbon_upZ;
        downZ = m_carbon_downZ;
    }
    else {
        PrintInfo(Form("CC1P1PiAnalysis::VertIsIn : Could not determine target name: %s. Please check, it may not be implemented.", targetRegion.Data()), m_print_cut_verbose);
        return false;
    }
    
    bool fidVertex = m_coordSysTool->inFiducial( vertex->position().x(), vertex->position().y(), vertex->position().z(), apothem, upZ, downZ );
    
    if(fidVertex){
        PrintInfo("Vertex is in fiducial volume", m_print_cut_verbose);
    }
    
    return fidVertex;
}

bool CC1P1PiAnalysis::FindParticles(Minerva::PhysicsEvent* event) const
{
    //This code is currently a little messy and needs cleaning up. There are currently lots of cross checks in the code here.
    
    PrintInfo("CC1P1PiAnalysis::FindParticles", m_print_cut_verbose);
    //Determine which track is most proton like and pion like:
    // 1) Get particle scores and compare which track is has the highest score for the given hypothosis.
    // 2) Look for Michel features.
    //
    //Check that they are contianed in det FV and they are not minos matched.
   
    Minerva::ProngVect prongs = event->primaryProngs();
    Minerva::ProngVect::iterator prong;
    
    std::vector<double> protonScore;
    std::vector<double> pionScore;
    std::vector<double> trackChi2NDF;
    
    SmartRef<Minerva::Particle> Part_from_Prong1;
    SmartRef<Minerva::Particle> Part_from_Prong2;
    
    int prong_count = 0;
    int hadron_counter = 0;
    //std::vector<Minerva::Particle::ID> hyp_order;
    
    std::vector<double> pr_score_N_Vec;
    std::vector<double> pi_score_N_Vec;
    
    int Prong1_PDG = -999;
    int Prong2_PDG = -999;
    
    for( prong = prongs.begin(); prong != prongs.end(); prong++ ){
        
        prong_count++;
        PrintInfo(Form("Checking Prong: %d", prong_count), m_print_cut_verbose);
        
        //Check prong isn't muon prong:
        if(m_MuonProng == (*prong)){
            PrintInfo("Prong already determined as muon.", m_print_cut_verbose);
            continue;
        }
        
        //Check if track is fully contained in detector FV and there are no minos matched tracks:
        Minerva::TrackVect tracks = (*prong)->minervaTracks();
        
        if( tracks.empty() ) {
            PrintInfo("  This prong contains an empty vector of tracks, skipping!", m_print_cut_verbose);
            continue;
            //Return false statement if found to minos match.
        }
        else if( (*prong)->MinosTrack() || (*prong)->MinosStub() ) {
            PrintInfo("  This is a MINOS matched prong, skipping!", m_print_cut_verbose);
            continue;
            //Return false statement if found to minos match.
        }
        
        hadron_counter++;
     
        //Why are we taking the last track? What is the size of the track?
        SmartRef<Minerva::Track> track = tracks[tracks.size() - 1];
        Gaudi::XYZPoint endpoint = track->lastState().position();
        
        if(!m_coordSysTool->inFiducial(endpoint.x(), endpoint.y(), endpoint.z(), m_det_apothem, m_det_upZ, m_det_downZ)){
            PrintInfo("Track not contained in detector fiducial volume.", m_print_cut_verbose);
            return false;
        }
        
        //The following code is based on code in ProtonUtils.
        std::vector<Minerva::Particle::ID> hypotheses;
        hypotheses.push_back(Minerva::Particle::Pion);
        hypotheses.push_back(Minerva::Particle::Proton);
        IParticleMakerTool::NameAliasListType toolsToUse;
        toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );
        
        bool found_particle = m_particleMaker->makeParticles((*prong), hypotheses, toolsToUse);
        
        if(found_particle){
            PrintInfo(Form("This prong has %d particle hypotheses attached.", (int)((*prong)->particles().size())), m_print_cut_verbose);
        }
        else{
            PrintInfo("Failed to produce particles", m_print_cut_verbose);
        }
        
        if((*prong)->particles().size() == 2){
            
            trackChi2NDF.push_back(track->chi2PerDoF());
            
            Minerva::ParticleVect partHypVec = (*prong)->particles();
            Minerva::ParticleVect::iterator part;
            
            double pr_score  = 0.0;
            double pi_score  = 0.0;
            double score_den = 0.0;
            
            for(part = partHypVec.begin(); part != partHypVec.end(); part++){
                //PrintInfo(Form("Testing %d hypothesis with signature: %d and score: %f", (*part)->idcode(), (*part)->methodSignature(), (*part)->score()), m_print_cut_verbose);
                
                std::string part_name;
                double minPartScore = -999.0;
                double maxPartChi2 = -999.0;
                
                //For now let's just compare the hyp of which is more proton/pion like
                //This simply sets the part_name, scores and chi2 are not used yet
                if((*part)->idcode() == Minerva::Particle::Proton){
                    part_name = "Proton";
                    minPartScore = m_minProtonScore;
                    maxPartChi2 = m_maxProtonChi2;
                    PrintInfo("        Running checks on Proton Hypothesis.", m_print_cut_verbose);
                }
                else if((*part)->idcode() == Minerva::Particle::Pion){
                    part_name = "Pion";
                    minPartScore = m_minPionScore;
                    maxPartChi2 = m_maxPionChi2;
                    PrintInfo("        Running checks on Pion Hypothesis.", m_print_cut_verbose);
                }
                
                //Actual PID bit: Does the particle
                if( (*part)->methodSignature().find("dEdX") != std::string::npos ) {
                    
                    if((*part)->idcode() == Minerva::Particle::Proton){
                        protonScore.push_back((*part)->score());
                        pr_score = (*part)->score();
                    }
                    else  if((*part)->idcode() == Minerva::Particle::Pion){
                        pionScore.push_back((*part)->score());
                        pi_score = (*part)->score();
                    }
            
                    score_den += (*part)->score();
                }
                
            }
        
            //Temp particle:
            //Add the Norms here:
            double pr_score_N = pr_score/score_den;
            double pi_score_N = pi_score/score_den;
            
            pr_score_N_Vec.push_back(pr_score_N);
            pi_score_N_Vec.push_back(pi_score_N);
            
            PrintInfo(Form("Prong %d:", hadron_counter), m_print_cut_verbose);
            PrintInfo(Form("        Proton Score %f (%f), Pion Score %f, (%f) Chi2NDF %f", pr_score, pr_score_N, pi_score, pi_score_N, track->chi2PerDoF()), m_print_cut_verbose);
            
            int found_p_no = 0;
            Minerva::Particle::ID part_name_check = Minerva::Particle::Unknown;
            int PDGCode = -999;
            
            if(pr_score_N < pi_score_N){
                
                for(int hyp_counter = 0; hyp_counter < (int)partHypVec.size(); hyp_counter++){
                    if( partHypVec[hyp_counter]->idcode() == Minerva::Particle::Proton){
                        found_p_no = hyp_counter;
                    }
                }
                
                part_name_check = Minerva::Particle::Proton;
                PDGCode = m_Proton_PDG;
                
            }
            else{
                for(int hyp_counter = 0; hyp_counter < (int)partHypVec.size(); hyp_counter++){
                    if( partHypVec[hyp_counter]->idcode() == Minerva::Particle::Pion){
                        found_p_no = hyp_counter;
                    }
                }
                
                part_name_check = Minerva::Particle::Pion;
                PDGCode = m_Pion_PDG;
            }
            
            //PrintInfo(Form("        Prong %d believed to be %d and has found_p_no = %d", hadron_counter, part_name_check, found_p_no), m_print_cut_verbose);

            
            if(hadron_counter == 1){
                Part_from_Prong1 = partHypVec[found_p_no];
                Prong1_PDG = PDGCode;
                
                PrintInfo(Form("        Part_from_Prong1 :: Consistent with %d Hyp?", part_name_check), m_print_cut_verbose);
                if(Part_from_Prong1->idcode() == part_name_check){
                    PrintInfo("YES!!!!", m_print_cut_verbose);
                }
                else{
                    PrintInfo("NO ********************** ?!", m_print_cut_verbose);
                }
                PrintInfo(Form("        IDCode: %d Score: %f", Part_from_Prong1->idcode(), Part_from_Prong1->score()), m_print_cut_verbose);
                
            }
            
            if(hadron_counter == 2){
                Part_from_Prong2 = partHypVec[found_p_no];
                Prong2_PDG = PDGCode;
                
                PrintInfo(Form("        Part_from_Prong2 :: Consistent with %d Hyp?", part_name_check), m_print_cut_verbose);
                if(Part_from_Prong2->idcode() == part_name_check){
                    PrintInfo("YES!!!!", m_print_cut_verbose);
                }
                else{
                    PrintInfo("NO ********************** ?!", m_print_cut_verbose);
                }
                
                PrintInfo(Form("        IDCode: %d Score: %f", Part_from_Prong2->idcode(), Part_from_Prong2->score()), m_print_cut_verbose);
            }
        }
        
        //Look for michels at end of the prong
    }
    
    //Given the particle hypotheses, set the candidate tracks:
    int pr_prong_no = -999;
    int pi_prong_no = -999;
    
    PrintInfo("******************************** Summary ********************************", m_print_cut_verbose);
    PrintInfo(Form("Vector Sizes Consistent: trackChi2NDF N = %d protonScore N = %d pionScore N = %d", (int)(trackChi2NDF.size()), (int)(protonScore.size()), (int)(pionScore.size())), m_print_cut_verbose);
    if(trackChi2NDF.size() == 2 && protonScore.size()  == 2 && pionScore.size() == 2){
        PrintInfo(" YES.", m_print_cut_verbose);
        
        double Prong1_Proton = protonScore[0]/(protonScore[0] + pionScore[0]);
        double Prong1_Pion = pionScore[0]/(protonScore[0] + pionScore[0]);
        
        double Prong2_Proton = protonScore[1]/(protonScore[1] + pionScore[1]);
        double Prong2_Pion = pionScore[1]/(protonScore[1] + pionScore[1]);
    
        PrintInfo("Prong 1:", m_print_cut_verbose);
        PrintInfo(Form("        Proton Score %f (%f) Pion Score %f (%f) Chi2NDF %f", protonScore[0], Prong1_Proton, pionScore[0], Prong1_Pion, trackChi2NDF[0]), m_print_cut_verbose);
        PrintInfo(Form("        PreCal Pr Sc N %f PreCal Pi Sc N %f", pr_score_N_Vec[0], pi_score_N_Vec[0]), m_print_cut_verbose);
        PrintInfo("Prong 2:", m_print_cut_verbose);
        PrintInfo(Form("        Proton Score %f (%f) Pion Score %f (%f) Chi2NDF %f", protonScore[1], Prong2_Proton, pionScore[1], Prong2_Pion, trackChi2NDF[1]), m_print_cut_verbose);
        PrintInfo(Form("        PreCal Pr Sc N %f PreCal Pi Sc N %f", pr_score_N_Vec[1], pi_score_N_Vec[1]), m_print_cut_verbose);
        PrintInfo("*************************************************************************", m_print_cut_verbose);

        int Prong1a_PDG = -999;
        if(Prong1_Proton > Prong1_Pion){
            Prong1a_PDG = m_Proton_PDG;
            //debug() << "Prong 1 is thought to be a Proton" << endmsg;
        }
        else{
            Prong1a_PDG = m_Pion_PDG;
            //debug() << "Prong 1 is thought to be a Pion" << endmsg;
        }
        
        int Prong2a_PDG = -999;
        if(Prong2_Proton > Prong2_Pion){
            Prong2a_PDG = m_Proton_PDG;
            //debug() << "Prong 2 is thought to be a Proton" << endmsg;
        }
        else{
            Prong2a_PDG = m_Pion_PDG;
            //debug() << "Prong 2 is thought to be a Pion" << endmsg;
        }
        
        PrintInfo("Checking PDG Codes:", m_print_cut_verbose);
        PrintInfo(Form("Prong 1: Pre: %d, Post %d", Prong1_PDG, Prong1a_PDG), m_print_cut_verbose);
        if(Prong1_PDG == Prong1a_PDG) PrintInfo(". They are the same!!!", m_print_cut_verbose);
        else PrintInfo(". Close but no cigar... :-(", m_print_cut_verbose);
        
        PrintInfo(Form("Prong 2: Pre: %d Post %d", Prong2_PDG, Prong2a_PDG), m_print_cut_verbose);
        if(Prong2_PDG == Prong2a_PDG) PrintInfo(". They are the same!!!", m_print_cut_verbose);
        else PrintInfo("Close but no cigar... :-(", m_print_cut_verbose);
        
        if(Prong1_PDG != Prong2_PDG){
            if(Prong1_PDG == m_Proton_PDG){
                pr_prong_no = 0;
                pi_prong_no = 1;
                
                m_ProtonParticle = Part_from_Prong1;
                m_PionParticle = Part_from_Prong2;
            }
            else{
                pr_prong_no = 1;
                pi_prong_no = 0;
                
                m_ProtonParticle = Part_from_Prong2;
                m_PionParticle = Part_from_Prong1;
            }
        }
        else if(Prong1_PDG == m_Proton_PDG){
            PrintInfo("Found two protons...", m_print_cut_verbose);
            return false;
        }
        else{
            PrintInfo("Found two pions...", m_print_cut_verbose);
            return false;
        }
    }
    else{
        PrintInfo(" No... Check this out!!!!", m_print_cut_verbose);
        return false;
    }

    m_ProtonProng = prongs[pr_prong_no];
    m_PionProng = prongs[pi_prong_no];
    
    m_ProtonScore[0] = protonScore[pr_prong_no];
    m_ProtonScore[1] = pionScore[pr_prong_no];
    
    m_PionScore[0] = protonScore[pi_prong_no];
    m_PionScore[1] = pionScore[pi_prong_no];
    
    m_Chi2NDF[0] = trackChi2NDF[pr_prong_no];
    m_Chi2NDF[1] = trackChi2NDF[pi_prong_no];
    
    bool pr_is_correct = false;
    bool pi_is_correct = false;
    PrintInfo("Final Check that the tracks are what we think they are", m_print_cut_verbose);
    PrintInfo(Form("ProtonParticle: %d", m_ProtonParticle->idcode()), m_print_cut_verbose);
    if(m_ProtonParticle->idcode() == Minerva::Particle::Proton){
        PrintInfo(" YES", m_print_cut_verbose);
        pr_is_correct = true;
    }
    PrintInfo(Form("PionParticle: %d", m_PionParticle->idcode()), m_print_cut_verbose);
    if(m_PionParticle->idcode() == Minerva::Particle::Pion){
        PrintInfo(" YES", m_print_cut_verbose);
        pi_is_correct = true;
    }
    
    if(!(pr_is_correct || pi_is_correct)){
        PrintInfo("Particles not correct... check code", m_print_cut_verbose);
        return false;
    }
    
    PrintInfo("Found Proton and Pion Tracks", m_print_cut_verbose);
    
    return true;
}

void CC1P1PiAnalysis::ResetParticles() const
{
    m_MuonProng = NULL;
    m_MuonParticle = NULL;

    m_ProtonProng = NULL;
    m_ProtonParticle = NULL;
    
    m_PionProng = NULL;
    m_PionParticle = NULL;
    
    m_ProtonScore[0] = -999.;
    m_ProtonScore[1] = -999.;
    
    m_PionScore[0] = -999.;
    m_PionScore[1] = -999.;
    
    m_Chi2NDF[0] = -999.;
    m_Chi2NDF[1] = -999.;
    
    for(int i = 0; i < 4; i++){
        m_Pion4Mom[i] = -999.;//TXYZ
        m_Piontrue4Mom[i] = -999.;
        
        m_Proton4Mom[i] = -999.;
        m_Protontrue4Mom[i] = -999.;
        
        m_Muon4Mom[i] = -999.;
        m_Muontrue4Mom[i] = -999.;
    }
    
}

//Generic Particle information builder:
void CC1P1PiAnalysis::SetPartInfo(std::string name)
{
    
    if(name == "pr" || name == "pi"){
        declareDoubleBranch(m_hypMeths, (name + "_prscore").c_str() , -999.);
        declareDoubleBranch(m_hypMeths, (name + "_piscore").c_str() , -999.);
    }
    else{
        declareDoubleBranch(m_hypMeths, (name + "_score").c_str() , -999.);
    }
    
    declareDoubleBranch(m_hypMeths, (name + "_det_frac").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_det_otherE").c_str(), -999.);
    
    declareDoubleBranch(m_hypMeths, (name + "_chi2ndf").c_str(), -999.);

    declareDoubleBranch(m_hypMeths, (name + "_E").c_str() , -999.);
    declareDoubleBranch(m_hypMeths, (name + "_trueE").c_str(), -999.);

    declareDoubleBranch(m_hypMeths, (name + "_mom").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_truemom").c_str(), -999.);
    
    declareContainerDoubleBranch(m_hypMeths, (name + "_4mom").c_str(), 4, -999.);
    declareContainerDoubleBranch(m_hypMeths, (name + "_true4mom").c_str(), 4, -999.);
    
    declareDoubleBranch(m_hypMeths, (name + "_pTMag").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_truepTMag").c_str(), -999.);
    
    declareContainerDoubleBranch(m_hypMeths, (name + "_pT").c_str(), 3, -999.);
    declareContainerDoubleBranch(m_hypMeths, (name + "_truepT").c_str(), 3, -999.);
    
    declareDoubleBranch(m_hypMeths, (name + "_pTT").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_truepTT").c_str(), -999.);
    
    declareIntBranch(m_hypMeths, (name + "_PDG").c_str(), -999);
    
    declareDoubleBranch(m_hypMeths, (name + "_Phi").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_truePhi").c_str(), -999.);
    
    declareDoubleBranch(m_hypMeths, (name + "_Theta").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_trueTheta").c_str(), -999.);
    
    declareDoubleBranch(m_hypMeths, (name + "_KE").c_str(), -999.);
    declareDoubleBranch(m_hypMeths, (name + "_trueKE").c_str(), -999.);
    
    declareContainerDoubleBranch(m_hypMeths, (name + "_startpos_xyz").c_str(), 3, -999.);
    declareContainerDoubleBranch(m_hypMeths, (name + "_truestartpos_xyz").c_str(), 3, -999.);
    declareContainerDoubleBranch(m_hypMeths, (name + "_endpos_xyz").c_str(), 3, -999.);
    declareContainerDoubleBranch(m_hypMeths, (name + "_trueendpos_xyz").c_str(), 3, -999.);
    
    if(name == "pr" || name == "pi"){
        declareIntBranch(m_hypMeths, (name + "_FSI").c_str(), -999);
    }
    
}

void CC1P1PiAnalysis::SetCommonBranches()
{
    declareDoubleBranch(m_hypMeths, "Enu", -999.);
    declareDoubleBranch(m_hypMeths, "trueEnu", -999.);
    
    declareDoubleBranch(m_hypMeths, "Q2", -999.0);
    declareDoubleBranch(m_hypMeths, "trueQ2", -999.0);
    
    declareDoubleBranch(m_hypMeths, "dpTT", -999.0);
    declareDoubleBranch(m_hypMeths, "truedpTT", -999.0);
    
    declareDoubleBranch(m_hypMeths, "dpTT_pi", -999.0);
    declareDoubleBranch(m_hypMeths, "truedpTT_pi", -999.0);
    
    declareDoubleBranch(m_hypMeths, "dpTT_pr", -999.0);
    declareDoubleBranch(m_hypMeths, "truedpTT_pr", -999.0);

    declareDoubleBranch(m_hypMeths, "dpT", -999.0);
    declareDoubleBranch(m_hypMeths, "truedpT", -999.0);
    
    declareContainerDoubleBranch(m_hypMeths, "dpT_vec", 3, -999.0);
    declareContainerDoubleBranch(m_hypMeths, "truedpT_vec", 3, -999.0);

    declareDoubleBranch(m_hypMeths, "dalphaT", -999.0);
    declareDoubleBranch(m_hypMeths, "truedalphaT", -999.0);

    declareDoubleBranch(m_hypMeths, "dphiT", -999.0);
    declareDoubleBranch(m_hypMeths, "truedphiT", -999.0);

}

void CC1P1PiAnalysis::FillCommonBranches(const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, Minerva::NeutrinoInt* cc1p1piHyp) const
{
    
    double Enu = m_Muon4Mom[0] + m_Proton4Mom[0] + m_Pion4Mom[0];
    cc1p1piHyp->setDoubleData("Enu", Enu);
    
    double Q2 = -999.;
    cc1p1piHyp->setDoubleData("Q2", Q2);
    
    const TVector3 * mu_p = new TVector3(m_Muon4Mom[1], m_Muon4Mom[2], m_Muon4Mom[3]);
    const TVector3 * pr_p = new TVector3(m_Proton4Mom[1], m_Proton4Mom[2], m_Proton4Mom[3]);
    const TVector3 * pi_p = new TVector3(m_Pion4Mom[1], m_Pion4Mom[2], m_Pion4Mom[3]);
    
    SmartRef<Minerva::Vertex> int_vert = event->interactionVertex();
    const Gaudi::XYZPoint vert_3vec = int_vert->position();
    double vertex[3] = {vert_3vec.x(), vert_3vec.y(), vert_3vec.z()};//{0.};
    
    double dpTT = -999.;
    double dpTMag = -999.;
    double dalphaT = -999.;
    double dphiT = -999.;
    
    TVector3 * dpT_3mom = GetTransverseVars(vertex, mu_p, pr_p, pi_p, dpTT, dpTMag, dalphaT, dphiT);
    std::vector<double> vec_dpT_3mom;
    vec_dpT_3mom.push_back(dpT_3mom->X());
    vec_dpT_3mom.push_back(dpT_3mom->Y());
    vec_dpT_3mom.push_back(dpT_3mom->Z());
    
    cc1p1piHyp->setDoubleData("dpTT", dpTT);
    cc1p1piHyp->setDoubleData("dpT", dpTMag);
    cc1p1piHyp->setDoubleData("dalphaT", dalphaT);
    cc1p1piHyp->setDoubleData("dphiT", dphiT);
   
    double dpTT_pi = GetDPTT(vertex, pi_p, mu_p, pr_p);
    
    cc1p1piHyp->setDoubleData("dpTT_pi", dpTT_pi);
    
    double dpTT_pr = GetDPTT(vertex, pr_p, pi_p, mu_p);
    cc1p1piHyp->setDoubleData("dpTT_pr", dpTT_pr);
    
    if(truth){
        
        double trueEnu = -999.;
        cc1p1piHyp->setDoubleData("trueEnu", trueEnu);
        
        double trueQ2 = -999.;
        cc1p1piHyp->setDoubleData("trueQ2", trueQ2);
        
        const TVector3 * truemu_p = new TVector3(m_Muontrue4Mom[1], m_Muontrue4Mom[2], m_Muontrue4Mom[3]);
        const TVector3 * truepr_p = new TVector3(m_Protontrue4Mom[1], m_Protontrue4Mom[2], m_Protontrue4Mom[3]);
        const TVector3 * truepi_p = new TVector3(m_Piontrue4Mom[1], m_Piontrue4Mom[2], m_Piontrue4Mom[3]);
        
        const Gaudi::LorentzVector nu_4vec = truth->IncomingPartVec();
        double nu_3vec_mag = sqrt(nu_4vec.px()*nu_4vec.px() + nu_4vec.py()*nu_4vec.py() + nu_4vec.pz()*nu_4vec.pz());
        double vertex_true[3] = { nu_4vec.px()/nu_3vec_mag, nu_4vec.py()/nu_3vec_mag, nu_4vec.pz()/nu_3vec_mag };
        
        double truedpTT = -999.;
        double truedpT = -999.;
        double truedalphaT = -999.;
        double truedphiT = -999.;
        
        TVector3 * dpT_3mom_true = GetTransverseVars(vertex_true, truemu_p, truepr_p, truepi_p, truedpTT, truedpT, truedalphaT, truedphiT, true);
        std::vector<double> vec_dpT_3mom_true;
        vec_dpT_3mom_true.push_back(dpT_3mom_true->X());
        vec_dpT_3mom_true.push_back(dpT_3mom_true->Y());
        vec_dpT_3mom_true.push_back(dpT_3mom_true->Z());
        
        cc1p1piHyp->setDoubleData("truedpTT", truedpTT);
        cc1p1piHyp->setDoubleData("truedpT", truedpT);
        cc1p1piHyp->setDoubleData("truedalphaT", truedalphaT);
        cc1p1piHyp->setDoubleData("truedphiT", truedphiT);
        
        double truedpTT_pi = GetDPTT(vertex_true, truepi_p, truemu_p, truepr_p, true);
        cc1p1piHyp->setDoubleData("truedpTT_pi", truedpTT_pi);
        
        double truedpTT_pr = GetDPTT(vertex_true, truepr_p, truepi_p, truemu_p, true);
        cc1p1piHyp->setDoubleData("truedpTT_pr", truedpTT_pr);
        
    }
}

void CC1P1PiAnalysis::FillPartInfo(std::string name, const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, Minerva::NeutrinoInt* cc1p1piHyp) const
{
    
    SmartRef<Minerva::Prong> prong;
    SmartRef<Minerva::Particle> particle;
    double mass = 0.0;
    int ch2_vec_no = 0;
    double tmp_scores[2] = {-999.};
    
    if(name == "mu"){
        prong = m_MuonProng;
        particle = m_MuonParticle;
        mass = MinervaUnits::M_mu;
    }
    else if(name == "pr"){
        prong = m_ProtonProng;
        particle = m_ProtonParticle;
        mass = MinervaUnits::M_p;
        ch2_vec_no = 0;
        tmp_scores[0] = m_ProtonScore[0];
        tmp_scores[1] = m_ProtonScore[1];
    }
    else if(name == "pi"){
        prong = m_PionProng;
        particle = m_PionParticle;
        mass = MinervaUnits::M_pion;
        ch2_vec_no = 1;
        tmp_scores[0] = m_PionScore[0];
        tmp_scores[1] = m_PionScore[1];
    }
    else{
        error() << "CC1P1PiAnalysis::FillPartInfo :: Could not find determine name \"" << name << "\". Please check";
    }
    
    if(mass == 0.0){
        warning() << "CC1P1PiAnalysis::FillPartInfo :: " << name << " mass is zero!!!" << endmsg;
    }

    double ch2ndf = -999.;
    
    if(name == "pr" || name == "pi"){
        cc1p1piHyp->setDoubleData( (name + "_prscore").c_str(), tmp_scores[0]);
        cc1p1piHyp->setDoubleData( (name + "_piscore").c_str(), tmp_scores[1]);
        
        double hasFSI = -999;
        cc1p1piHyp->setIntData( (name + "_FSI").c_str(), hasFSI);

        ch2ndf = m_Chi2NDF[ch2_vec_no];

    }
    else{
        double score = particle->score();
        cc1p1piHyp->setDoubleData( (name + "_score").c_str(), score );
        ch2ndf = -999.;
    }
    
    cc1p1piHyp->setDoubleData( (name + "_chi2ndf").c_str(), ch2ndf);
    
    //Reco vars:
    
    Gaudi::LorentzVector four_vec = particle->momentumVec();
    
    double Energy = four_vec.E();
    cc1p1piHyp->setDoubleData( (name + "_E").c_str(), Energy);
    
    double mom = sqrt( pow( four_vec.E(), 2 ) - pow( mass, 2 ) );
    cc1p1piHyp->setDoubleData( (name + "_mom").c_str(), mom);
    
    std::vector<double> sel4mom;
    sel4mom.push_back(four_vec.E());
    sel4mom.push_back(four_vec.px());
    sel4mom.push_back(four_vec.py());
    sel4mom.push_back(four_vec.pz());
    
    Rotate2BeamCoords(sel4mom);
    
    SetGlobal4Vec(name, sel4mom);

    cc1p1piHyp->setContainerDoubleData( (name + "_4mom").c_str(), sel4mom);
    
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    const Gaudi::XYZPoint vert_3vec = vertex->position();
    double vtx[3] = {vert_3vec.x(), vert_3vec.y(), vert_3vec.z()};
    
    const TVector3 * mom_vec = new TVector3(sel4mom[1], sel4mom[2], sel4mom[3]);
    const TVector3 * pT = GetPT(vtx, mom_vec);
    
    cc1p1piHyp->setDoubleData( (name + "_pTMag").c_str(), pT->Mag());

    std::vector<double> selpT;
    selpT.push_back(pT->X());
    selpT.push_back(pT->Y());
    selpT.push_back(pT->Z());
    
    cc1p1piHyp->setContainerDoubleData( (name + "_pT").c_str(), selpT);

    double pTT = -999.;
    cc1p1piHyp->setDoubleData( (name + "_pTT").c_str(), pTT);
    
    double Phi = m_coordSysTool->phiWRTBeam( four_vec );
    cc1p1piHyp->setDoubleData( (name + "_Phi").c_str(), Phi);
    
    double Theta = m_coordSysTool->thetaWRTBeam( four_vec );
    cc1p1piHyp->setDoubleData( (name + "_Theta").c_str(), Theta);
    
    double KE = four_vec.E() - mass;
    cc1p1piHyp->setDoubleData( (name + "_KE").c_str(), KE);
    
    Gaudi::XYZPoint upstream = (*prong->minervaTracks().front() ).upstreamState().position();
    std::vector<double> sel_start_xyz;
    sel_start_xyz.push_back(upstream.x());
    sel_start_xyz.push_back(upstream.y());
    sel_start_xyz.push_back(upstream.z());
    cc1p1piHyp->setContainerDoubleData( (name + "_startpos_xyz").c_str(), sel_start_xyz);
    
    Gaudi::XYZPoint downstream = (*prong->minervaTracks().front() ).downstreamState().position();
    std::vector<double> sel_end_xyz;
    sel_end_xyz.push_back(downstream.x());
    sel_end_xyz.push_back(downstream.y());
    sel_end_xyz.push_back(downstream.z());
    cc1p1piHyp->setContainerDoubleData( (name + "_endpos_xyz").c_str(), sel_end_xyz);
    
    //True vars:
    if(truth){
        //********************** Old Truth Infromation **********************//
       /* std::vector<const Minerva::TG4Trajectory*> trajectories;
        
        const Minerva::TG4Trajectory* tj = NULL;
        double other_energy = 0.0;
        std::map<const Minerva::TG4Trajectory*,double>::iterator it;
        std::map<const Minerva::TG4Trajectory*,double> trajMap = m_truthMatcher->getTG4Trajectories(prong, other_energy);
        if( !(trajMap.empty()) ){
            for(it = trajMap.begin(); it != trajMap.end(); it++) {
                tj = (*it).first;
                trajectories.push_back(tj);
            }
        }
        
        if(trajectories.size() == 1){
            debug() << "Prong has one tragectory!" << endmsg;
        }
        else{
            debug() << "Prong has " << trajectories.size() << " trajectories!" << endmsg;
        }*/
        
        //const Minerva::TG4Trajectory* traj = trajectories[0];
        //****************************** END ******************************//
        
        const Minerva::TG4Trajectory * traj = NULL;
        double fraction = -999.;
        double other_energy = -999.;
        
        StatusCode found_trag = m_truthMatcher->getTG4Trajectory(prong, traj, fraction, other_energy);
        
        if(found_trag && traj){
            
            cc1p1piHyp->setDoubleData( (name + "_det_frac").c_str(), fraction);
            cc1p1piHyp->setDoubleData( (name + "_det_otherE").c_str(), other_energy);
        
            Gaudi::LorentzVector traj_4p = traj->GetInitialMomentum();
            
            double trueEnergy = traj_4p.E();
            cc1p1piHyp->setDoubleData( (name + "_trueE").c_str(), trueEnergy);
            
            double truemom = sqrt( traj_4p.px()*traj_4p.px() + traj_4p.py()*traj_4p.py() + traj_4p.pz()*traj_4p.pz() );
            cc1p1piHyp->setDoubleData( (name + "_truemom").c_str(), truemom);
            
            std::vector<double> true4mom;
            true4mom.push_back(traj_4p.E());
            true4mom.push_back(traj_4p.px());
            true4mom.push_back(traj_4p.py());
            true4mom.push_back(traj_4p.pz());
            
            SetGlobal4Vec(name, true4mom, true);
            
            cc1p1piHyp->setContainerDoubleData( (name + "_true4mom").c_str(), true4mom);
            
            const Gaudi::LorentzVector nu_4vec = truth->IncomingPartVec();
            
            double nu_3vec_mag = sqrt(nu_4vec.px()*nu_4vec.px() + nu_4vec.py()*nu_4vec.py() + nu_4vec.pz()*nu_4vec.pz());
            double nu_3vec[3] = { nu_4vec.px()/nu_3vec_mag, nu_4vec.py()/nu_3vec_mag, nu_4vec.pz()/nu_3vec_mag };
            
            const TVector3 * true_mom_vec = new TVector3(true4mom[1], true4mom[2], true4mom[3]);
            const TVector3 * truepT_3vec = GetPT(nu_3vec, true_mom_vec);

            cc1p1piHyp->setDoubleData( (name + "_truepTMag").c_str(), truepT_3vec->Mag());
            
            std::vector<double> truepT;
            truepT.push_back(truepT_3vec->X());
            truepT.push_back(truepT_3vec->Y());
            truepT.push_back(truepT_3vec->Z());
            
            cc1p1piHyp->setContainerDoubleData( (name + "_truepT").c_str(), truepT);
            
            double truepTT = -999.;
            cc1p1piHyp->setDoubleData( (name + "_truepTT").c_str(), truepTT);
            
            int PDG = traj->GetPDGCode();
            cc1p1piHyp->setIntData( (name + "_PDG").c_str(), PDG);
            
            double truePhi = traj->GetInitialMomentum().Phi();
            cc1p1piHyp->setDoubleData( (name + "_truePhi").c_str(), truePhi);
            
            double trueTheta = traj->GetInitialMomentum().Theta();
            cc1p1piHyp->setDoubleData( (name + "_trueTheta").c_str(), trueTheta);
            
            double trueMass = 0.;
            switch (TMath::Abs(PDG)) {
                case 11:    trueMass = MinervaUnits::M_e;       break;
                case 13:    trueMass = MinervaUnits::M_mu;      break;
                case 15:    trueMass = MinervaUnits::M_tau;     break;
                case 2212:  trueMass = MinervaUnits::M_p;       break;
                case 2112:  trueMass = MinervaUnits::M_n;       break;
                case 211:   trueMass = MinervaUnits::M_pion;    break;
                case 111:   trueMass = MinervaUnits::M_pi0;     break;
                case 321:   trueMass = MinervaUnits::M_kaon;    break;
                case 311:   trueMass = MinervaUnits::M_k0;      break;
                default:    trueMass = 0.;                      break;
            }
            
            if(trueMass == 0. && PDG != 22){//Photons are massless so that would be okay but if it's not... Well oh noooo.
                warning() << "CC1P1PiAnalysis::FillPartInfo :: True Mass is Zero." << endmsg;
            }
            
            double trueKE = traj_4p.E() - trueMass;
            cc1p1piHyp->setDoubleData( (name + "_trueKE").c_str(), trueKE);
            
            Gaudi::LorentzVector inipos = traj->GetInitialPosition();
            std::vector<double> tru_start_xyz;
            tru_start_xyz.push_back(inipos.x());
            tru_start_xyz.push_back(inipos.y());
            tru_start_xyz.push_back(inipos.z());
            cc1p1piHyp->setContainerDoubleData( (name + "_truestartpos_xyz").c_str(), tru_start_xyz);
            
            Gaudi::LorentzVector finpos = traj->GetFinalPosition();
            std::vector<double> tru_end_xyz;
            tru_end_xyz.push_back(finpos.x());
            tru_end_xyz.push_back(finpos.y());
            tru_end_xyz.push_back(finpos.z());
            cc1p1piHyp->setContainerDoubleData( (name + "_trueendpos_xyz").c_str(), tru_end_xyz);
        }
    }
    
}

void CC1P1PiAnalysis::SetAccumLevel() const
{
    //m_accum_level[ cut - 1 ] = cut;
    m_accum_level++;
    PrintInfo(Form("***** Accum. Level %d ***** ", m_accum_level), m_print_acc_level);
}

void CC1P1PiAnalysis::ResetAccumLevel() const
{
    m_accum_level = 0;
    //for(int i = 0; i < m_ncuts; i++){
    //    m_accum_level[i] = 0;
    //}
}

void CC1P1PiAnalysis::SaveAccumLevel(Minerva::PhysicsEvent * event, Minerva::GenMinInteraction* truth) const
{
    //debug() << " " << endmsg;
    //debug() << " " << endmsg;
    //debug() << "Call to save accum_level" << endmsg;
    if(m_accum_level >= m_accum_level_to_save){
        //debug() << "Passed save requirement" << endmsg;
        event->setIntData("accum_level", m_accum_level);
        truth->setIntData("accum_level", m_accum_level);
        markEvent(event);
        
        PrintInfo(Form("++++ Saving Accum. Level %d ++++", m_accum_level), m_print_acc_level);
        if(m_accum_level < m_ncuts){
          //  debug() << "Event beleived to be below cut threshold." << endmsg;
            Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( m_anaSignature );
            
            std::vector<Minerva::NeutrinoInt*> nuInts;
            nuInts.push_back(nuInt);
            
            StatusCode sc = addInteractionHyp( event, nuInts );
            if(sc){
                PrintInfo("Added Int. Hyp", m_print_acc_level);
            }
            else {
                PrintInfo("Failed to add Int. Hyp", m_print_acc_level);
            }
        }
        else{
            PrintInfo("Event beleived to have passed all cuts.", m_print_other);
        }
        
        debug() << "Should have saved event" << endmsg;
    }
    else PrintInfo(Form("Failed to reach accum. level %d. Selection stopped at %d.", m_accum_level_to_save, m_accum_level), m_print_acc_level);

    //debug() << " " << endmsg;
    //debug() << " " << endmsg;
    //std::vector<int> tmp_vec;
    //for(int i = 0; i < m_ncuts; i++){
    //    tmp_vec.push_back(m_accum_level[i]);
    //}
    //event->setContainerIntData("accum_level",tmp_vec);
    
}

void CC1P1PiAnalysis::Rotate2BeamCoords(std::vector<double> val) const
{
    //Determine size of 4 vec at some point...
    PrintInfo("CC1P1PiAnalysis::Rotate2BeamCoords", m_print_other);

    if(!((int)val.size() == 4)){
        PrintInfo(Form("Warning : Not a 4 vector! Vector has dimension %d", (int)val.size()), m_print_other);
    }
    
    PrintInfo(Form("Initial 4Vec: P_E %f P_X %f P_Y %f P_Z %f", val[0], val[1], val[2], val[3]), m_print_other);
    
    double py = val[2];
    double pz = val[3];
    //! momentum rotated to beam coordinate system
    double py_prime = -1.0 *sin( MinervaUnits::numi_beam_angle_rad )*pz + cos( MinervaUnits::numi_beam_angle_rad )*py;
    double pz_prime = cos( MinervaUnits::numi_beam_angle_rad )*pz + sin( MinervaUnits::numi_beam_angle_rad )*py;
    
    val[2] = py_prime;
    val[3] = pz_prime;
    
    PrintInfo(Form("Rotated 4Vec: P_E %f P_X %f P_Y %f P_Z %f", val[0], val[1], val[2], val[3]), m_print_other);
    
}

void CC1P1PiAnalysis::SetGlobal4Vec(std::string name, Gaudi::LorentzVector vec, bool truth) const
{
    if(truth){
        if(name == "mu"){
            m_Muontrue4Mom[0] = vec.E();
            m_Muontrue4Mom[1] = vec.px();
            m_Muontrue4Mom[2] = vec.py();
            m_Muontrue4Mom[3] = vec.pz();
        }
        else if(name == "pr"){
            m_Protontrue4Mom[0] = vec.E();
            m_Protontrue4Mom[1] = vec.px();
            m_Protontrue4Mom[2] = vec.py();
            m_Protontrue4Mom[3] = vec.pz();
        }
        else if(name == "pi"){
            m_Piontrue4Mom[0] = vec.E();
            m_Piontrue4Mom[1] = vec.px();
            m_Piontrue4Mom[2] = vec.py();
            m_Piontrue4Mom[3] = vec.pz();
        }
        else{
            error() << "SetGlobal4Vec::FillPartInfo :: Could not determine particle hyp. Please check" << endmsg;
        }
    }
    else{
        if(name == "mu"){
            m_Muon4Mom[0] = vec.E();
            m_Muon4Mom[1] = vec.px();
            m_Muon4Mom[2] = vec.py();
            m_Muon4Mom[3] = vec.pz();
        }
        else if(name == "pr"){
            m_Proton4Mom[0] = vec.E();
            m_Proton4Mom[1] = vec.px();
            m_Proton4Mom[2] = vec.py();
            m_Proton4Mom[3] = vec.pz();
        }
        else if(name == "pi"){
            m_Pion4Mom[0] = vec.E();
            m_Pion4Mom[1] = vec.px();
            m_Pion4Mom[2] = vec.py();
            m_Pion4Mom[3] = vec.pz();
        }
        else{
            error() << "SetGlobal4Vec::FillPartInfo :: Could not determine particle hyp. Please check" << endmsg;
        }

    }
}

void CC1P1PiAnalysis::SetGlobal4Vec(std::string name, std::vector<double> vec, bool truth) const
{
    debug() << "CC1P1PiAnalysis::SetGlobal4Vec(std::string name, std::vector<double> vec, bool truth)" << endmsg;
    debug() << "         Note: Vector must have form: vev[4] = {E, PX, PY, PZ}" << endmsg;
    
    Gaudi::LorentzVector lvec(vec[1], vec[2], vec[3], vec[0]);

    SetGlobal4Vec(name, lvec, truth);
}

//Transverse variables:
//If is_truth the vtx becomes the true neutrino direction.

TVector3 * CC1P1PiAnalysis::GetTransverseVars(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT, bool is_truth) const
{
    TVector3 * nudir = new TVector3();
    
    if(is_truth){
        nudir->SetXYZ(vtx[0],vtx[1],vtx[2]);
    }
    else{
        const TVector3 * tmp_vec = GetNuDirRec(vtx);
        nudir->SetXYZ(tmp_vec->X(),tmp_vec->Y(),tmp_vec->Z());
    }
    
    const TVector3 * neutrino_dir = new TVector3(nudir->X(),nudir->Y(), nudir->Z());
    
    const TVector3 * mupT = GetVecT(neutrino_dir, mumom);
    const TVector3 * prpT = GetVecT(neutrino_dir, prmom);
    const TVector3 * pipT = GetVecT(neutrino_dir, pimom);
    
    TVector3 * deltapt = new TVector3();
    SetDPT(deltapt, mupT, prpT, pipT);
    
    dpTMag  = deltapt->Mag();
    dalphaT = (deltapt->Theta())*TMath::RadToDeg();
    dphiT   = (deltapt->Phi())*TMath::RadToDeg();
    dpTT    = GetDPTT(vtx, mumom, prmom, pimom, is_truth);
    
    return deltapt;
}


TVector3 * CC1P1PiAnalysis::GetPT(double vtx[], const TVector3 *& mom, bool is_truth) const
{
    TVector3 * nudir = new TVector3();
    
    if(is_truth){
        nudir->SetXYZ(vtx[0],vtx[1],vtx[2]);
    }
    else{
        const TVector3 * tmp_vec = GetNuDirRec(vtx);
        nudir->SetXYZ(tmp_vec->X(),tmp_vec->Y(),tmp_vec->Z());
    }
    
    const TVector3 * neutrino_dir = new TVector3(nudir->X(),nudir->Y(), nudir->Z());
    
    TVector3 * pT = GetVecT(neutrino_dir, mom);
    
    return pT;
}

void CC1P1PiAnalysis::SetDPT(TVector3 * deltapt, const TVector3 *& ptmuon, const TVector3 *& ptproton, const TVector3 *& ptpion) const
{
    //ptmuon and ptproton already in the same plain which is perpendicular to the neutrino and already in a near back-to-back configuration
    TVector3 tmpd = (*ptmuon) + (*ptproton) + (*ptpion);
    TVector3 tmp_had = (*ptproton) + (*ptpion);

    double phi = TMath::ACos( ptmuon->Dot(tmp_had)*(-1)/(ptmuon->Mag()*tmp_had.Mag()) );
    
    double theta = TMath::ACos( tmpd.Dot(*ptmuon)*(-1)/tmpd.Mag()/ptmuon->Mag()  );
    
    deltapt->SetMagThetaPhi(tmpd.Mag(),theta, phi);
}

TVector3 * CC1P1PiAnalysis::GetVecT(const TVector3 *& refdir, const TVector3 *& mom) const
{
    //
    //w.r.t. beam direction
    //
    if(!refdir){
        error() << "CC1P1PiAnalysis::GetVecT refdir null" << endmsg;
        exit(1);
    }
    
    
    TVector3 vRotated(*mom);
    vRotated.Rotate(TMath::Pi(), *refdir);
    
    TVector3 *vt = new TVector3( (*mom - vRotated)*0.5 );
    
    return vt;
}

double CC1P1PiAnalysis::GetDPTT(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, bool is_truth) const
{
    TVector3 * nudir = new TVector3();
    
    if(is_truth){
        nudir->SetXYZ(vtx[0],vtx[1],vtx[2]);
    }
    else{
        const TVector3 * tmp_vec = GetNuDirRec(vtx);
        nudir->SetXYZ(tmp_vec->X(),tmp_vec->Y(),tmp_vec->Z());
    }
    
    TVector3 tmp1_vec = nudir->Cross(*mumom);
    tmp1_vec *= 1/tmp1_vec.Mag();
    
    TVector3 sum_vec = *prmom + *pimom;
    
    return sum_vec.Dot(tmp1_vec);
}

TVector3 * CC1P1PiAnalysis::GetNuDirRec(double vtx[]) const
{
    
    TVector3 * nup1local = new TVector3(vtx[0], vtx[1], vtx[2]);
    (*nup1local) *= 0.001;
    
    if( m_PDP->Mag() < EPSILON || nup1local->Mag() < EPSILON ){
        debug() << "CC1P1PiAnalysis::CalcNuDir bad input " << m_PDP->Mag() << " " << nup1local->Mag() << endmsg;
        return 0x0;
    }
    
    TVector3 *nuDirCalc = new TVector3( (*nup1local) - (*m_PDP) );
    (*nuDirCalc) *= 1./nuDirCalc->Mag();
    
    return nuDirCalc;
}

void CC1P1PiAnalysis::PrintInfo(std::string var, bool print) const
{
    if(print){
        debug() << var << endmsg;
    }
}

void CC1P1PiAnalysis::EventFinished() const
{
    if(m_print_acc_level){
        PrintInfo("-------------------------------------------------------------------------------------------------");
        PrintInfo("----------------------------            Finished             ------------------------------------");
        PrintInfo("-------------------------------------------------------------------------------------------------");
    }
}




