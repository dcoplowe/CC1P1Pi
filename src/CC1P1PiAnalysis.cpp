//include the header file...
#include "CC1P1PiAnalysis.h"

#include "Event/TG4PrimaryTrajectory.h"

//Forward declared headers:
#include "AnaUtils/IMuonUtils.h"
#include "GeoUtils/IMinervaCoordSysTool.h"
#include "DetDesc/Material.h"
#include "GeoUtils/INuclearTargetTool.h"
#include "ParticleMaker/IParticleMakerTool.h"


//Root headers:
#include <TString.h>

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
    
    //int m_ncuts;// = 5;
    m_accum_level = new int [ m_ncuts ];
    
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

    SetPartInfo("mu");
    SetPartInfo("pi");
    SetPartInfo("pr");

    //---------------------------------------------------------------------
    // Declare the Truth block branches.
    // Truth branches contain information matched to a GenMinInteraction
    //---------------------------------------------------------------------
    
    declareBoolTruthBranch("reco_isMinosMatch");
    declareIntTruthBranch("should_be_accepted", 0); // Inherited from Template
    
    return sc;
}

StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{
    debug() << "CC1P1PiAnalysis::reconstructEvent" << endmsg;
    
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
    debug()<< "1) Find vertex" << endmsg;
    
    if( !event->hasInteractionVertex() ){
        debug() << "No event vertex. Quitting..." << endmsg;
        //event->setIntData("vert_exists", 0);
        return StatusCode::SUCCESS;
    }
    //else{
        debug() << "Found vertex!" << endmsg;
        event->setIntData("vert_exists", 1);
    
    counter("c_vertex")++;
    SetAccumLevel(1);
    
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
    
    //Check if there are the same number of prongs as tracks:
    Minerva::ProngVect n_prong_check = event->primaryProngs();
    unsigned int n_prongs = n_prong_check.size();
    
    debug() << "n_tracks = " << ntot_tracks << " n_prongs = " << n_prongs;
    if(ntot_tracks == n_prongs){
        debug() << " !! EQUAL !!";
    }
    debug() << " " << endmsg;
    
    if(!(ntot_tracks == nout_tracks && ntot_tracks == 3)){
        debug() << "Event doesn't contain extactly three tracks." << endmsg;
        return StatusCode::SUCCESS;
    }
    
    counter("c_3tracks")++;
    SetAccumLevel(2);
    
    debug()<< "Has 3 tracks!" << endmsg;
    
    event->setIntData("n_tracks3", 3);
    event->setIntData("n_prongs", n_prongs);
    
    //----------- 3 : Muon track coming from common vertex -----------//
    debug()<< "3) Muon Track" << endmsg;
   // SmartRef<Minerva::Prong>    muonProng = (Minerva::Prong*)NULL;
   // SmartRef<Minerva::Particle> muonPart = (Minerva::Particle*)NULL;
    
    if(!FindMuon(event, truth, m_MuonProng, m_MuonParticle)){
        debug() << "Muon not found..." << endmsg;
        return StatusCode::SUCCESS;
    }
    debug()<< "Muon track found!" << endmsg;
    
    counter("c_muon_trk")++;
    SetAccumLevel(3);
    
    //----------- 4 : Vertex in active tracker or carbon target -----------//
    //Not a cut but an action to determine the location of the vertex.
    debug() << "4) Vertex in Carbon or Scintillator" << endmsg;
    if(VertIsIn("Scint", event)){
        debug() << "Yes in SCINTILLATOR" << endmsg;
        event->setIntData("target_region", 1);
        counter("c_tar_scint")++;

    }
    else if (VertIsIn("Carbon", event)){
        debug() << "Yes in CARBON TARGET" << endmsg;
        event->setIntData("target_region", 2);
        counter("c_tar_carbon")++;
    }
    else{
        debug() << "Event not in either..." << endmsg;
        event->setIntData("target_region", 3);//Probably don't need this...
        counter("c_tar_other")++;
        return StatusCode::SUCCESS;
    }
    SetAccumLevel(4);
    
    
    //----------- 5 : PID on p/pi+ -----------//
    debug() << "5) PID: p/pi+" << endmsg;
    
    //HadronSystem hadrons;
    
    bool tFinPar = FindParticles(event);
    
    if(tFinPar){
        debug() << "Finished Selection Successfully. Pheeewwww ;)" << endmsg;
    }
    
    SetAccumLevel(5);
    
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
    //debug() << "CC1P1PiAnalysis::interpretEvent" << endmsg;
    
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
    //debug() << "CC1P1PiAnalysis::tagTruth" << endmsg;
    
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
        
        if(materialZ == 6) return false;
        apothem = m_carbon_apothem;
        upZ = m_carbon_upZ;
        downZ = m_carbon_downZ;
    }
    else {
        debug() << "CC1P1PiAnalysis::VertIsIn : Could not determine target name: "
                << targetRegion.Data() << ". Please check, it may not be implemented." << endmsg;
        return false;
    }
    
    bool fidVertex = m_coordSysTool->inFiducial( vertex->position().x(), vertex->position().y(), vertex->position().z(), apothem, upZ, downZ );
    
    if(fidVertex){
        debug() << "Vertex is in fiducial volume" << endmsg;
    }
    
    return fidVertex;
}

bool CC1P1PiAnalysis::FindParticles(Minerva::PhysicsEvent* event) const
{
    debug() << "CC1P1PiAnalysis::FindParticles" << endmsg;
    //Determine which track is most proton like and pion like:
    // 1) Get particle scores and compare which track is has the highest score for the given hypothosis.
    // 2) Look for Michel features.
    //
    //Check that they are contianed in det FV and they are not minos matched.
   
    Minerva::ProngVect prongs = event->primaryProngs();
    Minerva::ProngVect::iterator prong;
    
    std::vector<double> protonScore;
    std::vector<double> protonChi2NDF;//May not need.
    std::vector<double> pionScore;
    std::vector<double> pionChi2NDF;//May not need.
    
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
        debug() << "Checking Prong: " << prong_count << "." << endmsg;
        
        //Check prong isn't that of the muon:
        if(m_MuonProng == (*prong)){
            debug() << "Prong already determined as muon." << endmsg;
            continue;
        }
        
        //Check if track is fully contained in detector FV and there are no minos matched tracks:
        Minerva::TrackVect tracks = (*prong)->minervaTracks();
        
        if( tracks.empty() ) {
            debug() << "  This prong contains an empty vector of tracks, skipping!" << endmsg;
            continue;
            //Return false statement if found to minos match.
        }
        else if( (*prong)->MinosTrack() || (*prong)->MinosStub() ) {
            debug() << "  This is a MINOS matched prong, skipping!" << endmsg;
            continue;
            //Return false statement if found to minos match.
        }
        
        hadron_counter++;
     
        SmartRef<Minerva::Track> track = tracks[tracks.size() -1];
        Gaudi::XYZPoint endpoint = track->lastState().position();
        
        if(!m_coordSysTool->inFiducial(endpoint.x(), endpoint.y(), endpoint.z(), m_det_apothem, m_det_upZ, m_det_downZ)){
            debug() << "Track not contained in detector fiducial volume." << endmsg;
            return false;
        }
        
        //The following code is based on that of the ProtonUtils.
        std::vector<Minerva::Particle::ID> hypotheses;
        hypotheses.push_back(Minerva::Particle::Pion);
        hypotheses.push_back(Minerva::Particle::Proton);
        IParticleMakerTool::NameAliasListType toolsToUse;
        toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );
        
        bool found_particle = m_particleMaker->makeParticles((*prong), hypotheses, toolsToUse);
        
        if(found_particle){
            debug() << "This prong has " << (*prong)->particles().size() << " particle hypotheses attached." << endmsg;
        }
        else{
            debug() << "Failed to produce particles" << endmsg;
        }
        
        if((*prong)->particles().size() == 2){
            
            trackChi2NDF.push_back(track->chi2PerDoF());
            
            Minerva::ParticleVect partHypVec = (*prong)->particles();
            Minerva::ParticleVect::iterator part;
            
            double pr_score  = 0.0;
            double pi_score  = 0.0;
            double score_den = 0.0;
            
            for(part = partHypVec.begin(); part != partHypVec.end(); part++){
                debug() << "Testing " << (*part)->idcode() << " hypothesis with signature: " << (*part)->methodSignature() << " and score: " << (*part)->score() << endmsg;
                
                std::string part_name;
                double minPartScore = -999.0;
                double maxPartChi2 = -999.0;
                
                //For now let's just compare the hyp of which is more proton/pion like
                if((*part)->idcode() == Minerva::Particle::Proton){
                    part_name = "Proton";
                    minPartScore = m_minProtonScore;
                    maxPartChi2 = m_maxProtonChi2;
                    debug() << "        Running checks on Proton Hypothesis." << endmsg;
                }
                else if((*part)->idcode() == Minerva::Particle::Pion){
                    part_name = "Pion";
                    minPartScore = m_minPionScore;
                    maxPartChi2 = m_maxPionChi2;
                    debug() << "        Running checks on Pion Hypothesis." << endmsg;
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
            
            debug() << "Prong " << hadron_counter <<":" << endmsg;
            debug() << "        Proton Score " << pr_score << " (" << pr_score_N <<") Pion Score " << pi_score << " (" << pi_score_N <<") Chi2NDF " << track->chi2PerDoF() << endmsg;
            
            int found_p_no = 0;
            Minerva::Particle::ID part_name_check;
            int PDGCode = -999;
            
            if(pr_score_N > pi_score_N){
                
                for(int hyp_counter = 0; hyp_counter < (int)partHypVec.size(); hyp_counter++){
                    if( partHypVec[hyp_counter]->idcode() == Minerva::Particle::Proton){
                        found_p_no = hyp_counter;
                    }
                }
                
                part_name_check = Minerva::Particle::Proton;
                PDGCode = m_Proton_PDG;
                
            }
            else if(pr_score_N < pi_score_N){
                for(int hyp_counter = 0; hyp_counter < (int)partHypVec.size(); hyp_counter++){
                    if( partHypVec[hyp_counter]->idcode() == Minerva::Particle::Pion){
                        found_p_no = hyp_counter;
                    }
                }
                
                part_name_check = Minerva::Particle::Pion;
                PDGCode = m_Pion_PDG;
            }
            
            debug() << "        Prong " << hadron_counter << " believed to be " << part_name_check << " and has found_p_no = " << found_p_no << endmsg;

            
            if(hadron_counter == 1){
                Part_from_Prong1 = partHypVec[found_p_no];
                Prong1_PDG = PDGCode;
                
                debug() << "        Part_from_Prong1 :: Consistent with " << part_name_check << " Hyp? ";
                if(Part_from_Prong1->idcode() == part_name_check){
                    debug() << "YES!!!!";
                }
                else{
                    debug() << "NO ********************** ?!";
                }
                debug() << " " << endmsg;
                debug() << "        IDCode: " << Part_from_Prong1->idcode() <<" Score: " << Part_from_Prong1->score() << endmsg;
                
            }
            
            if(hadron_counter == 2){
                Part_from_Prong2 = partHypVec[found_p_no];
                Prong2_PDG = PDGCode;
                
                debug() << "        Part_from_Prong2 :: Consistent with " << part_name_check << " Hyp?";
                if(Part_from_Prong2->idcode() == part_name_check){
                    debug() << "YES!!!!";
                }
                else{
                    debug() << "NO ********************** ?!";
                }
                debug() << " " << endmsg;
                debug() << "        IDCode: " << Part_from_Prong2->idcode() <<" Score: " << Part_from_Prong2->score() << endmsg;
            }
        }
        
        //Look for michels at end of the prong
    }
    
    //Given the particle hypotheses, set the candidate tracks:
    int pr_prong_no = -999;
    int pi_prong_no = -999;
    
    debug() << "******************************** Summary ********************************" << endmsg;
    debug() << "Vector Sizes Consistent: trackChi2NDF N = " << trackChi2NDF.size() << " protonScore N = " << protonScore.size() << " pionScore N = " << pionScore.size();
    if(trackChi2NDF.size() == 2 && protonScore.size()  == 2 && pionScore.size() == 2){
        debug() << " YES." << endmsg;
        
        double Prong1_Proton = protonScore[0]/(protonScore[0] + pionScore[0]);
        double Prong1_Pion = pionScore[0]/(protonScore[0] + pionScore[0]);
        
        double Prong2_Proton = protonScore[1]/(protonScore[1] + pionScore[1]);
        double Prong2_Pion = pionScore[1]/(protonScore[1] + pionScore[1]);
    
        debug() << "Prong 1:" << endmsg;
        debug() << "        Proton Score " << protonScore[0] << " (" << Prong1_Proton <<") Pion Score " << pionScore[0] << " (" << Prong1_Pion <<") Chi2NDF " << trackChi2NDF[0] << endmsg;
        debug() << "        PreCal Pr Sc N " << pr_score_N_Vec[0] << " PreCal Pi Sc N " << pi_score_N_Vec[0] << endmsg;
        debug() << "Prong 2:" << endmsg;
        debug() << "        Proton Score " << protonScore[1] << " (" << Prong2_Proton <<") Pion Score " << pionScore[1] << " (" << Prong2_Pion <<") Chi2NDF " << trackChi2NDF[1] << endmsg;
        debug() << "        PreCal Pr Sc N " << pr_score_N_Vec[1] << " PreCal Pi Sc N " << pi_score_N_Vec[1] << endmsg;
        
        debug() << "*************************************************************************" << endmsg;
        
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
        
        debug() << "Checking PDG Codes:" << endmsg;
        debug() << "Prong 1: Pre:" << Prong1_PDG << " Post " << Prong1a_PDG;
        if(Prong1_PDG == Prong1a_PDG) debug() << ". They are the same!!!";
        else debug() << ". Close but no cigar... :-(";
        debug() << " " << endmsg;
        
        debug() << "Prong 2: Pre:" << Prong2_PDG << " Post " << Prong2a_PDG;
        if(Prong2_PDG == Prong2a_PDG) debug() << ". They are the same!!!";
        else debug() << ". Close but no cigar... :-(";
        debug() << " " << endmsg;
        
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
            debug() << "Found two protons..." << endmsg;
            return false;
        }
        else{
            debug() << "Found two pions..." << endmsg;
            return false;
        }
    }
    else{
        debug() << " No... Check this out!!!!" << endmsg;
        return false;
    }

    m_ProtonProng = prongs[pr_prong_no];
    m_PionProng = prongs[pi_prong_no];
    
    bool pr_is_correct = false;
    bool pi_is_correct = false;
    debug() << "Final Check that the tracks are what we think they are" << endmsg;
    debug() << "ProtonParticle: " << m_ProtonParticle->idcode();
    if(m_ProtonParticle->idcode() == Minerva::Particle::Proton){
        debug() << " YES" << endmsg;
        pr_is_correct = true;
    }
    debug() << "PionParticle: " << m_PionParticle->idcode();
    if(m_PionParticle->idcode() == Minerva::Particle::Pion){
        debug() << " YES" << endmsg;
        pi_is_correct = true;
    }
    
    if(!(pr_is_correct || pi_is_correct)){
        debug() << "Particles not correct... check code" << endmsg;
        return false;
    }
    
    debug() << "Found Proton and Pion Tracks" << endmsg;
    
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
}

//Generic Particle information builder:
void CC1P1PiAnalysis::SetPartInfo(std::string name)
{
    declareDoubleBranch(m_hypMeths, (name + "_score").c_str() , -999.);
    declareDoubleBranch(m_hypMeths, (name + "_chi2ndf").c_str(), -999.);
    
    //mom., energy, sel. pdg, true pdg.
    
    
}


void CC1P1PiAnalysis::FillPartInfo(SmartRef<Minerva::Prong> prong, SmartRef<Minerva::Particle> particle)
{
    
}

mutable void CC1P1PiAnalysis::SetAccumLevel(int cut)
{
    m_accum_level[ cut - 1 ] = 1;
}

mutable void CC1P1PiAnalysis::ResetAccumLevel()
{
    for(int i = 0; i < m_ncuts; i++){
        m_accum_level[i] = 0;
    }
}

