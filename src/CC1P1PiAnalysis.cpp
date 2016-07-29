//include the header file...
#include "CC1P1PiAnalysis.h"

#include "Event/TG4PrimaryTrajectory.h"

//Forward declared headers:
#include "AnaUtils/IMuonUtils.h"
#include "GeoUtils/IMinervaCoordSysTool.h"
#include "DetDesc/Material.h"
#include "GeoUtils/INuclearTargetTool.h"
#include "AnaUtils/IProtonUtils.h"
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
    
    //For Proton PID:
    declareProperty( "ProtonUtilsAlias", m_protonUtilsAlias = "CC1P1PiProtonUtils");
    declareProperty( "ProtonScoreThreshold", m_protonScoreThreshold = 0.0 );
    
    //For hadron PID:
    declareProperty("det_apothem", m_det_apothem = 1200.0*CLHEP::mm);//Same as Proton utils:
    declareProperty("det_upZ", m_det_upZ = 4000.0*CLHEP::mm);//Same as Proton utils:
    declareProperty("det_downZ", m_det_downZ = 10000.0*CLHEP::mm);//Same as Proton utils:
    
    declareProperty("ParticleMakerAlias", m_particleMakerAlias = "CC1P1PiParticleMaker");
    
    //These values are taken from ProtonUtils:
    declareProperty("minProtonScore", m_minProtonScore = 50.);
    declareProperty("maxProtonChi2",  m_maxProtonChi2 = 0.05);
    
    //Need to determine better scores (currently same as those in ProtonUtils for Proton Hyp:
    declareProperty("minPionScore", m_minPionScore = 50.);
    declareProperty("maxPionChi2",  m_maxPionChi2 = 0.05);

    
    
    
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
    
    try{ m_protonUtils = tool<IProtonUtils>("ProtonUtils", m_protonUtilsAlias); }
    catch( GaudiException& e ){
        error() << "Could not obtain ProtonUtils: " << m_protonUtilsAlias << endmsg;
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
    
    
    //----------- 5 : PID on p/pi+ -----------//
    debug() << "5) PID: p/pi+" << endmsg;
    
    //HadronSystem hadrons;
    
    bool tFinPar = FindParticles(event);
    
    if(tFinPar){
        debug() << "Okay" << endmsg;
    }
    
    
    
    
    
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

//May not need this:
bool CC1P1PiAnalysis::getProton( const Minerva::ProngVect& primaryProngs, SmartRef<Minerva::Prong>& protonProng, SmartRef<Minerva::Particle>& protonPart ) const
{
    
    if( m_protonUtils->findProtonProng( primaryProngs, protonProng, protonPart ) ) {
        if( !protonProng ) {
            warning() << "Identified a proton Prong, but it is NULL in CC1P1PiAnalysis::reconstructEvent!" << endmsg;
            return false;
        }
        //! Check that the proton particle is well identified. Tag the proton prong
        //! if there is one and only one Proton & it is of sufficient score.
        debug() << " Proton Particle Score: " << protonPart.data()->score() << endmsg;
        if( m_protonScoreThreshold < protonPart.data()->score() ) {
            debug() << "  Tagging the proton Prong" <<  endmsg;
            protonProng.data()->filtertaglist()->addFilterTag( "PrimaryProton", true );
        }
        else {
            debug() << "Proton Particle Score is below threshold - not tagging prong, proton score :  " << protonPart.data()->score() << endmsg;
        }
    }
    else {
        debug() << "Did not find a proton prong.  This cannot be a CCQE event." << endmsg;
        return false;
    } 
    
    return true;
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
    
    //SmartRef<Minerva::Prong>
    
    int prong_count = 0;
    
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
            
            for(part = partHypVec.begin(); part != partHypVec.end(); part++){
                debug() << "Found a " << (*part)->idcode() << " with signature: " << (*part)->methodSignature() << " and score: " << (*part)->score() << endmsg;
                
                std::string part_name;
                double minPartScore = -999.0;
                double maxPartChi2 = -999.0;
                
                //For now let's just compare the hyp of which is more proton/pion like
                if((*part)->idcode() == Minerva::Particle::Proton){
                    part_name = "Proton";
                    minPartScore = m_minProtonScore;
                    maxPartChi2 = m_maxProtonChi2;
                    debug() << "Found a proton particle hypothesis!" << endmsg;
                }
                else if((*part)->idcode() == Minerva::Particle::Pion){
                    part_name = "Pion";
                    minPartScore = m_minPionScore;
                    maxPartChi2 = m_maxPionChi2;
                    debug() << "Found a pion particle hypothesis!" << endmsg;
                }
                
               /* if(minPartScore == -999.0 || maxPartChi2 == -999.0){
                    error() << part_name << " particle requirements not changed from default (-999.0): ";
                    if(minPartScore == -999.0){
                        error() << "minPartScore ";
                    }
                    else if(maxPartChi2 == -999.0){
                        if(minPartScore == -999.0){
                            error() << "and ";
                        }
                        error() << "maxPartChi2 ";
                    }
                    error() << "unchanged..." << endmsg;
                }*/
                
                //Actual PID bit: Does the particle
                if( (*part)->methodSignature().find("dEdX") != std::string::npos ) {
                    
                    if((*part)->idcode() == Minerva::Particle::Proton){
                        protonScore.push_back((*part)->score());
                    }
                    else  if((*part)->idcode() == Minerva::Particle::Pion){
                        pionScore.push_back((*part)->score());
                    }
                    
                    /*if ( track->chi2PerDoF() < maxPartChi2 && (*part)->score() > minPartScore ) {
                        passProtonCut = true;
                    }
                    else {
                        debug() << "Chi2/DoF (" << track->chi2PerDoF() << ") and/or score (" << (*part)->score() << ") is not consistent with a " << part_name << endmsg;
                        debug() << "Max Chi2/DoF = " << maxPartChi2 << " Min Score = " << minPartScore << endmsg;
                    }*/
                }
                /*else {
                    debug() << "  Method signature is not consistent with a " << part_name << endmsg;
                }*/
                
                
            }
        }
    
        //Determine particle scores:
        /*double tmp_pr_sc = -999.0;
        double tmp_pi_sc = -999.0;
        bool pscores = m_protonUtils->getParticleScores( (*prong), tmp_pr_sc, tmp_pi_sc);
        if(pscores){
            debug() << "Particle scores found to be: Proton = " << tmp_pr_sc << " Pion = " << tmp_pi_sc << endmsg;
            protonScore.push_back(tmp_pr_sc);
            pionScore.push_back(tmp_pi_sc);
        }*/
        
        //Look for michels at end of the prong
    }
    
    //Given the particle hypotheses, set the candidate tracks:
    
    debug() << "**** Summary ****" << endmsg;
    debug() << "Vector Sizes Consistent: trackChi2NDF N = " << trackChi2NDF.size() << " protonScore N = " << protonScore.size() << " pionScore N = " << pionScore.size();
    if(trackChi2NDF.size() == 2 && protonScore.size()  == 2 && pionScore.size() == 2){
        debug() << " YES." << endmsg;
    
        debug() <<"Prong 1:" << endmsg;
        debug() <<"        Proton Score " << protonScore[0] << " Pion Score " << pionScore[0] << " Chi2NDF " << trackChi2NDF[0] << endmsg;
        debug() <<"Prong 2:" << endmsg;
        debug() <<"        Proton Score " << protonScore[1] << " Pion Score " << pionScore[1] << " Chi2NDF " << trackChi2NDF[1] << endmsg;
    }
    // m_ProtonParticle;
    // m_ProtonProng;
    
    // m_PionProng;
    //m_PionParticle;

    
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
