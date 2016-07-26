//include the header file...
#include "CC1P1PiAnalysis.h"

#include "Event/TG4PrimaryTrajectory.h"
#include <AnaUtils/IMuonUtils.h>

//this command allows other parts of Gaudi to use the tool
DECLARE_TOOL_FACTORY( CC1P1PiAnalysis );

CC1P1PiAnalysis::CC1P1PiAnalysis(const std::string& type, const std::string& name, const IInterface* parent ) : MinervaAnalysisTool( type, name, parent )
{
    // Declare the interface so other tools/algs can get this tool via an IInteractionHypothesis
    declareInterface<IInteractionHypothesis>(this);
    
    // Mandatory declaration of analysis signature:
    m_anaSignature = "CC1P1Pi";
    
    /*
     The m_hypMeths variable is from MinervaAnalysisTool and makes many things easier.
     The interpretations correspond to the methodSignatures of NeutrinoInts you will allow to be created.
     Push back the default interpretations you will allow.
     */
    /*
     I will include hypothesis methods by default:
     <ul>
     <li> InterpretationA
     <li> InterpretationB
     </ul>
     */
    m_hypMeths.push_back( "CC1P1Pi" );
    //m_hypMeths.push_back( "InterpretationB" );
    declareProperty("HypothesisMethods", m_hypMeths);
    
    // Declare other properties you can set from an options file.
    declareProperty( "SomeProperty", m_someProperty = 1 );
}

//! Initialize
StatusCode CC1P1PiAnalysis::initialize()
{
    info() << "DAVID: CC1P1PiAnalysis::initialize()" << endmsg;
    
    // Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() )
        return Error( "Failed to initialize!", sc );
    
    /*
     Start declaring branches you want to get in your analysis DST.
     These branches go into the AnaTuple associated with the PhysicsEvents that this tool \b reconstructs.
     See <a href="http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=6855">docdb6855</a> for more detail.
     */
    
    //---------------------------------------------------------------------
    // Declare the PhysicsEvent block branches
    //---------------------------------------------------------------------
    // If you add extra int data to the PhysicsEvent it will show up in your analysis DST.
    
    // This creates a double branch called time_width
    declareDoubleEventBranch( "time_width", -1 );
    //double branch linked with PhysicsEvent.
    // branch will be "time_width/D"
    
    // This branch is called n_orig_prongs and has a default value of -1.
    declareIntEventBranch( "n_orig_prongs", -1 );
    //int branch with default linked with PhysicsEvent ( "n_orig_prongs/I" )
    
    // You can also add container branches.
    // Here is a containerDouble branch called shower_momentum of constant size 4 with default value -999 for each element
    //  declareContainerDoubleEventBranch( "shower_momentum", 4, -999. );
    //constant-sized container double branch linked with PhysicsEvent ( "shower_momentum[4]/D" )
    
    // Here is a containerInt branch called hit_module with variable size.
    // Size will be set by a branch named hit_module_sz (created for you).
    //  declareContainerIntEventBranch( "hit_module" );
    //variable-sized container int branch linked with PhysicsEvent.
    // branches will be "hit_module_sz/I" and "hit_module[hit_module_sz]/I"
    
    //---------------------------------------------------------------------
    // Declare the Interpretations block branches
    //---------------------------------------------------------------------
    // If you add extra containerDouble data to a NeutrinoInt it will show up in your analysis DST
    
    // The size of the iso_blob_energy branch is controlled by the n_iso_blobs branch, which you should not manually create.
    declareContainerDoubleBranch( m_hypMeths, "iso_blob_energy", "n_iso_blobs" );
    // variable-sized containerDouble branches for each hypMeth.
    // branches will be "InterpretationA_iso_blob_energy[n_iso_blobs]/D", "InterpretationB_..." )
    // a branch "n_iso_blobs/I" will be created automatically if it does not already exist.
    
    // Another example of variable length container branch using the same sizer: n_iso_blobs
    declareContainerIntBranch( m_hypMeths, "iso_blob_nclusters", "n_iso_blobs" );
    // another variable-sized container double branch ( "<hypmeth>_iso_blob_nclusters[n_iso_blobs]/I" )
    
    //---------------------------------------------------------------------
    // Declare the Truth block branches.
    // Truth branches contain information matched to a GenMinInteraction
    //---------------------------------------------------------------------
    // If a GenMinInteraction has 'should_be_accepted' in the extra int data it will show up in you analysis DST
    declareIntTruthBranch( "should_be_accepted", 0 );
    // int branch linked to GenMinInteraction ( "truth_should_be_accepted/I" )
    
    return sc;
}

//! Reconstruct a PhysicsEvent
/*!
 You can:
 Break down unwanted prongs, tracks, blobs etc.
 Create new prongs, tracks, etc.
 Promote prongs to be primary prongs.
 */
StatusCode CC1P1PiAnalysis::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth /* = NULL */ ) const
{
    debug() << "CC1P1PiAnalysis::reconstructEvent( PhysicsEvent *event, GenMinInteraction* truth )" << endmsg;
    
    // Add your own extra data to the PhysicsEvent.  Remember that you declared branches so that extra data could be written to an analysis DST...
    // Do the same to Prongs and Particles
    event->setIntData( "n_orig_prongs", event->primaryProngs().size() );
    
    // You can also tag the GenMinInteraction with any special truth matching stuff here
    std::vector<int> intData;
    SmartRefVector<Minerva::TG4Trajectory> truePrimaries = truth->trajectories();
    for( SmartRefVector<Minerva::TG4Trajectory>::iterator i = truePrimaries.begin(); i != truePrimaries.end(); ++i )
    {
        Minerva::TG4PrimaryTrajectory* traj = dynamic_cast<Minerva::TG4PrimaryTrajectory*>( static_cast<Minerva::TG4Trajectory*>(*i) );
        if( traj )
        {
            //if you can find a reconstructed track for this...
            intData.push_back( 1 );
            //else
            //intData.push_back(0);
        }
    }
    truth->setContainerIntData("tracked_FSPart",intData);
    intData.clear();
    
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



//! Attach interpretation to the event if it passed the filters
/*!
 Create NeutrinoInts and attach analysis data to them.
 Here are the NeutrinoInt signatures I use:
 <table border="1">
 <tr>  <th>signature</th> <th>who gets it</th>  </tr>
 <tr> <td>InterpretationA</td>  <td>Some explanation of what this interpretation means</td> </tr>
 <tr> <td>InterpretationB</td>  <td>Some explanation of what this interpretation means</td> </tr>
 </table>
 
 */
StatusCode CC1P1PiAnalysis::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *interaction, std::vector<Minerva::NeutrinoInt*>& nuInts ) const
{
    debug() << "CC1P1PiAnalysis::interpretEvent( const PhysicsEvent*, const GenMinInteraction*, std::vector<NeutrinoInt*>& )const" << endmsg;
    
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

