#ifndef TEMPLATEANALYSISTOOL_H 
#define TEMPLATEANALYSISTOOL_H 1

//I inherit thee
#include "AnaUtils/MinervaAnalysisTool.h"

//! This class is an example implementation of MinervaAnalysisTool
/*!
  It doesn't do much.  It show you the bare minumum you must do to create an analysis tool.
  If you want to see what you tool could do, then look at existing analysis tool.
  There are comments about why things are needed.
  @sa NukeCCReconstructor
  @sa CCInclusiveReco
  */
class TemplateAnalysisTool : public MinervaAnalysisTool
{
 public:

  //! Standard constructor
  TemplateAnalysisTool( const std::string& type, const std::string& name, const IInterface* parent );

  //! Destructor (mandatory for inheritance)
  ~TemplateAnalysisTool(){};

  //! Initialize (mandatory for inheritance)
  StatusCode initialize();
  //! Finalize (mandatory for inheritance)
  StatusCode finalize();

  //! Reconstruct the event (mandatory for inheritance)
  StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const;

  //! Attach an interpretations to the event (mandatory for inheritance)
  StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, std::vector<Minerva::NeutrinoInt*>& nuInts ) const;

  //! Tag truth via GenMinInteraction
  StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;

  bool truthIsPlausible( const Minerva::PhysicsEvent * event ) const;

 protected:


 private:
  //===== BEGIN Tool Properties =====//
  /** @name Properties
    Properties that you can set from an options file.
    @{*/
  /*! @brief An example of something you can set from the options file
    @property{SomeProperty,"theDefaultValue"}

    The only reason I put this in is to encourage doxygen comments.
   */
  int m_someProperty;

  //!@}
  //===== END Tool Properties =====//
};

#endif

