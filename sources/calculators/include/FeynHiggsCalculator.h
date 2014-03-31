/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        FeynHiggsCalculator.h                                        *
*                                                                              *
* Description Wrapper class for FeynHiggs                                      *
*                                                                              *
* Authors     Bjoern Sarrazin  <sarrazin@physik.uni-bonn.de>                   *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_FEYNHIGGSCALCULATOR_H
#define FITTINO_FEYNHIGGSCALCULATOR_H

#include "FeynHiggsTypes.h"
#include "CalculatorBase.h"
#include "Collection.h"
#include "PtreeForwardDeclaration.h"


/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  class SimplePrediction;
    class FeynHiggsChannel;
    class SLHADataStorageBase;

  /*!
   *  \ingroup calculators
   *  \brief Wrapper class for FeynHiggs.
   */
  class FeynHiggsCalculator : public CalculatorBase {

    public:
      /*!
       *  Standard constructor.
       */
      FeynHiggsCalculator( const PhysicsModel* model, const boost::property_tree::ptree& ptree, std::string inputMethod );
      /*!
       *  Standard destructor.
       */
      virtual ~FeynHiggsCalculator();

      void CalculatePredictions();

    protected:
      int         _error;
      std::string _fileName;

      /*! \cond UML */

      double _mass_h0;
      double _mass_H0;
      double _mass_A0;
      double _mass_Hp;

  private:
      void AddChannel( std::string higgsName, std::string channelName, int channelNumber, bool fermionic, bool SM );
      void AddChannels_H0SfSf  ( int iHiggs, std::string higgsName, unsigned int type, std::string* names );
      void AddChannels_HpSfSf  ( unsigned int type, std::string* names1, std::string* names2 );
      void AddChannels_HpFF    ( unsigned int type, std::string* names1, std::string* names2 );
      void AddChannels_H0FF    ( unsigned int iHiggs, std::string higgsName, unsigned type, std::string* names );
      void AddChannels_H0HH    ( unsigned int iHiggs, std::string higgsName );
      void AddChannels_H0VV    ( unsigned int iHiggs, std::string higgsName );
      void AddChannels_H0HV    ( unsigned int iHiggs, std::string higgsName );
      void AddChannels_H0NeuNeu( unsigned int iHiggs, std::string higgsName );
      void AddChannels_H0ChaCha( unsigned int iHiggs, std::string higgsName );
      void AddChannels_HpHV();
      void AddChannels_HpNeuCha();

      std::vector< FeynHiggsChannel* > _channels;
      FHRealType*    _gammas;
      FHRealType*    _gammasms;
      FHComplexType* _couplings;
      FHComplexType* _couplingsms;

      std::string _inputMethod;

      double _normSM_sigma_ggh;
      double _normSM_sigma_ggh_2;
      double _normSM_sigma_bbh;
      double _normSM_sigma_qqh;
      double _normSM_sigma_tth;
      double _normSM_sigma_Wh;
      double _normSM_sigma_Zh;

      std::string _higgs [5];
      std::string _nu    [4];
      std::string _lepton[4];
      std::string _up    [4];
      std::string _down  [4];
      std::string _neu   [5];
      std::string _cha   [3];

      Collection<SimplePrediction*> _input;

      SLHADataStorageBase* _slhadatastorageSPheno;
      SLHADataStorageBase* _slhadatastorageFeynHiggs;

      bool _smallObsSet;


      /*! \endcond UML */

  };


}

#endif
