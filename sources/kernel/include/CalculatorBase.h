/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        CalculatorBase.h                                                 *
*                                                                              *
* Description Base class for calculators                                       *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_CALCULATORBASE_H
#define FITTINO_CALCULATORBASE_H

#include <boost/property_tree/ptree_fwd.hpp>

#include "Collection.h"
#include "VariableBase.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  class FormulaQuantity;
  class ModelBase;
  class Quantity;

  /*!
   *  \ingroup kernel
   *  \brief Base class for calculators.
   */
  class CalculatorBase {

    public:
      /*!
       *  Standard constructor.
       */
      CalculatorBase( const ModelBase* model, const boost::property_tree::ptree* ptree = 0 );
      /*!
       *  Standard destructor.
       */
      const std::string&                      GetName() const;
      const Collection<Quantity*>&            GetCollectionOfQuantities() const;
      const Collection<const VariableBase<std::string>*>& GetCollectionOfStringVariables() const;

    public:
      virtual                                 ~CalculatorBase();
      bool                                    Calculate();
      virtual void                            Initialize();
      /*!
       *  \todo Remove when no longer used by derived classes (Matthias).
       */
      virtual void                            SetupMeasuredValues();

    protected:
      bool _requirementsFulfilled;
      std::string                             _name;
      std::string                             _tag;
      const ModelBase*                        _model;

    protected:
      const double&                           GetInput( std::string name ) const;
      const double&                           GetOutput( std::string name ) const;
      void                                    AddInput( std::string name );
      void                                    AddInput( std::string name, std::string value );
      void                                    AddOutput( std::string name, const double& value );
      void                                    AddOutput( std::string name );
      void                                    AddQuantity( Quantity* prediction );
      void                                    AddStringVariable( const std::string& name, const std::string& value);

      void                                    PrintInput() const;

      void                                    SetOutput( std::string name, const double& value );
      void                                    SetName( std::string name );
    
      void                                    UpdateInput();
      const boost::property_tree::ptree*      GetConfiguration() const;
      void                                    AddFile( std::string name );
      void                                    RemoveFiles();
      virtual void                            CalculatePredictions() = 0;


      /*!; \cond UML */
    private:
      std::string                             _className;
      std::map<std::string, FormulaQuantity*> _input;
      std::map<std::string, Quantity* >       _settableOutput;
      const boost::property_tree::ptree*      _ptree;
      Collection<Quantity*>                   _collectionOfQuantities;
      Collection<const VariableBase<std::string>*>        _collectionOfStringVariables;
      std::vector<std::string>                _outFiles;


      template<class T>
      void ModifyVariableName(VariableBase<T>* );
      void                                    AddStringVariable( VariableBase<std::string>* variable );
      /*! \endcond UML */

  };

}

#endif // FITTINO_CALCULATORBASE_H
