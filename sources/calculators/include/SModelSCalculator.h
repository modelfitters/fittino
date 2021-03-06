#ifndef FITTINO_SMODELSCALCULATOR_H
#define FITTINO_SMODELSCALCULATOR_H

#include "boost/python.hpp"
#include "CalculatorBase.h"
#include <set>

namespace Fittino {

    class Executor;

  class SModelSCalculator : public CalculatorBase {
    
  public:
    
    SModelSCalculator( const ModelBase* model, const boost::property_tree::ptree& ptree);
    ~SModelSCalculator();

  public:
    
    virtual void CalculatePredictions();
    virtual void Initialize();

  private:

      Executor* _crossSections_LO;
      Executor* _crossSections_NLL;

      std::string _xmlFile;
      std::string _fileName;
      std::string _databaseVersion;
      std::string _parameterFile;
      double _rValue;
      boost::python::object _parser;
      boost::python::object _inputFiles;
      boost::python::object _listOfExpRes;
      boost::python::object _testPoints;
      double _decompositionStatus;
      double _fileStatus;


      std::set<std::string> _txNamesWithResults;

      std::vector<std::string> _missingModels_TxNames;
      std::vector<std::string> _missingModels_Brackets;
      std::vector<double> _missingModels_Weights_Total;
      std::vector<double> _missingModels_Fractions_OutsideGrid;
      std::vector<double> _missingModels_Fractions_InsideGrid;
      unsigned int _numberOfMissingModelsConsidered;
      double _numberOfMissingModelsDetermined;
      
      std::vector<std::string> _constraintsMissing_Brackets;
      std::vector<double> _constraintsMissing_Weights;
      unsigned int _constraintsMissing_NumberConsidered;
      double _constraintsMissing_NumberDetermined;
      
      std::vector<std::string> _constraintsOutsideGrid_Brackets;
      std::vector<double> _constraintsOutsideGrid_Weights;
      unsigned int _constraintsOutsideGrid_NumberConsidered;
      double _constraintsOutsideGrid_NumberDetermined;

      void ReadXML();
      void ReadMissingConstraints();
      void ReadConstraintsOutsideGrid();


  };

}

#endif
