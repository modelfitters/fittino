#ifndef FITTINO_SMODELSCALCULATOR_H
#define FITTINO_SMODELSCALCULATOR_H

#include "CalculatorBase.h"
#include "boost/python.hpp"

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
      boost::python::object _fileList;
      boost::python::object _listOfExpRes;
      boost::python::object _testPoints;

      void ReadXML();

  };

}

#endif
