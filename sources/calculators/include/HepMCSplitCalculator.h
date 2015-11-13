#ifndef FITTINO_HEPMCSPLITCALCULATOR_H
#define FITTINO_HEPMCSPLITCALCULATOR_H

#include <set>
#include <vector>
#include <string>
#include "boost/filesystem.hpp"

#include "CalculatorBase.h"

namespace Fittino {
  class HepMCSplitCalculator : public CalculatorBase {
    
  public:

    HepMCSplitCalculator( const ModelBase*, const boost::property_tree::ptree& ptree );
    ~HepMCSplitCalculator();

  public:

    virtual void CalculatePredictions();

 private:
    std::set<int> _squarks;
    boost::filesystem::path _hepMCFile;
    std::string _rootFile;
    std::string _tree;
    std::string _eventNumber;
    std::vector<std::string> _pdgIDs;

  };
}

#endif 
