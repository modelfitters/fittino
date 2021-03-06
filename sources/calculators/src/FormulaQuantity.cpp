/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        FormulaQuantity.cpp                                              *
*                                                                              *
* Description                                                                  *
*                                                                              *
* Authors     Bjoern Sarrazin  <sarrazin@physik.uni-bonn.de>                   *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include<set>

#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "TFormula.h"

#include "FormulaQuantity.h"
#include "ModelBase.h"


Fittino::FormulaQuantity::FormulaQuantity(std::string name, std::string formula, const Fittino::ModelBase *model, std::map<std::string, const Quantity*> map)
: Quantity( name ){

    for ( const auto& kv : map ) {

        if ( model->GetCollectionOfQuantities().Count( kv.first ) ) {

            throw ConfigurationException( "Reserved name " + kv.first + " used for a model quantity." );

        }

    }

    int begin = formula.find( "[" );

    std::set<std::string> set;

    while ( begin != std::string::npos ) {

        begin+=1;
        int end = formula.find( "]", begin );
        std::string substring = formula.substr( begin, end - begin  );
        set.insert( substring );
        begin = formula.find( "[", end );

    }

    std::set<std::string>::iterator it;

    for ( it = set.begin(); it != set.end(); ++it ) {

        std::string index = boost::lexical_cast<std::string>( _parameters.size() );
        boost::replace_all( formula, "[" + *it + "]", "[" + index + "]" );

        if ( map.count(*it) ) {

            _parameters.push_back( map.at( *it ) );

        }
        else {

            _parameters.push_back( model->GetCollectionOfQuantities().At(*it) );

        }

    }

    _formula = new TFormula( "", formula.c_str() );

    if ( ( _formula->GetNdim() != 1 && _formula->GetNdim()!= 0 ) || _formula->GetNpar() != set.size() ) {

        throw ConfigurationException( "Invalid formula." );

    }

    Update();

}


Fittino::FormulaQuantity::~FormulaQuantity() {

    delete _formula;

}

void Fittino::FormulaQuantity::Update() {

    for ( unsigned int i = 0; i < _parameters.size(); i++ ) {

        _formula->SetParameter( i, _parameters[i]->GetValue() );

    }

    _value = _formula->Eval( 0 );

}
