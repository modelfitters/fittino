
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include "TMath.h"
#include "ClsCalculator_2.h"
#include "ConfigurationException.h"
#include "PhysicsModel.h"
#include "ModelParameter.h"
#include "Factory.h"
#include "SimplePrediction.h"
#include <iostream>
#include <fstream>
#include <CalculatorException.h>
#include <boost/regex.hpp>
#include "CalculatorBase.h"
#include <cstdlib>
#include <set>
#include <algorithm>
#include <functional>

std::map<std::string, double> s0;
std::map<std::string, double> ds;
std::map<std::string, double> ds_sys;
std::map<std::string, double> s95exp;
std::map<std::string, double> rvalue;
std::map<std::string, double> b0;
std::map<std::string, double> db;
std::map<std::string, double> n0;
std::map<std::string, double> clsexp;
std::vector<std::pair<std::string, double>> clsexp_vec;
std::vector<std::pair<std::string, double>> clsexp_sel;
std::vector<std::pair<std::string, double>> s0_sel;
std::vector<std::pair<std::string, double>> ds_sel;
std::vector<std::pair<std::string, double>> ds_sys_sel;
std::vector<std::pair<std::string, double>> b0_sel;
std::vector<std::pair<std::string, double>> db_sel;
std::vector<std::pair<std::string, double>> n0_sel;
std::vector<std::pair<std::string, double>> s95exp_sel;
std::vector<std::pair<std::string, double>> rvalue_sel;
std::vector<std::pair<std::string, double>> rvalue_vec;

int zero_signal = 0;

Fittino::ClsCalculator_2::ClsCalculator_2( const ModelBase* model, const boost::property_tree::ptree& ptree ) 
  : SLHACalculatorBase( model, &ptree ),
  
    _nr_SR (ptree.get<std::string>("number_SR")),    
    _energy (ptree.get<std::string>("energy")),    
    _nr_toys (ptree.get<std::string>("number_toys")),    
    _CM_tot_result (ptree.get<std::string>( "results" ) ),
    _cls_script ( ptree.get<std::string> ("python_script") ),    
    _cls_python_result ( ptree.get<std::string> ("python_result") ),  
    _ORG ( ptree.get<std::string> ("orthogonal_group") )  

{
    AddOutput("Cls");
    AddOutput("Chi_2");
    AddOutput("Mu_hat");
}

Fittino::ClsCalculator_2::~ClsCalculator_2() {
  
}

void Fittino::ClsCalculator_2::CalculatePredictions() {
 std::cout << "Beginn" << std::endl;

 ReadFile();
 if (_energy== "8"){
    SetupMeasuredValues();}
 if (_energy== "13"){
    NewMethod();}
 if(zero_signal ==0){
 ReadPythonResult();
}
}
void Fittino::ClsCalculator_2::ReadFile() {

 if (!boost::filesystem::exists (_CM_tot_result)){
        throw CalculatorException(_name, "Cls_value_Error");
    }

 else{

        std::string dummy, temp1, temp2;
        std::string line;
        double temp3 = 0, temp4=0, temp5=0, temp6=0, temp7=0, temp8=0, temp9=0, temp10=0, temp11=0, temp12=0, temp13=0, temp14=0, temp15=0, temp16=0, temp17=0, temp18=0;
        std::vector<std::string> analysis;
        std::vector<std::string> sr;
        std::vector<double> si_norm_ev;
        std::vector<double> si_er_stat;
        std::vector<double> si_er_sys;
        std::vector<double> cls_exp;
        std::vector<double> s95_exp;
        std::vector<double> cls_exp_err;
        std::vector<double> cls_obs;
        std::vector<double> cls_obs_err;
        std::vector<double> bkg;
        std::vector<double> bkg_err;
        std::vector<double> obs;
        std::vector<double> r_value;
        std::ifstream infile(_CM_tot_result);
        int n = 0;
        while(std::getline(infile, line)){
         if(n>=1){
         if(line.empty())  break ;

 //        double temp3 = 0, temp4=0, temp5=0, temp6=0, temp7=0, temp8=0, temp9=0, temp10=0, temp11=0, temp12=0, temp13=0, temp14=0, temp15=0, temp16=0, temp17=0, temp18=0, temp19=0, temp20=0, temp21=0, temp22=0, temp23=0,temp24=0, temp25=0, temp26=0, temp27=0, temp28=0, temp29=0, temp30=0, temp31=0, temp32=0, temp33=0;

         std::istringstream ss(line);
            ss>>temp1>>temp2>>temp3>>temp4>>dummy>>dummy>>dummy>>dummy>>temp5>>temp6>>temp7>>temp8>>temp11>>temp12>>temp13>>dummy>>dummy>>dummy>>dummy>>dummy>>temp18>>dummy>>temp10>>dummy>>dummy>>dummy>>dummy>>temp14>>temp15>>temp16>>temp17>>temp9;
            analysis.push_back(temp1);
            sr.push_back(temp2);
            si_norm_ev.push_back(temp5);
            si_er_stat.push_back(temp6);
            si_er_sys.push_back(temp7);
//            r_value.push_back(temp10);
            r_value.push_back(temp5/temp18);
            bkg.push_back(temp12);
            bkg_err.push_back(temp13);
            obs.push_back(temp11); 
            cls_obs.push_back(temp14);
            cls_obs_err.push_back(temp15);
            cls_exp.push_back(temp16);
            s95_exp.push_back(temp18);
            cls_exp_err.push_back(temp17);     
            
           }
            n+=1;
}
        infile.close();


      std::string output3;
      std::string output4;



int p = 0;
for(const auto& elem : analysis ) {
    output3 = elem + "_" + sr[p];
    boost::replace_all(output3, ".", "_" );
    boost::replace_all(output3, "-", "_" );
    boost::replace_all(output3, "+", "_" );
    boost::replace_all(output3, "lox", "low" );
    s0[output3] = si_norm_ev[p] ;
    b0[output3] = bkg[p] ;
    ds[output3] = si_er_stat[p] ;
    ds_sys[output3] = si_er_sys[p] ;
    db[output3] = bkg_err[p] ;
    n0[output3] = obs[p] ;
    clsexp[output3] = cls_exp[p] ;
    s95exp[output3] = s95_exp[p] ;
    rvalue[output3] = r_value[p] ;
    clsexp_vec.push_back(std::make_pair(output3, cls_exp[p]));
    rvalue_vec.push_back(std::make_pair(output3, r_value[p]));
    p+=1 ;
    }
  }
}

void Fittino::ClsCalculator_2::SetupMeasuredValues() {
int nr = 0;


nr = std::stoi(_nr_SR) ;

std::sort(clsexp_vec.begin(), clsexp_vec.end(), [](const std::pair<std::string, double> &left, const std::pair<std::string, double> &right) {
   return left.second < right.second;
});




int k = 0;
for (int i = 0; i< clsexp_vec.size() ; i++)
{
   std::cout << clsexp_vec[i].first << "::" << clsexp_vec[i].second << std::endl;
   if (k < nr){
   
   s0_sel.push_back(std::make_pair(clsexp_vec[i].first, s0[ clsexp_vec[i].first]));
   b0_sel.push_back(std::make_pair(clsexp_vec[i].first, b0[ clsexp_vec[i].first]));
   ds_sel.push_back(std::make_pair(clsexp_vec[i].first, ds[ clsexp_vec[i].first]));
   db_sel.push_back(std::make_pair(clsexp_vec[i].first, db[ clsexp_vec[i].first]));
   n0_sel.push_back(std::make_pair(clsexp_vec[i].first, n0[ clsexp_vec[i].first]));
   }
   k +=1;
}
 std::cout << "SIGNAL" << std::endl; 
for (int i = 0; i< s0_sel.size() ; i++)
{
  std::cout << s0_sel[i].first << "::" << s0_sel[i].second << std::endl;
  std::cout << b0_sel[i].first << "::" << b0_sel[i].second << std::endl;
  std::cout << ds_sel[i].first << "::" << ds_sel[i].second << std::endl;
  std::cout << db_sel[i].first << "::" << db_sel[i].second << std::endl;
  std::cout << n0_sel[i].first << "::" << n0_sel[i].second << std::endl;
}

std::vector<double> s0_sel_num;
std::vector<double> b0_sel_num;
std::vector<double> ds_sel_num;
std::vector<double> db_sel_num;
std::vector<double> n0_sel_num;

for (int i = 0; i< s0_sel.size() ; i++)
{  
  s0_sel_num.push_back(s0_sel[i].second);
  b0_sel_num.push_back(b0_sel[i].second);
  ds_sel_num.push_back(ds_sel[i].second);
  db_sel_num.push_back(db_sel[i].second);
  n0_sel_num.push_back(n0_sel[i].second);
}

std::string s0_string;
std::string b0_string;
std::string ds_string;
std::string db_string;
std::string n0_string;

using boost::algorithm::join;
using boost::adaptors::transformed;
auto tostr = static_cast<std::string(*)(double)>(std::to_string);

std::cout << join(s0_sel_num | transformed(tostr), ",") << std::endl;
s0_string = "[" + join(s0_sel_num | transformed(tostr), ",") + "]";
b0_string = "[" + join(b0_sel_num | transformed(tostr), ",") + "]";
ds_string = "[" + join(ds_sel_num | transformed(tostr), ",") + "]";
db_string = "[" + join(db_sel_num | transformed(tostr), ",") + "]";
n0_string = "[" + join(n0_sel_num | transformed(tostr), ",") + "]";
std::cout << s0_string << std::endl;
std::cout << b0_string << std::endl;
std::cout << ds_string << std::endl;
std::cout << db_string << std::endl;
std::cout << n0_string << std::endl;




std:: string execute;
std::string cls_script_ohnepy;
cls_script_ohnepy = _cls_script.substr(0, _cls_script.size()-3);
execute = "python -c \"import " + cls_script_ohnepy + " as bla; print bla.calc_CLs(" + n0_string + ", 0.0, " + b0_string + "," + db_string + ", 0.0, " + s0_string + "," + ds_string + ", " + _nr_toys +"," + _nr_SR +", 1, 1 ,1) \" ";
const char *python = execute.c_str();
std::cout << execute ;


for (int i=0; i<s0_sel_num.size() ; i++){
if(s0_sel_num[i]< 0.01){
zero_signal =1;
}
}

if(zero_signal ==0){
 std:: system( python );
}
//std:: system( python );

if(zero_signal ==1){
 
  SetOutput( "Cls", 100  );
  SetOutput( "Chi_2", 100  );
}



}

void Fittino::ClsCalculator_2::NewMethod() {

std::vector <std::string> LISTE;
// LISTE

LISTE ={"atlas_1602_09058-SR0b3j", "atlas_1602_09058-SR0b5j",
       "atlas_1602_09058-SR1b", "atlas_1602_09058-SR3b",
       "atlas_1604_01306-SR1", "atlas_1604_07773-EM1",
       "atlas_1604_07773-EM2", "atlas_1604_07773-EM3",
       "atlas_1604_07773-EM4", "atlas_1604_07773-EM5",
       "atlas_1604_07773-EM6", "atlas_1604_07773-IM1",
       "atlas_1604_07773-IM2", "atlas_1604_07773-IM3",
       "atlas_1604_07773-IM4", "atlas_1604_07773-IM5",
       "atlas_1604_07773-IM6", "atlas_1604_07773-IM7",
       "atlas_1605_03814-2jl", "atlas_1605_03814-2jm",
       "atlas_1605_03814-2jt", "atlas_1605_03814-4jt",
       "atlas_1605_03814-5j", "atlas_1605_03814-6jm",
       "atlas_1605_03814-6jt", "atlas_1605_04285-2jSoftLep",
       "atlas_1605_04285-4jHighX", "atlas_1605_04285-4jLowX",
       "atlas_1605_04285-5j", "atlas_1605_04285-5jSoftLep",
       "atlas_1605_04285-6j", "atlas_1605_09318-SR_Gbb_A",
       "atlas_1605_09318-SR_Gbb_B", "atlas_1605_09318-SR_Gbb_C",
       "atlas_1605_09318-SR_Gtt_0L_A", "atlas_1605_09318-SR_Gtt_0L_B",
       "atlas_1605_09318-SR_Gtt_0L_C", "atlas_1605_09318-SR_Gtt_1L_A",
       "atlas_1605_09318-SR_Gtt_1L_B", "atlas_1606_03903-SR1",
       "atlas_1606_03903-SR2", "atlas_1606_03903-SR3",
       "atlas_conf_2015_082-SRZ", "atlas_conf_2016_013-0J_6j_2b",
       "atlas_conf_2016_013-0J_6j_3b", "atlas_conf_2016_013-0J_6j_m4b",
       "atlas_conf_2016_013-1J_6j_3b_HM",
       "atlas_conf_2016_013-1J_6j_3b_LM",
       "atlas_conf_2016_013-1J_6j_m4b_HM",
       "atlas_conf_2016_013-1J_6j_m4b_LM", "atlas_conf_2016_013-m2J_6j_2b",
       "atlas_conf_2016_013-m2J_6j_3b", "atlas_conf_2016_013-m2J_6j_m4b",
       "atlas_conf_2016_050-SR1", "atlas_conf_2016_050-tN_high",
       "atlas_conf_2016_050-bC2x_diag", "atlas_conf_2016_050-bC2x_med",
       "atlas_conf_2016_050-bCbv", "atlas_conf_2016_054-GG2J",
       "atlas_conf_2016_054-GG6J_bulk",
       "atlas_conf_2016_054-GG6J_high_mass",
       "atlas_conf_2016_054-GG4J_low_x",
       "atlas_conf_2016_054-GG4J_low_x_b_veto",
       "atlas_conf_2016_054-GG4J_high_x", "atlas_conf_2016_054-SS4J",
       "atlas_conf_2016_054-SS5J", "atlas_conf_2016_054-SS4J_low_x",
       "atlas_conf_2016_054-SS5J_high_x", "atlas_conf_2016_076-SR_low",
       "atlas_conf_2016_076-SR_high", "atlas_conf_2016_076-SR_w_sf",
       "atlas_conf_2016_076-SR_t_sf", "atlas_conf_2016_076-SR_w_df",
       "atlas_conf_2016_076-SR_t_df", "atlas_conf_2016_078-Meff_2j_1200",
       "atlas_conf_2016_078-Meff_2j_1600",
       "atlas_conf_2016_078-Meff_2j_2000",
       "atlas_conf_2016_078-Meff_2j_800",
       "atlas_conf_2016_078-Meff_3j_1200",
       "atlas_conf_2016_078-Meff_4j_1000",
       "atlas_conf_2016_078-Meff_4j_1400",
       "atlas_conf_2016_078-Meff_4j_1800",
       "atlas_conf_2016_078-Meff_4j_2200",
       "atlas_conf_2016_078-Meff_4j_2600",
       "atlas_conf_2016_078-Meff_5j_1400",
       "atlas_conf_2016_078-Meff_6j_1800",
       "atlas_conf_2016_078-Meff_6j_2200", "atlas_conf_2016_096-2LADF",
       "atlas_conf_2016_096-2LASF", "atlas_conf_2016_096-2LBDF",
       "atlas_conf_2016_096-2LBSF", "atlas_conf_2016_096-2LCDF",
       "atlas_conf_2016_096-2LCSF", "atlas_conf_2016_096-3LH",
       "atlas_conf_2016_096-3LI", "cms_pas_sus_15_011-SRA_b0_a.lowMET",
       "cms_pas_sus_15_011-SRA_b0_b.mediumMET",
       "cms_pas_sus_15_011-SRA_b0_c.highMET",
       "cms_pas_sus_15_011-SRA_b0_d.vhighMET",
       "cms_pas_sus_15_011-SRA_b1_a.lowMET",
       "cms_pas_sus_15_011-SRA_b1_b.mediumMET",
       "cms_pas_sus_15_011-SRA_b1_c.highMET",
       "cms_pas_sus_15_011-SRA_b1_d.vhighMET",
       "cms_pas_sus_15_011-SRB_b0_a.lowMET",
       "cms_pas_sus_15_011-SRB_b0_b.mediumMET",
       "cms_pas_sus_15_011-SRB_b0_c.highMET",
       "cms_pas_sus_15_011-SRB_b0_d.vhighMET",
       "cms_pas_sus_15_011-SRB_b1_a.lowMET",
       "cms_pas_sus_15_011-SRB_b1_b.mediumMET",
       "cms_pas_sus_15_011-SRB_b1_c.highMET",
       "cms_pas_sus_15_011-SRB_b1_d.vhighMET",
       "cms_pas_sus_15_011-SRC_ATLAS",
       "cms_pas_sus_15_011-SRD_cen_a.vlowMll_b0",
       "cms_pas_sus_15_011-SRD_cen_a.vlowMll_b1",
       "cms_pas_sus_15_011-SRD_cen_a.vlowMll_b_all",
       "cms_pas_sus_15_011-SRD_cen_b.lowMll_b0",
       "cms_pas_sus_15_011-SRD_cen_b.lowMll_b1",
       "cms_pas_sus_15_011-SRD_cen_b.lowMll_b_all",
       "cms_pas_sus_15_011-SRD_cen_c.mediumMll_b0",
       "cms_pas_sus_15_011-SRD_cen_c.mediumMll_b1",
       "cms_pas_sus_15_011-SRD_cen_c.mediumMll_b_all",
       "cms_pas_sus_15_011-SRD_cen_d.highMll_b0",
       "cms_pas_sus_15_011-SRD_cen_d.highMll_b1",
       "cms_pas_sus_15_011-SRD_cen_d.highMll_b_all",
       "cms_pas_sus_15_011-SRD_cen_e.vhighMll_b0",
       "cms_pas_sus_15_011-SRD_cen_e.vhighMll_b1",
       "cms_pas_sus_15_011-SRD_cen_e.vhighMll_b_all",
       "cms_pas_sus_15_011-SRE_for_a.vlowMll_b0",
       "cms_pas_sus_15_011-SRE_for_a.vlowMll_b1",
       "cms_pas_sus_15_011-SRE_for_a.vlowMll_b_all",
       "cms_pas_sus_15_011-SRE_for_b.lowMll_b0",
       "cms_pas_sus_15_011-SRE_for_b.lowMll_b1",
       "cms_pas_sus_15_011-SRE_for_b.lowMll_b_all",
       "cms_pas_sus_15_011-SRE_for_c.mediumMll_b0",
       "cms_pas_sus_15_011-SRE_for_c.mediumMll_b1",
       "cms_pas_sus_15_011-SRE_for_c.mediumMll_b_all",
       "cms_pas_sus_15_011-SRE_for_d.highMll_b0",
       "cms_pas_sus_15_011-SRE_for_d.highMll_b1",
       "cms_pas_sus_15_011-SRE_for_d.highMll_b_all",
       "cms_pas_sus_15_011-SRE_for_e.vhighMll_b0",
       "cms_pas_sus_15_011-SRE_for_e.vhighMll_b1",
       "cms_pas_sus_15_011-SRE_for_e.vhighMll_b_all"};


//ORG

std::ifstream infile(_ORG);
int zeile = 0;
int ORG_num[142][142];
double blub;
std::string bla;
while(std::getline(infile, bla)){
  if(zeile>=0){
    if(bla.empty())  break ;
     std::istringstream ss(bla);
     for(int i =0; i<142; i++){
        ss>>blub;
        ORG_num[zeile][i] = (int)blub;
        }
    }
   zeile+=1;
}

infile.close();


double Cls_sel[142];
std::fill_n(Cls_sel, 142, 10000000);
int k = 0;

// oders the SR respecting the OG and cls_exp, new method
int nr = 0;
nr = std::stoi(_nr_SR) ;

std::sort(rvalue_vec.begin(), rvalue_vec.end(), [](const std::pair<std::string, double> &left, const std::pair<std::string, double> &right) {
   return left.second > right.second;
});

for (int i = 0; i< 142 ; i++){
   boost::replace_all(LISTE[i], ".", "_" );
   boost::replace_all(LISTE[i], "-", "_" );
}
int number = 0;
int t = 0;
int index_1[142];
std::fill_n(index_1, 142, -1.0);
std::ptrdiff_t pos_0 = 0 ;    
for (int i = 0; i< 100 ; i++){
   std::ptrdiff_t pos = std::distance(LISTE.begin(), find(LISTE.begin(), LISTE.end(), rvalue_vec[i].first));    
   bool test = false;

   if(i < 142){
     if(i==0){
       s0_sel.push_back(std::make_pair(rvalue_vec[i].first, s0[rvalue_vec[i].first]));
       b0_sel.push_back(std::make_pair(rvalue_vec[i].first, b0[rvalue_vec[i].first]));
       ds_sel.push_back(std::make_pair(rvalue_vec[i].first, ds[rvalue_vec[i].first]));
       ds_sys_sel.push_back(std::make_pair(rvalue_vec[i].first, ds_sys[rvalue_vec[i].first]));
       db_sel.push_back(std::make_pair(rvalue_vec[i].first, db[rvalue_vec[i].first]));
       n0_sel.push_back(std::make_pair(rvalue_vec[i].first, n0[rvalue_vec[i].first]));
       rvalue_sel.push_back(std::make_pair(rvalue_vec[i].first, rvalue[rvalue_vec[i].first]));
       pos_0 = std::distance(LISTE.begin(), std::find(LISTE.begin(), LISTE.end(), rvalue_vec[i].first));    
     }
     if(i>0){
       number =  std::count(index_1, index_1+142, 1);
       if(ORG_num[pos_0][pos]==1){ 
          test = true;
          for(int l =0; l<142; l++){
   	         if(number ==0) test = true;
                 else if (number - l >= 1 && index_1[l]==1 && ORG_num[index_1[l]][pos]==0) test = false;	
           }
        }
      }         
     if(test){
       s0_sel.push_back(std::make_pair(rvalue_vec[i].first, s0[rvalue_vec[i].first]));
       b0_sel.push_back(std::make_pair(rvalue_vec[i].first, b0[rvalue_vec[i].first]));
       ds_sel.push_back(std::make_pair(rvalue_vec[i].first, ds[rvalue_vec[i].first]));
       ds_sys_sel.push_back(std::make_pair(rvalue_vec[i].first, ds_sys[rvalue_vec[i].first]));
       db_sel.push_back(std::make_pair(rvalue_vec[i].first, db[rvalue_vec[i].first]));
       n0_sel.push_back(std::make_pair(rvalue_vec[i].first, n0[rvalue_vec[i].first]));
       rvalue_sel.push_back(std::make_pair(rvalue_vec[i].first, rvalue[rvalue_vec[i].first]));
       index_1[pos] = 1;
     } 
    } 
    k +=1;
}
 
std::vector<std::pair<std::string, double>> s0_sel_sel;
std::vector<std::pair<std::string, double>> ds_sel_sel;
std::vector<std::pair<std::string, double>> ds_sys_sel_sel;
std::vector<std::pair<std::string, double>> b0_sel_sel;
std::vector<std::pair<std::string, double>> db_sel_sel;
std::vector<std::pair<std::string, double>> n0_sel_sel;
std::vector<std::pair<std::string, double>> s95exp_sel_sel;
std::vector<std::pair<std::string, double>> rvalue_sel_sel;


std::sort(rvalue_sel.begin(), rvalue_sel.end(), [](const std::pair<std::string, double> &left, const std::pair<std::string, double> &right) {
   return left.second > right.second;
});
for (int i =0; i < rvalue_sel.size(); i++){
  if (t < nr){
    s0_sel_sel.push_back(std::make_pair(rvalue_sel[i].first, s0[ rvalue_sel[i].first]));
    b0_sel_sel.push_back(std::make_pair(rvalue_sel[i].first, b0[ rvalue_sel[i].first]));
    ds_sel_sel.push_back(std::make_pair(rvalue_sel[i].first, ds[ rvalue_sel[i].first]));
    ds_sys_sel_sel.push_back(std::make_pair(rvalue_sel[i].first, ds_sys[ rvalue_sel[i].first]));
    db_sel_sel.push_back(std::make_pair(rvalue_sel[i].first, db[ rvalue_sel[i].first]));
    n0_sel_sel.push_back(std::make_pair(rvalue_sel[i].first, n0[ rvalue_sel[i].first]));
    
    t +=1;
  }
}


std::vector<double> s0_sel_num;
std::vector<double> b0_sel_num;
std::vector<double> ds_sel_num;
std::vector<double> ds_sys_sel_num;
std::vector<double> db_sel_num;
std::vector<double> n0_sel_num;

for (int i = 0; i< s0_sel_sel.size() ; i++){  
  s0_sel_num.push_back(s0_sel_sel[i].second);
  b0_sel_num.push_back(b0_sel_sel[i].second);
  ds_sel_num.push_back(ds_sel_sel[i].second);
  ds_sys_sel_num.push_back(ds_sys_sel_sel[i].second);
  db_sel_num.push_back(db_sel_sel[i].second);
  n0_sel_num.push_back(n0_sel_sel[i].second);
  
}

std::string s0_string;
std::string b0_string;
std::string ds_string;
std::string ds_sys_string;
std::string db_string;
std::string n0_string;

using boost::algorithm::join;
using boost::adaptors::transformed;
auto tostr = static_cast<std::string(*)(double)>(std::to_string);

s0_string = "[" + join(s0_sel_num | transformed(tostr), ",") + "]";
b0_string = "[" + join(b0_sel_num | transformed(tostr), ",") + "]";
ds_string = "[" + join(ds_sel_num | transformed(tostr), ",") + "]";
ds_sys_string = "[" + join(ds_sys_sel_num | transformed(tostr), ",") + "]";
db_string = "[" + join(db_sel_num | transformed(tostr), ",") + "]";
n0_string = "[" + join(n0_sel_num | transformed(tostr), ",") + "]";


std::cout << s0_string << std::endl;
std::cout << b0_string << std::endl;
std::cout << ds_string << std::endl;
std::cout << ds_sys_string << std::endl;
std::cout << db_string << std::endl;
std::cout << n0_string << std::endl;


std:: string execute;
std::string cls_script_ohnepy;
cls_script_ohnepy = _cls_script.substr(0, _cls_script.size()-3);
execute = "python -c \"import " + cls_script_ohnepy + " as bla; print bla.main(" + n0_string+ "," + b0_string+ "," + db_string+ "," + s0_string+ "," + ds_string+ "," +  ds_sys_string+ "," + ") \" ";
const char *python = execute.c_str();


if(zero_signal ==0){
 std:: system( python );
}
 

}
void Fittino::ClsCalculator_2::SetupMeasuredValues_13TeV() {

std::vector< std::vector  <std::string> > ORG;

// Orthogonal group

ORG ={{"atlas_1602_09058-SR0b3j", "atlas_1602_09058-SR0b5j"},
{"atlas_1602_09058-SR1b", "atlas_1602_09058-SR3b"},
{"atlas_1605_04285-2jSoftLep", "atlas_1605_04285-5jSoftLep",
 "atlas_conf_2016_054-GG2J", "atlas_conf_2016_054-GG6J_bulk",
 "atlas_1605_04285-4jHighX", "atlas_1605_04285-4jLowX", "atlas_1605_04285-5j",
 "atlas_1605_04285-6j", "atlas_conf_2016_054-GG6J_high_mass",
 "atlas_conf_2016_054-GG4J_low_x", "atlas_conf_2016_054-GG4J_low_x_b_veto",
 "atlas_conf_2016_050-bCbv", "atlas_conf_2016_054-SS4J",
 "atlas_conf_2016_054-SS5J", "atlas_conf_2016_054-SS5J_high_x",
 "atlas_conf_2016_054-SS4J_low_x", "atlas_conf_2016_054-GG4J_high_x"},
{"atlas_1605_09318-SR_Gbb_A", "atlas_1605_09318-SR_Gbb_B",
 "atlas_1605_09318-SR_Gbb_C", "atlas_1605_09318-SR_Gtt_0L_A",
 "atlas_1605_09318-SR_Gtt_0L_B", "atlas_1605_09318-SR_Gtt_0L_C"},
{"atlas_1605_09318-SR_Gtt_1L_A", "atlas_1605_09318-SR_Gtt_1L_B",
 "atlas_1606_03903-SR1", "atlas_1606_03903-SR2", "atlas_1606_03903-SR3",
 "atlas_conf_2016_013-0J_6j_2b", "atlas_conf_2016_050-SR1",
 "atlas_conf_2016_013-0J_6j_3b", "atlas_conf_2016_050-tN_high",
 "atlas_conf_2016_013-0J_6j_m4b", "atlas_conf_2016_050-bC2x_diag",
 "atlas_conf_2016_013-1J_6j_3b_HM", "atlas_conf_2016_050-bC2x_med",
 "atlas_conf_2016_013-1J_6j_3b_LM", "atlas_conf_2016_013-1J_6j_m4b_HM",
 "atlas_conf_2016_013-1J_6j_m4b_LM", "atlas_conf_2016_013-m2J_6j_2b",
 "atlas_conf_2016_013-m2J_6j_3b", "atlas_conf_2016_013-m2J_6j_m4b"},
{"atlas_conf_2015_082-SRZ", "atlas_conf_2016_096-3LH"},
{"atlas_conf_2016_076-SR_low", "atlas_conf_2016_076-SR_high",
 "atlas_conf_2016_096-2LADF", "atlas_conf_2016_076-SR_w_df",
 "atlas_conf_2016_096-2LBDF", "atlas_conf_2016_076-SR_t_df",
 "atlas_conf_2016_096-2LCDF", "atlas_conf_2016_096-2LASF",
 "atlas_conf_2016_076-SR_w_sf", "atlas_conf_2016_096-2LBSF",
 "atlas_conf_2016_076-SR_t_sf", "atlas_conf_2016_096-2LCSF"},
{"atlas_conf_2016_078-Meff_5j_1400"},
{"atlas_conf_2016_078-Meff_6j_1800", "atlas_conf_2016_078-Meff_6j_2200"},
{"atlas_conf_2016_096-3LI"},
{"cms_pas_sus_15_011-SRA_b0_a.lowMET"},
{"cms_pas_sus_15_011-SRA_b0_b.mediumMET"},
{"cms_pas_sus_15_011-SRA_b0_c.highMET"},
{"cms_pas_sus_15_011-SRA_b0_d.vhighMET"},
{"cms_pas_sus_15_011-SRA_b1_a.lowMET"},
{"cms_pas_sus_15_011-SRA_b1_b.mediumMET"},
{"cms_pas_sus_15_011-SRA_b1_c.highMET"},
{"cms_pas_sus_15_011-SRA_b1_d.vhighMET"},
{"cms_pas_sus_15_011-SRB_b0_a.lowMET"},
{"cms_pas_sus_15_011-SRB_b0_b.mediumMET"},
{"cms_pas_sus_15_011-SRB_b0_c.highMET"},
{"cms_pas_sus_15_011-SRB_b0_d.vhighMET"},
{"cms_pas_sus_15_011-SRB_b1_a.lowMET"},
{"cms_pas_sus_15_011-SRB_b1_b.mediumMET"},
{"cms_pas_sus_15_011-SRB_b1_c.highMET"},
{"cms_pas_sus_15_011-SRB_b1_d.vhighMET"},
{"cms_pas_sus_15_011-SRC_ATLAS"},
{"cms_pas_sus_15_011-SRD_cen_a.vlowMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_a.vlowMll_b1"},
{"cms_pas_sus_15_011-SRD_cen_b.lowMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_b.lowMll_b1"},
{"cms_pas_sus_15_011-SRD_cen_c.mediumMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_c.mediumMll_b1"},
{"cms_pas_sus_15_011-SRD_cen_d.highMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_d.highMll_b1"},
{"cms_pas_sus_15_011-SRD_cen_e.vhighMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_e.vhighMll_b1"},
{"cms_pas_sus_15_011-SRE_for_a.vlowMll_b0"},
{"cms_pas_sus_15_011-SRE_for_a.vlowMll_b1"},
{"cms_pas_sus_15_011-SRE_for_b.lowMll_b0"},
{"cms_pas_sus_15_011-SRE_for_b.lowMll_b1"},
{"cms_pas_sus_15_011-SRE_for_c.mediumMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_c.mediumMll_b1"},
{"cms_pas_sus_15_011-SRD_cen_d.highMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_d.highMll_b1"},
{"cms_pas_sus_15_011-SRD_cen_e.vhighMll_b0"},
{"cms_pas_sus_15_011-SRD_cen_e.vhighMll_b1"},
{"cms_pas_sus_15_011-SRE_for_a.vlowMll_b0"},
{"cms_pas_sus_15_011-SRE_for_a.vlowMll_b1"},
{"cms_pas_sus_15_011-SRE_for_b.lowMll_b0"},
{"cms_pas_sus_15_011-SRE_for_b.lowMll_b1"},
{"cms_pas_sus_15_011-SRE_for_c.mediumMll_b0"},
{"cms_pas_sus_15_011-SRE_for_c.mediumMll_b1"},
{"cms_pas_sus_15_011-SRE_for_d.highMll_b0"},
{"cms_pas_sus_15_011-SRE_for_d.highMll_b1"},
{"cms_pas_sus_15_011-SRE_for_e.vhighMll_b0"},
{"cms_pas_sus_15_011-SRE_for_e.vhighMll_b1"}};

double Cls_sel[ORG.size()];
std::fill_n(Cls_sel, ORG.size(), 10000000);
int k = 0;

// oders the SR respecting the OG and cls_exp

for ( std::vector<std::vector<std::string>>::size_type i = 0; i < ORG.size(); i++ )
{
 k=0;
   for ( std::vector<std::string>::size_type j = 0; j < ORG[i].size(); j++ )
   {
      boost::replace_all(ORG[i][j], ".", "_" );
      boost::replace_all(ORG[i][j], "-", "_" );
      std::cout << clsexp[ORG[i][j]] << std::endl;
      if (clsexp[ORG[i][j]] <= Cls_sel[i]){
           if (j >=1 && clsexp_sel.size()!=0 && k!=0){
             try{
                clsexp_sel.pop_back();
                }
             catch(...){}
              }
             std::cout << "Test";
             Cls_sel[i] = clsexp[ORG[i][j]];
             clsexp_sel.push_back(std::make_pair(ORG[i][j], Cls_sel[i]));
           k+=1; 
     }
   }
}



/*
for (int i = 0; i< clsexp_sel.size() ; i++)
{
   std::cout << clsexp_sel[i].first << "::" << clsexp_sel[i].second << std::endl;
   
}
*/

int nr = 0;

// gets SR, BKG etc. for the selected number of SR

nr = std::stoi(_nr_SR) ;

std::sort(clsexp_sel.begin(), clsexp_sel.end(), [](const std::pair<std::string, double> &left, const std::pair<std::string, double> &right) {
   return left.second < right.second;
});




int t = 0;
for (int i = 0; i< clsexp_sel.size() ; i++)
{
   std::cout << clsexp_sel[i].first << "::" << clsexp_sel[i].second << std::endl;
   if (t < nr){
   
   s0_sel.push_back(std::make_pair(clsexp_sel[i].first, s0[ clsexp_sel[i].first]));
   b0_sel.push_back(std::make_pair(clsexp_sel[i].first, b0[ clsexp_sel[i].first]));
   ds_sel.push_back(std::make_pair(clsexp_sel[i].first, ds[ clsexp_sel[i].first]));
   db_sel.push_back(std::make_pair(clsexp_sel[i].first, db[ clsexp_sel[i].first]));
   n0_sel.push_back(std::make_pair(clsexp_sel[i].first, n0[ clsexp_sel[i].first]));
   }
   t +=1;
}
 std::cout << "SIGNAL" << std::endl; 
for (int i = 0; i< s0_sel.size() ; i++)
{
  std::cout << s0_sel[i].first << "::" << s0_sel[i].second << std::endl;
  std::cout << b0_sel[i].first << "::" << b0_sel[i].second << std::endl;
  std::cout << ds_sel[i].first << "::" << ds_sel[i].second << std::endl;
  std::cout << db_sel[i].first << "::" << db_sel[i].second << std::endl;
  std::cout << n0_sel[i].first << "::" << n0_sel[i].second << std::endl;
}

std::vector<double> s0_sel_num;
std::vector<double> b0_sel_num;
std::vector<double> ds_sel_num;
std::vector<double> db_sel_num;
std::vector<double> n0_sel_num;

for (int i = 0; i< s0_sel.size() ; i++)
{  
  s0_sel_num.push_back(s0_sel[i].second);
  b0_sel_num.push_back(b0_sel[i].second);
  ds_sel_num.push_back(ds_sel[i].second);
  db_sel_num.push_back(db_sel[i].second);
  n0_sel_num.push_back(n0_sel[i].second);
  
}

std::string s0_string;
std::string b0_string;
std::string ds_string;
std::string db_string;
std::string n0_string;

using boost::algorithm::join;
using boost::adaptors::transformed;
auto tostr = static_cast<std::string(*)(double)>(std::to_string);

std::cout << join(s0_sel_num | transformed(tostr), ",") << std::endl;
s0_string = "[" + join(s0_sel_num | transformed(tostr), ",") + "]";
b0_string = "[" + join(b0_sel_num | transformed(tostr), ",") + "]";
ds_string = "[" + join(ds_sel_num | transformed(tostr), ",") + "]";
db_string = "[" + join(db_sel_num | transformed(tostr), ",") + "]";
n0_string = "[" + join(n0_sel_num | transformed(tostr), ",") + "]";
/*
using boost::algorithm::join;
using boost::adaptors::transformed;
auto tostr = static_cast<std::string(*)(double)>(std::to_string);

std::cout << join(s0_sel_num | transformed(tostr), ",") << std::endl;
s0_string = "[" + join(s0_sel_num | transformed(tostr), ",") + "]";
b0_string = "[" + join(b0_sel_num | transformed(tostr), ",") + "]";
ds_string = "[" + join(ds_sel_num | transformed(tostr), ",") + "]";
db_string = "[" + join(db_sel_num | transformed(tostr), ",") + "]";
n0_string = "[" + join(n0_sel_num | transformed(tostr), ",") + "]";
*/
std::cout << s0_string << std::endl;
std::cout << b0_string << std::endl;
std::cout << ds_string << std::endl;
std::cout << db_string << std::endl;
std::cout << n0_string << std::endl;



std:: string execute;
std::string cls_script_ohnepy;
cls_script_ohnepy = _cls_script.substr(0, _cls_script.size()-3);
execute = "python -c \"import " + cls_script_ohnepy + " as bla; print bla.calc_CLs(" + n0_string + ", 0.0, " + b0_string + "," + db_string + ", 0.0, " + s0_string + "," + ds_string + ", " + _nr_toys +"," + _nr_SR +", 1, 1 ,1) \" ";
const char *python = execute.c_str();
std::cout << execute ;



for (int i=0; i<s0_sel_num.size() ; i++){
if(s0_sel_num[i]< 0.01){
zero_signal =1;
}
}

if(zero_signal ==0){
 std:: system( python );
}
if(zero_signal ==1){
 
  SetOutput( "Cls", 100  );
  SetOutput( "Chi_2", 100  );
}
 
}

void Fittino::ClsCalculator_2::ReadPythonResult() {
  
 if (!boost::filesystem::exists (_cls_python_result)){
        throw CalculatorException(_name, "Cls_value_Error");}
 else{
        double cls = 0;
        double chi_2= 0;
        double mu_hat = 0;
        std::string line;
        std::ifstream infile(_cls_python_result);
        while(std::getline(infile, line)){
         if(line.empty())  break ;


         std::istringstream ss(line);
            ss>>chi_2 >> cls >> mu_hat;
         }
        infile.close();
    if(_nr_toys=="None"){   
    cls = -500.0; 
    SetOutput( "Cls", cls  );
    SetOutput( "Chi_2", chi_2  );
    SetOutput( "Mu_hat", mu_hat  );
    }
    else {    
    SetOutput( "Cls", cls  );
    SetOutput( "Chi_2", chi_2  );
    SetOutput( "Mu_hat", mu_hat  );
    }
   }
}

//void Fittino::ClsCalculator_2::Initialize() {
  
//}

