
///////////////////////////////////////////////////////////
// MONITOR OF THE OPTIMIZATION PROCESS OF THE MARKOV CHAIN
///////////////////////////////////////////////////////////
//
// provides with plots to check the optimization reliability:
// [1] chi2 vs. variable -> where are the minima ?
// [2] var vs. steps vs. chi2 -> which region of the parameter space are covered ?
// [3] success rate vs. the chain numbers -> how is the convergence ?
// [4] proposal width vs. the chain numbers -> how is the convergence ?
//
// How to use it ?
// Two arguments to the function "monitor"
// path = "/myPath/MarkovChain.root"
// var = "P_A0"

// This macro handles so far only the individual optimization
// (more to come)
//
// 18th February 2011
// Xavier Prudent
// prudent[at]physik.tu-dresden.de


void monitor( TString path = "", TString var = "" ){

  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");

  TFile* file = TFile::Open( path );
  TTree* tree = (TTree*) file->Get( "widthOptimization" );
  TTree* tree2 = (TTree*) file->Get( "widthVectorNtuple" );

  TCanvas* c1 = new TCanvas ( "c1", "", 800, 800 );
  c1->Divide(2,2);


  c1->cd( 1 );
  tree->Draw( "chi2:" + var, "chi2 < 1000" );
  TH1* histo1 = (TH1*) gPad->GetPrimitive("htemp"); 
  histo1->SetTitle( "#chi2 vs. " + var + ", #chi2 < 10^{3}");
  histo1->SetYTitle("#chi2");
  histo1->SetXTitle( var );


  c1->cd( 2 );
  tree->Draw( var + ":steps:chi2", "chi2 < 1000", "colz" );
  TH2* histo2 = (TH2*) gPad->GetPrimitive("htemp"); 
  histo2->SetTitle( var + " vs. steps, #chi2 < 10^{3}");
  histo2->SetYTitle( var );


  c1->cd( 3 );
  Int_t nbChain = tree2->GetEntries();
  tree2->Draw( "success:numChain","","LP" );
  histo1 = (TH1*) gPad->GetPrimitive("htemp");
  histo1->GetXaxis()->SetNdivisions( nbChain );
  histo1->SetTitle( "Success rate for each chain" );
  histo1->SetYTitle("Success rate");
  histo1->SetXTitle("Chain number");


  c1->cd( 4 );
  tree2->Draw( "width:numChain","","L" );
  histo1 = (TH1*) gPad->GetPrimitive("htemp");
  histo1->GetXaxis()->SetNdivisions( nbChain );
  histo1->SetTitle( "Proposal width for each chain" );
  histo1->SetYTitle("Proposal width");
  histo1->SetXTitle("Chain number");

  cout << " >>>> Evolution of proposal widths..." << endl;
  tree2->Scan("numChain:width");
  cout << " >>>> Pick the last value as optimised" << endl;


  




  return;
}
