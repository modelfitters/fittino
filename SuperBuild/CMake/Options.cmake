include( CMakeDependentOption )

option( BUILD_OFFLINE "" OFF )

cmake_dependent_option( INSTALL_Python2 "Install Python 2 environment fulfilling the requirements of CheckMATE, SCYNet, and SModelS." OFF "NOT INSTALL_Python3" OFF )
add_feature_info( INSTALL_Python2 INSTALL_Python2 "the Python 2 evironment installed fulfills the requirements of CheckMATE, SCYNet, and SModelS. This option depends on the option INSTALL_Python3 to be OFF." )

cmake_dependent_option( INSTALL_Python3 "Install Python 3 environment including flavio and its requirements" OFF "NOT INSTALL_Python2" OFF )
add_feature_info( INSTALL_Python3 INSTALL_Python3 "the Python 3 installation includes flavio and its requirements. Flavio is used by the FlavioCalculator of Fittino." )

if( INSTALL_Python2 OR INSTALL_Python3 ) 
set( installPython ON) 
else()
set (installPython OFF)
endif()
cmake_dependent_option( Conda_REQUIREMENTSFILE  "Create Conda environment as specified in requiremts file." ON "installPython" OFF )

cmake_dependent_option( Python_REQUIREMENTSFILE  "Create Python environment as specified in requiremts file." ON "installPython" OFF )

option( INSTALL_MadGraph5 "Install MadGraph5" OFF )
add_feature_info( INSTALL_MadGraph5 INSTALL_MadGraph5 "MadGraph5 is required for installation of CheckMATE." )

option( INSTALL_Prospino "Install Prospino" OFF )
add_feature_info( INSTALL_Prospino INSTALL_Prospino "Prospino is required for compiling helper program used by ProspinoCalculator." )

option( INSTALL_HepMC2 "Install HepMC2" OFF )
add_feature_info( INSTALL_HepMC2 INSTALL_HepMC2 "HepMC2 is required for installation of Pythia." )

option( INSTALL_Pythia8 "Install Pythia8" OFF )
add_feature_info( INSTALL_Pythia8 INSTALL_Pythia8 "Pythia8 is required for installation of CheckMATE." )

option( INSTALL_Delphes "Install Delphes" OFF )
add_feature_info( INSTALL_Delphes INSTALL_Delphes "Delphes is required for installation of CheckMATE." )

cmake_dependent_option( INSTALL_CheckMATE "Install CheckMATE" OFF "INSTALL_Delphes;INSTALL_Pythia8;INSTALL_HepMC2;INSTALL_MadGraph5;INSTALL_Python2" OFF )
add_feature_info( INSTALL_CheckMATE INSTALL_CheckMATE "This option depends on the options INSTALL_Delphes, INSTALL_Pythia8, INSTALL_HepMC2, INSTALL_MadGraph5 to be ON." )

option( INSTALL_SLHAea "Install SLHAea" OFF )
add_feature_info( INSTALL_SLHAea INSTALL_SLHAea "SLHAea is required for using the SPhenoCalculator of Fittino." )

option( INSTALL_SPheno "Install SPheno" OFF )
add_feature_info( INSTALL_SPheno INSTALL_SPheno "SPheno is used by the SPhenoCalculator of Fittino." )

cmake_dependent_option( INSTALL_MSSMTriLnV "Install MSSMTriLnV" OFF "INSTALL_SPheno" OFF )
add_feature_info( INSTALL_MSSMTriLnV INSTALL_MSSMTriLnV "This option depends on the option INSTALL_SPheno to be ON." )

option( USE_WCxf "Use WCxf" ON )
add_feature_info( USE_WCxf USE_WCxf "SPhenoMSSMTriLnV and flavio communicate via WCxf." )

option( USE_SPhenoMSSMTriLnV_LOOPS "Use loops in SPhenoMSSMTriLnV." ON )
add_feature_info( USE_SPhenoMSSMTriLnV_LOOPS USE_SPhenoMSSMTriLnV_LOOPS "If OFF the epsTree parameter of SPhenoMSSMTriLnV is set to 0 such that loop corrections to wilson coefficients are never calculated. If ON it stays at 0.1 such that loops are calculated if the tree level contribution is less than 0.1 (in units of GeV)." )

cmake_dependent_option( INSTALL_SModelS "Install SModelS" OFF "INSTALL_Python2" OFF )
add_feature_info( INSTALL_SModelS INSTALL_SModelS "SModelS is used by the SModelsCalculator of Fittino. This option depends on the option INSTALL_Python2 to be ON. It requires the package UnixCommands to be found in order to untar the database of SModelS." )

cmake_dependent_option( INSTALL_SCYNet "Install SCYNet" OFF "INSTALL_Python2" OFF )
add_feature_info( INSTALL_SCYNet INSTALL_SCYNet "SCYNet is used by the SCYNetCalculator of Fittino. This option depends on the option INSTALL_Python2 to be ON.")

cmake_dependent_option( INSTALL_SCYNet2 "Install SCYNet2" OFF "INSTALL_Python2" OFF )
add_feature_info( INSTALL_SCYNet2 INSTALL_SCYNet2 "SCYNet is used by the SCYNetCalculator of Fittino. This option depends on the option INSTALL_Python2 to be ON.")


option( INSTALL_Fittino "Install Fittino" OFF )
add_feature_info( INSTALL_Fittino INSTALL_Fittino "" )
