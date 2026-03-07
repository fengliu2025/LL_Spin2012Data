#include <iostream>
#include <fstream>
#include <string>
#include "ntp_Lambda_Analyzer.h"
#include "ntp_Lambda_Calculator.h"
#include "ntp_Lambda_Histogram.h"
#include "ntp_Lambda_Reader.h"
void main_Analysis_SE(){ 
	std::vector<std::string> FullInputFiles;
	std::vector<std::string> MixEventInputFiles;
	std::vector<std::string> SameEventInputFiles;
	//std::string OutputFile=Form("/gpfs/mnt/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012result/Density_Matrix_SE_%2f.root",1.0);
	
	//std::ifstream filelist("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/Density_Matrix/Inputfilelist.txt");
	std::string OutputFile=Form("Density_Matrix_SE_%2f_LittleTree.root",1.0);
        //std::ifstream filelist("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/Density_Matrix/Inputfilelist.txt");
        std::ifstream filelist("/gpfs/mnt/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012data/LittleTree/FileList.txt");
	if (!filelist.is_open()) {
		std::cerr <<"can't open the file" <<std::endl;
	}
	std::string directory="";
	std::string line;
	
	while (std::getline(filelist,line)){
		FullInputFiles.push_back(line);
	}


	ntp_Lambda_Reader *mySameEventReader  = new ntp_Lambda_Reader(FullInputFiles);
	ntp_Lambda_Reader *myMixEventReader   = new ntp_Lambda_Reader(FullInputFiles);

	ntp_Lambda_Calculator *myCalculator   = new ntp_Lambda_Calculator();
	ntp_Lambda_Histogram *myHistogram     = new ntp_Lambda_Histogram(mySameEventReader,myCalculator,OutputFile);

	myHistogram->InitHitogram();

    ntp_Lambda_Analyzer *myAnalyzer       = new ntp_Lambda_Analyzer(mySameEventReader, myMixEventReader, myCalculator, myHistogram );

    myAnalyzer->Analysis_SameEvent();
    //myAnalyzer->Analysis_MixEvent();
	myHistogram->WriteAll();

	return 0;
}
