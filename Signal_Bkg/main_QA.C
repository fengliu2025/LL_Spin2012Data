#include <iostream>
#include <fstream>
#include <string>

#include "ntp_Lambda_Reader.h"
#include "ntp_Lambda_Histogram.h"
#include "ntp_Lambda_Analyzer.h"


void main_QA(){
	std::vector<std::string> InputFiles;
	std::string OutputFile="QA.root";
	std::ifstream filelist("Inputfilelist.txt");
	if (!filelist.is_open()) {
		std::cerr <<"can't open the file" <<std::endl;
	}
	std::string directory="/star/u/jjiastar/pwg/Spin/production/output/";
	std::string line;
	while (std::getline(filelist,line)){
		InputFiles.push_back(directory+line);

	}
	
	//InputFiles.push_back("25130036_0.root.picoLambdaAnaMaker.root");
	//create the ntp_Lambda_Reader
	ntp_Lambda_Reader *myReader = new ntp_Lambda_Reader(InputFiles);
	//create the ntp_Lambda_Histogram
	ntp_Lambda_Histogram *myHistogram = new ntp_Lambda_Histogram(myReader,OutputFile);
	myHistogram->InitHistogram();
	//create the ntp_Lambda_Analyzer
	ntp_Lambda_Analyzer *myAnalyzer = new ntp_Lambda_Analyzer(myReader,myHistogram);
	myAnalyzer->Analysis_MassSpectrum();

	myHistogram->WriteAll();
	

}
