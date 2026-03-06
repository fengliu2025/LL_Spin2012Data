#include <iostream>
#include <fstream>
#include <string>
#include "ntp_Lambda_Analyzer.h"
#include "ntp_Lambda_Calculator.h"
#include "ntp_Lambda_Histogram.h"
#include "ntp_Lambda_Reader.h"
//void main_Analysis(){
int main(int argc, char* argv[]){
	if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <number>\n";
        return 1;
    }
	int i_file = std::stoi(argv[1]);  
	std::vector<std::string> FullInputFiles;
	std::vector<std::string> MixEventInputFiles;
	std::vector<std::string> SameEventInputFiles;
	std::string OutputFile=Form("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012result/Density_Matrix_ME_%d.root",i_file);
	
	std::ifstream filelist("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/Density_Matrix/Inputfilelist.txt");
	//std::ifstream filelist("Inputfilelist.txt");
	if (!filelist.is_open()) {
		std::cerr <<"can't open the file" <<std::endl;
	}
	//std::string directory="/star/u/jjiastar/pwg/Spin/production/output/";
	std::string line;
	
	while (std::getline(filelist,line)){
		FullInputFiles.push_back(line);
	}
	std::cout<<FullInputFiles.size()<<std::endl;
	MixEventInputFiles.push_back(FullInputFiles[i_file]);
	SameEventInputFiles.push_back(FullInputFiles[i_file]);

	ntp_Lambda_Reader *mySameEventReader  = new ntp_Lambda_Reader(SameEventInputFiles);
	ntp_Lambda_Reader *myMixEventReader   = new ntp_Lambda_Reader(MixEventInputFiles);

	ntp_Lambda_Calculator *myCalculator   = new ntp_Lambda_Calculator();
	ntp_Lambda_Histogram *myHistogram     = new ntp_Lambda_Histogram(mySameEventReader,myCalculator,OutputFile);

	myHistogram->InitHitogram();

    ntp_Lambda_Analyzer *myAnalyzer       = new ntp_Lambda_Analyzer(mySameEventReader, myMixEventReader, myCalculator, myHistogram );

    //myAnalyzer->Analysis_SameEvent();
    myAnalyzer->Analysis_MixEvent();
	myHistogram->WriteAll();

	return 0;
}
