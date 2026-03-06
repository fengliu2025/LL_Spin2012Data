#include <iostream>
#include <fstream>
#include <string>
#include "ntp_Lambda.C"

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

	//Create the ntp_Lambda class 
	//ntp_Lambda t = ntp_Lambda(InputFiles,OutputFile);
	//t.Loop();
	std::cout<<InputFiles.size()<<std::endl;

}
