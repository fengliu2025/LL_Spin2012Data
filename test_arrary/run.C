#include "ntp_Lambda.C"


void run(){

	ntp_Lambda t;



	for(int i=0;i<6;i++){
		t.InPutFileList.push_back(Form("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012data/BigTree/Output_%d.root",i+1));
	}
	t.OutPutFile="My_QA_Plot.root";


	t.Loop();






}