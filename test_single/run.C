#include "ntp_Lambda.C"


void run(){

	ntp_Lambda t;
	for(int i=0;i<6;i++){
		t.InPutFileList.push_back(Form("/star/u/vanekjan/pwg/vanekjan/CodeQA/L_spin_spin_correlations/codes/Analysis/input/data/output_%d.root",i+1));
	}
	t.OutPutFile("Jan_QA_Plot.root");


	t.Loop();






}