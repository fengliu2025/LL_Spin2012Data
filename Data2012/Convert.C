#include "ntp_Lambda.C"


void Convert(){
	ntp_Lambda myClass;
	for(int i=1;i<7;i++){
                std::cout<<i<<std::endl;
		myClass.InPutFileList.push_back( Form("/star/u/vanekjan/pwg/vanekjan/CodeQA/L_spin_spin_correlations/codes/Analysis/input/data/output_%d.root",i)  );
                myClass.OutPutFileList.push_back(Form("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012data/BigTree/Output_%d.root",i)  );
     }

	myClass.Loop();

}
