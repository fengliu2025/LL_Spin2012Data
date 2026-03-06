#include "ntp_Lambda.C"


void Convert(){
	ntp_Lambda myClass;
	for(int i=1;i<7;i++){
                std::cout<<i<<std::endl;
		myClass.InPutFileList.push_back( Form("output_%d.root",i)  );
                myClass.OutPutFileList.push_back(Form("/gpfs01/star/scratch/fengliu/LL_Spin_2012/Output_%d.root",i)  );
     }

	myClass.Loop();

}
