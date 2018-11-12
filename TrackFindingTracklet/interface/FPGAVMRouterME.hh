//This class implementes the VM router
#ifndef FPGAVMROUTERME_H
#define FPGAVMROUTERME_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAVMRouterME:public FPGAProcessBase{

public:

  FPGAVMRouterME(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    
    layer_=0;
    disk_=0;
    
    if (name_[6]=='L') layer_=name_[7]-'0';    
    if (name_[6]=='D') disk_=name_[7]-'0';    

    assert((layer_!=0)||(disk_!=0));

    if (layer_!=0) {
      nbitsfinebintable_=8;
      unsigned int nbins=1<<nbitsfinebintable_;
      
      for(unsigned int i=0;i<nbins;i++) {
	finebintable_.push_back(-1);
      }
      
      for(unsigned int i=0;i<nbins;i++) {

	
	int ibin=(i>>(nbitsfinebintable_-3));
	
	int zfine=(i>>(nbitsfinebintable_-6))-(ibin<<3);
	
	//awkward bit manipulations since the index is from a signed number...
	int index=i+(1<<(nbitsfinebintable_-1));
	
	if (index>=(1<<nbitsfinebintable_)){
	  index-=(1<<nbitsfinebintable_);
	}
	
	finebintable_[index]=zfine;
	
      }
    }

    if (disk_!=0) {

      nbitsfinebintable_=8;
      unsigned int nbins=1<<nbitsfinebintable_;
      
      for(unsigned int i=0;i<nbins;i++) {

	double rstub=0.0;
	
	if (i<10) {
	  if (disk_<=2) {
	    rstub=rDSSinner[i];
	  } else {
	    rstub=rDSSouter[i];
	  }
	} else {
	  rstub=kr*(i<<(nrbitsdisk-nbitsfinebintable_));
	}

	if (rstub<rmindiskvm) {
	  finebintable_.push_back(-1);
	} else {	
	  int bin=8.0*(rstub-rmindiskvm)/(rmaxdisk-rmindiskvm);
	  assert(bin>=0);
	  assert(bin<MEBinsDisks);	  
	  int rfine=64*((rstub-rmindiskvm)-bin*(rmaxdisk-rmindiskvm)/8.0)/(rmaxdisk-rmindiskvm);
	  finebintable_.push_back(rfine);
	}
      }
    }
    
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="allstubout"||
	output=="allstuboutn1"||
	output=="allstuboutn2"||
	output=="allstuboutn3"||
	output=="allstuboutn4"||
	output=="allstuboutn5"||
	output=="allstuboutn6"
	){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      allstubs_.push_back(tmp);
      return;
    }

    bool print=getName()=="VMRME_L3PHI2";
    print=false;
    
    if (print) cout << "FPGAVMRouterME in "<<getName()<<endl;
    
    for (int i=0;i<32;i++) {
      std::ostringstream oss;
      oss<<(i+1);
      string s=oss.str();

      if (print) {
	cout << "FPGAVMRouterME strings to match "<<output<<" "<<"vmstuboutPHI"+s+"n1"<<endl;
      }
	
      if (output=="vmstuboutPHI"+s+"n1"||
	  output=="vmstuboutPHI"+s+"n2"||
	  output=="vmstuboutPHIA"+s+"n1"||
	  output=="vmstuboutPHIB"+s+"n1"||
	  output=="vmstuboutPHIC"+s+"n1"||
	  output=="vmstuboutPHID"+s+"n1"||
	  output=="vmstuboutPHIE"+s+"n1"||
	  output=="vmstuboutPHIF"+s+"n1"||
	  output=="vmstuboutPHIG"+s+"n1"||
	  output=="vmstuboutPHIH"+s+"n1"||
	  output=="vmstuboutPHI"+s+"n3"
	  ){
	FPGAVMStubsME* tmp=dynamic_cast<FPGAVMStubsME*>(memory);
	assert(tmp!=0);
	if (print){
	  cout << "FPGAVMRouter found match "<<getName() << " " << i << endl;
	}
	vmstubsPHI_[i].push_back(tmp);
	return;
      }
    }

    cout << "Could not find : "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="stubin"){
      FPGAInputLink* tmp1=dynamic_cast<FPGAInputLink*>(memory);
      assert(tmp1!=0);
      if (tmp1!=0){
	stubinputs_.push_back(tmp1);
      }
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute(){

    //cout << "In FPGAVMRouterME "<<getName()<<" "<<stubinputs_.size()<<endl;

    
    unsigned int count=0;

    if (stubinputs_.size()!=0){
      for(unsigned int j=0;j<stubinputs_.size();j++){
	for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	  count++;
	  if (count>MAXVMROUTER) continue;
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);

	  int iphiRaw=stub.first->iphivmRaw();
	  int iphiRawPlus=stub.first->iphivmRawPlus();
	  int iphiRawMinus=stub.first->iphivmRawMinus();

	  int iphistub=iphiRaw;

	  if (hourglass) {

	    assert(0);
	    
	    int layer=stub.first->layer().value();
	    int disk=abs(stub.first->disk().value());

	    //cout << getName()<<" layer,disk : "<<layer<<" "<<disk<<endl;
	    
	    int nvm=-1;
	    if (layer!=-1) {
	      nvm=nallstubslayers[layer]*nvmtelayers[layer];
	    }
	    if (disk!=0){
	      nvm=nallstubsdisks[disk-1]*nvmtedisks[disk-1];
	    }
	    assert(nvm>0&&nvm<=32);
	    iphiRaw=iphiRaw/(32/nvm);
	    iphiRawPlus=iphiRawPlus/(32/nvm);
	    iphiRawMinus=iphiRawMinus/(32/nvm);
	    if (iphiRawPlus<0) iphiRawPlus=0;
	    if (iphiRawPlus>=nvm) iphiRawPlus=nvm-1;
	    if (iphiRawMinus<0) iphiRawMinus=0;
	    if (iphiRawMinus>=nvm) iphiRawMinus=nvm-1;
	    
	  } else {

	    assert(iphiRaw>=4 and iphiRaw<=27);
	  
	  
	    iphiRaw-=4;
	    iphiRaw=iphiRaw>>1;
	    assert(iphiRaw>=0);
	    assert(iphiRaw<=11);
	    
	    iphiRawPlus-=4;
	    iphiRawPlus=iphiRawPlus>>1;
	    if (iphiRawPlus<0) iphiRawPlus=0;
	    if (iphiRawPlus>11) iphiRawPlus=11;
	    
	    iphiRawMinus-=4;
	    iphiRawMinus=iphiRawMinus>>1;
	    if (iphiRawMinus<0) iphiRawMinus=0;
	    if (iphiRawMinus>11) iphiRawMinus=11;

	  }

	  if (disk_!=0) {

	    int index=stub.first->r().value();
	    if (stub.first->isPSmodule()){
	      index=stub.first->r().value()>>(stub.first->r().nbits()-nbitsfinebintable_);
	    }

	    int rfine=finebintable_[index];

	    stub.first->setfiner(rfine);

	  }

	  if (layer_!=0) {

	    int index=(stub.first->z().value()>>(stub.first->z().nbits()-nbitsfinebintable_))&((1<<nbitsfinebintable_)-1);
	    
	    int zfine=finebintable_[index];

	    stub.first->setfinez(zfine);

	  }
	  
	  bool insert=false;
	  
	  for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
	    if (debug1) {
	      cout << "FPGAVMRouterME "<<getName()<<" add stub ( r = "<<stub.second->r()<<" phi = "<<stub.second->phi()<<" ) in : "<<vmstubsPHI_[iphiRaw][l]->getName()<<" iphistub = " << iphistub << endl;
	    }
	    vmstubsPHI_[iphiRaw][l]->addStub(stub);
	    insert=true;
	  }

	  if (iphiRaw!=iphiRawPlus) {
	    for (unsigned int l=0;l<vmstubsPHI_[iphiRawPlus].size();l++){
	      vmstubsPHI_[iphiRawPlus][l]->addStub(stub);
	    }
	  }
	  if (iphiRaw!=iphiRawMinus) {
	    for (unsigned int l=0;l<vmstubsPHI_[iphiRawMinus].size();l++){
	      vmstubsPHI_[iphiRawMinus][l]->addStub(stub);
	    }
	  }
	 
	  
	  
	  stub.first->setAllStubIndex(allstubs_[0]->nStubs());
	  stub.second->setAllStubIndex(allstubs_[0]->nStubs());

	  for (unsigned int l=0;l<allstubs_.size();l++){
	    allstubs_[l]->addStub(stub);
	  }

	  if (!insert){
	    cout << "In "<<getName()<<" did not insert stub from input "<<stubinputs_[j]->getName()<<endl;
	  }
	  assert(insert);

	}
      }
      
      if (writeVMOccupancyME) {
	static ofstream out("vmoccupancyme.txt");

	for (int i=0;i<24;i++) {
	  if (vmstubsPHI_[i].size()!=0) {
	    out<<vmstubsPHI_[i][0]->getName()<<" "<<vmstubsPHI_[i][0]->nStubs()<<endl;
	  }
	}

      }

    }

    if (writeAllStubs) {
      static ofstream out("allstubsme.txt");
      out<<allstubs_[0]->getName()<<" "<<allstubs_[0]->nStubs()<<endl;
    }
 

  }



private:

  int layer_;
  int disk_;

  int nbitsfinebintable_;
  vector<int> finebintable_;

  
  vector<FPGAInputLink*> stubinputs_;
  vector<FPGAAllStubs*> allstubs_;

  vector<FPGAVMStubsME*> vmstubsPHI_[32];

};

#endif

