//This class implementes the VM router
#ifndef FPGAVMROUTER_H
#define FPGAVMROUTER_H

#include "FPGAProcessBase.hh"
#include "FPGATETableOuter.hh"
#include "FPGATETableInner.hh"
#include "FPGATETableOuterDisk.hh"
#include "FPGATETableInnerDisk.hh"
#include "FPGATETableInnerOverlap.hh"

using namespace std;

class FPGAVMRouter:public FPGAProcessBase{

public:

  FPGAVMRouter(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){

    layer_=0;
    disk_=0;
    
    if (name_[4]=='L') layer_=name_[5]-'0';    
    if (name_[4]=='D') disk_=name_[5]-'0';    

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
	output=="allstuboutn6"||
	output=="allstuboutn7"||
	output=="allstuboutn8"||
	output=="allstuboutn9"||
	output=="allstuboutn10"||
	output=="allstuboutn11"||
	output=="allstuboutn12"||
	output=="allstuboutn13"||
	output=="allstuboutn14"||
	output=="allstuboutn15"||
	output=="allstuboutn16"||
	output=="allstuboutn17"||
	output=="allstuboutn18"
	){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      allstubs_.push_back(tmp);
      return;
    }
    for (int i=0;i<32;i++) {
      std::ostringstream oss;
      oss<<(i+1);
      string s=oss.str();
      for (int n=0;n<10;n++) {
	std::ostringstream oss1;
	oss1<<(n+1);
	string ns=oss1.str();	
	if (output=="vmstuboutPHIA"+s+"n"+ns||
	    output=="vmstuboutPHIB"+s+"n"+ns||
	    output=="vmstuboutPHIC"+s+"n"+ns||
	    output=="vmstuboutPHID"+s+"n"+ns||
	    output=="vmstuboutPHIE"+s+"n"+ns||
	    output=="vmstuboutPHIF"+s+"n"+ns||
	    output=="vmstuboutPHIG"+s+"n"+ns||
	    output=="vmstuboutPHIH"+s+"n"+ns||
	    output=="vmstuboutPHII"+s+"n"+ns||
	    output=="vmstuboutPHIJ"+s+"n"+ns||
	    output=="vmstuboutPHIK"+s+"n"+ns||
	    output=="vmstuboutPHIL"+s+"n"+ns||
	    output=="vmstuboutPHIX"+s+"n"+ns||
	    output=="vmstuboutPHIY"+s+"n"+ns||
	    output=="vmstuboutPHIZ"+s+"n"+ns||
	    output=="vmstuboutPHIW"+s+"n"+ns||
	    output=="vmstuboutPHIQ"+s+"n"+ns||
	    output=="vmstuboutPHIR"+s+"n"+ns||
	    output=="vmstuboutPHIS"+s+"n"+ns||
	    output=="vmstuboutPHIT"+s+"n"+ns
	  ){
	  //cout << "memory name : "<<memory->getName()<<" "
	  //     <<memory->getName().substr(3,2)<<" "
	  //     <<memory->getName().substr(11,1)
	  //     <<endl;

	  if (memory->getName().substr(3,2)=="TE") {
	    FPGAVMStubsTE* tmp=dynamic_cast<FPGAVMStubsTE*>(memory);
	    assert(tmp!=0);
	    if (memory->getName().substr(11,1)[0]<'I') {
	      vmstubsTEPHI_[i].push_back(tmp);
	    } else if (memory->getName().substr(11,1)[0]<'M') {
	      vmstubsTEExtraPHI_[i].push_back(tmp);
	    } else {
	      vmstubsTEOverlapPHI_[i].push_back(tmp);
	    }
	  } else if (memory->getName().substr(3,2)=="ME") {
	    FPGAVMStubsME* tmp=dynamic_cast<FPGAVMStubsME*>(memory);
	    assert(tmp!=0);
	    vmstubsMEPHI_[i].push_back(tmp);
	  } else {
	    assert(0);
	  }

	  
	  return;
	}
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

    assert(allstubs_.size()!=0);

    unsigned int count=0;
    for(unsigned int j=0;j<stubinputs_.size();j++){
      for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	if (count>MAXVMROUTER) continue;
	std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);
	
	stub.first->setAllStubIndex(count);
	stub.second->setAllStubIndex(count);

	stub.first->setAllStubAddressTE(count);

	for (unsigned int l=0;l<allstubs_.size();l++){
	  allstubs_[l]->addStub(stub);
	}
	
	count++;
      }
    }	


    if (writeAllStubs) {
      static ofstream out("allstubs.txt");
      out<<allstubs_[0]->getName()<<" "<<allstubs_[0]->nStubs()<<endl;
    }

    
    if (disk_!=5) {
      executeTE(false);
    }
    if (layer_==1||layer_==2||disk_==1) {
      executeTE(true);
    }
    executeME();
  }


  
  
  void executeTE(bool overlap){

    
    //cout << "In FPGAVMRouterTE"<<endl;
    
    assert(stubinputs_.size()!=0);

    unsigned int count=0;

 
    if (layer_!=0){  //First handle layer stubs

      for(unsigned int j=0;j<stubinputs_.size();j++){
	for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	  count++;
	  if (count>MAXVMROUTER) continue;
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);


	  int iphiRaw=stub.first->iphivmRaw();

	  bool insert=false;


	  int binlookup=-1;
	  int binlookupextra=-1;
	  if (overlap) {
	    assert(layer_==1||layer_==2);
	    binlookup=lookupInnerOverlapLayer(stub.first);
	  } else {
	    switch (layer_) {
	    case 2 : binlookup=lookupOuterLayer(stub.first);
	      binlookupextra=lookupInnerLayer(stub.first);
	      break;
	    case 4 : binlookup=lookupOuterLayer(stub.first);
	      break;
	    case 6 : binlookup=lookupOuterLayer(stub.first);
	      break;
	    case 1 : binlookup=lookupInnerLayer(stub.first);
	      break;
	    case 3 : binlookup=lookupInnerLayer(stub.first);
	      binlookupextra=lookupOuterLayer(stub.first);
	      break;
	    case 5 : binlookup=lookupInnerLayer(stub.first);
	      break;
	    default : assert(0);
	    }
	  }
	  if ((layer_==2 or layer_==3) && binlookupextra!=-1) {
	    stub.first->setVMBitsExtra(binlookupextra);
	  }

	  if (binlookup!=-1) {
	    if (overlap) {
	      stub.first->setVMBitsOverlap(binlookup);
	    } else {
	      stub.first->setVMBits(binlookup);
	    }
	  }
	  
	  unsigned int layer=stub.first->layer().value();
	  if ((layer==1 || layer== 2) && binlookupextra!=-1 ) {
	    int iphiRawTmp=iphiRaw/(32/(nallstubslayers[layer]*nvmteextralayers[layer]));
	    for (unsigned int l=0;l<vmstubsTEExtraPHI_[iphiRawTmp].size();l++){
	      if (debug1) {
		cout << getName()<<" try adding extra stub to "<<vmstubsTEExtraPHI_[iphiRawTmp][l]->getName()<<endl;
	      }
	      vmstubsTEExtraPHI_[iphiRawTmp][l]->addStub(stub);
	      insert=true;
	    }
	  }
	  
	  if (binlookup!=-1) {
	    if (overlap) {
	      int iphiRawTmp=iphiRaw/(32/(nallstubsoverlaplayers[layer]*nvmteoverlaplayers[layer]));
	      for (unsigned int l=0;l<vmstubsTEOverlapPHI_[iphiRawTmp].size();l++){
		if (debug1) {
		  cout << getName()<<" try adding overlap stub to "<<vmstubsTEOverlapPHI_[iphiRawTmp][l]->getName()<<endl;
		}
		vmstubsTEOverlapPHI_[iphiRawTmp][l]->addStub(stub);
		insert=true;
	      }
	    } else {
	      int iphiRawTmp=iphiRaw/(32/(nallstubslayers[layer]*nvmtelayers[layer]));
	      for (unsigned int l=0;l<vmstubsTEPHI_[iphiRawTmp].size();l++){
		if (debug1) {
		  cout << getName()<<" try adding stub to "<<vmstubsTEPHI_[iphiRawTmp][l]->getName()<<endl;
		}
		vmstubsTEPHI_[iphiRawTmp][l]->addStub(stub);
		insert=true;
	      }
	    }
	  }
	  
	  if (false) {
	    if (!insert) {
	      cout << getName()<<" did not insert stub"<<endl;
	    }
	    assert(insert);
	  }
	}
      }

    }

    if (disk_!=0) {
      for(unsigned int j=0;j<stubinputs_.size();j++){
	for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);

	  if (!stub.second->isPSmodule()) {
	    if (debug1) {
	      cout << getName() <<" stub at r = "<<stub.second->r()<<" is 2S module"<<endl;
	    }
	    continue;
	  }
	  
	  int iphiRaw=stub.first->iphivmRaw();

	  bool insert=false;


	  if (overlap) {
	    
	    int binlookup=lookupOuterOverlapD1(stub.first);
	    assert(binlookup>=0);
	    stub.first->setVMBitsOverlap(binlookup);

	    iphiRaw=iphiRaw/(32/(nallstubsoverlapdisks[0]*nvmteoverlapdisks[0]));

	    for (unsigned int l=0;l<vmstubsTEOverlapPHI_[iphiRaw].size();l++){
	      if (debug1) {
		cout << getName()<<" added stub to : "<<vmstubsTEOverlapPHI_[iphiRaw][l]->getName()<<endl;
	      }
	      vmstubsTEOverlapPHI_[iphiRaw][l]->addStub(stub);
	      insert=true;
	    }
	      
	  } else {

	    int binlookup=-1;

	    if (!overlap) {
	      switch (disk_) {
	      case 2 : binlookup=lookupOuterDisk(stub.first);
		break;
	      case 4 : binlookup=lookupOuterDisk(stub.first);
		break;
	      case 1 : binlookup=lookupInnerDisk(stub.first);
		break;
	      case 3 : binlookup=lookupInnerDisk(stub.first);
		break;
	      default : assert(0);  
	      }
	    } else {
	      assert(disk_==1);
	      binlookup=lookupOuterOverlapD1(stub.first);
	    }
            
	    if (binlookup==-1) continue;
	    stub.first->setVMBits(binlookup);

              
	    unsigned int disk=abs(stub.first->disk().value());
	    iphiRaw=iphiRaw/(32/(nallstubsdisks[disk-1]*nvmtedisks[disk-1]));

	    for (unsigned int l=0;l<vmstubsTEPHI_[iphiRaw].size();l++){
	      if (debug1) {
		cout << getName()<<" added stub to : "<<vmstubsTEPHI_[iphiRaw][l]->getName()<<endl;
	      }
	      vmstubsTEPHI_[iphiRaw][l]->addStub(stub);
	      insert=true;
	    }
	    
	  }

	  if (!insert) {
	    cout << getName() << " did not insert stub"<<endl;
	  }
	  assert(insert);

	}
      }
    }


    if (writeVMOccupancyTE) {
      static ofstream out("vmoccupancyte.txt");
      
      for (int i=0;i<24;i++) {
	if (vmstubsTEPHI_[i].size()!=0) {
	  out<<vmstubsTEPHI_[i][0]->getName()<<" "<<vmstubsTEPHI_[i][0]->nStubs()<<endl;
	}
      }
    }
    
    //cout << "Done in FPGAVMRouterTE"<<endl;
    
  }


  void executeME(){

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

	    int layer=stub.first->layer().value();
	    int disk=abs(stub.first->disk().value());

	    //cout << getName()<<" layer,disk : "<<layer<<" "<<disk<<endl;
	    
	    int nvm=-1;
	    if (layer!=-1) {
	      nvm=nallstubslayers[layer]*nvmmelayers[layer];
	    }
	    if (disk!=0){
	      nvm=nallstubsdisks[disk-1]*nvmmedisks[disk-1];
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

	    assert(rfine>=0);
	    
	    stub.first->setfiner(rfine);

	  }

	  if (layer_!=0) {

	    int index=(stub.first->z().value()>>(stub.first->z().nbits()-nbitsfinebintable_))&((1<<nbitsfinebintable_)-1);
	    
	    int zfine=finebintable_[index];

	    stub.first->setfinez(zfine);

	  }
	  
	  bool insert=false;
	  
	  
	  for (unsigned int l=0;l<vmstubsMEPHI_[iphiRaw].size();l++){
	    if (debug1) {
	      cout << "FPGAVMRouterME "<<getName()<<" add stub ( r = "<<stub.second->r()<<" phi = "<<stub.second->phi()<<" ) in : "<<vmstubsMEPHI_[iphiRaw][l]->getName()<<" iphistub = " << iphistub << " iphivmRaw Minus Plus "<<stub.first->iphivmRaw()<<" "<<stub.first->iphivmRawMinus()<<" "<<stub.first->iphivmRawPlus()<<" bins "
		   <<iphiRawMinus<<" "<<iphiRawPlus<<endl;
	    }
	    vmstubsMEPHI_[iphiRaw][l]->addStub(stub);
	    insert=true;
	  }

	  if (iphiRaw!=iphiRawPlus) {
	    for (unsigned int l=0;l<vmstubsMEPHI_[iphiRawPlus].size();l++){
	      vmstubsMEPHI_[iphiRawPlus][l]->addStub(stub);
	    }
	  }
	  if (iphiRaw!=iphiRawMinus) {
	    for (unsigned int l=0;l<vmstubsMEPHI_[iphiRawMinus].size();l++){
	      vmstubsMEPHI_[iphiRawMinus][l]->addStub(stub);
	    }
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
	  if (vmstubsMEPHI_[i].size()!=0) {
	    out<<vmstubsMEPHI_[i][0]->getName()<<" "<<vmstubsMEPHI_[i][0]->nStubs();
	    for (unsigned int ibin=0;ibin<MEBinsDisks*2;ibin++){
	      out <<" "<<vmstubsMEPHI_[i][0]->nStubsBin(ibin);
	    }
	    out<<endl;
	  }
	}

      }

    }


  }



  

  
  int lookupOuterOverlapD1(FPGAStub* stub){

    assert(disk_==1);
    
    static FPGATETableOuterDisk outerTableOverlapD1;
    static bool first=true;

    if (first) {
      outerTableOverlapD1.init(1,7,3);
      first=false;
    }
    
    FPGAWord r=stub->r();
    FPGAWord z=stub->z();
    int rbin=(r.value())>>(r.nbits()-7);
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-3);
    bool negdisk=stub->disk().value()<0;
    if (negdisk) zbin=7-zbin; //Should this be separate table?
    return outerTableOverlapD1.lookup(rbin,zbin);

  }


  int lookupOuterDisk(FPGAStub* stub){

    assert(disk_==2||disk_==4);
    
    static FPGATETableOuterDisk outerTableD2;
    static FPGATETableOuterDisk outerTableD4;
    static bool first=true;

    if (first) {
      outerTableD2.init(2,7,3);
      outerTableD4.init(4,7,3);
      first=false;
    }
    
    FPGAWord r=stub->r();
    FPGAWord z=stub->z();
    int rbin=(r.value())>>(r.nbits()-7);
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-3);
    bool negdisk=stub->disk().value()<0;
    if (negdisk) zbin=7-zbin; //Should this be separate table?
    switch (disk_){
    case 2: return outerTableD2.lookup(rbin,zbin);
      break;
    case 4: return outerTableD4.lookup(rbin,zbin);
      break;
    }
    assert(0);
  }


  int lookupInnerDisk(FPGAStub* stub){

    assert(disk_==1||disk_==3);
    
    static FPGATETableInnerDisk innerTableD1;
    static FPGATETableInnerDisk innerTableD3;
    static bool first=true;

    if (first) {
      innerTableD1.init(1,2,7,3);
      innerTableD3.init(3,4,7,3);
      first=false;
    }
    
    FPGAWord r=stub->r();
    FPGAWord z=stub->z();
    int rbin=(r.value())>>(r.nbits()-7);
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-3);
    bool negdisk=stub->disk().value()<0;
    if (negdisk) zbin=7-zbin; //Should this be separate table?
    switch (disk_){
    case 1: return innerTableD1.lookup(rbin,zbin);
      break;
    case 3: return innerTableD3.lookup(rbin,zbin);
      break;
    }
    assert(0);
  }

  int lookupOuterLayer(FPGAStub* stub){

    assert(layer_==2||layer_==3||layer_==4||layer_==6);
    
    static FPGATETableOuter outerTableL2;
    static FPGATETableOuter outerTableL3;
    static FPGATETableOuter outerTableL4;
    static FPGATETableOuter outerTableL6;
    static bool first=true;

    if (first) {
      outerTableL2.init(2,7,4);
      outerTableL3.init(3,7,4);
      outerTableL4.init(4,7,4);
      outerTableL6.init(6,7,4);
      first=false;
    }
    
    FPGAWord r=stub->r();
    FPGAWord z=stub->z();
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-7);
    int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-4);
    switch (layer_){
    case 2: return outerTableL2.lookup(zbin,rbin);
      break;
    case 3: return outerTableL3.lookup(zbin,rbin);
      break;
    case 4: return outerTableL4.lookup(zbin,rbin);
      break;
    case 6: return outerTableL6.lookup(zbin,rbin);
      break;
    }
    assert(0);
  }


  int lookupInnerLayer(FPGAStub* stub){

    assert(layer_==1||layer_==2||layer_==3||layer_==5);
    
    static FPGATETableInner innerTableL1;
    static FPGATETableInner innerTableL2;
    static FPGATETableInner innerTableL3;
    static FPGATETableInner innerTableL5;
    static bool first=true;

    if (first) {
      innerTableL1.init(1,2,7,4);
      innerTableL2.init(2,3,7,4);
      innerTableL3.init(3,4,7,4);
      innerTableL5.init(5,6,7,4);
      first=false;
    }
    
    FPGAWord r=stub->r();
    FPGAWord z=stub->z();
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-7);
    int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-4);
    switch (layer_){
    case 1: return innerTableL1.lookup(zbin,rbin);
      break;
    case 2: return innerTableL2.lookup(zbin,rbin);
      break;
    case 3: return innerTableL3.lookup(zbin,rbin);
      break;
    case 5: return innerTableL5.lookup(zbin,rbin);
      break;
    }
    assert(0);
  }


  int lookupInnerOverlapLayer(FPGAStub* stub){

    assert(layer_==1||layer_==2);
    
    static FPGATETableInnerOverlap innerTableOverlapL1;
    static FPGATETableInnerOverlap innerTableOverlapL2;
    static bool first=true;

    if (first) {
      innerTableOverlapL1.init(1,1,7,3);
      innerTableOverlapL2.init(2,1,7,3);
      first=false;
    }
    
    FPGAWord r=stub->r();
    FPGAWord z=stub->z();
    int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-7);
    int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-3);
    switch (layer_){
    case 1: return innerTableOverlapL1.lookup(zbin,rbin);
      break;
    case 2: return innerTableOverlapL2.lookup(zbin,rbin);
      break;
    }
    assert(0);
  }


  

private:

  int layer_;
  int disk_;

  int nbitsfinebintable_;
  vector<int> finebintable_;

  vector<FPGAInputLink*> stubinputs_;
  vector<FPGAAllStubs*> allstubs_;

  vector<FPGAVMStubsTE*> vmstubsTEPHI_[32];
  vector<FPGAVMStubsTE*> vmstubsTEExtraPHI_[32];
  vector<FPGAVMStubsTE*> vmstubsTEOverlapPHI_[32];
  vector<FPGAVMStubsME*> vmstubsMEPHI_[32];


};

#endif

