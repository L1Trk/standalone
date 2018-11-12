//This class implementes the VM router
#ifndef FPGAVMROUTERTE_H
#define FPGAVMROUTERTE_H

#include "FPGAProcessBase.hh"
#include "FPGATETableOuter.hh"
#include "FPGATETableInner.hh"
#include "FPGATETableOuterDisk.hh"
#include "FPGATETableInnerDisk.hh"
#include "FPGATETableInnerOverlap.hh"

using namespace std;

class FPGAVMRouterTE:public FPGAProcessBase{

public:

  FPGAVMRouterTE(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){

    layer_=0;
    disk_=0;
    
    if (name_[6]=='L') layer_=name_[7]-'0';    
    if (name_[6]=='D') disk_=name_[7]-'0';    
    if (name_[6]=='F') disk_=name_[7]-'0';    
    if (name_[6]=='B') disk_=-(name_[7]-'0');

    overlap_=false;
    
    overlap_=(name_[11]=='X'||name_[11]=='Y'||name_[11]=='W'||name_[11]=='Q');

    assert((layer_!=0)||(disk_!=0));

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
    for (int i=0;i<24;i++) {
      std::ostringstream oss;
      oss<<(i+1);
      string s=oss.str();
      if (output=="vmstuboutA"+s+"n1"||
	  output=="vmstuboutB"+s+"n1"||
	  output=="vmstuboutC"+s+"n1"||
	  output=="vmstuboutD"+s+"n1"||
	  output=="vmstuboutE"+s+"n1"||
	  output=="vmstuboutF"+s+"n1"||
	  output=="vmstuboutX"+s+"n1"||
	  output=="vmstuboutY"+s+"n1"||
	  output=="vmstuboutW"+s+"n1"||
	  output=="vmstuboutQ"+s+"n1"||
	  output=="vmstuboutA"+s+"n2"||
	  output=="vmstuboutB"+s+"n2"||
	  output=="vmstuboutC"+s+"n2"||
	  output=="vmstuboutD"+s+"n2"||
	  output=="vmstuboutE"+s+"n2"||
	  output=="vmstuboutF"+s+"n2"||
	  output=="vmstuboutX"+s+"n2"||
	  output=="vmstuboutY"+s+"n2"||
	  output=="vmstuboutW"+s+"n2"||
	  output=="vmstuboutQ"+s+"n2"||
	  output=="vmstuboutA"+s+"n3"||
	  output=="vmstuboutB"+s+"n3"||
	  output=="vmstuboutC"+s+"n3"||
	  output=="vmstuboutD"+s+"n3"||
	  output=="vmstuboutE"+s+"n3"||
	  output=="vmstuboutF"+s+"n3"||
	  output=="vmstuboutX"+s+"n3"||
	  output=="vmstuboutY"+s+"n3"||
	  output=="vmstuboutW"+s+"n3"||
	  output=="vmstuboutQ"+s+"n3"||
	  output=="vmstuboutA"+s+"n4"||
	  output=="vmstuboutB"+s+"n4"||
	  output=="vmstuboutC"+s+"n4"||
	  output=="vmstuboutD"+s+"n4"||
	  output=="vmstuboutE"+s+"n4"||
	  output=="vmstuboutF"+s+"n4"||
	  output=="vmstuboutA"+s+"n5"||
	  output=="vmstuboutB"+s+"n5"||
	  output=="vmstuboutC"+s+"n5"||
	  output=="vmstuboutD"+s+"n5"||
	  output=="vmstuboutE"+s+"n5"||
	  output=="vmstuboutF"+s+"n5"||
	  output=="vmstuboutA"+s+"n6"||
	  output=="vmstuboutB"+s+"n6"||
	  output=="vmstuboutC"+s+"n6"||
	  output=="vmstuboutD"+s+"n6"||
	  output=="vmstuboutE"+s+"n6"||
	  output=="vmstuboutF"+s+"n6"||
	  output=="vmstuboutA"+s+"n7"||
	  output=="vmstuboutB"+s+"n7"||
	  output=="vmstuboutC"+s+"n7"||
	  output=="vmstuboutD"+s+"n7"||
	  output=="vmstuboutE"+s+"n7"||
	  output=="vmstuboutF"+s+"n7"||
	  output=="vmstuboutA"+s+"n8"||
	  output=="vmstuboutB"+s+"n8"||
	  output=="vmstuboutC"+s+"n8"||
	  output=="vmstuboutD"+s+"n8"||
	  output=="vmstuboutE"+s+"n8"||
	  output=="vmstuboutF"+s+"n8"||
	  output=="vmstuboutA"+s+"n9"||
	  output=="vmstuboutB"+s+"n9"||
	  output=="vmstuboutC"+s+"n9"||
	  output=="vmstuboutD"+s+"n9"||
	  output=="vmstuboutE"+s+"n9"||
	  output=="vmstuboutF"+s+"n9"||
	  output=="vmstuboutA"+s+"n10"||
	  output=="vmstuboutB"+s+"n10"||
	  output=="vmstuboutC"+s+"n10"||
	  output=="vmstuboutD"+s+"n10"||
	  output=="vmstuboutE"+s+"n10"||
	  output=="vmstuboutF"+s+"n10"
	  ){
	FPGAVMStubsTE* tmp=dynamic_cast<FPGAVMStubsTE*>(memory);
	assert(tmp!=0);
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
  
    assert(stubinputs_.size()!=0);

    unsigned int count=0;

    vector<unsigned int> asindex_count = {0,0,0,0,0};
 
    if (layer_!=0){  //First handle layer stubs
      for(unsigned int j=0;j<stubinputs_.size();j++){
	for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	  count++;
	  if (count>MAXVMROUTER) continue;
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);
          
	  int iphiRaw=stub.first->iphivmRaw();

	  bool insert=false;
      
      // get preliminary stub index set in InputLink
      unsigned int asindex_pre = stub.first->stubindex().value();
      unsigned int asindex = asindex_pre;
      // correct the index
      // necessary if the InputLink is not the first one read by VMRouterTE
      unsigned int phiregion = stub.first->phiregion().value();
      if (phiregion==0) asindex = asindex_pre + asindex_count[0];
      else if (phiregion==1) asindex = asindex_pre + asindex_count[1];
      else if (phiregion==2) asindex = asindex_pre + asindex_count[2];
      else if (phiregion==3) asindex = asindex_pre + asindex_count[3];
      else {
        // Stub iphiRaw<4 or iphiRaw>27. It is used in TE but not ME, and is not stored in the AS memory filled by VMRouterME. Use an extra index to keep track of them.
        asindex = asindex_pre + asindex_count[4];
      }

      stub.first->setAllStubIndex(asindex);
      stub.second->setAllStubIndex(asindex);
      
      stub.first->setAllStubAddressTE(allstubs_[0]->nStubs());
		  
	  for (unsigned int l=0;l<allstubs_.size();l++){
	    allstubs_[l]->addStub(stub);
	  }
      
	  if (getName()=="VMRTE_L2PHIW"||getName()=="VMRTE_L2PHIQ") {
	    //special case where even even layer is treated as an odd (inner layer)
            assert(layer_==2);
            int binlookup=lookupInnerOverlapLayer(stub.first);
            if (binlookup==-1) continue;
            stub.first->setVMBitsOverlap(binlookup);

	    iphiRaw-=4;
	    assert(iphiRaw>=0);
	    assert(iphiRaw<24);
	    if (overlap_) iphiRaw>>=2;

	    for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
	      if (debug1) {
		cout << "FPGAVMRouterTE added stub to : "<<vmstubsPHI_[iphiRaw][l]->getName()<<endl;
	      }
	      vmstubsPHI_[iphiRaw][l]->addStub(stub);
	      insert=true;
	    }
	  }  else {
        int binlookup=-1;
        if (overlap_) {
          assert(layer_==1);
          binlookup=lookupInnerOverlapLayer(stub.first);
        } else {
          switch (layer_) {
          case 2 : binlookup=lookupOuterLayer(stub.first);
            break;
          case 4 : binlookup=lookupOuterLayer(stub.first);
            break;
          case 6 : binlookup=lookupOuterLayer(stub.first);
            break;
          case 1 : binlookup=lookupInnerLayer(stub.first);
            break;
          case 3 : binlookup=lookupInnerLayer(stub.first);
            break;
          case 5 : binlookup=lookupInnerLayer(stub.first);
            break;
          default : assert(0);
          }
        }
        if (binlookup==-1) continue;
	    if (overlap_) {
	      stub.first->setVMBitsOverlap(binlookup);
	    } else {
	      stub.first->setVMBits(binlookup);
	    }
	    if (stub.first->layer().value()%2==0) { //odd layers
	      iphiRaw-=4;
	      assert(iphiRaw>=0);
	      assert(iphiRaw<24);
	      if (overlap_) {
            iphiRaw>>=2;
          }
	      for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
		vmstubsPHI_[iphiRaw][l]->addStub(stub);
		if (debug1) {
		  cout << getName()<<" adding stub to "<<vmstubsPHI_[iphiRaw][l]->getName()<<endl;
		}
		insert=true;
	      }
	    } else {  //even layers
	      iphiRaw/=2;
	      assert(iphiRaw>=0);
	      assert(iphiRaw<16);
	      for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
		vmstubsPHI_[iphiRaw][l]->addStub(stub);
		if (debug1) {
		  cout << getName()<<" adding stub to "<<vmstubsPHI_[iphiRaw][l]->getName()<<endl;
		}
		insert=true;
	      }
	    }
	  }

	  assert(insert);
	} // end of stubs loop in each InputLink

    auto IL_indices = stubinputs_[j]->getASPhiIndices();
    asindex_count[0] += IL_indices[0];
    asindex_count[1] += IL_indices[1];
    asindex_count[2] += IL_indices[2];
    asindex_count[3] += IL_indices[3];
    asindex_count[4] += IL_indices[4];
    
      } // end of InputLinks loop

    }
    if (disk_!=0) {
      for(unsigned int j=0;j<stubinputs_.size();j++){
	//cout << "Input link "<<stubinputs_[j]->getName()<<endl;
	for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);

	  
	  int iphiRaw=stub.first->iphivmRaw();

	  bool insert=false;

      // get preliminary stub index set in InputLink
      unsigned int asindex_pre = stub.first->stubindex().value();
      unsigned int asindex = asindex_pre;
      // correct the index
      // necessary if the InputLink is not the first one read by VMRouterTE
      unsigned int phiregion = stub.first->phiregion().value();
      if (phiregion==0) asindex = asindex_pre + asindex_count[0];
      else if (phiregion==1) asindex = asindex_pre + asindex_count[1];
      else if (phiregion==2) asindex = asindex_pre + asindex_count[2];
      else if (phiregion==3) asindex = asindex_pre + asindex_count[3];
      else {
        // Stub iphiRaw<4 or iphiRaw>27. It is used in TE but not ME, and is not stored in the AS memory filled by VMRouterME. Use an extra index to keep track of them.
        asindex = asindex_pre + asindex_count[4];
      }

      stub.first->setAllStubIndex(asindex);
      stub.second->setAllStubIndex(asindex);

	  stub.first->setAllStubAddressTE(allstubs_[0]->nStubs());
      
	  if (!stub.second->isPSmodule()) continue;
	  
	  for (unsigned int l=0;l<allstubs_.size();l++){
        //if (allstubs_[l]->getName()=="AS_D1PHIQn1") {
        //  cout << "Adding stub to memory "<<iSector_<<" "<<allstubs_[l]->getName()<<" r= "<<stub.second->r()<<endl;
        //}
        allstubs_[l]->addStub(stub);
	  }
	  
	  if (getName()=="VMRTE_D1PHIW"||getName()=="VMRTE_D1PHIQ") {
        //special case where odd disk is treated as outer disk
        
        if (overlap_){
	  assert(stub.first->layer().value()<0); //can not handle stubs in layers here
          int binlookup=lookupOuterOverlapD1(stub.first);
          assert(binlookup>=0);
          stub.first->setVMBitsOverlap(binlookup);
        } else {
          assert(0);
          
        }
            
        iphiRaw/=2;
        assert(iphiRaw>=0);
        assert(iphiRaw<16);
        iphiRaw>>=1; //only 8 VMS in even disks
        for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
          if (debug1) { 
            cout << "FPGAVMRouterTE added stub to : "<<vmstubsPHI_[iphiRaw][l]->getName()<<endl;
          }
          vmstubsPHI_[iphiRaw][l]->addStub(stub);
          insert=true;
        }
	  } else {
        if (getName()=="VMRTE_D1PHIA"||
            getName()=="VMRTE_D1PHIB"||
            getName()=="VMRTE_D1PHIC"||
            getName()=="VMRTE_D1PHID"||
            getName()=="VMRTE_D2PHIA"||
            getName()=="VMRTE_D2PHIB"||
            getName()=="VMRTE_D2PHIC"||
            getName()=="VMRTE_D2PHID"
            ) {
          
          int binlookup=-1;
          
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
          if (binlookup==-1) continue;
          stub.first->setVMBits(binlookup);

          
          
          if (abs(stub.first->disk().value())%2==1) {
            //odd disks here
            if (overlap_) {
              iphiRaw>>=1; //only 12 VMS in odd disks
            } else {
              iphiRaw-=4;
              assert(iphiRaw>=0);
              assert(iphiRaw<24);
              iphiRaw>>=1; //only 12 VMS in odd disks
            }
			
            for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
              vmstubsPHI_[iphiRaw][l]->addStub(stub);
              insert=true;
            }
          }
          else {
            //even disks here
            iphiRaw/=2;
            assert(iphiRaw>=0);
            assert(iphiRaw<16);
            for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
              vmstubsPHI_[iphiRaw][l]->addStub(stub);
              insert=true;
            }
          }
        } else {
          
          int binlookup=-1;
          
          if (!overlap_) {
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
	  if (overlap_) {
	    stub.first->setVMBitsOverlap(binlookup);
	  } else {
	    stub.first->setVMBits(binlookup);
          }
          
	      if (abs(stub.first->disk().value())%2==1) {
		//odd disks here
		if (overlap_) {
		  iphiRaw>>=2; //only 6 VMS in odd disks
		} else {
		  iphiRaw-=4;
		  assert(iphiRaw>=0);
		  assert(iphiRaw<24);
		  iphiRaw>>=2; //only 6 VMS in odd disks
		}
		for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
		  vmstubsPHI_[iphiRaw][l]->addStub(stub);
		  insert=true;
		}
	      }
	      else {
		//even disks here
		iphiRaw/=2;
		assert(iphiRaw>=0);
		assert(iphiRaw<16);
		iphiRaw>>=1; //only 8 VMS in even disks
		for (unsigned int l=0;l<vmstubsPHI_[iphiRaw].size();l++){
		  vmstubsPHI_[iphiRaw][l]->addStub(stub);
		  insert=true;
		}
	      }

	    }
	  }

	  if (!insert) {
	    cout << getName() << " did not insert stub"<<endl;
	  }
	  assert(insert);

	} // end of stubs loop in each InputLink

    auto IL_indices = stubinputs_[j]->getASPhiIndices();
    asindex_count[0] += IL_indices[0];
    asindex_count[1] += IL_indices[1];
    asindex_count[2] += IL_indices[2];
    asindex_count[3] += IL_indices[3];
    asindex_count[4] += IL_indices[4];
    
      } // end of InputLinks loop
    }


    if (writeAllStubs) {
      static ofstream out("allstubste.txt");
      out<<allstubs_[0]->getName()<<" "<<allstubs_[0]->nStubs()<<endl;
      //if (allstubs_[0]->getName()=="AS_D1PHIQn1") {
      //	cout << "Number of stubs in : "<<allstubs_[0]->getName()<<" "<<allstubs_[0]->nStubs()<<endl;
      //}
    }


    if (writeVMOccupancyTE) {
      static ofstream out("vmoccupancyte.txt");
      
      for (int i=0;i<24;i++) {
	if (vmstubsPHI_[i].size()!=0) {
	  out<<vmstubsPHI_[i][0]->getName()<<" "<<vmstubsPHI_[i][0]->nStubs()<<endl;
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

    assert(layer_==2||layer_==4||layer_==6);
    
    static FPGATETableOuter outerTableL2;
    static FPGATETableOuter outerTableL4;
    static FPGATETableOuter outerTableL6;
    static bool first=true;

    if (first) {
      outerTableL2.init(2,7,4);
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
    case 4: return outerTableL4.lookup(zbin,rbin);
      break;
    case 6: return outerTableL6.lookup(zbin,rbin);
      break;
    }
    assert(0);
  }


  int lookupInnerLayer(FPGAStub* stub){

    assert(layer_==1||layer_==3||layer_==5);
    
    static FPGATETableInner innerTableL1;
    static FPGATETableInner innerTableL3;
    static FPGATETableInner innerTableL5;
    static bool first=true;

    if (first) {
      innerTableL1.init(1,2,7,4);
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

  bool overlap_;	     
	     
  vector<FPGAInputLink*> stubinputs_;
  vector<FPGAAllStubs*> allstubs_;

  vector<FPGAVMStubsTE*> vmstubsPHI_[24];


};

#endif

