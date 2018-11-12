//This class implementes the tracklet engine
#ifndef FPGAMATCHENGINE_H
#define FPGAMATCHENGINE_H

#include "FPGAProcessBase.hh"
#include "FPGATrackletCalculator.hh"

using namespace std;

class FPGAMatchEngine:public FPGAProcessBase{

public:

  FPGAMatchEngine(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    layer_=0;
    disk_=0;
    string subname=name.substr(3,2);
    if (!hourglass) {
      subname=name.substr(8,2);
    }
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    if (subname=="D1") disk_=1;
    if (subname=="D2") disk_=2;
    if (subname=="D3") disk_=3;
    if (subname=="D4") disk_=4;
    if (subname=="D5") disk_=5;
    if (layer_==0&&disk_==0) {
      cout << name<<" subname = "<<subname<<" "<<layer_<<" "<<disk_<<endl;
    }
    assert((layer_!=0)||(disk_!=0));

    if (layer_>0) {

      unsigned int nbits=3;
      if (layer_>=4) nbits=4;
      
      for(unsigned int irinv=0;irinv<32;irinv++){
	double rinv=(irinv-15.5)*(1<<(nbitsrinv-5))*krinvpars;
	double projbend=bend(rmean[layer_-1],rinv);
	for(unsigned int ibend=0;ibend<(unsigned int)(1<<nbits);ibend++){
	  double stubbend=FPGAStub::benddecode(ibend,layer_<=3);
	  bool pass=fabs(stubbend-projbend)<2.0;
	  table_.push_back(pass);
	}
      }
      
    }

    if (disk_>0) {

      for(unsigned int iprojbend=0;iprojbend<32;iprojbend++){
	double projbend=0.5*(iprojbend-15.0);
	for(unsigned int ibend=0;ibend<8;ibend++){
	  double stubbend=FPGAStub::benddecode(ibend,true);
	  bool pass=fabs(stubbend-projbend)<1.5;
	  tablePS_.push_back(pass);
	}
	for(unsigned int ibend=0;ibend<16;ibend++){
	  double stubbend=FPGAStub::benddecode(ibend,false);
	  bool pass=fabs(stubbend-projbend)<1.5;
	  table2S_.push_back(pass);
	}
      }
      
    }

    
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="matchout") {
      FPGACandidateMatch* tmp=dynamic_cast<FPGACandidateMatch*>(memory);
      assert(tmp!=0);
      candmatches_=tmp;
      return;
    }
    assert(0);

  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="vmstubin") {
      FPGAVMStubsME* tmp=dynamic_cast<FPGAVMStubsME*>(memory);
      assert(tmp!=0);
      vmstubs_=tmp;
      return;
    }
    if (input=="vmprojin") {
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      vmprojs_=tmp;
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute() {

    unsigned int countall=0;
    unsigned int countpass=0;

    for(unsigned int j=0;j<vmprojs_->nTracklets();j++){
      FPGATracklet* proj=vmprojs_->getFPGATracklet(j);

      int nmatches=0;
      
      if (debug1) {
	cout << "Found projection in "<<getName()<<endl;
      }
     	
      if (layer_>0){

	unsigned int zbin1 =proj->zbin1projvm(layer_);
	unsigned int zbin2 = zbin1; 

	unsigned int finez =proj->finezvm(layer_);

	if (proj->zbin2projvm(layer_)==1) zbin2 += 1;
	
	for (unsigned int ibin=zbin1;ibin<=zbin2;ibin++) {

	  unsigned int nstub=vmstubs_->nStubsBin(ibin);
	  
	  for(unsigned int i=0;i<nstub;i++){
	    if (debug1) {
	      cout << "Found stub in "<<getName()<<endl;
	    }
	    std::pair<FPGAStub*,L1TStub*> stub=vmstubs_->getStubBin(ibin,i);
	    countall++;

	    if (finephiME) {
	      FPGAWord projphi=proj->fpgaphiproj(layer_);
	      FPGAWord stubphicorr=stub.first->phicorr();

	      int ifinephiproj=(projphi.value()>>(projphi.nbits()-8));
	      int ifinephistubcorr=(stubphicorr.value()>>(stubphicorr.nbits()-8));
	      
	      if (abs(ifinephiproj-ifinephistubcorr)>1) {
		continue;
	      }
	    }
							   
	    int z=stub.first->finez().value();

	    if (ibin!=zbin1) z+=8;
	    
	    int idz=z-finez;
	    	    
	    double dz=proj->zproj(layer_)-stub.second->z();
	    
	    if (proj->layer()==1) {
	      //cout << getName()<<" dz idz : "<<dz<<" "<<idz
	      //   <<" proj: "<<proj->zproj(layer_)<<" "<<proj->fpgazproj(layer_).value()*16*kz<<" "
	      //   <<proj->fpgazproj(layer_).value()<<" "
	      //   <<((proj->fpgazproj(layer_).value()>>(proj->fpgazproj(layer_).nbits()-6))&7)
	      //   <<" "<<finez
	      //   <<" stub: "<<stub.second->z()<<" "<<stub.first->z().value()*16*kz<<" "
	      //   <<stub.first->z().value()<<" "
	      //   <<((stub.first->z().value()>>(stub.first->z().nbits()-6))&7)
	      //   <<" "<<z<<endl;
	      if (abs(idz)>2) {
		if (debug1) {
		  cout << getName()<<" Match rejected for L1L2 seed with dz = "
		       <<dz<<" idz = "<<idz<<endl;
		}
		continue;
	      }
	    } else {
	      if (abs(idz)>5) continue; //needed for L5L6 seeds...
	    }

	    int irinvvm=16+(proj->fpgarinv().value()>>(proj->fpgarinv().nbits()-5));

	    int nbits=3;
	    if (layer_>=4) nbits=4;
	    assert(nbits==stub.first->bend().nbits());
	    
	    unsigned int index=(irinvvm<<nbits)+stub.first->bend().value();

	    assert(index<table_.size());
	    bool pass=table_[index];

	    if (!pass) {
	      if (debug1) {
		cout << "Match rejected with bend lookup index = "
		     <<index<<endl; 
	      }
	      continue;
	    }
	    if (debug1) {
	      cout << "Adding match in "<<getName()<<endl;
	    }
	    
	    countpass++;
	    if (nmatches<1000) {
	      candmatches_->addMatch(proj,stub);
	    }
	    nmatches++;
	    if (countall>=MAXME) break;
	  }
	}
      } // if (layer_>0)      
      else if (disk_!=0) {

	unsigned int rbin1 =proj->rbin1projvm(disk_);
	unsigned int rbin2 = rbin1; 

	unsigned int finer =proj->finervm(disk_);

	if (proj->rbin2projvm(disk_)==1) rbin2 += 1;
	
	for (unsigned int ibin=rbin1;ibin<=rbin2;ibin++) {

	  unsigned int nstub=vmstubs_->nStubsBin(ibin);

	  for(unsigned int i=0;i<nstub;i++){
	    if (debug1) {
	      cout << "Found stub in "<<getName()<<endl;
	    }
	    std::pair<FPGAStub*,L1TStub*> stub=vmstubs_->getStubBin(ibin,i);
	    countall++;

	    if (finephiME) {
	      FPGAWord projphi=proj->fpgaphiprojdisk(disk_);
	      FPGAWord stubphicorr=stub.first->phicorr();
	      
	      int ifinephiproj=(projphi.value()>>(projphi.nbits()-8));
	      int ifinephistubcorr=(stubphicorr.value()>>(stubphicorr.nbits()-8));
	    
	      if (abs(ifinephiproj-ifinephistubcorr)>1) {
		continue;
	      }
	    }
	    
	    
	    int r=stub.first->finer().value();
	    
	    if (ibin!=rbin1) r+=8;
	    
	    int idr=r-finer;

	    if (stub.first->isPSmodule()){
	      if (abs(idr)>1) {
		if (debug1) {
		  cout << getName() << "PS stub rejected with idr = "<<idr<<endl;
		}
		continue;
	      }
	    } else {
	      if (abs(idr)>5) {
		if (debug1) {
		  cout << getName() << "2S stub rejected with idr = "<<idr<<endl;
		}
		continue;
	      }
	    }

	    int ibendproj=proj->getBendIndex(disk_).value();

	    int nbits=3;
	    if (!stub.first->isPSmodule()) nbits=4;

	    assert(nbits==stub.first->bend().nbits());
	    
	    unsigned int index=(ibendproj<<nbits)+stub.first->bend().value();

	    bool pass=false;

	    if (stub.first->isPSmodule()) {
	      assert(index<tablePS_.size());
	      pass=tablePS_[index];
	    } else {
	      assert(index<table2S_.size());
	      pass=table2S_[index];
	    }
	    
	    if (!pass) {
	      if (debug1) {
		cout << getName() << "stub rejected with index = "<<index<<endl;
	      }
	      continue;
	    }
	    
	    
	    
	    countpass++;
	    if (nmatches<1000) {
	      if (debug1) {
		cout << getName() << " adding match "<<endl;
	      }
	      candmatches_->addMatch(proj,stub);
	    }
	    nmatches++;
	    if (countall>=MAXME) break;
	    
	  }
	}
	
      } else {
	//neither disk nor layer?
	assert(0);
      } 
      
	
    } // outer for loop
     
    
    if (writeME) {
      static ofstream out("matchengine.txt");
      out << getName()<<" "<<countall<<" "<<countpass<<endl;
    }

  } // execute()

  double bend(double r, double rinv) {

    double dr=0.18;
    
    double delta=r*dr*0.5*rinv;

    double bend=-delta/0.009;
    if (r<55.0) bend=-delta/0.01;

    return bend;
    
  }

  
private:

  FPGAVMStubsME* vmstubs_;
  FPGAVMProjections* vmprojs_;

  FPGACandidateMatch* candmatches_;

  int layer_;
  int disk_;

  //used in the layers
  vector<bool> table_;

  //used in the disks
  vector<bool> tablePS_;
  vector<bool> table2S_;

};

#endif
