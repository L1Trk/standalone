//This class implementes the projection router
#ifndef FPGAPROJECTIONROUTER_H
#define FPGAPROJECTIONROUTER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAProjectionRouter:public FPGAProcessBase{

public:

  FPGAProjectionRouter(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    string subname=name.substr(8,2);
    if (hourglass) {
      subname=name.substr(3,2);
    }
    //cout << "name subname : "<<name<<" "<<subname<<endl;
    layer_=0;
    disk_=0;
    vmprojPHI1_=0;
    vmprojPHI2_=0;
    vmprojPHI3_=0;
    vmprojPHI4_=0;
    vmprojPHI5_=0;
    vmprojPHI6_=0;
    vmprojPHI7_=0;
    vmprojPHI8_=0;
    
    
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
    assert(disk_!=0||layer_!=0);
    allproj_=0;

    nrbits_=5;
    nphiderbits_=6;

    
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="allprojout"){
      FPGAAllProjections* tmp=dynamic_cast<FPGAAllProjections*>(memory);
      assert(tmp!=0);
      allproj_=tmp;
      return;
    }
    if (hourglass) {

      int nproj=-1;
      int nprojvm=-1;
      if (layer_>0) {
	nproj=nallprojlayers[layer_-1];
	nprojvm=nvmmelayers[layer_-1];
      }
      if (disk_>0) {
	nproj=nallprojdisks[disk_-1];
	nprojvm=nvmmedisks[disk_-1];
      }
      assert(nproj!=-1);

      for (int iproj=0;iproj<nproj;iproj++) {
	for (int iprojvm=0;iprojvm<nprojvm;iprojvm++) {
	  ostringstream oss;
	  oss << "vmprojoutPHI"<<char(iproj+'A')<<iproj*nprojvm+iprojvm+1;
	  string name=oss.str();
	  if (output==name) {
	    FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	    assert(tmp!=0);
	    if (iprojvm==0) vmprojPHI1_=tmp;
	    if (iprojvm==1) vmprojPHI2_=tmp;
	    if (iprojvm==2) vmprojPHI3_=tmp;
	    if (iprojvm==3) vmprojPHI4_=tmp;
	    if (iprojvm==4) vmprojPHI5_=tmp;
	    if (iprojvm==5) vmprojPHI6_=tmp;
	    if (iprojvm==6) vmprojPHI7_=tmp;
	    if (iprojvm==7) vmprojPHI8_=tmp;
	    return;
	  }
	}
      }

      /*
      if (output=="vmprojoutPHIA1"||output=="vmprojoutPHIB5"||output=="vmprojoutPHIC9"||output=="vmprojoutPHID13"||output=="vmprojoutPHIE17"||output=="vmprojoutPHIF21"||output=="vmprojoutPHIG25"||output=="vmprojoutPHIH29"){
	FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	assert(tmp!=0);
	vmprojPHI1_=tmp;
	return;
      }
      if (output=="vmprojoutPHIA2"||output=="vmprojoutPHIB6"||output=="vmprojoutPHIC10"||output=="vmprojoutPHID14"||output=="vmprojoutPHIE18"||output=="vmprojoutPHIF22"||output=="vmprojoutPHIG26"||output=="vmprojoutPHIH30"){
	FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	assert(tmp!=0);
	vmprojPHI2_=tmp;
	return;
      }
      if (output=="vmprojoutPHIA3"||output=="vmprojoutPHIB7"||output=="vmprojoutPHIC11"||output=="vmprojoutPHID15"||output=="vmprojoutPHIE19"||output=="vmprojoutPHIF23"||output=="vmprojoutPHIG27"||output=="vmprojoutPHIH31"){
	FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	assert(tmp!=0);
	vmprojPHI3_=tmp;
	return;
      }
      if (output=="vmprojoutPHIA4"||output=="vmprojoutPHIB8"||output=="vmprojoutPHIC12"||output=="vmprojoutPHID16"||output=="vmprojoutPHIE20"||output=="vmprojoutPHIF24"||output=="vmprojoutPHIG28"||output=="vmprojoutPHIH32"){
	FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	assert(tmp!=0);
	vmprojPHI4_=tmp;
	return;
	}
      */
    } else {
      if (layer_==2||layer_==4||layer_==6) {
	if (output=="vmprojoutPHI1"||output=="vmprojoutPHI5"||output=="vmprojoutPHI9"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI1_=tmp;
	  return;
	}
	if (output=="vmprojoutPHI2"||output=="vmprojoutPHI6"||output=="vmprojoutPHI10"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI2_=tmp;
	  return;
	}
	if (output=="vmprojoutPHI3"||output=="vmprojoutPHI7"||output=="vmprojoutPHI11"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI3_=tmp;
	  return;
	}
	if (output=="vmprojoutPHI4"||output=="vmprojoutPHI8"||output=="vmprojoutPHI12"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI4_=tmp;
	  return;
	}      
      } else {
	if (output=="vmprojoutPHI1"||output=="vmprojoutPHI5"||output=="vmprojoutPHI9"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	vmprojPHI1_=tmp;
	return;
	}
	if (output=="vmprojoutPHI2"||output=="vmprojoutPHI6"||output=="vmprojoutPHI10"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI2_=tmp;
	  return;
	}
	if (output=="vmprojoutPHI3"||output=="vmprojoutPHI7"||output=="vmprojoutPHI11"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI3_=tmp;
	  return;
	}
	if (output=="vmprojoutPHI4"||output=="vmprojoutPHI8"||output=="vmprojoutPHI12"){
	  FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
	  assert(tmp!=0);
	  vmprojPHI4_=tmp;
	return;
	}
      }
    }


    
    cout << "Did not find output : "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="proj1in"||input=="proj2in"||
	input=="proj3in"||input=="proj4in"||
	input=="proj5in"||input=="proj6in"||
	input=="proj7in"||input=="proj8in"||
	input=="proj9in"||input=="proj10in"||
	input=="proj11in"||input=="proj12in"||
	input=="proj13in"||input=="proj14in"||
	input=="proj15in"||input=="proj16in"||
	input=="proj17in"||input=="proj18in"||
	input=="proj19in"||input=="proj20in"||
	input=="proj21in"||input=="proj22in"||
	input=="proj23in"||input=="proj24in"||
	input=="proj25in"||input=="proj26in"||
	input=="proj27in"||input=="proj28in"||
	input=="proj29in"||input=="proj30in"||
	input=="proj31in"||input=="proj32in"||
	input=="proj33in"||input=="proj34in"||
	input=="proj35in"||input=="proj36in"||
	input=="proj37in"||input=="proj38in"||
	input=="proj39in"||input=="proj40in"||
	input=="proj41in"||input=="proj42in"||
	input=="proj43in"||input=="proj44in"||
	input=="proj45in"||input=="proj46in"||
	input=="proj47in"||input=="proj48in"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputproj_.push_back(tmp);
      return;
    }
    cout << "Could not find input : "<<input<<" in "<<getName()<<endl;
    assert(0);
  }

  void execute() {

    //cout << "FPGAProjectionRouter::execute : "<<getName()<<" "<<inputproj_.size()<<endl;

    unsigned int count=0;

    int lastTCID=-1;
    int lastTCIDplus=-1;
    int lastTCIDminus=-1;
    
    if (layer_!=0) {
      for (unsigned int j=0;j<inputproj_.size();j++){

	for (unsigned int i=0;i<inputproj_[j]->nTracklets();i++){

	  count++;
	  if (count>MAXPROJROUTER) continue;

	  FPGAWord fpgaphi=inputproj_[j]->getFPGATracklet(i)->fpgaphiproj(layer_);
	  FPGAWord fpgaz=inputproj_[j]->getFPGATracklet(i)->fpgazproj(layer_);

	  
	  //skip if projection is out of range!
	  if (fpgaz.atExtreme()) continue;
	  if (fpgaphi.atExtreme()) continue;


	  int iphitmp=fpgaphi.value();
	  int iphi=iphitmp>>(fpgaphi.nbits()-5);

	  if (hourglass) {
	    int nvm=-1;
	    int nbins=-1;
	    if (layer_!=0) {
	      nvm=nvmmelayers[layer_-1]*nallstubslayers[layer_-1];
	      nbins=nvmmelayers[layer_-1];
	    }
	    if (disk_!=0){
	      assert(0); //should never get here
	      nvm=nvmmedisks[disk_-1]*nallstubsdisks[disk_-1];
	      nbins=nvmmedisks[disk_-1];
	    }
	    assert(nvm>0);

	    iphi=(iphi/(32/nvm))&(nbins-1);	    
	  } else {
	    assert(iphi>=4);
	    assert(iphi<=27);
	    iphi-=4;
	    iphi=(iphi>>1);
	    iphi=iphi&3;
	    assert(iphi>=0);
	    assert(iphi<=3);
	  }
	  

	  assert(allproj_!=0);

	  unsigned int index=allproj_->nTracklets();

	  FPGATracklet* tracklet=inputproj_[j]->getFPGATracklet(i);

	  //This block of code just checks that the configuration is consistent
	  if (tracklet->minusNeighbor(layer_)) {
	    //cout << "For minus: tracklet TCID "<<tracklet<<" "<<tracklet->TCID()<<" "<<inputproj_[j]->getName()<<endl;
	    if (lastTCIDminus>=tracklet->TCID()) {
	      cout << "Wrong TCID ordering for Minus projections in "<<getName()<<" last "<<lastTCIDminus<<" "<<tracklet->TCID()<<endl;
	    } else {
	      lastTCIDminus=tracklet->TCID();
	    }
	  } else if (tracklet->plusNeighbor(layer_)) {
	    if (lastTCIDplus>=tracklet->TCID()) {
	      cout << "Wrong TCID ordering for Plus projections in "<<getName()<<" last "<<lastTCIDplus<<" "<<tracklet->TCID()<<endl;
	    } else {
	      lastTCIDplus=tracklet->TCID();
	    }
	  } else {
	    if (lastTCID>=tracklet->TCID()) {
	      cout << "Wrong TCID ordering for projections in "<<getName()<<endl;
	    } else {
	      lastTCID=tracklet->TCID();
	    }
	  }

	  
	  if (!(inputproj_[j]->getFPGATracklet(i)->minusNeighbor(layer_)||
		inputproj_[j]->getFPGATracklet(i)->plusNeighbor(layer_))){
	    //cout << "layer minusNeighbor plusNeighbor homeSector iSector :"
	    //	 <<layer_<<" "
	    //	 <<inputproj_[j]->getFPGATracklet(i)->minusNeighbor(layer_)<<" "
	    //	 <<inputproj_[j]->getFPGATracklet(i)->plusNeighbor(layer_)<<" "
	    // 	 <<inputproj_[j]->getFPGATracklet(i)->homeSector()<<" "
	    // 	 <<iSector_<<" "
	    //	 <<inputproj_[j]->getName()
	    //	 <<endl;
	    //assert(inputproj_[j]->getFPGATracklet(i)->homeSector()==iSector_);
	  }
	  allproj_->addTracklet(inputproj_[j]->getFPGATracklet(i));

	  //cout << "index iphi : "<<index<<" "<<iphi<<endl;
	  
	  if (iphi==0) {
	    assert(vmprojPHI1_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI1_->getName()<<endl;
	    }
	    vmprojPHI1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==1) {
	    assert(vmprojPHI2_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI2_->getName()<<endl;
	    }
	    vmprojPHI2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==2) {
	    assert(vmprojPHI3_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI3_->getName()<<endl;
	    }
	    vmprojPHI3_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==3) {
	    assert(vmprojPHI4_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI4_->getName()<<endl;
	    }
	    vmprojPHI4_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==4) {
	    assert(vmprojPHI5_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI5_->getName()<<endl;
	    }
	    vmprojPHI5_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==5) {
	    assert(vmprojPHI6_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI6_->getName()<<endl;
	    }
	    vmprojPHI6_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==6) {
	    assert(vmprojPHI7_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI7_->getName()<<endl;
	    }
	    vmprojPHI7_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==7) {
	    assert(vmprojPHI8_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI8_->getName()<<endl;
	    }
	    vmprojPHI8_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	}
      }
    } else {  //do the disk now
      for (unsigned int j=0;j<inputproj_.size();j++){

	for (unsigned int i=0;i<inputproj_[j]-> nTracklets();i++){
	  count++;
	  if (count>MAXPROJROUTER) continue;

	  int disk=disk_;
	  if (inputproj_[j]->getFPGATracklet(i)->t()<0.0) disk=-disk_;

	  
	  FPGAWord fpgaphi=inputproj_[j]->getFPGATracklet(i)->fpgaphiprojdisk(disk);
	  FPGAWord fpgar=inputproj_[j]->getFPGATracklet(i)->fpgarprojdisk(disk);

	  //skip if projection is out of range!
	  if (fpgar.atExtreme()) continue;
	  if (fpgaphi.atExtreme()) continue;

	  int iphitmp=fpgaphi.value();
	  int iphi=iphitmp>>(fpgaphi.nbits()-5);

	  if (hourglass) {
	    int nvm=-1;
	    int nbins=-1;
	    if (disk_!=0){
	      nvm=nvmmedisks[disk_-1]*nallstubsdisks[disk_-1];
	      nbins=nvmmedisks[disk_-1];
	    }
	    assert(nvm>0);
	    iphi=(iphi/(32/nvm))&(nbins-1);
	  } else {
	    assert(iphi>=4);
	    assert(iphi<=27);
	    iphi-=4;
	    iphi=(iphi>>1);
	    iphi=iphi&3;
	    assert(iphi>=0);
	    assert(iphi<=3);
	  }
	    
	  assert(allproj_!=0);

	  unsigned int index=allproj_->nTracklets();
	  allproj_->addTracklet(inputproj_[j]->getFPGATracklet(i));

	  FPGATracklet* tracklet=inputproj_[j]->getFPGATracklet(i);

	  int rindex=(tracklet->fpgarprojdisk(disk_).value()>>(tracklet->fpgarprojdisk(disk_).nbits()-nrbits_))&((1<<nrbits_)-1);

	  int phiderindex=(tracklet->fpgaphiprojderdisk(disk_).value()>>(tracklet->fpgaphiprojderdisk(disk_).nbits()-nphiderbits_))&((1<<nphiderbits_)-1);

	  int signindex=(tracklet->fpgarprojderdisk(disk_).value()<0);

	  int bendindex=(signindex<<(nphiderbits_+nrbits_))+
	    (rindex<<(nphiderbits_))+
	    phiderindex;
	  
	  int ibendproj=bendTable(abs(disk_)-1,bendindex);

	  tracklet->setBendIndex(ibendproj,disk_);
	  
	  if (iphi==0) {
	    assert(vmprojPHI1_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI1_->getName()<<endl;
	    }
	    vmprojPHI1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==1) {
	    assert(vmprojPHI2_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI2_->getName()<<endl;
	    }
	    vmprojPHI2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==2) {
	    assert(vmprojPHI3_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI3_->getName()<<endl;
	    }
	    vmprojPHI3_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }

	  if (iphi==3) {
	    assert(vmprojPHI4_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI4_->getName()<<endl;
	    }
	    vmprojPHI4_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }	  

	  if (iphi==4) {
	    assert(vmprojPHI5_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI5_->getName()<<endl;
	    }
	    vmprojPHI5_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }	  

	  if (iphi==5) {
	    assert(vmprojPHI6_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI6_->getName()<<endl;
	    }
	    vmprojPHI6_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }	  

	  if (iphi==6) {
	    assert(vmprojPHI7_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI7_->getName()<<endl;
	    }
	    vmprojPHI7_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }	  

	  if (iphi==7) {
	    assert(vmprojPHI8_!=0);
	    if (debug1){
	      cout << "FPGAProjectionRouter "<<getName()<<" add projection to : "<<vmprojPHI8_->getName()<<endl;
	    }
	    vmprojPHI8_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }	  
	}
      }
    }

    //cout << "Done in "<<getName()<<endl;
    
    if (writeAllProjections) {
      static ofstream out("allprojections.txt"); 
      out << getName() << " " << allproj_->nTracklets() << endl;
    } 
   

    if (writeVMProjections) {
      static ofstream out("vmprojections.txt"); 
      if (vmprojPHI1_!=0) out << vmprojPHI1_->getName() << " " << vmprojPHI1_->nTracklets() << endl;
      if (vmprojPHI2_!=0) out << vmprojPHI2_->getName() << " " << vmprojPHI2_->nTracklets() << endl;
      if (vmprojPHI3_!=0) out << vmprojPHI3_->getName() << " " << vmprojPHI3_->nTracklets() << endl;
      if (vmprojPHI4_!=0) out << vmprojPHI4_->getName() << " " << vmprojPHI4_->nTracklets() << endl;
      if (vmprojPHI5_!=0) out << vmprojPHI5_->getName() << " " << vmprojPHI5_->nTracklets() << endl;
      if (vmprojPHI6_!=0) out << vmprojPHI6_->getName() << " " << vmprojPHI6_->nTracklets() << endl;
      if (vmprojPHI7_!=0) out << vmprojPHI7_->getName() << " " << vmprojPHI7_->nTracklets() << endl;
      if (vmprojPHI8_!=0) out << vmprojPHI8_->getName() << " " << vmprojPHI8_->nTracklets() << endl;
    }
  }

  double bend(double r, double rinv) {
    
    double dr=0.18;
    
    double delta=r*dr*0.5*rinv;
    
    double bend=-delta/0.009;
    if (r<55.0) bend=-delta/0.01;
    
    return bend;
    
  }
  
  int bendTable(int diskindex,int bendindex) {

    static vector<int> bendtable[5];

    static bool first=true;

    if (first) {
      first=false;
    
      for (unsigned int idisk=0;idisk<5;idisk++) {

	unsigned int nsignbins=2;
	unsigned int nrbins=1<<(nrbits_);
	unsigned int nphiderbins=1<<(nphiderbits_);
      
	for(unsigned int isignbin=0;isignbin<nsignbins;isignbin++) {
	  for(unsigned int irbin=0;irbin<nrbins;irbin++) {
	    int ir=irbin;
	    if (ir>(1<<(nrbits_-1))) ir-=(1<<nrbits_);
	    ir=ir<<(nrbitsprojdisk-nrbits_);
	    for(unsigned int iphiderbin=0;iphiderbin<nphiderbins;iphiderbin++) {
	      int iphider=iphiderbin;
	      if (iphider>(1<<(nphiderbits_-1))) iphider-=(1<<nphiderbits_);
	      iphider=iphider<<(nbitsphiprojderL123-nphiderbits_);
	      
	      double rproj=ir*krprojshiftdisk;
	      double phider=iphider*FPGATrackletCalculator::ITC_L1L2.der_phiD_final.get_K();
	      double t=zmean[idisk]/rproj;
	      
	      if (isignbin) t=-t;
	  
	      double rinv=-phider*(2.0*t);

	      double bendproj=0.5*bend(rproj,rinv);

	    
	      int ibendproj=2.0*bendproj+15.5;
	      if (ibendproj<0) ibendproj=0;
	      if (ibendproj>31) ibendproj=31;
	      
	      bendtable[idisk].push_back(ibendproj);

	    }
	  }
	}
      }
    }

    

    return bendtable[diskindex][bendindex];

  }
  
  
private:

  int layer_; 
  int disk_; 

  int nrbits_;
  int nphiderbits_;

  
  vector<FPGATrackletProjections*> inputproj_;

  FPGAAllProjections* allproj_;
  FPGAVMProjections* vmprojPHI1_;
  FPGAVMProjections* vmprojPHI2_;
  FPGAVMProjections* vmprojPHI3_;
  FPGAVMProjections* vmprojPHI4_;
  FPGAVMProjections* vmprojPHI5_;
  FPGAVMProjections* vmprojPHI6_;
  FPGAVMProjections* vmprojPHI7_;
  FPGAVMProjections* vmprojPHI8_;


};

#endif
