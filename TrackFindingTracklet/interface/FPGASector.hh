//This class holds functional blocks of a sector

 
#ifndef FPGASECTOR_H
#define FPGASECTOR_H

#include "FPGAInputLink.hh"
#include "FPGAAllStubs.hh"
#include "FPGAVMStubsTE.hh"
#include "FPGAVMStubsME.hh"
#include "FPGAStubPairs.hh"
#include "FPGATrackletParameters.hh"
#include "FPGATrackletProjections.hh"
#include "FPGAAllProjections.hh"
#include "FPGAVMProjections.hh"
#include "FPGACandidateMatch.hh"
#include "FPGAFullMatch.hh"
#include "FPGATrackFit.hh"
#include "FPGACleanTrack.hh"

#include "FPGAVMRouter.hh"
#include "FPGAVMRouterME.hh"
#include "FPGAVMRouterTE.hh"
#include "FPGATrackletEngine.hh"
#include "FPGATrackletCalculator.hh"
#include "FPGAProjectionRouter.hh"
#include "FPGAProjectionTransceiver.hh"
#include "FPGAMatchEngine.hh"
#include "FPGAMatchCalculator.hh"
#include "FPGAMatchProcessor.hh"
#include "FPGAMatchTransceiver.hh"
#include "FPGAFitTrack.hh"
#include "FPGAPurgeDuplicate.hh"

using namespace std;

class FPGASector{

public:

  FPGASector(unsigned int i){
    isector_=i;
    double dphi=two_pi/NSector;
    double dphiHG=0.0;
    if (hourglass) {
      dphiHG=0.5*(dphisectorHG-two_pi/NSector);
    }
    phimin_=isector_*dphi-dphiHG;
    phimax_=phimin_+dphi+2*dphiHG;
    if (hourglass) {
      phimin_-=0.5*two_pi/NSector;
      phimax_-=0.5*two_pi/NSector;
    }
    if (phimin_>0.5*two_pi) phimin_-=two_pi;
    if (phimax_>0.5*two_pi) phimax_-=two_pi;
    if (phimin_>phimax_)  phimin_-=two_pi;
  }

  void addStub(L1TStub stub, string dtc) {
    
    if (hourglass) {
      double phi=stub.phi();
      //cout << "FPGASector::addStub layer phi phimin_ phimax_ : "<<stub.layer()+1<<" "<<" "<<dtc<<" "<<phi<<" "<<phimin_<<" "<<phimax_<<endl;
      double dphi=two_pi/NSector/6.0;
      if (hourglass) {
	dphi=0.5*(dphisectorHG-two_pi/NSector);
      }

      static std::map<string,std::vector<int> > ILindex;
      std::vector<int>& tmp=ILindex[dtc];
      if (tmp.size()==0){
	//cout << "Adding entries for dtc : "<<dtc;
	for (unsigned int i=0;i<IL_.size();i++){
	  if (IL_[i]->getName().find("_"+dtc)!=string::npos){
	    //cout << " Adding link "<<IL_[i]->getName()<<endl;
	    tmp.push_back(i);
	  }
	}
	//cout << endl;
      }
      
      if (((phi>phimin_-dphi)&&(phi<phimax_+dphi))||
	  ((phi>two_pi+phimin_-dphi)&&(phi<two_pi+phimax_+dphi))) {
	FPGAStub fpgastub(stub,phimin_,phimax_);
	std::vector<int>& tmp=ILindex[dtc];
	assert(tmp.size()!=0);
	for (unsigned int i=0;i<tmp.size();i++){
	  //cout << "Add stub to link"<<IL_[tmp[i]]->getName()<<endl;
	  IL_[tmp[i]]->addStub(stub,fpgastub,dtc);
	}
      }
    }  else {
      
      double phi=stub.phi();
      int layer=stub.layer()+1;
      //cout << "FPGASector::addStub phi phimin_ phimax_ : "<<phi<<" "<<phimin_<<" "<<phimax_<<endl;
      double dphi=two_pi/NSector/6.0;
      if (layer<999) {
	if (((phi>phimin_-dphi)&&(phi<phimax_+dphi))||
	    ((phi>two_pi+phimin_-dphi)&&(phi<two_pi+phimax_+dphi))) {
	  FPGAStub fpgastub(stub,phimin_,phimax_);
	  //cout << "Trying to add stub in sector : "<<isector_<<" layer = "<<layer<<endl;
	  for (unsigned int i=0;i<IL_.size();i++){
	    //cout << i<<" "<<IL_[i]->getName()<<" "<<isector_<<endl;
	    IL_[i]->addStub(stub,fpgastub);
	  }
	}
      } else {
	//int disk=stub.disk();
	if (((phi>phimin_-dphi)&&(phi<phimax_+dphi))||
	    ((phi>two_pi+phimin_-dphi)&&(phi<two_pi+phimax_+dphi))) {
	  //cout << "Trying to add stub in sector : "<<isector_<<" disk = "<<disk<<endl;
	  for (unsigned int i=0;i<IL_.size();i++){
	    FPGAStub fpgastub(stub,phimin_,phimax_);
	    IL_[i]->addStub(stub,fpgastub);
	  }      
	}
      }
    }
    
  }
    
  void addMem(string memType,string memName){
    if (memType=="InputLink:") {
      IL_.push_back(new FPGAInputLink(memName,isector_,phimin_,phimax_));
      Memories_[memName]=IL_.back();
      MemoriesV_.push_back(IL_.back());
    } else if (memType=="AllStubs:") {
      AS_.push_back(new FPGAAllStubs(memName,isector_,phimin_,phimax_));
      Memories_[memName]=AS_.back();
      MemoriesV_.push_back(AS_.back());
    } else if (memType=="VMStubsTE:") {
      VMSTE_.push_back(new FPGAVMStubsTE(memName,isector_,phimin_,phimax_));
      Memories_[memName]=VMSTE_.back();
      MemoriesV_.push_back(VMSTE_.back());
    } else if (memType=="VMStubsME:") {
      VMSME_.push_back(new FPGAVMStubsME(memName,isector_,phimin_,phimax_));
      Memories_[memName]=VMSME_.back();
      MemoriesV_.push_back(VMSME_.back());
    } else if (memType=="StubPairs:") {
      SP_.push_back(new FPGAStubPairs(memName,isector_,phimin_,phimax_));
      Memories_[memName]=SP_.back();
      MemoriesV_.push_back(SP_.back());
    } else if (memType=="TrackletParameters:") {
      TPAR_.push_back(new FPGATrackletParameters(memName,isector_,phimin_,phimax_));
      Memories_[memName]=TPAR_.back();
      MemoriesV_.push_back(TPAR_.back());
    } else if (memType=="TrackletProjections:") {
      TPROJ_.push_back(new FPGATrackletProjections(memName,isector_,phimin_,phimax_));
      Memories_[memName]=TPROJ_.back();
      MemoriesV_.push_back(TPROJ_.back());
    } else if (memType=="AllProj:") {
      AP_.push_back(new FPGAAllProjections(memName,isector_,phimin_,phimax_));
      Memories_[memName]=AP_.back();
      MemoriesV_.push_back(AP_.back());
    } else if (memType=="VMProjections:") {
      VMPROJ_.push_back(new FPGAVMProjections(memName,isector_,phimin_,phimax_));
      Memories_[memName]=VMPROJ_.back();
      MemoriesV_.push_back(VMPROJ_.back());
    } else if (memType=="CandidateMatch:") {
      CM_.push_back(new FPGACandidateMatch(memName,isector_,phimin_,phimax_));
      Memories_[memName]=CM_.back();
      MemoriesV_.push_back(CM_.back());
    } else if (memType=="FullMatch:") {
      FM_.push_back(new FPGAFullMatch(memName,isector_,phimin_,phimax_));
      Memories_[memName]=FM_.back();
      MemoriesV_.push_back(FM_.back());
    } else if (memType=="TrackFit:") {
      TF_.push_back(new FPGATrackFit(memName,isector_,phimin_,phimax_));
      Memories_[memName]=TF_.back();
      MemoriesV_.push_back(TF_.back());
    } else if (memType=="CleanTrack:") {
      CT_.push_back(new FPGACleanTrack(memName,isector_,phimin_,phimax_));
      Memories_[memName]=CT_.back();
      MemoriesV_.push_back(CT_.back());
    } else {
      cout << "Don't know of memory type: "<<memType<<endl;
      exit(0);
    }
    
  }

  void addProc(string procType,string procName){
    if (procType=="VMRouter:") {
      VMR_.push_back(new FPGAVMRouter(procName,isector_));
      Processes_[procName]=VMR_.back();
    } else if (procType=="VMRouterTE:") {
      VMRTE_.push_back(new FPGAVMRouterTE(procName,isector_));
      Processes_[procName]=VMRTE_.back();
    } else if (procType=="VMRouterME:") {
      VMRME_.push_back(new FPGAVMRouterME(procName,isector_));
      Processes_[procName]=VMRME_.back();
    } else if (procType=="TrackletEngine:") {
      TE_.push_back(new FPGATrackletEngine(procName,isector_));
      Processes_[procName]=TE_.back();
    } else if (procType=="TrackletCalculator:"||
	       procType=="TrackletDiskCalculator:") {
      TC_.push_back(new FPGATrackletCalculator(procName,isector_));
      Processes_[procName]=TC_.back();
    } else if (procType=="ProjectionRouter:") {
      PR_.push_back(new FPGAProjectionRouter(procName,isector_));
      Processes_[procName]=PR_.back();
    } else if (procType=="ProjectionTransceiver:") {
      PT_.push_back(new FPGAProjectionTransceiver(procName,isector_));
      Processes_[procName]=PT_.back();
    } else if (procType=="MatchEngine:") {
      ME_.push_back(new FPGAMatchEngine(procName,isector_));
      Processes_[procName]=ME_.back();
    } else if (procType=="MatchCalculator:"||
	       procType=="DiskMatchCalculator:") {
      MC_.push_back(new FPGAMatchCalculator(procName,isector_));
      Processes_[procName]=MC_.back();
    } else if (procType=="MatchProcessor:") {
      MP_.push_back(new FPGAMatchProcessor(procName,isector_));
      Processes_[procName]=MP_.back();
    } else if (procType=="MatchTransceiver:") {
      MT_.push_back(new FPGAMatchTransceiver(procName,isector_));
      Processes_[procName]=MT_.back();
    } else if (procType=="FitTrack:") {
      FT_.push_back(new FPGAFitTrack(procName,isector_));
      Processes_[procName]=FT_.back();
    } else if (procType=="PurgeDuplicate:") {
      PD_.push_back(new FPGAPurgeDuplicate(procName,isector_));
      Processes_[procName]=PD_.back();
    } else {
      cout << "Don't know of processing type: "<<procType<<endl;
      exit(0);      
    }
  }

  void addWire(string mem,string procinfull,string procoutfull){

    //cout << "Mem : "<<mem<<" input from "<<procinfull
    //	 << " output to "<<procoutfull<<endl;

    stringstream ss1(procinfull);
    string procin, output;
    getline(ss1,procin,'.');
    getline(ss1,output);

    stringstream ss2(procoutfull);
    string procout, input;
    getline(ss2,procout,'.');
    getline(ss2,input);

    //cout << "Procin  : "<<procin<<" "<<output<<endl;
    //cout << "Procout : "<<procout<<" "<<input<<endl;

    FPGAMemoryBase* memory=getMem(mem);

    if (procin!="") {
      FPGAProcessBase* inProc=getProc(procin);
      inProc->addOutput(memory,output);
      }

    if (procout!="") {
      FPGAProcessBase* outProc=getProc(procout);
      outProc->addInput(memory,input);
    }



  }

  FPGAProcessBase* getProc(string procName){

    map<string, FPGAProcessBase*>::iterator it=Processes_.find(procName);

    if (it!=Processes_.end()) {
      return it->second;
    }
    cout << "Could not find process with name : "<<procName<<endl;
    assert(0);
    return 0;
  }

  FPGAMemoryBase* getMem(string memName){

    map<string, FPGAMemoryBase*>::iterator it=Memories_.find(memName);

    if (it!=Memories_.end()) {
      return it->second;
    }
    cout << "Could not find memory with name : "<<memName<<endl;
    assert(0);
    return 0;
  }

  void writeInputStubs(bool first) {
    for (unsigned int i=0;i<IL_.size();i++){
      IL_[i]->writeStubs(first);
    }
  }

  void writeInputStubs_in2(bool first) {
    for (unsigned int i=0;i<IL_.size();i++){
      IL_[i]->writeInputStubs(first, writestubs_in2,padding);
    }
  }

  void writeVMSTE(bool first) {
    for (unsigned int i=0;i<VMSTE_.size();i++){
      VMSTE_[i]->writeStubs(first);
    }
  }
  
  void writeVMSME(bool first) {
    for (unsigned int i=0;i<VMSME_.size();i++){
      VMSME_[i]->writeStubs(first);
    }
  }

  void writeAS(bool first) {
    for (unsigned int i=0;i<AS_.size();i++){
      AS_[i]->writeStubs(first);
    }
  }

  void writeSP(bool first) {
    for (unsigned int i=0;i<SP_.size();i++){
      SP_[i]->writeSP(first);
    }
  }

  void writeTPAR(bool first) {
    for (unsigned int i=0;i<TPAR_.size();i++){
      TPAR_[i]->writeTPAR(first);
    }
  }

  void writeTPROJ(bool first) {
    for (unsigned int i=0;i<TPROJ_.size();i++){
      TPROJ_[i]->writeTPROJ(first);
    }
  }

  void writeAP(bool first) {
    for (unsigned int i=0;i<AP_.size();i++){
      AP_[i]->writeAP(first);
    }
  }

  void writeVMPROJ(bool first) {
    for (unsigned int i=0;i<VMPROJ_.size();i++){
      VMPROJ_[i]->writeVMPROJ(first);
    }
  }

  void writeCM(bool first) {
    for (unsigned int i=0;i<CM_.size();i++){
     CM_[i]->writeCM(first);
    }
  }

  void writeMC(bool first) {
    for (unsigned int i=0;i<FM_.size();i++){
      FM_[i]->writeMC(first);
    }
  }

  void writeTF(bool first){
    for(unsigned int i=0; i<TF_.size(); ++i){
      TF_[i]->writeTF(first);
    }
  }

  void writeCT(bool first) {
    for(unsigned int i=0; i<CT_.size(); ++i){
      CT_[i]->writeCT(first);
    }
  }

  void clean() {
    if (writeNMatches) {
      int matchesL1=0;
      int matchesL3=0;
      int matchesL5=0;
      for(unsigned int i=0;i<TPAR_.size();i++) {
	TPAR_[i]->writeMatches(matchesL1,matchesL3,matchesL5);
      }
      static ofstream out("nmatchessector.txt");
      out <<matchesL1<<" "<<matchesL3<<" "<<matchesL5<<endl;
    }
    
    
    for(unsigned int i=0;i<MemoriesV_.size();i++) {
      MemoriesV_[i]->clean();
    }
  }

  void executeVMR(){

    if (writeIL) {
      static ofstream out("inputlink.txt");
      for (unsigned int i=0;i<IL_.size();i++){
	out<<IL_[i]->getName()<<" "<<IL_[i]->nStubs()<<endl;
      } 
    }
    
    for (unsigned int i=0;i<VMR_.size();i++){
      VMR_[i]->execute();
    }
    for (unsigned int i=0;i<VMRTE_.size();i++){
      VMRTE_[i]->execute();
    }
    for (unsigned int i=0;i<VMRME_.size();i++){
      VMRME_[i]->execute();
    }
  }

  void executeTE(){
    for (unsigned int i=0;i<TE_.size();i++){
      TE_[i]->execute();
    }
  }

  void executeTC(){
    for (unsigned int i=0;i<TC_.size();i++){
      TC_[i]->execute();
    }
  }

  void executePR(){
    for (unsigned int i=0;i<PR_.size();i++){
      PR_[i]->execute();
    }
  }

  void executeME(){
    for (unsigned int i=0;i<ME_.size();i++){
      ME_[i]->execute();
    }
  }

  void executeMC(){
    for (unsigned int i=0;i<MC_.size();i++){
      MC_[i]->execute();
    }
  }

  void executeMP(){
    for (unsigned int i=0;i<MP_.size();i++){
      MP_[i]->execute();
    }
  }

  void executeFT(){
    fpgatracks_.clear();
    for (unsigned int i=0;i<FT_.size();i++){
      FT_[i]->execute(fpgatracks_);
    }
  }

  void executePD(std::vector<FPGATrack*>& tracks){
    for (unsigned int i=0;i<PD_.size();i++){
      PD_[i]->execute(tracks);
    }
  }

  void executePT(FPGASector* sectorPlus,FPGASector* sectorMinus){
    //For now the order is assumed
    for (unsigned int i=0;i<PT_.size();i++){
      string name=PT_[i]->getName();
      //cout << "FPGASector:executePT "<<name<<endl;
      //cout << "name.find(\"Minus\") : "<<name.find("Minus")<<endl;
      if (name.find("Minus")!=std::string::npos) {
	name.replace(name.find("Minus"),5,"Plus");
	//cout << "New name : "<<name<<endl;
	for (unsigned int j=0;j<sectorMinus->PT_.size();j++){
	  if (sectorMinus->PT_[j]->getName()==name) {
	    PT_[i]->execute(sectorMinus->PT_[j]);
	  }
	}
      } else if (name.find("Plus")!=std::string::npos) {
	name.replace(name.find("Plus"),4,"Minus");
	//cout << "New name : "<<name<<endl;
	for (unsigned int j=0;j<sectorPlus->PT_.size();j++){
	  if (sectorPlus->PT_[j]->getName()==name) {
	    PT_[i]->execute(sectorPlus->PT_[j]);
	  }
	}
      } else {
	assert(0);
      }
      
    }

    if (writeTrackProjOcc) {
      static ofstream out("trackprojocc.txt");
      for (unsigned int i=0; i<TPROJ_.size();i++){
	out << TPROJ_[i]->getName()<<" "<<TPROJ_[i]->nTracklets()<<endl;
      }
    }
    

  }


  void executeMT(FPGASector* sectorPlus,FPGASector* sectorMinus){
    //For now the order is assumed
    for (unsigned int i=0;i<MT_.size();i++){
      string name=MT_[i]->getName();
      //cout << "FPGASector:executePT "<<name<<endl;
      if (name.find("Minus")!=std::string::npos) {
	name.replace(6,5,"Plus");
	//cout << "New name : "<<name<<endl;
	for (unsigned int j=0;j<sectorMinus->MT_.size();j++){
	  if (sectorMinus->MT_[j]->getName()==name) {
	    MT_[i]->execute(sectorMinus->MT_[j]);
	  }
	}
      } else if (name.find("Plus")!=std::string::npos) {
	name.replace(6,4,"Minus");
	//cout << "New name : "<<name<<endl;
	for (unsigned int j=0;j<sectorPlus->MT_.size();j++){
	  if (sectorPlus->MT_[j]->getName()==name) {
	    MT_[i]->execute(sectorPlus->MT_[j]);
	  }
	}
      } else {
	assert(0);
      }
    }
  }


  bool foundTrack(ofstream& outres, L1SimTrack simtrk){
    bool match=false;
    for (unsigned int i=0;i<TF_.size();i++){
      if (TF_[i]->foundTrack(outres,simtrk)) match=true;
    }
    return match;
  }

  std::vector<FPGATracklet*> getAllTracklets() {
    std::vector<FPGATracklet*> tmp;
    for(unsigned int i=0;i<TPAR_.size();i++){
      for(unsigned int j=0;j<TPAR_[i]->nTracklets();j++){
	tmp.push_back(TPAR_[i]->getFPGATracklet(j));
      }
    }
    return tmp;
  }

  std::vector<std::pair<FPGAStub*,L1TStub*> > getStubs() const {
    std::vector<std::pair<FPGAStub*,L1TStub*> > tmp;

    for(unsigned int imem=0;imem<IL_.size();imem++){
      for(unsigned int istub=0;istub<IL_[imem]->nStubs();istub++){
	tmp.push_back(IL_[imem]->getStub(istub));
      }
    }
    
    return tmp;
  }
    

  double phimin() const {return phimin_;}
  double phimax() const {return phimax_;}

private:

  int isector_;
  double phimin_;
  double phimax_;

  std::vector<FPGATrack*> fpgatracks_;


  std::map<string, FPGAMemoryBase*> Memories_;
  std::vector<FPGAMemoryBase*> MemoriesV_;
  std::vector<FPGAInputLink*> IL_;
  std::vector<FPGAAllStubs*> AS_;
  std::vector<FPGAVMStubsTE*> VMSTE_;
  std::vector<FPGAVMStubsME*> VMSME_;
  std::vector<FPGAStubPairs*> SP_;
  std::vector<FPGATrackletParameters*> TPAR_;
  std::vector<FPGATrackletProjections*> TPROJ_;
  std::vector<FPGAAllProjections*> AP_;
  std::vector<FPGAVMProjections*> VMPROJ_;
  std::vector<FPGACandidateMatch*> CM_;
  std::vector<FPGAFullMatch*> FM_;
  std::vector<FPGATrackFit*> TF_;
  std::vector<FPGACleanTrack*> CT_;
  
  std::map<string, FPGAProcessBase*> Processes_;
  std::vector<FPGAVMRouter*> VMR_;
  std::vector<FPGAVMRouterTE*> VMRTE_;
  std::vector<FPGAVMRouterME*> VMRME_;
  std::vector<FPGATrackletEngine*> TE_;
  std::vector<FPGATrackletCalculator*> TC_;
  std::vector<FPGAProjectionRouter*> PR_;
  std::vector<FPGAProjectionTransceiver*> PT_;
  std::vector<FPGAMatchEngine*> ME_;
  std::vector<FPGAMatchCalculator*> MC_;
  std::vector<FPGAMatchProcessor*> MP_;
  std::vector<FPGAMatchTransceiver*> MT_;
  std::vector<FPGAFitTrack*> FT_;
  std::vector<FPGAPurgeDuplicate*> PD_;


  
};

#endif
