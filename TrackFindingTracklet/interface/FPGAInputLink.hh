// This class holds a list of stubs for an input link.
//This modules 'owns' the pointers to the stubs. All
//subsequent modules that handles stubs uses a pointer
//to the original stored here

 
#ifndef FPGAINPUTLINK_H
#define FPGAINPUTLINK_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGAMemoryBase.hh"
#include "FPGAVMRouterPhiCorrTable.hh"
#include <math.h>
#include <sstream>
#include <ctype.h>

using namespace std;

class FPGAInputLink:public FPGAMemoryBase{

public:

  FPGAInputLink(string name, unsigned int iSector, 
		double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;

    indexphi_ = {0,0,0,0,0};

      string subname=name.substr(5,7);

      phiregion_=subname[3]-'A';
      phiregionoverlap_=-1;

      if (phiregion_>8) {
	if (subname[3]=='X') phiregionoverlap_=0;
	if (subname[3]=='Y') phiregionoverlap_=1;
	if (subname[3]=='Z') phiregionoverlap_=2;
	if (subname[3]=='W') phiregionoverlap_=3;
	if (subname[3]=='Q') phiregionoverlap_=4;
	if (subname[3]=='R') phiregionoverlap_=5;
	if (subname[3]=='S') phiregionoverlap_=6;
	if (subname[3]=='T') phiregionoverlap_=7;
	if (phiregionoverlap_!=-1) phiregion_=-1;
      }

      string subnamelayer=getName().substr(3,2);
      layerdisk_=subnamelayer[1]-'0';

      isLayer_=subnamelayer[0]=='L';
      isDisk_=subnamelayer[0]=='D';

    
  }

  void addStub(L1TStub& al1stub, FPGAStub& stub, string dtc="") {

    static bool first=true;
    static FPGAVMRouterPhiCorrTable phiCorrLayers[6];

    if (first) {
      for (int l=0;l<6;l++){
	int nbits=3;
	if (l>=3) nbits=4;
        phiCorrLayers[l].init(l+1,nbits,3);
      }
      first=false;
    }
    
    //cout << getName()<<" addStub "<<stub.layer().value()+1<<" "<<al1stub.phi()<<" "<<al1stub.z()<<endl;
    
    if (stub.layer().value()!=-1) {

      FPGAWord r=stub.r();
      FPGAWord bend=stub.bend();
      int bendbin=bend.value();
      //cout << "bend:"<<bend.value()<<" "<<bend.nbits()<<" "<<bendbin<<endl;
      int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-3);
      
      int iphicorr=phiCorrLayers[stub.layer().value()].getphiCorrValue(bendbin,rbin);

      stub.setPhiCorr(iphicorr);

    }
    
    bool add=false;
    unsigned int asindex = 0;
    int iphivmRaw=stub.iphivmRaw();
    
    
    if (hourglass) {

      //cout << getName()<<" layer = "<<stub.layer().value()+1<<endl;
      
      if (stub.layer().value()==-1 && isLayer() ) return; 
      if (stub.layer().value()!=-1 && isDisk() ) return;
      if (stub.layer().value()!=-1){
	if (stub.layer().value()+1!=layerdisk()) return;
      } else {
	if (abs(stub.disk().value())!=layerdisk()) return;
      }

      int phibin=-1;
      if (stub.layer().value()!=-1) {
	if (phiregionoverlap()==-1) {
	  phibin=iphivmRaw/(32/nallstubslayers[layerdisk()-1]);
	} else {
	  phibin=iphivmRaw/(32/nallstubsoverlaplayers[layerdisk()-1]);
	}
      } else {
	phibin=iphivmRaw/(32/nallstubsdisks[layerdisk()-1]); 
      }

      //cout << getName()<<" "<<iSector_<<" phibin phiregion : "<<phibin<<" "<<phiregion()<<endl;
      
      if (phibin!=phiregion() && phibin!=phiregionoverlap()) return;


      
 
      //cout << getName()<<" "<<getName().substr(10,dtc.size())<<" "<<dtc<<endl;
      
      if (getName().substr(10,dtc.size())!=dtc) return;

      string half=getName().substr(getName().size()-3,3);
      if (half[1]!='n') {
	half=getName().substr(getName().size()-1,1);
      }
      //cout << "iSector half iphivmRaw: "<<iSector_<<" "<<half<<" "<<iphivmRaw<<endl;
      assert(half[0]=='A' || half[0]=='B');


      double sectorphi=al1stub.phi();
      if (sectorphi<0.0) sectorphi+=two_pi;

      while (sectorphi>3*two_pi/(2*NSector)) {
	sectorphi-=(two_pi/NSector);
      }
      
      
      if (half[0]=='B' && iphivmRaw<=15) return;
      if (half[0]=='A' && iphivmRaw>15) return;
      

      
      //cout << getName() << " adding stub to phibin : "<<phibin<<" iphivmRaw "<<iphivmRaw<<endl;
      add=true;


    } else {
      if (stub.layer().value()!=-1) {
	string subname=getName().substr(5,7);
	string subnamelayer=getName().substr(3,2);
	
	int layer=stub.layer().value()+1;

	if (!((layer==1&&subnamelayer=="L1")||
	      (layer==2&&subnamelayer=="L2")||
	      (layer==3&&subnamelayer=="L3")||
	      (layer==4&&subnamelayer=="L4")||
	      (layer==5&&subnamelayer=="L5")||
	      (layer==6&&subnamelayer=="L6"))){
	  return;
	}
	
	if ((subname.substr(5,2)=="ZP"&&stub.z().value()>0)||
	    (subname.substr(5,2)=="ZM"&&stub.z().value()<=0)) {
	  if (stub.phiregion().value()==0) asindex = indexphi_[0]++;
	  else if (stub.phiregion().value()==1) asindex = indexphi_[1]++;
	  else if (stub.phiregion().value()==2) asindex = indexphi_[2]++;
	  else if (stub.phiregion().value()==3) asindex = indexphi_[3]++;
	  else {
	    // Stub iphiRaw<4 or iphiRaw>27. It is used in TE but not ME, and is not stored in the AS memory filled by VMRouterME. Use an extra index to keep track of them.
	    asindex = indexphi_[4]++;
	  }
	}

	if (layer==2&&subnamelayer=="L2"&&(subname=="PHIW_ZP"||subname=="PHIQ_ZP")) {
	  if (al1stub.z()>25.0) {
	  if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIW_ZP") add=true;     //overlap
	  if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIQ_ZP") add=true;    //overlap
	  }
	}
	if (layer==2&&subnamelayer=="L2"&&(subname=="PHIW_ZM"||subname=="PHIQ_ZM")) {
	  if (al1stub.z()<-25.0) {
	    if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIW_ZM") add=true;     //overlap
	  if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIQ_ZM") add=true;    //overlap
	  }
	}
	
	//Special case --- L3 is grouped with D1 for overlap seeding
	
	//if (layer==3&&subnamelayer=="L3"&&(subname=="PHIW_ZP"||subname=="PHIQ_ZP")) {
	//	if (al1stub.z()>40.0) {
	//	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW_ZP") add=true;     //overlap
	//	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ_ZP") add=true;    //overlap
	//	  }
	//}
	
	//if (layer==3&&subnamelayer=="L3"&&(subname=="PHIW_ZM"||subname=="PHIQ_ZM")) {
	//	if (al1stub.z()<-40.0) {
	//	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW_ZM") add=true;     //overlap
	//	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ_ZM") add=true;    //overlap
	//	  }
	//}
	
	
	if (!((layer==1&&subnamelayer=="L1")||
	      (layer==2&&subnamelayer=="L2")||
	      (layer==3&&subnamelayer=="L3")||
	      (layer==4&&subnamelayer=="L4")||
	      (layer==5&&subnamelayer=="L5")||
	      (layer==6&&subnamelayer=="L6"))){
	  return;
	}
      
	//cout << "Stub candidate in "<<getName()<<" "<<subnamelayer<<" "<<subname<<" "<<iphivmRaw<<" "
	//	 <<al1stub.phi()<<" "
	// 	 <<al1stub.z()<<endl;
	
	int ASPHI = stub.phiregion().value()+1;  // AllStub phi region
	
	if (subnamelayer=="L1"||subnamelayer=="L3"||subnamelayer=="L5"){
	  if (stub.z().value()>0) {
	    if (ASPHI==1 && subname=="PHI1_ZP") add=true;
	    if (ASPHI==2 && subname=="PHI2_ZP") add=true;
	    if (ASPHI==3 && subname=="PHI3_ZP") add=true;
	    /*
	      if (subnamelayer=="L1") {
	      if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZP") add=true;
	      if (iphivmRaw>=8 && iphivmRaw<=11 && subname=="PHI2_ZP") add=true;
	      if (iphivmRaw>=12 && iphivmRaw<=15 && subname=="PHI3_ZP") add=true;
	      if (iphivmRaw>=16 && iphivmRaw<=19 && subname=="PHI4_ZP") add=true;
	      if (iphivmRaw>=20 && iphivmRaw<=23 && subname=="PHI5_ZP") add=true;
	      if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI6_ZP") add=true;
	      } else {
	      if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1_ZP") add=true;
	      if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2_ZP") add=true;
	      if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3_ZP") add=true;
	      }
	    */
	    //these are for TE
	    if ((subnamelayer=="L1"&&fabs(al1stub.z())<87.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	      if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHIA_ZP") add=true;
	      if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHIB_ZP") add=true;
	      if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHIC_ZP") add=true;
	    }
	    if ((subnamelayer=="L1"&&fabs(al1stub.z())>80.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	      if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIX_ZP") add=true;     //overlap
	      if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIY_ZP") add=true;    //overlap
	    }
	  } else {
	    if (ASPHI==1 && subname=="PHI1_ZM") add=true;
	    if (ASPHI==2 && subname=="PHI2_ZM") add=true;
	    if (ASPHI==3 && subname=="PHI3_ZM") add=true;
	    /*
	      if (subnamelayer=="L1") {
	      if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZM") add=true;
	      if (iphivmRaw>=8 && iphivmRaw<=11 && subname=="PHI2_ZM") add=true;
	      if (iphivmRaw>=12 && iphivmRaw<=15 && subname=="PHI3_ZM") add=true;
	      if (iphivmRaw>=16 && iphivmRaw<=19 && subname=="PHI4_ZM") add=true;
	      if (iphivmRaw>=20 && iphivmRaw<=23 && subname=="PHI5_ZM") add=true;
	      if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI6_ZM") add=true;
	      } else {
	      if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1_ZM") add=true;
	      if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2_ZM") add=true;
	      if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3_ZM") add=true;
	      }
	    */
	    //these are for TE
	    if ((subnamelayer=="L1"&&fabs(al1stub.z())<87.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	      if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHIA_ZM") add=true;
	      if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHIB_ZM") add=true;
	      if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHIC_ZM") add=true;
	    }
	    if ((subnamelayer=="L1"&&fabs(al1stub.z())>80.0)||subnamelayer=="L3"||subnamelayer=="L5"){
	      if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIX_ZM") add=true;   //overlap
	      if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIY_ZM") add=true;	//overlap  
	    }
	  }
	}

      
	if (subnamelayer=="L2"||subnamelayer=="L4"||subnamelayer=="L6"){
	  if (stub.z().value()>0) {
	    //remember that these are for ME
	    if (ASPHI==1 && subname=="PHI1_ZP") add=true;
	    if (ASPHI==2 && subname=="PHI2_ZP") add=true;
	    if (ASPHI==3 && subname=="PHI3_ZP") add=true;
	    if (ASPHI==4 && subname=="PHI4_ZP") add=true;
	    /*
	      if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZP") add=true;
	      if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHI2_ZP") add=true;
	      if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHI3_ZP") add=true;
	      if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI4_ZP") add=true;
	    */
	    //these are for TE
	    if (iphivmRaw>=0 && iphivmRaw<=7 && subname=="PHIA_ZP") add=true;
	    if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHIB_ZP") add=true;
	    if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHIC_ZP") add=true;
	    if (iphivmRaw>=24 && iphivmRaw<=31 && subname=="PHID_ZP") add=true;
	  } else {
	    //remember that these are for ME
	    if (ASPHI==1 && subname=="PHI1_ZM") add=true;
	    if (ASPHI==2 && subname=="PHI2_ZM") add=true;
	    if (ASPHI==3 && subname=="PHI3_ZM") add=true;
	    if (ASPHI==4 && subname=="PHI4_ZM") add=true;
	    /*
	      if (iphivmRaw>=4 && iphivmRaw<=7 && subname=="PHI1_ZM") add=true;
	      if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHI2_ZM") add=true;
	      if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHI3_ZM") add=true;
	      if (iphivmRaw>=24 && iphivmRaw<=27 && subname=="PHI4_ZM") add=true;
	    */
	    //these are for TE
	    if (iphivmRaw>=0 && iphivmRaw<=7 && subname=="PHIA_ZM") add=true;
	    if (iphivmRaw>=8 && iphivmRaw<=15 && subname=="PHIB_ZM") add=true;
	    if (iphivmRaw>=16 && iphivmRaw<=23 && subname=="PHIC_ZM") add=true;
	    if (iphivmRaw>=24 && iphivmRaw<=31 && subname=="PHID_ZM") add=true;
	  }
	}
      }
      else if (stub.disk().value()!=0) {
	//cout << "FPGAInputLink::addStub in disk "<<stub.disk().value()<<endl;
	
	string subname=getName().substr(5,4);
	string subnamelayer=getName().substr(3,2);
      
      int disk=stub.disk().value();

      if (!((disk==1&&subnamelayer=="F1")||
	    (disk==2&&subnamelayer=="F2")||
	    (disk==3&&subnamelayer=="F3")||
	    (disk==4&&subnamelayer=="F4")||
	    (disk==5&&subnamelayer=="F5")||
	    (disk==-1&&subnamelayer=="B1")||
	    (disk==-2&&subnamelayer=="B2")||
	    (disk==-3&&subnamelayer=="B3")||
	    (disk==-4&&subnamelayer=="B4")||
	    (disk==-5&&subnamelayer=="B5"))){
	return;
      }    

      if (stub.phiregion().value()==0) asindex = indexphi_[0]++;
      else if (stub.phiregion().value()==1) asindex = indexphi_[1]++;
      else if (stub.phiregion().value()==2) asindex = indexphi_[2]++;
      else if (stub.phiregion().value()==3) asindex = indexphi_[3]++;
      else {
        // Stub iphiRaw<4 or iphiRaw>27. It is used in TE but not ME, and is not stored in the AS memory filled by VMRouterME. Use an extra index to keep track of them.
        asindex = indexphi_[4]++;
      }     
      
      if (abs(disk)==1&&(subnamelayer=="F1"||subnamelayer=="B1")&&(subname=="PHIW"||subname=="PHIQ")) {
	if (al1stub.isPSmodule()&&al1stub.r()>40.0) {
	  if (al1stub.z()>0.0&&subnamelayer=="F1") {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW") add=true;     //overlap
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ") add=true;    //overlap
	  }
	  if (al1stub.z()<0.0&&subnamelayer=="B1") {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIW") add=true;     //overlap
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIQ") add=true;    //overlap
	  }
	}
      }
      
      //cout << "Stub candidate in "<<getName()<<" "<<subnamelayer<<" "<<subname<<" "<<iphivmRaw<<" "
      //	 <<al1stub.phi()<<" "
      // 	 <<al1stub.z()<<endl;

      int ASPHI = stub.phiregion().value()+1;  // AllStub phi region
      
      if (subnamelayer=="F1"||subnamelayer=="F3"||subnamelayer=="F5"||
	  subnamelayer=="B1"||subnamelayer=="B3"||subnamelayer=="B5"){
	  if (ASPHI==1 && subname=="PHI1") add=true;
	  if (ASPHI==2 && subname=="PHI2") add=true;
	  if (ASPHI==3 && subname=="PHI3") add=true;
	  /*
	    if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3") add=true;
	  */
	  //these are for TE
	  if ((subnamelayer=="F1"||subnamelayer=="B1")&&(subname=="PHIA"||subname=="PHIB"||subname=="PHIC"||subname=="PHID")){
	    if (iphivmRaw>=4 && iphivmRaw<=17 && subname=="PHIA") add=true;
	    if (iphivmRaw>=4 && iphivmRaw<=17 && subname=="PHIB") add=true;
	    if (iphivmRaw>=14 && iphivmRaw<=27 && subname=="PHIC") add=true;
	    if (iphivmRaw>=14 && iphivmRaw<=27 && subname=="PHID") add=true;
	  }else{
	    if (iphivmRaw>=4 && iphivmRaw<=15 && subname=="PHIA") add=true;
	    if (iphivmRaw>=16 && iphivmRaw<=27 && subname=="PHIB") add=true;
	    if (al1stub.r()<47.0) {
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIX") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIY") add=true;
	    }
	  }
	}
	
	
	if (subnamelayer=="F2"||subnamelayer=="F4"||
	    subnamelayer=="B2"||subnamelayer=="B4"){
	  //remember that these are for ME
	  if (ASPHI==1 && subname=="PHI1") add=true;
	  if (ASPHI==2 && subname=="PHI2") add=true;
	  if (ASPHI==3 && subname=="PHI3") add=true;
	  /*
	    if (iphivmRaw>=4 && iphivmRaw<=11 && subname=="PHI1") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=19 && subname=="PHI2") add=true;
	    if (iphivmRaw>=20 && iphivmRaw<=27 && subname=="PHI3") add=true;
	  */
	  //these are for TE
	  if ((subnamelayer=="F2"||subnamelayer=="B2")&&(subname=="PHIA"||subname=="PHIB"||subname=="PHIC"||subname=="PHID")){
	    if (iphivmRaw>=0 && iphivmRaw<=13 && subname=="PHIA") add=true;
	    if (iphivmRaw>=4 && iphivmRaw<=19 && subname=="PHIB") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=27 && subname=="PHIC") add=true;
	    if (iphivmRaw>=18 && iphivmRaw<=31 && subname=="PHID") add=true;
	  }else{
	    if (iphivmRaw>=0 && iphivmRaw<=19 && subname=="PHIA") add=true;
	    if (iphivmRaw>=12 && iphivmRaw<=31 && subname=="PHIB") add=true;
	  }
	}
      }
    }
      
    if (!add) {
      //cout << "Will not add stub" << endl;
      return;
    }
    if (debug1) {
      cout << "Will add stub in "<<getName()<<" phimin_ phimax_ "<<phimin_<<" "<<phimax_<<" "<<"iphiwmRaw = "<<iphivmRaw<<" phi="<<al1stub.phi()<<" z="<<al1stub.z()<<" r="<<al1stub.r()<<endl;
    }
    if (stubs_.size()<MAXSTUBSLINK) {
      L1TStub* l1stub=new L1TStub(al1stub);
      //FPGAStub* stub=new FPGAStub(*l1stub,phimin_,phimax_);
      FPGAStub* stubptr=new FPGAStub(stub);

      // set stub index only for those sent to VMRouterTE
      if (isalpha(getName()[8])) {
        l1stub->setAllStubIndex(asindex);
        stubptr->setAllStubIndex(asindex);
      }
      
      std::pair<FPGAStub*,L1TStub*> tmp(stubptr,l1stub);
      stubs_.push_back(tmp);
    }
  }

  unsigned int nStubs() const {return stubs_.size();}

  vector<unsigned int> getASPhiIndices() const {return indexphi_;}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  void writeStubs(bool first)
  {
    string fname="MemPrints/InputStubs/InputStubs_";
    fname+=getName();
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    fname+=".dat";
    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    out_ << "BX = "<<(bitset<3>)bx_ << " Event : " << event_ << endl;

    for (unsigned int j=0;j<stubs_.size();j++){
      assert(stubs_[j].first->isBarrel() or stubs_[j].first->isDisk());
      string stub = stubs_[j].first->isBarrel() ? stubs_[j].first->str()
        : stubs_[j].first->strdisk();
      if (j<16) out_ <<"0";
      out_ << hex << j << dec;
      out_ << " " << stub << endl;
    }
    out_.close();

    bx_++;
    event_++;
    if (bx_>7) bx_=0;
  }
	
  void writeInputStubs(bool first, bool w2, bool padded) {
    
    //Barrel
    std::string fname="MemPrints/InputStubs/InputStubs_";
    fname+=getName();
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    if(w2)
      fname+="_in2.dat";
    else
      fname+=".dat";

    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    int nlay[6];
    for (unsigned int i=0;i<6;i++) {
      nlay[i]=0;
    }

    unsigned int maxstubslink=0;
    if (TMUX==4) maxstubslink=21;
    else if (TMUX==6) maxstubslink=33;
    else if (TMUX==8) maxstubslink=45;
    else {
      cout << "ERROR! Only TMUX=4/6/8 are supported! Exiting..." << endl;
      return;
    }

    for (unsigned int i=0;i<stubs_.size();i++){
      if(stubs_[i].first->isBarrel()){
	int lay=stubs_[i].first->layer().value();
	assert(lay>=0);
	assert(lay<6);
	nlay[lay]++;
      }
    }
    unsigned long int nlay2[6];
    nlay2[0] = nlay[0];
    for(unsigned int i=1; i<6; ++i) {
      nlay2[i] = nlay2[i-1]+nlay[i];
      if (nlay2[i]>maxstubslink) {
	nlay2[i]=maxstubslink;
      }
      assert(nlay2[i]<=maxstubslink);
    }

    string marker="";

    if(!w2) {
      //Write out the header
      out_ << "11111111111"<<marker;  
      out_ << (bitset<3>)bx_; //Dummy BX
      for (unsigned int i=0;i<3;i++) {
	out_ <<marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      out_ << "00000000000"<<marker<<(bitset<3>)bx_;  
      for (unsigned int i=3;i<6;i++) {
	out_ << marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      
      //now write the stubs in layer order
      string leftover="";

      
      for (unsigned int i=0;i<6;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isBarrel()){
	    unsigned int lay=stubs_[j].first->layer().value();
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      string stub=leftover+tmp.substr(0,32-leftover.size());
	      leftover=tmp.substr(32-leftover.size(),leftover.size()+4);
	      out_ << stub << endl;
	      if (leftover.size()==32){
		out_ << leftover << endl;
		leftover="";
	      }
	    }
	  }
	}
      }
      if (leftover.size()!=0) {
	string tmp="00000000000000000000000000000000";
	out_ << leftover << tmp.substr(0,32-leftover.size())<<endl;
      }
    }
    else{ 
      //format in2 for IPBUS
      unsigned int xline = 10;
      unsigned long int xword;
      //Header
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25)|(0x1ffffff);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //Stub counts
      xline = (xline+1)&15;
      xword = (nlay2[5])|(nlay2[4]<<6)|(nlay2[3]<<12)|(nlay2[2]<<18)
	|(nlay2[1]<<24)|(nlay2[0]<<30); 
      out_<< std::hex << xline<<"000";
      out_.fill('0');
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //the stubs      
      int scount = 0;
      for (unsigned int i=0;i<6;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isBarrel()){
	    unsigned int lay=stubs_[j].first->layer().value();
	    if (lay==i) {
	      scount++;
	      if (scount>(int)maxstubslink) continue;
	      string tmp=stubs_[j].first->strbare();
	      bitset<36> btmp(tmp);
	      xword = btmp.to_ulong();
	      xline = (xline+1)&15;
	      out_<< std::hex << xline<<"000";
	      out_.width(4);
	      out_<<(xword>>20);
	      out_<< " "<< xline<<"00";
	      out_.width(5);
	      out_<<(xword & 0xfffff);
	      out_<<" 51000003 51000007\n";
	    }
	  }
	}
      }
      if (padded) {
    for (int k=0;k<(int)maxstubslink-scount;k++) {
	  xline = (xline+1)&15;
	  out_<< std::hex << xline<<"0000000 "<< std::hex << xline<<"0000001";
	  out_<<" 51000003 51000007\n";
	}
      }
      //Trailer
      xline = (xline+1)&15;
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
    }    
  
    out_.close();
      
    //Now the disks

    //forward

    fname="MemPrints/InputStubs/InputStubsDiskF_";
    fname+=getName();
    fname+="_";
    ostringstream ossF;
    ossF << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=ossF.str();
    if(w2)
      fname+="_in2.dat";
    else
      fname+=".dat";

    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    for (unsigned int i=0;i<5;i++) {
      nlay[i]=0;
      nlay2[i]=0;
    }

    for (unsigned int i=0;i<stubs_.size();i++){
      if(stubs_[i].first->isDisk()&&stubs_[i].first->disk().value()>0){
	int lay=stubs_[i].first->disk().value()-1;
	assert(lay>=0);
	assert(lay<5);
	nlay[lay]++;
      }
    }
    nlay[5] = 0;
    nlay2[5] = 0;
    nlay2[0] = nlay[0];
    for(unsigned int i=1; i<5; ++i)
      nlay2[i] = nlay2[i-1]+nlay[i];

    if(!w2) {
      //Write out the header
      out_ << "11111111111"<<marker;  
      out_ << (bitset<3>)bx_; //Dummy BX
      for (unsigned int i=0;i<3;i++) {
	out_ <<marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      out_ << "00000000000"<<marker<<(bitset<3>)bx_;  
      for (unsigned int i=3;i<6;i++) {
	out_ << marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      
      //now write the stubs in disk order
      string leftover="";
      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()>0){
	    unsigned int lay=stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      string stub=leftover+tmp.substr(0,32-leftover.size());
	      leftover=tmp.substr(32-leftover.size(),leftover.size()+4);
	      out_ << stub << endl;
	      if (leftover.size()==32){
		out_ << leftover << endl;
		leftover="";
	      }
	    }
	  }
	}
      }
      if (leftover.size()!=0) {
	string tmp="00000000000000000000000000000000";
	out_ << leftover << tmp.substr(0,32-leftover.size())<<endl;
      }
    }
    else{ 
      //format in2 for IPBUS
      unsigned int xline = 10;
      unsigned long int xword;
      //Header
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25)|(0x1ffffff);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //Stub counts
      xline = (xline+1)&15;
      xword = (nlay2[5])|(nlay2[4]<<6)|(nlay2[3]<<12)|(nlay2[2]<<18)
	|(nlay2[1]<<24)|(nlay2[0]<<30); 
      out_<< std::hex << xline<<"000";
      out_.fill('0');
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //the stubs      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()>0){
	    unsigned int lay=stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      bitset<36> btmp(tmp);
	      xword = btmp.to_ulong();
	      xline = (xline+1)&15;
	      out_<< std::hex << xline<<"000";
	      out_.width(4);
	      out_<<(xword>>20);
	      out_<< " "<< xline<<"00";
	      out_.width(5);
	      out_<<(xword & 0xfffff);
	      out_<<" 51000003 51000007\n";
	    }
	  }
	}
      }

      //Trailer
      xline = (xline+1)&15;
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";

    }    
  
    out_.close();

    //back

    fname="MemPrints/InputStubs/InputStubsDiskB_";
    fname+=getName();
    fname+="_";
    ostringstream ossB;
    ossB << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=ossB.str();
    if(w2)
      fname+="_in2.dat";
    else
      fname+=".dat";

    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else 
      out_.open(fname.c_str(),std::ofstream::app);

    for (unsigned int i=0;i<5;i++) {
      nlay[i]=0;
      nlay2[i]=0;
    }

    for (unsigned int i=0;i<stubs_.size();i++){
      if(stubs_[i].first->isDisk()&&stubs_[i].first->disk().value()<0){
	int lay=-stubs_[i].first->disk().value()-1;
	assert(lay>=0);
	assert(lay<5);
	nlay[lay]++;
      }
    }
    nlay[5] = 0;
    nlay2[5] = 0;
    nlay2[0] = nlay[0];
    for(unsigned int i=1; i<5; ++i)
      nlay2[i] = nlay2[i-1]+nlay[i];

    if(!w2) {
      //Write out the header
      out_ << "11111111111"<<marker;  
      out_ << (bitset<3>)bx_; //Dummy BX
      for (unsigned int i=0;i<3;i++) {
	out_ <<marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      out_ << "00000000000"<<marker<<(bitset<3>)bx_;  
      for (unsigned int i=3;i<6;i++) {
	out_ << marker<<(bitset<6>)nlay[i];
      }
      out_ << endl;
      
      //now write the stubs in disk order
      string leftover="";
      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()<0){
	    unsigned int lay=-stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      string stub=leftover+tmp.substr(0,32-leftover.size());
	      leftover=tmp.substr(32-leftover.size(),leftover.size()+4);
	      out_ << stub << endl;
	      if (leftover.size()==32){
		out_ << leftover << endl;
		leftover="";
	      }
	    }
	  }
	}
      }
      if (leftover.size()!=0) {
	string tmp="00000000000000000000000000000000";
	out_ << leftover << tmp.substr(0,32-leftover.size())<<endl;
      }
    }
    else{ 
      //format in2 for IPBUS
      unsigned int xline = 10;
      unsigned long int xword;
      //Header
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25)|(0x1ffffff);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //Stub counts
      xline = (xline+1)&15;
      xword = (nlay2[5])|(nlay2[4]<<6)|(nlay2[3]<<12)|(nlay2[2]<<18)
	|(nlay2[1]<<24)|(nlay2[0]<<30); 
      out_<< std::hex << xline<<"000";
      out_.fill('0');
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";
      //the stubs      
      for (unsigned int i=0;i<5;i++) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  if(stubs_[j].first->isDisk()&&stubs_[j].first->disk().value()<0){
	    unsigned int lay=-stubs_[j].first->disk().value()-1;
	    if (lay==i) {
	      string tmp=stubs_[j].first->strbare();
	      bitset<36> btmp(tmp);
	      xword = btmp.to_ulong();
	      xline = (xline+1)&15;
	      out_<< std::hex << xline<<"000";
	      out_.width(4);
	      out_<<(xword>>20);
	      out_<< " "<< xline<<"00";
	      out_.width(5);
	      out_<<(xword & 0xfffff);
	      out_<<" 51000003 51000007\n";
	    }
	  }
	}
      }

      //Trailer
      xline = (xline+1)&15;
      xword = ((unsigned long int)7<<33)|(((unsigned long int) bx_)<<25);
      out_<< std::hex << xline<<"000";
      out_.width(4);
      out_<<(xword>>20);
      out_<< " "<< xline<<"00";
      out_.width(5);
      out_<<(xword & 0xfffff);
      out_<<" 51000003 51000007\n";

    }    
  
    out_.close();

    //increment BX
    bx_++;
    if (bx_>7) bx_=0;

  }
  
  

  void clean() {
    for(unsigned int i=0;i<stubs_.size();i++){
      delete stubs_[i].first;
      delete stubs_[i].second;
    }
    stubs_.clear();
    indexphi_.clear();
    indexphi_.resize(5,0);
  }

  double phimin() const {return phimin_;}
  double phimax() const {return phimax_;}

  int phiregion() const {return phiregion_;}
  int phiregionoverlap() const {return phiregionoverlap_;}

  int layerdisk() const {return layerdisk_;}

  bool isLayer() const {return isLayer_;}
  bool isDisk() const {return isDisk_;}
  
  
private:

  double phimin_;
  double phimax_;
  vector<std::pair<FPGAStub*,L1TStub*> > stubs_;

  // index array for counting the number of stubs in each phi region for AS memories
  vector<unsigned int> indexphi_;

  int phiregion_;
  int phiregionoverlap_;
  
  int layerdisk_;

  int isLayer_;
  int isDisk_;

  
};

#endif
