#ifndef FPGATETABLEINNERDISK_H
#define FPGATETABLEINNERDISK_H

#include "FPGATETableBase.hh"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>


using namespace std;

class FPGATETableInnerDisk:public FPGATETableBase{

public:

  FPGATETableInnerDisk() {
    nbits_ = 9;
  }

  ~FPGATETableInnerDisk() {

  }


  void init(int disk1,
	    int disk2,
	    int rbits,
	    int zbits
	    ) {

    disk1_=disk1;
    disk2_=disk2;
    rbits_=rbits;
    zbits_=zbits;

    rbins_=(1<<rbits);
    rmind1_=0.0;
    rmaxd1_=rmaxdisk;
    dr_=rmaxdisk/rbins_;

    zbins_=(1<<zbits);
    zmind1_=zmean[disk1-1]-dzmax;
    zmaxd1_=zmean[disk1-1]+dzmax;
    dz_=2*dzmax/zbins_;

    zmeand2_=zmean[disk2-1];

    for (int irbin=0;irbin<rbins_;irbin++) {
      for (int izbin=0;izbin<zbins_;izbin++) {
	//int ibin=irbin+izbin*rbins_;
	int value=getLookupValue(irbin,izbin);
	//cout << "table "<<table_.size()<<" "<<value<<" "<<rmeanl2_<<endl;
	table_.push_back(value);
      }
    }
    if (writeVMTables) {
      writeVMTable("VMTableInnerD"+std::to_string(disk1_)+"D"+std::to_string(disk2_)+".txt");
    }
  }

  // negative return means that seed can not be formed
  int getLookupValue(int irbin, int izbin){

    double r1=rmind1_+irbin*dr_;
    double r2=rmind1_+(irbin+1)*dr_;

    double z1=zmind1_+izbin*dz_;
    double z2=zmind1_+(izbin+1)*dz_;


    double rmaxd2=-2*rmaxdisk;
    double rmind2=2*rmaxdisk;

    findr(r1,z1,rmind2,rmaxd2);
    findr(r1,z2,rmind2,rmaxd2);
    findr(r2,z1,rmind2,rmaxd2);
    findr(r2,z2,rmind2,rmaxd2);

    //cout << "rmind2 rmaxd2 "<<rmind2<<" "<<rmaxd2<<endl;
    
    assert(rmind2<rmaxd2);

    if (rmind2>rmaxdiskvm) return -1;
    if (rmaxd2<rmindiskvm) return -1;

    int NBINS=NLONGVMBINS*NLONGVMBINS/2; //divide by two for + and - z
    
    int rbinmin=NBINS*(rmind2-rmindiskvm)/(rmaxdiskvm-rmindiskvm);
    int rbinmax=NBINS*(rmaxd2-rmindiskvm)/(rmaxdiskvm-rmindiskvm);

    //cout << "zbinmin zminl2 "<<zbinmin<<" "<<zminl2<<endl;
    //cout << "zbinmax zmaxl2 "<<zbinmax<<" "<<zmaxl2<<endl;
    
    if (rbinmin<0) rbinmin=0;
    if (rbinmax>=NBINS) rbinmax=NBINS-1;

    //cout <<"rbminmin rbinmax "<<rbinmin<<" "<<rbinmax<<endl;
    
    assert(rbinmin<=rbinmax);
    assert(rbinmax-rbinmin<=(int)NLONGVMBINS);

    int value=rbinmin/8;
    value*=2;
    if (rbinmax/8-rbinmin/8>0) value+=1;
    value*=8;
    value+=(rbinmin&7);
    //cout << "zbinmax/8 zbinmin/8 value "<<zbinmax/8<<" "<<zbinmin/8<<" "<<value<<endl;
    assert(value/8<7);
    int deltar=rbinmax-rbinmin;
    assert(deltar<8);
    value+=(deltar<<6);
    assert(value<(1<<9));
    return value;
    
  }


  void findr(double r, double z, double& rmind2, double& rmaxd2){

    double rd2=rintercept(z0cut,r,z);

    //cout << "rd2 : "<<r<<" "<<z<<" "<<rd2<<endl;
    
    if (rd2<rmind2) rmind2=rd2;
    if (rd2>rmaxd2) rmaxd2=rd2;
    
    rd2=rintercept(-z0cut,r,z);

    //cout << "rd2 : "<<rd2<<endl;

    if (rd2<rmind2) rmind2=rd2;
    if (rd2>rmaxd2) rmaxd2=rd2;

  }

  double rintercept(double zcut, double r, double z) {

    return (zmeand2_-zcut)*r/(z-zcut);
    
  }

  int lookup(int rbin, int zbin) {

    int index=rbin*zbins_+zbin;
    assert(index<(int)table_.size());
    return table_[index];
    
  }
    
  /*
  
  void writephi(std::string fname) {

    ofstream out(fname.c_str());

    //cout << "writephi 2 phitableentries_ : "<<phitableentries_<<endl;

    for (int i=0;i<phitableentries_;i++){
      FPGAWord entry;
      //cout << "phitablebits_ : "<<phitablebits_<<endl;
      entry.set(i,phitablebits_);
      //out << entry.str()<<" "<<tablephi_[i]<<endl;
      out <<tablephi_[i]<<endl;
    }
    out.close();
  
  }

  */


private:

  int disk1_;
  int disk2_;
  int zbits_;
  int rbits_;
  
  int rbins_;
  double rmind1_;
  double rmaxd1_;
  double dr_;

  int zbins_;
  double zmind1_;
  double zmaxd1_;
  double dz_;
  
  double zmeand2_;
  

  
};



#endif



