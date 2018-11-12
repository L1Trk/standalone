#ifndef FPGATETABLEINNER_H
#define FPGATETABLEINNER_H

#include "FPGATETableBase.hh"
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>


using namespace std;

class FPGATETableInner:public FPGATETableBase{

public:

  FPGATETableInner() {
    nbits_ = 10;
  }

  ~FPGATETableInner() {

  }


  void init(int layer1,
	    int layer2,
	    int zbits,
	    int rbits
	    ) {

    layer1_=layer1;
    layer2_=layer2;
    zbits_=zbits;
    rbits_=rbits;

    bool extra=layer1==2&&layer2==3;
    
    rbins_=(1<<rbits);
    rminl1_=rmean[layer1-1]-drmax;
    rmaxl1_=rmean[layer1-1]+drmax;
    dr_=2*drmax/rbins_;

    zbins_=(1<<zbits);
    zminl1_=-zlength;
    zminl2_=zlength;
    dz_=2*zlength/zbins_;

    rmeanl2_=rmean[layer2-1];

    for (int izbin=0;izbin<zbins_;izbin++) {
      for (int irbin=0;irbin<rbins_;irbin++) {
	//int ibin=irbin+izbin*rbins_;
	int value=getLookupValue(izbin,irbin,extra);
	//cout << "table "<<table_.size()<<" "<<value<<" "<<rmeanl2_<<endl;
	table_.push_back(value);
      }
    }

    if (writeVMTables) {
      writeVMTable("VMTableInnerL"+std::to_string(layer1_)+"L"+std::to_string(layer2_)+".txt");
    }
    
  }

  // negative return means that seed can not be formed
  int getLookupValue(int izbin, int irbin, bool extra){

    double z1=zminl1_+izbin*dz_;
    double z2=zminl1_+(izbin+1)*dz_;

    double r1=rminl1_+irbin*dr_;
    double r2=rminl1_+(irbin+1)*dr_;

    if (extra and fabs(0.5*(z1+z2))<52.0) {   //This seeding combinations should not be for central region of detector
      return -1; 
    }
    
    double zmaxl2=-2*zlength;
    double zminl2=2*zlength;

    findz(z1,r1,zminl2,zmaxl2);
    findz(z1,r2,zminl2,zmaxl2);
    findz(z2,r1,zminl2,zmaxl2);
    findz(z2,r2,zminl2,zmaxl2);

    assert(zminl2<zmaxl2);

    if (zminl2>zlength) return -1;
    if (zmaxl2<-zlength) return -1;

    int NBINS=NLONGVMBINS*NLONGVMBINS;
    
    int zbinmin=NBINS*(zminl2+zlength)/(2*zlength);
    int zbinmax=NBINS*(zmaxl2+zlength)/(2*zlength);

    //cout << "zbinmin zminl2 "<<zbinmin<<" "<<zminl2<<endl;
    //cout << "zbinmax zmaxl2 "<<zbinmax<<" "<<zmaxl2<<endl;
    
    if (zbinmin<0) zbinmin=0;
    if (zbinmax>=NBINS) zbinmax=NBINS-1;

    assert(zbinmin<=zbinmax);
    assert(zbinmax-zbinmin<=(int)NLONGVMBINS);

    int value=zbinmin/8;
    value*=2;
    if (zbinmax/8-zbinmin/8>0) value+=1;
    value*=8;
    value+=(zbinmin&7);
    //cout << "zbinmax/8 zbinmin/8 value "<<zbinmax/8<<" "<<zbinmin/8<<" "<<value<<endl;
    assert(value/8<15);
    int deltaz=zbinmax-zbinmin;
    if (deltaz>7) {
      //cout << "deltaz = "<<deltaz<<endl;
      deltaz=7;
    }
    assert(deltaz<8);
    value+=(deltaz<<7);
    return value;
    
  }


  void findz(double z, double r, double& zminl2, double& zmaxl2){

    double zl2=zintercept(z0cut,z,r);

    if (zl2<zminl2) zminl2=zl2;
    if (zl2>zmaxl2) zmaxl2=zl2;
    
    zl2=zintercept(-z0cut,z,r);

    if (zl2<zminl2) zminl2=zl2;
    if (zl2>zmaxl2) zmaxl2=zl2;

  }

  double zintercept(double zcut, double z, double r) {

    return zcut+(z-zcut)*rmeanl2_/r;

  }

  int lookup(int zbin, int rbin) {

    int index=zbin*rbins_+rbin;
    return table_[index];
    
  }


private:


  int layer1_;
  int layer2_;
  int zbits_;
  int rbits_;
  
  int rbins_;
  double rminl1_;
  double rmaxl1_;
  double dr_;

  int zbins_;
  double zminl1_;
  double zminl2_;
  double dz_;
  
  double rmeanl2_;
  

  
};



#endif



