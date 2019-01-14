#ifndef FPGATRACKLET_H
#define FPGATRACKLET_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGAWord.hh"
#include "FPGATrack.hh"
#include "FPGALayerProjection.hh"
#include "FPGADiskProjection.hh"
#include "FPGALayerResidual.hh"
#include "FPGADiskResidual.hh"

using namespace std;

class FPGATracklet{

public:

  FPGATracklet(L1TStub* innerStub, L1TStub* outerStub,
	       FPGAStub* innerFPGAStub, FPGAStub* outerFPGAStub,
	       unsigned int homeSector,
	       double phioffset,
	       double rinv, double phi0, double z0, double t,
	       double rinvapprox, double phi0approx, 
	       double z0approx, double tapprox,
	       int irinv, int iphi0, 
	       int iz0, int it,
	       bool validproj[4],
	       int iphiproj[4], int izproj[4],
	       int iphider[4], int izder[4],
	       bool minusNeighbor[4], bool plusNeighbor[4],
	       double phiproj[4], double zproj[4],
	       double phider[4], double zder[4],
	       double phiprojapprox[4], double zprojapprox[4],
	       double phiderapprox[4], double zderapprox[4],
	       bool validprojdisk[5],
	       int iphiprojDisk[5], int irprojDisk[5],
	       int iphiderDisk[5], int irderDisk[5],
	       bool minusNeighborDisk[5], bool plusNeighborDisk[5],
	       double phiprojDisk[5], double rprojDisk[5],
	       double phiderDisk[5], double rderDisk[5],
	       double phiprojapproxDisk[5], double rprojapproxDisk[5],
	       double phiderapproxDisk[5], double rderapproxDisk[5],
	       bool disk, bool overlap=false){


    //cout << "New FPGATracklet"<<endl;

    homeSector_=homeSector;
    
    overlap_=overlap;
    disk_=disk;
    assert(!(disk&&overlap));
    barrel_=(!disk)&&(!overlap);

    trackletIndex_=-1;
    TCIndex_=-1;


    phioffset_=phioffset;

    assert(disk_||barrel_||overlap_);

    if (barrel_) assert(innerStub->layer()<6);
  
    innerStub_=innerStub;
    outerStub_=outerStub;
    innerFPGAStub_=innerFPGAStub;
    outerFPGAStub_=outerFPGAStub;
    rinv_=rinv;
    phi0_=phi0;
    z0_=z0;
    t_=t;
    rinvapprox_=rinvapprox;
    phi0approx_=phi0approx;
    z0approx_=z0approx;
    tapprox_=tapprox;

    fpgarinv_.set(irinv,nbitsrinv,false,__LINE__,__FILE__); 
    fpgaphi0_.set(iphi0,nbitsphi0,false,__LINE__,__FILE__); 
    fpgaz0_.set(iz0,nbitsz0,false,__LINE__,__FILE__);
    fpgat_.set(it,nbitst,false,__LINE__,__FILE__);       

    assert(innerStub_->layer()<6||innerStub_->disk()<5);
    assert(outerStub_->layer()<6||outerStub_->disk()<5);


    if (innerStub_->layer()==0) {
      rproj_[0]=rmeanL3;
      rproj_[1]=rmeanL4;
      rproj_[2]=rmeanL5;
      rproj_[3]=rmeanL6;
      projlayer_[0]=3;
      projlayer_[1]=4;
      projlayer_[2]=5;
      projlayer_[3]=6;
    }

    if (innerStub_->layer()==1) {
      rproj_[0]=rmeanL1;
      rproj_[1]=rmeanL4;
      rproj_[2]=rmeanL5;
      rproj_[3]=rmeanL6;
      projlayer_[0]=1;
      projlayer_[1]=4;
      projlayer_[2]=5;
      projlayer_[3]=6;
    }
    
    if (innerStub_->layer()==2) {
      rproj_[0]=rmeanL1;
      rproj_[1]=rmeanL2;
      rproj_[2]=rmeanL5;
      rproj_[3]=rmeanL6;
      projlayer_[0]=1;
      projlayer_[1]=2;
      projlayer_[2]=5;
      projlayer_[3]=6;
    }

    if (innerStub_->layer()==4) {
      rproj_[0]=rmeanL1;
      rproj_[1]=rmeanL2;
      rproj_[2]=rmeanL3;
      rproj_[3]=rmeanL4;
      projlayer_[0]=1;
      projlayer_[1]=2;
      projlayer_[2]=3;
      projlayer_[3]=4;
    }

    if (innerStub_->layer()>999) {
      
      assert((innerStub->disk()==1)||(innerStub->disk()==3)||
	     (innerStub->disk()==-1)||(innerStub->disk()==-3));

      rproj_[0]=rmeanL1;
      rproj_[1]=rmeanL2;
      projlayer_[0]=1;
      projlayer_[1]=2;
      projlayer_[2]=-1;
      projlayer_[3]=-1;

      if (overlap_) {
	assert((innerStub->disk()==1)||(innerStub->disk()==-1));
	if (innerStub->disk()==1) {
	  zprojdisk_[0]=zmeanD1;
	  zprojdisk_[1]=zmeanD2;
	  zprojdisk_[2]=zmeanD3;
	  zprojdisk_[3]=zmeanD4;
	  zprojdisk_[4]=zmeanD5;
	} else {
	  zprojdisk_[0]=-zmeanD1;
	  zprojdisk_[1]=-zmeanD2;
	  zprojdisk_[2]=-zmeanD3;
	  zprojdisk_[3]=-zmeanD4;
	  zprojdisk_[4]=-zmeanD5;
	}
	projdisk_[0]=1;
	projdisk_[1]=2;
	projdisk_[2]=3;
	projdisk_[3]=4;
	projdisk_[4]=5;
      } else {
	
	if (innerStub->disk()==1) {
	  zprojdisk_[0]=zmeanD3;
	  zprojdisk_[1]=zmeanD4;
	  zprojdisk_[2]=zmeanD5;
	  projdisk_[0]=3;
	  projdisk_[1]=4;
	  projdisk_[2]=5;
	}
	if (innerStub->disk()==3) {
	  zprojdisk_[0]=zmeanD1;
	  zprojdisk_[1]=zmeanD2;
	  zprojdisk_[2]=zmeanD5;
	  projdisk_[0]=1;
	  projdisk_[1]=2;
	  projdisk_[2]=5;
	}
	
	if (innerStub->disk()==-1) {
	  zprojdisk_[0]=-zmeanD3;
	  zprojdisk_[1]=-zmeanD4;
	  zprojdisk_[2]=-zmeanD5;
	  projdisk_[0]=-3;
	  projdisk_[1]=-4;
	  projdisk_[2]=-5;
	}
	if (innerStub->disk()==-3) {
	  zprojdisk_[0]=-zmeanD1;
	  zprojdisk_[1]=-zmeanD2;
	  zprojdisk_[2]=-zmeanD5;
	  projdisk_[0]=-1;
	  projdisk_[1]=-2;
	  projdisk_[2]=-5;
	}	
      }

    } else {
      
      int sign=1;
      if (it<0) sign=-1;

      zprojdisk_[0]=sign*zmeanD1;
      zprojdisk_[1]=sign*zmeanD2;
      zprojdisk_[2]=sign*zmeanD3;
      zprojdisk_[3]=sign*zmeanD4;
      zprojdisk_[4]=sign*zmeanD5;
      projdisk_[0]=sign*1;
      projdisk_[1]=sign*2;
      projdisk_[2]=sign*3;
      projdisk_[3]=sign*4;
      projdisk_[4]=sign*5;

    }
    

    
      
    if (barrel_) {
      //barrel
      for (int i=0;i<4;i++) {

	if (!validproj[i]) continue;

	layerproj_[projlayer_[i]-1].init(projlayer_[i],
					 rproj_[i],
					 iphiproj[i],
					 izproj[i],
					 iphider[i],
					 izder[i],
					 minusNeighbor[i],
					 plusNeighbor[i],
					 phiproj[i],
					 zproj[i],
					 phider[i],
					 zder[i],
					 phiprojapprox[i],
					 zprojapprox[i],
					 phiderapprox[i],
					 zderapprox[i]);	  
	
      }
      //Now handle projections to the disks
      for (int i=0;i<5;i++) {
	
	if (!validprojdisk[i]) continue;

	diskproj_[i].init(i+1,
			  zprojdisk_[i],
			  iphiprojDisk[i],
			  irprojDisk[i],
			  iphiderDisk[i],
			  irderDisk[i],
			  minusNeighborDisk[i],
			  plusNeighborDisk[i],
			  phiprojDisk[i],
			  rprojDisk[i],
			  phiderDisk[i],
			  rderDisk[i],
			  phiprojapproxDisk[i],
			  rprojapproxDisk[i],
			  phiderapproxDisk[i],
			  rderapproxDisk[i]);	
      }
    } 

    if (disk_) {
      //cout << "FPGATracklet making disk tracklet" << endl;
      //disk stub 
      for (int i=0;i<3;i++) {

	if (!validprojdisk[i]) continue;

	diskproj_[abs(projdisk_[i])-1].init(projdisk_[i],
				     zprojdisk_[i],
				     iphiprojDisk[i],
				     irprojDisk[i],
				     iphiderDisk[i],
				     irderDisk[i],
				     minusNeighborDisk[i],
				     plusNeighborDisk[i],
				     phiprojDisk[i],
				     rprojDisk[i],
				     phiderDisk[i],
				     rderDisk[i],
				     phiprojapproxDisk[i],
				     rprojapproxDisk[i],
				     phiderapproxDisk[i],
				     rderapproxDisk[i]);	

      }
      for (int i=0;i<2;i++) {

	if (!validproj[i]) continue;
	
	layerproj_[i].init(i+1,
			   rproj_[i],
			   iphiproj[i],
			   izproj[i],
			   iphider[i],
			   izder[i],
			   minusNeighbor[i],
			   plusNeighbor[i],
			   phiproj[i],
			   zproj[i],
			   phider[i],
			   zder[i],
			   phiprojapprox[i],
			   zprojapprox[i],
			   phiderapprox[i],
			   zderapprox[i]);	  
      
      }
    }

    if (overlap_) {
      //projections to layers
      for (int i=0;i<1;i++) {
	assert(rproj_[i]<60.0);

	if (!validproj[i]) continue;

	layerproj_[i].init(i+1,
			   rproj_[i],
			   iphiproj[i],
			   izproj[i],
			   iphider[i],
			   izder[i],
			   minusNeighbor[i],
			   plusNeighbor[i],
			   phiproj[i],
			   zproj[i],
			   phider[i],
			   zder[i],
			   phiprojapprox[i],
			   zprojapprox[i],
			   phiderapprox[i],
			   zderapprox[i]);	  
      }
      //Now handle projections to the disks
      for (int i=0;i<4;i++) {
	if (!validprojdisk[i]) continue;
	int offset=1;
	if (outerStub->layer()+1==2&&innerStub->layer()+1==3) offset=0;
	diskproj_[i+offset].init(i+offset+1,
				 zprojdisk_[i],
				 iphiprojDisk[i],
				 irprojDisk[i],
				 iphiderDisk[i],
				 irderDisk[i],
				 minusNeighborDisk[i],
				 plusNeighborDisk[i],
				 phiprojDisk[i],
				 rprojDisk[i],
				 phiderDisk[i],
				 rderDisk[i],
				 phiprojapproxDisk[i],
				 rprojapproxDisk[i],
				 phiderapproxDisk[i],
				 rderapproxDisk[i]);	

      }
    } 


    ichisqfit_.set(-1,8,false);

  }




  ~FPGATracklet() {

  }



  L1TStub* innerStub() {return innerStub_;}

  L1TStub* outerStub() {return outerStub_;}

  std::string addressstr() {
    std::ostringstream oss;
    oss << innerFPGAStub_->phiregionaddressstr()<<"|" 
	<< outerFPGAStub_->phiregionaddressstr();

    return oss.str();

  }


  //Tracklet parameters print out
  std::string trackletparstr() {
    std::ostringstream oss;
    if(writeoutReal){
      oss << fpgarinv_.value()*krinvpars<<" "
        << fpgaphi0_.value()*kphi0pars<<" "
        << fpgaz0_.value()*kz<<" "
        << fpgat_.value()*ktpars;
    }

    //Binary Print out
      if(!writeoutReal){
	oss << innerFPGAStub_->stubindex().str()<<"|" 
            << outerFPGAStub_->stubindex().str()<<"|"
            << fpgarinv_.str()<<"|"
	    << fpgaphi0_.str()<<"|"
	    << fpgaz0_.str()<<"|"
	    << fpgat_.str();
      }
    return oss.str();
  }

  std::string vmstrlayer(int layer, unsigned int allstubindex) {
    std::ostringstream oss;
    FPGAWord index;
    if (allstubindex>=(1<<7)) {
      cout << "Warning projection number too large!"<<endl;	
      index.set((1<<7)-1,7,true,__LINE__,__FILE__);
    } else {
      index.set(allstubindex,7,true,__LINE__,__FILE__);
    }
    oss << index.str()<<"|"<<layerproj_[layer-1].fpgaphiprojvm().str()
	<<"|"<< layerproj_[layer-1].fpgazbin1projvm().str() 
        <<"|"<< layerproj_[layer-1].fpgazbin2projvm().str();
	//<<"|"<< layerproj_[layer-1].fpgazprojvm().str();
    return oss.str();

  }

  std::string vmstrdisk(int disk, unsigned int allstubindex) {
    std::ostringstream oss;
    FPGAWord index;
    if (allstubindex>=(1<<7)) {
      cout << "Warning projection number too large!"<<endl;	
      index.set((1<<7)-1,7,true,__LINE__,__FILE__);
    } else {
      index.set(allstubindex,7,true,__LINE__,__FILE__);
    } 
    oss << index.str()<<"|"<<diskproj_[disk-1].fpgaphiprojvm().str()
	<<"|"<< diskproj_[disk-1].fpgarprojvm().str();
    return oss.str();

  }

  
  std::string trackletprojstr(int layer) const {
    assert(layer>=1&&layer<=6);
    std::ostringstream oss;
    FPGAWord tmp;
    if (trackletIndex_<0||trackletIndex_>63) {
      cout << "trackletIndex_ = "<<trackletIndex_<<endl;
      assert(0);
    }
    tmp.set(trackletIndex_,6,true,__LINE__,__FILE__);
    FPGAWord tcid;
    tcid.set(TCIndex_,7,true,__LINE__,__FILE__);
    oss << layerproj_[layer-1].plusNeighbor()<<"|"
        << layerproj_[layer-1].minusNeighbor()<<"|"
        << tcid.str()<<"|"
	<< tmp.str()<<"|"
        << layerproj_[layer-1].fpgaphiproj().str()<<"|"
	<< layerproj_[layer-1].fpgazproj().str()<<"|"
	<< layerproj_[layer-1].fpgaphiprojder().str()<<"|"
	<< layerproj_[layer-1].fpgazprojder().str();

    return oss.str();

  }
  std::string trackletprojstrD(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    std::ostringstream oss;
    FPGAWord tmp;
    if (trackletIndex_<0||trackletIndex_>63) {
      cout << "trackletIndex_ = "<<trackletIndex_<<endl;
      assert(0);
    }
    tmp.set(trackletIndex_,6,true,__LINE__,__FILE__);
    FPGAWord tcid;
    tcid.set(TCIndex_,7,true,__LINE__,__FILE__);

    oss << diskproj_[abs(disk)-1].plusNeighbor()<<"|"
        << diskproj_[abs(disk)-1].minusNeighbor()<<"|" 
        << tcid.str()<<"|" 
        << tmp.str()<<"|"
	<< diskproj_[abs(disk)-1].fpgaphiproj().str()<<"|"
	<< diskproj_[abs(disk)-1].fpgarproj().str()<<"|"
	<< diskproj_[abs(disk)-1].fpgaphiprojder().str()<<"|"
	<< diskproj_[abs(disk)-1].fpgarprojder().str();

    return oss.str();

  }

  std::string trackletprojstrlayer(int layer) const {
    return trackletprojstr(layer);
  }
  std::string trackletprojstrdisk(int disk) const {
    std::ostringstream oss;
    return trackletprojstrD(disk);

  }

  bool validProj(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].valid();
  }

  FPGAWord fpgaphiprojder(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgaphiprojder();
  }

  FPGAWord fpgazproj(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgazproj();
  }

  FPGAWord fpgaphiproj(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgaphiproj();
  }

  FPGAWord fpgazprojder(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgazprojder();
  }

  int zbin1projvm(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgazbin1projvm().value();
  }

  int zbin2projvm(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgazbin2projvm().value();
  }

  int finezvm(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgafinezvm().value();
  }

  int rbin1projvm(int disk) const {
    assert(disk>=1&&disk<=5);
    return diskproj_[disk-1].fpgarbin1projvm().value();
  }

  int rbin2projvm(int disk) const {
    assert(disk>=1&&disk<=5);
    return diskproj_[disk-1].fpgarbin2projvm().value();
  }

  int finervm(int disk) const {
    assert(disk>=1&&disk<=5);
    return diskproj_[disk-1].fpgafinervm().value();
  }

  int phiprojvm(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgaphiprojvm().value();
  }

  int zprojvm(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].fpgazprojvm().value();
  }


  
  double phiproj(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].phiproj();
  }

  double phiprojder(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].phiprojder();
  }

  double zproj(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].zproj();
  }

  double zprojder(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].zprojder();
  }



  
  double zprojapprox(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].zprojapprox();
  }

  double zprojderapprox(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].zprojderapprox();
  }

  double phiprojapprox(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].phiprojapprox();
  }

  double phiprojderapprox(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].phiprojderapprox();
  }

  

  double rproj(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].rproj();
  }

  bool minusNeighbor(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].minusNeighbor();
  }

  bool plusNeighbor(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerproj_[layer-1].plusNeighbor();
  }



  double rstub(int layer) {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].rstub();
  }





  //Disks residuals

  bool validProjDisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].valid();
  }


  FPGAWord fpgaphiresiddisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].fpgaphiresid();
  }

  FPGAWord fpgarresiddisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].fpgarresid();
  }

  double phiresiddisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].phiresid();
  }

  double rresiddisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].rresid();
  }

  double phiresidapproxdisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].phiresidapprox();
  }

  double rresidapproxdisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].rresidapprox();
  }


  double zstubdisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].zstub();
  }


  double zprojdisk(int disk) const {
    if (!disk_) return zprojdisk_[abs(disk)-1];
    for (int i=0;i<3;i++) {
      if (projdisk_[i]==disk) {
	return zprojdisk_[i];
      }
    }
    assert(0);
    return 0.0;

  }

  void setBendIndex(int bendIndex,int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    diskproj_[abs(disk)-1].setBendIndex(bendIndex);
  }

  FPGAWord getBendIndex(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].getBendIndex();
  }
  

  double alphadisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].alpha();
  }

  FPGAWord ialphadisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].ialpha();
  }

  

  FPGAWord fpgaphiprojdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].fpgaphiproj();
  }

  FPGAWord fpgaphiprojderdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].fpgaphiprojder();
  }

  FPGAWord fpgarprojdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].fpgarproj();
  }
  
  FPGAWord fpgarprojderdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].fpgarprojder();
  }

  

  double phiprojapproxdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].phiprojapprox();
  }

  double phiprojderapproxdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].phiprojderapprox();
  }
  
  double rprojapproxdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].rprojapprox();
  }

  double rprojderapproxdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].rprojderapprox();
  }



  double phiprojdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].phiproj();
  }

  double phiprojderdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].phiprojder();
  }
  
  double rprojdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].rproj();
  }
  
  double rprojderdisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].rprojder();
  }



  bool minusNeighborDisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].minusNeighbor();
  }


  bool plusNeighborDisk(int disk) const {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskproj_[abs(disk)-1].plusNeighbor();
  }


  bool matchdisk(int disk) {
    assert(abs(disk)>=1&&abs(disk)<=5);
    return diskresid_[abs(disk)-1].valid();
  }

  void addMatch(int layer, int ideltaphi, int ideltaz, 
		double dphi, double dz, 
		double dphiapprox, double dzapprox, 
		int stubid,double rstub,
		std::pair<FPGAStub*,L1TStub*> stubptrs){

    //cout << "FPGATracklet::addMatch adding match in layer = "<<layer<<endl;

    assert(layer>=1&&layer<=6);

    layerresid_[layer-1].init(layer,ideltaphi,ideltaz,stubid,dphi,dz,dphiapprox,dzapprox,rstub,stubptrs);
    
  }



  void addMatchDisk(int disk, int ideltaphi, int ideltar, 
		    double dphi, double dr, 
		    double dphiapprox, double drapprox, double alpha,
		    int stubid,double zstub,
		    std::pair<FPGAStub*,L1TStub*> stubptrs){

    //cout << "addMatchDisk1 "<<disk<<" "<<ideltaphi<<endl; 

    assert(abs(disk)>=1&&abs(disk)<=5);

    diskresid_[abs(disk)-1].init(disk,ideltaphi,ideltar,stubid,dphi,dr,dphiapprox,drapprox,zstub,alpha,stubptrs.first->alphanew(),stubptrs);
      

  }

  std::string candmatchstr(int layer) {
    assert(layer>=1&&layer<=6);
    std::ostringstream oss;
    FPGAWord tmp;
    if (trackletIndex_<0||trackletIndex_>63) {
      cout << "trackletIndex_ = "<<trackletIndex_<<endl;
      assert(0);
    }
    tmp.set(trackletIndex_,6);
    oss << tmp.str()<<"|";

    return oss.str();
  }

  std::string candmatchdiskstr(int disk) {
    assert(disk>=1&&disk<=5);
    std::ostringstream oss;
    FPGAWord tmp;
    if (trackletIndex_<0||trackletIndex_>63) {
      cout << "trackletIndex_ = "<<trackletIndex_<<endl;
      assert(0);
    }
    tmp.set(trackletIndex_,6);
    oss << tmp.str()<<"|";

    return oss.str();
  }






  int nMatches() {

    int nmatches=0;

    for (int i=0;i<6;i++) {
      if (layerresid_[i].valid()) {
	nmatches++;
      }
    }
    
    return nmatches;
    
  }

  int nMatchesDisk(bool skipD5=false) {

    int nmatches=0;

    int lastdisk=5;
    if (skipD5) lastdisk--; 
    for (int i=0;i<lastdisk;i++) {
      if (diskresid_[i].valid()) {
	nmatches++;
      }
    }
    return nmatches;

  }


  int nMatchesOverlap() {

    int nmatches=0;

    assert(overlap_);

    for (int i=1;i<5;i++) {
      if (diskresid_[i].valid()) {
	nmatches++;
      }
    }

    for (int i=0;i<2;i++) {
      if (layerresid_[i].valid()) nmatches++;
    }
    return nmatches;
    
  }



  bool match(int layer) {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].valid();
  }


  std::string fullmatchstr(int layer) {
    assert(layer>=1&&layer<=6);
    std::ostringstream oss;

    FPGAWord tmp;
    if (trackletIndex_<0||trackletIndex_>63) {
      cout << "trackletIndex_ = "<<trackletIndex_<<endl;
      assert(0);
    }
    tmp.set(trackletIndex_,6,true,__LINE__,__FILE__);
    FPGAWord tcid;
    tcid.set(TCIndex_,7,true,__LINE__,__FILE__);
    oss << tcid.str()<<"|"
        << tmp.str()<<"|"
	<< layerresid_[layer-1].fpgastubid().str()<<"|"
	<< layerresid_[layer-1].fpgaphiresid().str()<<"|"
	<< layerresid_[layer-1].fpgazresid().str();

    return oss.str();

  }

  std::string fullmatchdiskstr(int disk) {
    assert(disk>=1&&disk<=5);
    std::ostringstream oss;

    FPGAWord tmp;
    if (trackletIndex_<0||trackletIndex_>63) {
      cout << "trackletIndex_ = "<<trackletIndex_<<endl;
      assert(0);
    }
    tmp.set(trackletIndex_,6,true,__LINE__,__FILE__);
    FPGAWord tcid;
    tcid.set(TCIndex_,7,true,__LINE__,__FILE__);
    oss << tcid.str()<<"|"
        << tmp.str()<<"|"
	<< diskresid_[disk-1].fpgastubid().str()<<"|"
	<< diskresid_[disk-1].fpgaphiresid().str()<<"|"
	<< diskresid_[disk-1].fpgarresid().str();

    return oss.str();

  }


  bool validResid(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].valid();
  }

  
  std::pair<FPGAStub*,L1TStub*> stubptrs(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].stubptrs();
  }

  
  double phiresid(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].phiresid();
  }

  double phiresidapprox(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].phiresidapprox();
  }

  double zresid(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].zresid();
  }

  double zresidapprox(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].zresidapprox();
  }





  FPGAWord fpgaphiresid(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].fpgaphiresid();
  }

  FPGAWord fpgazresid(int layer) const {
    assert(layer>=1&&layer<=6);
    return layerresid_[layer-1].fpgazresid();
  }


  std::vector<L1TStub*> getL1Stubs() {

    std::vector<L1TStub*> tmp;

    tmp.push_back(innerStub_);
    tmp.push_back(outerStub_);

    for (unsigned int i=0;i<6;i++) {
      if (layerresid_[i].valid()) {
	tmp.push_back(layerresid_[i].stubptrs().second);
      }
    }

    for (unsigned int i=0;i<5;i++) {
      if (diskresid_[i].valid()) tmp.push_back(diskresid_[i].stubptrs().second);
    }

    return tmp;

  }

  std::map<int, int> getStubIDs() {


    std::map<int, int> stubIDs;
    
    //printf("Track Type Barrel %i   Disk %i   Overlap %i  \n",barrel_,disk_,overlap_);
    //printf(" irinv %i,  iphi %i  it %i \n",irinvfit().value(),iphi0fit().value(),itfit().value());
    //printf(" Matches:  barrel  %i   disks %i \n",nMatches(),nMatchesDisk());
    //for(int i=0; i<16; i++)  stubIDs[i] = 63; //this is no stub code value

    // For future reference, *resid_[i] uses i as the absolute stub index. (0-5 for barrel, 0-4 for disk)
    // On the other hand, proj*_[i] uses i almost like *resid_[i], except the *seeding* layer indices are removed entirely.
    // E.g. An L3L4 track has 0=L1, 1=L2, 2=L4, 3=L5 for the barrels (for proj*_[i])

    assert(innerFPGAStub_->stubindex().nbits()==7);
    assert(innerFPGAStub_->phiregion().nbits()==3);
    assert(outerFPGAStub_->stubindex().nbits()==7);
    assert(outerFPGAStub_->phiregion().nbits()==3);
	
    if(barrel_) {
      for(int i=0; i<6; i++) {

        //check barrel
	if (layerresid_[i].valid()) {
		  // two extra bits to indicate if the matched stub is local or from neighbor
		  int location = 1;  // local
		  if (minusNeighbor(i+1)) location = 0; // phi-
		  if (plusNeighbor(i+1)) location = 2;  // phi+
		  location<<=layerresid_[i].fpgastubid().nbits();
		  
		  stubIDs[1+i] = layerresid_[i].fpgastubid().value()+location;
        }			  		  
		  
	//check disk
	if(diskresid_[i].valid()) {
	  // two extra bits to indicate if the matched stub is local or from neighbor
	  int location = 1;  // local
	  if (minusNeighborDisk(i+1)) location = 0; // phi-
	  if (plusNeighborDisk(i+1)) location = 2;  // phi+
	  location<<=diskresid_[i].fpgastubid().nbits();

	  if(itfit().value() < 0) {
		stubIDs[-11-i] = diskresid_[i].fpgastubid().value()+location;
	  } else {
	    stubIDs[11+i] = diskresid_[i].fpgastubid().value()+location;
	  }  
	}	 			  		  		  
      }
		  
		  
      //get stubs making up tracklet
      //printf(" inner %i  outer %i layers \n",innerFPGAStub_.layer().value(),outerFPGAStub_.layer().value());
      if (hourglass) {
	stubIDs[innerFPGAStub_->layer().value()+1] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<10);
	stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<10);
      } else {
	stubIDs[innerFPGAStub_->layer().value()+1] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<9);
	stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<9);
      }
		  		  
    } else if (disk_) {
      for(int i=0; i<5; i++) {
	
	//check barrel
	if(layerresid_[i].valid()) {
		  // two extra bits to indicate if the matched stub is local or from neighbor
		  int location = 1;  // local
		  if (minusNeighbor(i+1)) location = 0; // phi-
		  if (plusNeighbor(i+1)) location = 2;  // phi+
		  location<<=layerresid_[i].fpgastubid().nbits();
		  
		  stubIDs[1+i] = layerresid_[i].fpgastubid().value()+location;
        }
	
	//check disks
        if(i==4 && layerresid_[1].valid()) continue; // Don't add D5 if track has L1 stub
	if(diskresid_[i].valid()) {
	  // two extra bits to indicate if the matched stub is local or from neighbor
	  int location = 1;  // local
	  if (minusNeighborDisk(i+1)) location = 0; // phi-
	  if (plusNeighborDisk(i+1)) location = 2;  // phi+
	  location<<=diskresid_[i].fpgastubid().nbits();
	  
	  if(innerStub_->disk() < 0) {
	    stubIDs[-11-i] = diskresid_[i].fpgastubid().value()+location;
	  } else {
	    stubIDs[11+i] = diskresid_[i].fpgastubid().value()+location;
	  }
	}	 			  
      }
		  
      //get stubs making up tracklet
      //printf(" inner %i  outer %i disks \n",innerFPGAStub_.disk().value(),outerFPGAStub_.disk().value());
      if(innerFPGAStub_->disk().value() < 0) { //negative side runs 6-10
	if (hourglass) {
	  stubIDs[innerFPGAStub_->disk().value()-10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<10);
	  stubIDs[outerFPGAStub_->disk().value()-10] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<10);
	} else {
	  stubIDs[innerFPGAStub_->disk().value()-10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<9);
	  stubIDs[outerFPGAStub_->disk().value()-10] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<9);
	}
      } else { // positive side runs 11-15]
	if (hourglass) {
	  stubIDs[innerFPGAStub_->disk().value()+10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<10);
	  stubIDs[outerFPGAStub_->disk().value()+10] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<10);
	} else {
	  stubIDs[innerFPGAStub_->disk().value()+10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<9);
	  stubIDs[outerFPGAStub_->disk().value()+10] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<9);
	}
      }		  

    } else if (overlap_) {
      
      for(int i=0; i<5; i++) {
	
	//check barrel
	if(layerresid_[i].valid()) {
		  // two extra bits to indicate if the matched stub is local or from neighbor
		  int location = 1;  // local
		  if (minusNeighbor(i+1)) location = 0; // phi-
		  if (plusNeighbor(i+1)) location = 2;  // phi+
		  location<<=layerresid_[i].fpgastubid().nbits();
		  
		  stubIDs[1+i] = layerresid_[i].fpgastubid().value()+location;
	}

	//check disks
	if(diskresid_[i].valid()) {
	  // two extra bits to indicate if the matched stub is local or from neighbor
	  int location = 1;  // local
	  if (minusNeighborDisk(i+1)) location = 0; // phi-
	  if (plusNeighborDisk(i+1)) location = 2;  // phi+
	  location<<=diskresid_[i].fpgastubid().nbits();
	  
	  if(innerStub_->disk() < 0) { // if negative overlap
            if(innerFPGAStub_->layer().value()!=2 || !layerresid_[0].valid() || i!=3 ) { // Don't add D4 if this is an L3L2 track with an L1 stub
	      stubIDs[-11-i] = diskresid_[i].fpgastubid().value()+location;
            }
	  } else {
            if(innerFPGAStub_->layer().value()!=2 || !layerresid_[0].valid() || i!=3 ) {
	      stubIDs[11+i] = diskresid_[i].fpgastubid().value()+location;
            }
	  }
	}	 			  
      }
      
      //get stubs making up tracklet
      //printf(" inner %i  outer %i layers \n",innerFPGAStub_.layer().value(),outerFPGAStub_.layer().value());
      //printf(" inner %i  outer %i disks \n",innerFPGAStub_.disk().value(),outerFPGAStub_.disk().value());
      if(innerFPGAStub_->layer().value()==2) { // L3L2 track
	if (hourglass) {
	  stubIDs[innerFPGAStub_->layer().value()+1] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<10);
	  stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<10);
	} else {
	  stubIDs[innerFPGAStub_->layer().value()+1] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<9);
	  stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<9);
	}
      } else if(innerFPGAStub_->disk().value() < 0) { //negative side runs -11 - -15
	if (hourglass) {
	  stubIDs[innerFPGAStub_->disk().value()-10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<10);
	  stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<10);
	} else {
	  stubIDs[innerFPGAStub_->disk().value()-10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<9);
	  stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<9);
	}
      } else { // positive side runs 11-15]
	if (hourglass) {
	  stubIDs[innerFPGAStub_->disk().value()+10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<10);
	  stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<10);
	} else {
	  stubIDs[innerFPGAStub_->disk().value()+10] = ((innerFPGAStub_->phiregion().value())<<7)+innerFPGAStub_->stubindex().value()+(1<<9);
	  stubIDs[outerFPGAStub_->layer().value()+1] = ((outerFPGAStub_->phiregion().value())<<7)+outerFPGAStub_->stubindex().value()+(1<<9);
	}
      }		  
		  
    }


    //cout << "New track layer = "<<innerFPGAStub_->layer().value()+1<<endl;
    //for(std::map<int,int>::iterator i=stubIDs.begin(); i!=stubIDs.end(); i++) {
    //      printf("Layer/Disk %i  ID  %i \n",i->first, i->second);
    //}

	//for (const auto & id : stubIDs) {
	//	if (id.second<0) cout<<"i stubID : "<<id.first<<" " <<id.second<<endl;
	//	assert(id.second>=0);
	//}
	
    return stubIDs;
  }


  double rinv() const { return rinv_; }
  double phi0() const { return phi0_; }
  double t() const { return t_; }
  double z0() const { return z0_; }

  double rinvapprox() const { return rinvapprox_; }
  double phi0approx() const { return phi0approx_; }
  double tapprox() const { return tapprox_; }
  double z0approx() const { return z0approx_; }


  FPGAWord fpgarinv() const { return fpgarinv_; }
  FPGAWord fpgaphi0() const { return fpgaphi0_; }
  FPGAWord fpgat() const { return fpgat_; }
  FPGAWord fpgaz0() const { return fpgaz0_; }

  double rinvfit() const { return rinvfit_; }
  double phi0fit() const { return phi0fit_; }
  double tfit() const { return tfit_; }
  double z0fit() const { return z0fit_; }
  double chiSqfit() const { return chisqfit_; }

  double rinvfitexact() const { return rinvfitexact_; }
  double phi0fitexact() const { return phi0fitexact_; }
  double tfitexact() const { return tfitexact_; }
  double z0fitexact() const { return z0fitexact_; }

  FPGAWord irinvfit() const { return irinvfit_; }
  FPGAWord iphi0fit() const { return iphi0fit_; }
  FPGAWord itfit() const { return itfit_; }
  FPGAWord iz0fit() const { return iz0fit_; }
  FPGAWord ichiSqfit() const { return ichisqfit_; }


  void setFitPars(double rinvfit, double phi0fit, double tfit,
		  double z0fit, double chisqfit,
		  double rinvfitexact, double phi0fitexact, double tfitexact,
		  double z0fitexact, double chisqfitexact,
		  int irinvfit, int iphi0fit, int itfit,
		  int iz0fit, int ichisqfit,
                  const vector<L1TStub*>& l1stubs = vector<L1TStub*>()){

    rinvfit_=rinvfit;
    phi0fit_=phi0fit;
    tfit_=tfit;
    z0fit_=z0fit;
    chisqfit_=chisqfit;

    rinvfitexact_=rinvfitexact;
    phi0fitexact_=phi0fitexact;
    tfitexact_=tfitexact;
    z0fitexact_=z0fitexact;
    chisqfitexact_=chisqfitexact;
    
    if (irinvfit>(1<<14)) irinvfit=(1<<14);
    if (irinvfit<=-(1<<14)) irinvfit=-(1<<14)+1;
    irinvfit_.set(irinvfit,15,false,__LINE__,__FILE__);
    iphi0fit_.set(iphi0fit,19,false,__LINE__,__FILE__);
    itfit_.set(itfit,14,false,__LINE__,__FILE__);

    if (iz0fit>=(1<<(nbitsz0-1))) {
      iz0fit=(1<<(nbitsz0-1))-1; 
    }

    if (iz0fit<=-(1<<(nbitsz0-1))) {
      iz0fit=1-(1<<(nbitsz0-1)); 
    }

    iz0fit_.set(iz0fit,nbitsz0,false,__LINE__,__FILE__);
    ichisqfit_.set(ichisqfit,8,true,__LINE__,__FILE__);

    fpgatrack_=new FPGATrack(makeTrack(l1stubs));

  }


  std::string trackfitstr() {
    std::ostringstream oss;

    string stubid0="111111111";
    string stubid1="111111111";
    string stubid2="111111111";
    string stubid3="111111111";

    if (isBarrel()) {
      if (layer()==1) {
	if (layerresid_[2].valid()) {
	  stubid0=layerresid_[2].fpgastubid().str();
	}
	if (layerresid_[3].valid()) {
	  stubid1=layerresid_[3].fpgastubid().str();
	}
	if (layerresid_[4].valid()) {
	  stubid2=layerresid_[4].fpgastubid().str();
	}
	if (layerresid_[5].valid()) {
	  stubid3=layerresid_[5].fpgastubid().str();
	}
	if (diskresid_[0].valid()) {
	  stubid3=diskresid_[0].fpgastubid().str();
	}
if (diskresid_[1].valid()) {
	  stubid2=diskresid_[1].fpgastubid().str();
	}
	if (diskresid_[2].valid()) {
	  stubid1=diskresid_[2].fpgastubid().str();
	}
	if (diskresid_[3].valid()) {
	  stubid0=diskresid_[3].fpgastubid().str();
	}
      }

      if (layer()==3) {
	if (layerresid_[0].valid()) {
	  stubid0=layerresid_[0].fpgastubid().str();
	}
	if (layerresid_[1].valid()) {
	  stubid1=layerresid_[1].fpgastubid().str();
	}
	if (layerresid_[4].valid()) {
	  stubid2=layerresid_[4].fpgastubid().str();
	}
	if (layerresid_[5].valid()) {
	  stubid3=layerresid_[5].fpgastubid().str();
	}
	if (diskresid_[0].valid()) {
	  stubid3=diskresid_[0].fpgastubid().str();
	}
	if (diskresid_[1].valid()) {
	  stubid2=diskresid_[1].fpgastubid().str();
	}
      }

      if (layer()==5) {
	if (layerresid_[0].valid()) {
	  stubid0=layerresid_[0].fpgastubid().str();
	}
	if (layerresid_[1].valid()) {
	  stubid1=layerresid_[1].fpgastubid().str();
	}
	if (layerresid_[2].valid()) {
	  stubid2=layerresid_[2].fpgastubid().str();
	}
	if (layerresid_[3].valid()) {
	  stubid3=layerresid_[3].fpgastubid().str();
	}
      }
    }

    if (isDisk()) {
      if (disk()==1) {
	if (layerresid_[0].valid()) {
	  stubid0=layerresid_[0].fpgastubid().str();
	}
	if (diskresid_[2].valid()) {
	  stubid1=diskresid_[2].fpgastubid().str();
	}
	if (diskresid_[3].valid()) {
	  stubid2=diskresid_[3].fpgastubid().str();
	}
	if (diskresid_[4].valid()) {
	  stubid3=diskresid_[4].fpgastubid().str();
	} else 	if (layerresid_[1].valid()) {
	  stubid3=layerresid_[1].fpgastubid().str();
	}
      }

      if (disk()==3) {
	if (layerresid_[0].valid()) {
	  stubid0=layerresid_[0].fpgastubid().str();
	}
	if (diskresid_[0].valid()) {
	  stubid1=diskresid_[0].fpgastubid().str();
	}
	if (diskresid_[1].valid()) {
	  stubid2=diskresid_[1].fpgastubid().str();
	}
	if (diskresid_[4].valid()) {
	  stubid3=diskresid_[4].fpgastubid().str();
	} else 	if (layerresid_[1].valid()) {
	  stubid3=layerresid_[1].fpgastubid().str();
	}
      }
      
    }
      
    if (isOverlap()) {
      if (layer()==1) {
	if (diskresid_[1].valid()) {
	  stubid0=diskresid_[1].fpgastubid().str();
	}
	if (diskresid_[2].valid()) {
	  stubid1=diskresid_[2].fpgastubid().str();
	}
	if (diskresid_[3].valid()) {
	  stubid2=diskresid_[3].fpgastubid().str();
	}
	if (diskresid_[4].valid()) {
	  stubid3=diskresid_[4].fpgastubid().str();
	}

      }

    }
      


    
    
    
    // real Q print out for fitted tracks
    if(writeoutReal){
    oss << (irinvfit_.value())*krinvpars<<" "
        << (iphi0fit_.value())*kphi0pars<<" "
        << (itfit_.value())*ktpars<<" "
        << (iz0fit_.value())*kz<<" "
      //<< ichisqfit_.str()<< "|"                            
        << innerFPGAStub_->phiregionaddressstr()<<" "
        << outerFPGAStub_->phiregionaddressstr()<<" "
	<< stubid0<<"|"
	<< stubid1<<"|"
	<< stubid2<<"|"
	<< stubid3;
    }
    //Binary print out
    if(!writeoutReal){
      oss << irinvfit_.str()<<"|"
	<< iphi0fit_.str()<<"|"
      //<< "xxxxxxxxxxx|"
	<< itfit_.str()<<"|"
	<< iz0fit_.str()<<"|"
      //<< ichisqfit_.str()<< "|"
	<< innerFPGAStub_->phiregionaddressstr()<<"|"
	<< outerFPGAStub_->phiregionaddressstr()<<"|"
	<< stubid0<<"|"
	<< stubid1<<"|"
	<< stubid2<<"|"
	<< stubid3;
    }
    return oss.str();
  }


  FPGATrack makeTrack(vector<L1TStub*> l1stubs) {
    assert(fit());
    
    if (l1stubs.size() == 0) {
      l1stubs = getL1Stubs(); // If fitter produced no stub list, take it from original tracklet.
    };
  

    FPGATrack tmpTrack(irinvfit_.value(),
                       iphi0fit_.value(),
                       itfit_.value(),
                       iz0fit_.value(),
                       ichisqfit_.value(),
                       chisqfit_,
                       getStubIDs(), // This needs fixing in case fitter updated stub list.
                       l1stubs,
                       seed());

    return tmpTrack;
    
  }

  FPGATrack* getTrack() {
    assert(fpgatrack_!=0);
    return fpgatrack_;
  }
  

  bool fit() const { return ichisqfit_.value()!=-1; }

  int layer() const {
    if (isDisk()) return 0; 
    if (isOverlap()) {
      //assert(outerStub_->layer()+1<7);
      return outerStub_->layer()+1;
    }
    assert(isBarrel());
    // assert(innerStub_->layer()+1<7);
    return innerStub_->layer()+1; 
  }

  int disk() const {
    if (innerFPGAStub_->isDisk()) return innerStub_->disk(); 
    if (outerFPGAStub_->isDisk()) return outerStub_->disk(); 
    return 0;
  }

  int disk2() const {
    if (innerStub_->disk()>0) {
      return innerStub_->disk()+1;
    }
    return innerStub_->disk()-1;
  }

  int overlap() const {
    return innerStub_->layer()+21;
  }

  int seed() const {
    // Returns integer code for tracklet seed.
    // Barrel: L1L2=1, L3L4=3, L5L6=6
    // Disk: D1D2=11, D3D4=13 (+/- for forward/backward disk)
    // Overlap: L1D1=21, L2D1=22 (+/- for forward/backward disk), L3L2=2
    if(barrel_) return innerStub_->layer()+1;
    if(disk_) {
        if(innerStub_->disk() < 0) return innerStub_->disk()-10;
        else return innerStub_->disk()+10;
    }
    if(overlap_) {
        if(innerFPGAStub_->isBarrel()) return 2;
        else if(innerStub_->disk() < 0) return -21-outerStub_->layer();
        else return 21+outerStub_->layer();
    }
    return 0;
  }

  bool isBarrel() const { 
    return barrel_;
  }

  bool isOverlap() const { 
    return overlap_;
  }

  int isDisk() const { 
    return disk_;
  }

  bool foundTrack(L1SimTrack simtrk){

    double deta=simtrk.eta()-asinh(itfit().value()*ktpars);
    double dphi=simtrk.phi()-(iphi0fit().value()*kphi0pars+phioffset_);

    if (dphi>0.5*two_pi) dphi-=two_pi;
    if (dphi<-0.5*two_pi) dphi+=two_pi;

    //if (match(1)==false) {
    //  return false; 
    //}
    
    bool found=(fabs(deta)<0.06)&&(fabs(dphi)<0.01);
    
    //cout << "deta dphi "<<deta<<" "<<dphi<<" "<<found<<endl;
    
    return found;

  }

  double phioffset() const {return phioffset_;}

  void setTrackletIndex(int index) {
    trackletIndex_=index;
    if (hourglass) {
      assert(index<128);
    } else {
      assert(index<64);
    }
  }

  int trackletIndex() const {return trackletIndex_;}

  void setTCIndex(int index) {TCIndex_=index;}

  int TCIndex() const {return TCIndex_;}
	
  int TCID() const {
    if (hourglass) {
      return TCIndex_*(1<<7)+trackletIndex_;
    } else {
      return TCIndex_*(1<<6)+trackletIndex_;
    }
  }
	
  unsigned int homeSector() const {return homeSector_;}
  
private:

  //Three types of tracklets... Overly complicated
  bool barrel_;
  bool disk_;
  bool overlap_;

  unsigned int homeSector_;
  
  double phioffset_;

  FPGAStub* innerFPGAStub_;
  FPGAStub* outerFPGAStub_;


  L1TStub* innerStub_;
  L1TStub* outerStub_;

  int trackletIndex_;
  int TCIndex_;

  //Tracklet track parameters  

  FPGAWord fpgarinv_;
  FPGAWord fpgaphi0_;
  FPGAWord fpgaz0_;
  FPGAWord fpgat_;

  double rinv_;
  double phi0_;
  double z0_;
  double t_;
  
  double rinvapprox_;
  double phi0approx_;
  double z0approx_;
  double tapprox_;

  double rproj_[6];
  double zprojdisk_[5];

  int      projlayer_[4];
  int      projdisk_[5];



  //Track  parameters from track fit
  
  FPGAWord irinvfit_;
  FPGAWord iphi0fit_;
  FPGAWord itfit_;
  FPGAWord iz0fit_;
  FPGAWord ichisqfit_;

  double rinvfit_;
  double phi0fit_;
  double tfit_;
  double z0fit_;
  double chisqfit_;

  double rinvfitexact_;
  double phi0fitexact_;
  double tfitexact_;
  double z0fitexact_;
  double chisqfitexact_;

  FPGATrack *fpgatrack_;


  FPGALayerProjection layerproj_[6];
  FPGADiskProjection diskproj_[5];

  FPGALayerResidual layerresid_[6];
  FPGADiskResidual diskresid_[6];

  
};



#endif
