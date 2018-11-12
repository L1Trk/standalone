//This class implementes the projection tranceiver
#ifndef FPGAMATCHTRANSCEIVER_H
#define FPGAMATCHTRANSCEIVER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAMatchTransceiver:public FPGAProcessBase{

public:

  FPGAMatchTransceiver(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="matchout1"||
	output=="matchout2"||
	output=="matchout3"||
	output=="matchout4"||
	output=="matchout5"){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      outputmatches_.push_back(tmp);
      return;
    }
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="proj1in"||
	input=="proj2in"||
	input=="proj3in"||
	input=="proj4in"||
	input=="proj5in"||
	input=="proj6in"||
	input=="proj7in"||
	input=="proj8in"||
	input=="proj9in"||
	input=="proj10in"||
	input=="proj11in"||
	input=="proj12in"||
	input=="proj13in"||
	input=="proj14in"||
	input=="proj15in"||
	input=="proj16in"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      inputmatches_.push_back(tmp);
      return;
    }

    assert(0);
  }

  std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > orderedMatches(vector<FPGAFullMatch*>& fullmatch) {

    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > tmp;

    std::vector<unsigned int> indexArray;
    for (unsigned int i=0;i<fullmatch.size();i++) {
      indexArray.push_back(0);
    }

    /*
    for (unsigned int i=0;i<fullmatch.size();i++) {
      for (unsigned int j=0;j<fullmatch[i]->nMatches();j++) {
	cout << fullmatch[i]->getName()<<" "
	     << fullmatch[i]->getFPGATracklet(j)->TCID()<<endl;
      }
    }
    */
    
    int bestIndex=-1;
    do {
      int bestTCID=(1<<16);
      bestIndex=-1;
      for (unsigned int i=0;i<fullmatch.size();i++) {
	if (indexArray[i]>=fullmatch[i]->nMatches()) {
	  //skip as we were at the end
	  continue;
	}
	int TCID=fullmatch[i]->getFPGATracklet(indexArray[i])->TCID();
	if (TCID<bestTCID) {
	  bestTCID=TCID;
	  bestIndex=i;
	}
      }
      if (bestIndex!=-1) {
	tmp.push_back(fullmatch[bestIndex]->getMatch(indexArray[bestIndex]));
	indexArray[bestIndex]++;
      }
    } while (bestIndex!=-1);

    //cout << "In FPGAFitTrack::orderedMatches : "<<tmp.size()<<endl;
    for (unsigned int i=0;i<tmp.size();i++) {
      //cout <<" "<<tmp[i]->TCID()<<endl;
      if (i>0) {
	//This allows for equal TCIDs. This means that we can e.g. have a track seeded
	//in L1L2 that projects to both L3 and D4. The algorithm will pick up the first hit and
	//drop the second
	if (tmp[i-1].first->TCID()>tmp[i].first->TCID()) {
	  cout << "Wrong TCID order in "<<getName()<<" "<<tmp[i-1].first->TCID()
				     <<" "<<tmp[i].first->TCID()<<endl;
	  //assert(0);
	}
      }
    }
    //cout << endl;

    return tmp;
  }
  
  //this->inputmatches_ to other->outputmatches_ 

  void execute(FPGAMatchTransceiver* other){
    int count=0;
    //assert(inputmatches_.size()==3);
	auto inputmatchesordered = orderedMatches(inputmatches_);
      //cout << "InputMatch name : "<<inputmatches_[i]->getName() <<endl;
      //cout << getName()<<" output match size : "<<other->outputmatches_.size()<<endl;
	for(unsigned int l=0;l<inputmatchesordered.size();l++){
  FPGATracklet* tracklet=inputmatchesordered[l].first;
  //cout << "FPGAMatchTransceiver "<<getName()<<" seed layer and disk : "<<tracklet->layer()<<" "<<tracklet->disk()<<endl;
  int nMatches=0;
  for(unsigned int j=0;j<other->outputmatches_.size();j++){
	assert(other->outputmatches_[j]!=0);
	//cout << "OutputMatch name : "<<outputmatches_[j]->getName() <<endl;
	string subname=outputmatches_[j]->getName().substr(3,4);
	//cout << "FPGAMatchTransceiver "<<getName()<<" target subname = "<<subname<<endl;
	if ((subname=="D1L1"&&tracklet->layer()==1&&abs(tracklet->disk())==1)||
		(subname=="D1L2"&&tracklet->layer()==2)||  //dangerous to only check layer
		(subname=="D1D2"&&tracklet->layer()==0&&abs(tracklet->disk())==1)||
		(subname=="D3D4"&&tracklet->layer()==0&&abs(tracklet->disk())==3)||
		(subname=="L1L2"&&tracklet->layer()==1&&tracklet->disk()==0)||
		(subname=="L3L4"&&tracklet->layer()==3&&tracklet->disk()==0)||
		(subname=="L5L6"&&tracklet->layer()==5&&tracklet->disk()==0)
		) {
	  other->outputmatches_[j]->addMatch(inputmatchesordered[l]);
	  count++;
	  nMatches++;
	}
  }
  assert(nMatches==1);
	}

    if (writeMatchTransceiver) {
      static ofstream out("matchtransceiver.txt");
      out << getName() << " " 
	  << count << endl;
    }
  }
  

private:

  vector<FPGAFullMatch*> inputmatches_;

  vector<FPGAFullMatch*> outputmatches_;

};

#endif
