
#include "../include/SingleTopLHEAnalyzer.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>

using namespace std;


int main()
{

  int nSamples = 17;
  int modes = 1;
  string* suffix = new string[nSamples];
  string* prefix = new string[modes];

prefix[0]="madgraph_";
//prefix[1]="madspin_";

suffix[0] = "SM";
suffix[1] = "cbwi_m1";
suffix[2] = "cbwi_m2";
suffix[3] = "cbwi_p1";
suffix[4] = "cbwi_p2";
suffix[5] = "cptbi_m1";
suffix[6] = "cptbi_m2";
suffix[7] = "cptbi_p1";
suffix[8] = "cptbi_p2";
suffix[9] = "ctw_m1";
suffix[10] = "ctw_m2";
suffix[11] = "ctw_p1";
suffix[12] = "ctw_p2";
suffix[13] = "ctwi_m1";
suffix[14] = "ctwi_m2";
suffix[15] = "ctwi_p1";
suffix[16] = "ctwi_p2";



  TFile** fInput = new TFile*[nSamples];
  TTree** tInput = new TTree*[nSamples];
  SingleTopLHEAnalyzer** singleTopLHEAnalyzer = new SingleTopLHEAnalyzer*[nSamples];

  for (int j=0; j<modes; j++)
  {
    for (int i=0; i<nSamples; i++)
    {
      cout<<"Treating"+prefix[j] + suffix[i]+".root"<<endl;
      string inputName = "data/madgraph/"+ prefix[j] + suffix[i] + ".root";
      string outputName = "output_" +prefix[j] + suffix[i] +  ".root";

      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("LHEF");
      singleTopLHEAnalyzer[i] = new SingleTopLHEAnalyzer(tInput[i]);
      singleTopLHEAnalyzer[i]->Loop();
      string commandline = "mv output.root data/madgraph/output/" + outputName;
      system(commandline.c_str());
    }
  }



return 0;
}
