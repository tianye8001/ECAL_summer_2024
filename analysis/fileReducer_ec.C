#include <iostream> 
#include <fstream>
#include <cmath> 
#include "math.h" 
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TArc.h"
#include "TString.h"
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"

using namespace std;

#include "analysis_tree_solid_ec.C"
#define MAX_CLUSTERS_PER_PLANE 2000
#define MAX_CHANNEL 16
// some numbers to be hard coded 
// make sure they are correct while using this script
//################################################################################################################################################## 

const double DEG=180./3.1415926;   //rad to degree

//#####################################################################################################################################################

//input:
//    infile: the path of input root file from GEMC_SOLID
//    numberofFile: how many number of files in the infile, usually 1.0E4 events per file,
//                  for beam on target root tree, 1.0E9 corresponds to 80 triggers, use 80 for it
//    key:  the string used to tell what kind of run it is
//    evgen: event type, 0 is beam on target, 1 is eDIS, 2 is eAll, 3 is bggen, 4 is even file 
//int fileReducer_beamtest_gem_survey(string inputfile_name,int numberOfFile=1, double event_actual=1, int evgen=1){
int fileReducer_ec(string inputfile_name,int numberOfFile=1, double event_actual=1, int evgen=1){

	char the_filename[500];
	sprintf(the_filename, "%s",inputfile_name.substr(0,inputfile_name.rfind(".")).c_str());
	TFile* outFile = new TFile(Form("%s_reduce_tree_analysis.root",the_filename), "RECREATE");
	const int t=1;
	TFile *file[t];
	TTree *tree_generated[t];
	TTree *tree_flux[t];
	TTree *tree_header[t];
	TTree *tree_solid_ec[t];
	TTree *tree_solid_ec_ps[t];
	TTree* outTree = new TTree("T", "simulation tree");
	cout<<"numberOfFile="<<numberOfFile<<"event_actual="<<event_actual<<"evgen="<<evgen<<endl;
	Float_t npe_hgc_total=0,npe_hgc_total_trigged=0;
	float px_gen = 0;
	float py_gen = 0;
	float pz_gen = 0;
	float vx_gen = 0;
	float vy_gen = 0;
	float vz_gen = 0;
	int pid_gen = 0;

	float p_gen=0, theta_gen=0, phi_gen=0;

	float rate = 0;
	float Q2 = 0;
	float rateRad = 0;
	float Npesum=0;
	float Cer[MAX_CHANNEL];
	int ecN=0;
	float PreShP,PreShP_e,PreShPx,PreShPy,PreShPz,PreShSum,PreShE, PreSh_l, PreSh_r, PreSh_t,PreShtheta,GEM00theta;
	float ShowerSum, Shower_l, Shower_r, Shower_t;

	outTree->Branch("rate",       &rate,       "rate/F"      ); //vx at vertex
	outTree->Branch("Q2",       &Q2,       "Q2/F"      ); //vx at vertex
	outTree->Branch("vx",       &vx_gen,       "vx/F"      ); //vx at vertex
	outTree->Branch("vy",       &vy_gen,       "vy/F"      ); //vy at vertex
	outTree->Branch("vz",       &vz_gen,       "vz/F"      ); //vz at vertex
	outTree->Branch("px",       &px_gen,       "px/F"      ); //px at vertex
	outTree->Branch("py",       &py_gen,       "py/F"      ); //py at vertex
	outTree->Branch("pz",       &pz_gen,       "pz/F"      ); //pz at vertex
	outTree->Branch("p",        &p_gen,        "p/F"       ); //ptot at vertex
	outTree->Branch("pid",      &pid_gen,        "pid/I"       ); //ptot at vertex
	// preshower tree
	outTree->Branch("PreShP",      &PreShP,   "PreShP/F"    ); //momentum hit on the virtual plane of presh 
	outTree->Branch("PreShP_e",      &PreShP_e,   "PreShP_e/F"    ); //momentum hit on the virtual plane of presh 
	outTree->Branch("PreShPx",      &PreShPx,   "PreShPx/F"    ); //momentum px hit on the virtual plane of presh
	outTree->Branch("PreShPy",      &PreShPy,   "PreShPy/F"    ); //momentum py hit on the virtual plane of presh
	outTree->Branch("PreShPz",      &PreShPz,   "PreShPz/F"    ); //momentum pz hit on the virtual plane of presh
	outTree->Branch("PreShtheta",      &PreShtheta,   "PreShtheta/F"    ); //momentum pz hit on the virtual plane of presh
	outTree->Branch("PreShSum",      &PreShSum,   "PreShSum/F"    ); //Sum of the deposit energy in three preshower modules
	outTree->Branch("PreShE",      &PreShE,   "PreShE/F"    ); //Sum of the hit energy in three preshower modules
	outTree->Branch("PreSh_l",      &PreSh_l,   "PreSh_l/F"    ); //Deposit energy in the left preshower module
	outTree->Branch("PreSh_r",      &PreSh_r,   "PreSh_r/F"    ); //Deposit energy in the right preshower module
	outTree->Branch("PreSh_t",      &PreSh_t,   "PreSh_t/F"    ); //Deposit energy in the top preshower module
	// shower tree
	outTree->Branch("ShowerSum",      &ShowerSum,   "ShowerSum/F"    ); //Sum of the deposit energy in three shower modules
	outTree->Branch("Shower_l",      &Shower_l,   "Shower_l/F"    ); //Deposit energy in the left shower module
	outTree->Branch("Shower_r",      &Shower_r,   "Shower_r/F"    ); //Deposit energy in the right shower module
	outTree->Branch("Shower_t",      &Shower_t,   "Shower_t/F"    ); //Deposit energy in the top shower module
	file[0]=new TFile(inputfile_name.c_str());
	long int nentries[t];
	vector <double> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0, *var9=0,*var10;
	vector <int> *gen_pid=0;
	vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	vector<double> *flux_id=0,*flux_hitn=0;
	vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0, *flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;

	for(int n=0;n<t;n++){

		tree_header[n] = (TTree*) file[n]->Get("userHeader");
		tree_header[n]->SetBranchAddress("userVar001",&var1);     //1
		tree_header[n]->SetBranchAddress("userVar002",&var2);     //x
		tree_header[n]->SetBranchAddress("userVar003",&var3);     //y
		tree_header[n]->SetBranchAddress("userVar004",&var4);     //W
		tree_header[n]->SetBranchAddress("userVar005",&var5);     //Q2
		tree_header[n]->SetBranchAddress("userVar006",&var6);     //rate
		tree_header[n]->SetBranchAddress("userVar007",&var7);     //radrate
		tree_header[n]->SetBranchAddress("userVar008",&var8);     //Ei, incoming beam energy after energy loss????
		tree_header[n]->SetBranchAddress("userVar009",&var9);     //Abeam
		tree_header[n]->SetBranchAddress("userVar010",&var10);     //Abeam

		tree_generated[n] = (TTree*) file[n]->Get("generated");
		tree_generated[n]->SetBranchAddress("pid",&gen_pid);   //particle ID
		tree_generated[n]->SetBranchAddress("px",&gen_px);     //momentum of the generated particle at target
		tree_generated[n]->SetBranchAddress("py",&gen_py);
		tree_generated[n]->SetBranchAddress("pz",&gen_pz);
		tree_generated[n]->SetBranchAddress("vx",&gen_vx);     //vertex of the generated particle at target
		tree_generated[n]->SetBranchAddress("vy",&gen_vy);
		tree_generated[n]->SetBranchAddress("vz",&gen_vz);

		tree_flux[n] = (TTree*) file[n]->Get("flux");
		tree_flux[n]->SetBranchAddress("hitn",&flux_hitn);     // hit number
		tree_flux[n]->SetBranchAddress("id",&flux_id);         //hitting detector ID
		tree_flux[n]->SetBranchAddress("pid",&flux_pid);       //pid
		tree_flux[n]->SetBranchAddress("mpid",&flux_mpid);     // mother pid
		tree_flux[n]->SetBranchAddress("tid",&flux_tid);       // track id
		tree_flux[n]->SetBranchAddress("mtid",&flux_mtid);     // mother track id
		tree_flux[n]->SetBranchAddress("otid",&flux_otid);     // original track id
		tree_flux[n]->SetBranchAddress("trackE",&flux_trackE);   //track energy of 1st step,  track here is G4 track
		tree_flux[n]->SetBranchAddress("totEdep",&flux_totEdep); //totEdep in all steps, track here is G4 track
		tree_flux[n]->SetBranchAddress("avg_x",&flux_avg_x);     //average x, weighted by energy deposition in each step
		tree_flux[n]->SetBranchAddress("avg_y",&flux_avg_y);     //average y
		tree_flux[n]->SetBranchAddress("avg_z",&flux_avg_z);     //average z
		tree_flux[n]->SetBranchAddress("avg_lx",&flux_avg_lx);   // local average x
		tree_flux[n]->SetBranchAddress("avg_ly",&flux_avg_ly);   // local average y
		tree_flux[n]->SetBranchAddress("avg_lz",&flux_avg_lz);   // local average z
		tree_flux[n]->SetBranchAddress("px",&flux_px);          // px of 1st step
		tree_flux[n]->SetBranchAddress("py",&flux_py);          // px of 1st step
		tree_flux[n]->SetBranchAddress("pz",&flux_pz);          // px of 1st step
		tree_flux[n]->SetBranchAddress("vx",&flux_vx);          // x coordinate of 1st step
		tree_flux[n]->SetBranchAddress("vy",&flux_vy);          // y coordinate of 1st step
		tree_flux[n]->SetBranchAddress("vz",&flux_vz);          // z coordinate of 1st step
		//information recorded by ec
		tree_solid_ec[n]= (TTree*) file[n]->Get("solid_ec");
		tree_solid_ec_ps[n]= (TTree*) file[n]->Get("solid_ec_ps");
		setup_tree_solid_ec(tree_solid_ec[n]);	
		setup_tree_solid_ec_ps(tree_solid_ec_ps[n]);	
		TRandom3 rand;
		rand.SetSeed(0);

		int sensor_good=0;
		int event_good=0,event_trig_good=0;
		// 	long int N_events = (long int)tree_header->GetEntries();
		nentries[n] = (long int)tree_generated[n]->GetEntries();	
		printf("Entries = %li \n",nentries[n]);

		cout << "total number of events : " << nentries[n] << endl;	

		//----------------------------
		//      loop trees
		//---------------------------
		double Ek=0;
		double Ec=0;
		double trigger_ec=0;
		int hit_id=-1,pid_max=0,mpid_max=0;
		double px_max=0,py_max=0,pz_max=0,p_max=0,p_max_e=0,theta_max=0,theta_GEM00=0;
		double Eend_ec_Esum=0;
		double Eend_ec_ps_Esum=0;
		int N_preshower=0;
		int pid_gen1=0;
		double hit_pf=0;
		double edep_6p1_max= 0;
		int Is_prime=0;					  
		int Is_ECback=-1;			
		double dE=0;		  
		bool writeFlag = false;
		for(long int i=0;i<nentries[n];i++){	  		
			cout<<"event " << i << "\r";
			if (evgen==2){
				tree_header[n]->GetEntry(i);
				Q2 = var5->at(0);
				rate=var6->at(0)/numberOfFile; // new eDIS and eAll generator
			}
			else if (evgen==0) {
				rate=40e-6/1.6e-19/event_actual;  //beamOnTarget  file I= 5 uA	 	  
			}
			else if (evgen==4) rate=1; //even event file
			else if (evgen==3){         // bggen file
				tree_header[n]->GetEntry(i);
				rate=var10->at(0)/numberOfFile; 
			}
			else {
				cout << "Not right filemode" << endl;    
				return 0; 
			}
			tree_generated[n]->GetEntry(i);
			for (std::size_t j=0;j<gen_pid->size() ;j++) {
				pid_gen= gen_pid->at(j);
				px_gen=gen_px->at(j)/1e3;
				py_gen=gen_py->at(j)/1e3;
				pz_gen=gen_pz->at(j)/1e3;
				vx_gen=gen_vx->at(j)*0.1;
				vy_gen=gen_vy->at(j)*0.1;
				vz_gen=gen_vz->at(j)*0.1;
				p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
				theta_gen=TMath::ACos(pz_gen/p_gen)*DEG;
				//		cout<<"event="<<i<<"j="<<j<<"gen_p="<<p_gen<<endl; 
			}
			//---
			tree_flux[n]->GetEntry(i);		  
			int sec_hgc=0;		
			int Is_trig=0;					  
			Is_prime=0;					  
			Is_ECback=-1;					  
			hit_id=-1;
			double Eec=0,Eec_photonele=0,Eec_ele=0,EdepSC_D=0,EdepSC_C=0;
			px_max = 0;
			theta_max = 0;
			py_max = 0;
			pz_max = 0;
			p_max = 0;
			p_max_e = 0;
			pid_max = 0;
			mpid_max = 0;
			edep_6p1_max=0;
			Eend_ec_ps_Esum=0;
			//Cer_n=0;
			Eend_ec_Esum=0;
			N_preshower=0;
			//if(Q2>0.01){
			//                            cout<<"event="<<i<<"flux_size="<<flux_hitn->size()<<endl;
			for (std::size_t j=0;j<flux_hitn->size() ;j++) {
				if(flux_pz->at(j)>0){  
					Is_prime ++;
					//check tid tid==1 prime particle
					//	if(flux_tid->at(j) !=1) continue;
					hit_pf=sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j));  //MeV to GeV
					double hit_th=atan2(sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)),flux_pz->at(j))*DEG;
					//double E=flux_trackE->at(j)/1e3;		  
					double E=flux_trackE->at(j);		  
					if(flux_id->at(j)==1) hit_id=0; // front of prelead
					else if(flux_id->at(j)==2) hit_id=1; // EC back
					else cout << "wrong flux_id" << flux_id->at(j) << endl;
					if(hit_id==0 ){
						//if( flux_tid->at(j)==1 && hit_id==1){

						//Is_prime=1;
						N_preshower ++;
						if (E >= edep_6p1_max){
							edep_6p1_max=E;
							p_max=hit_pf;
							px_max=flux_px->at(j);
							py_max=flux_py->at(j);
							pz_max=flux_pz->at(j);
							//theta_max = hit_th;
							theta_max = theta_gen;
							pid_max = flux_pid->at(j);
							mpid_max = flux_mpid->at(j);
						}
						Eend_ec_ps_Esum += E;
					} // hit id cut 
					}//pz>0
				}	// end of flux		

				// process ec
				tree_solid_ec[n]->GetEntry(i);
				tree_solid_ec_ps[n]->GetEntry(i);
				dE=gRandom->Gaus(1.0,sqrt(pow(0.04149,2)+pow(0.07938,2)/7.0));
				//dE=1.0;
				double Eend_ec_sum=0;
				double Eend_ec_ps_sum=0;	
				double Eend_ec[4]={0};
				double Eend_ec_ps[4]={0};
				process_tree_solid_ec(tree_solid_ec[n],tree_solid_ec_ps[n],Eend_ec_sum,Eend_ec_ps_sum,Eend_ec,Eend_ec_ps);
				PreShSum=Eend_ec_ps_sum;
				PreShE=Eend_ec_ps_Esum;
				PreSh_l=Eend_ec_ps[2];
				PreSh_r=Eend_ec_ps[3];
				PreSh_t=Eend_ec_ps[1];
				ShowerSum=Eend_ec_sum*dE;
				Shower_l=Eend_ec[2]*dE;
				Shower_r=Eend_ec[3]*dE;
				Shower_t=Eend_ec[1]*dE;
				PreShP=edep_6p1_max;
				PreShPx=px_max;
				PreShPy=py_max;
				PreShPz=pz_max;
				PreShtheta=theta_max;
				outTree->Fill();
				Is_prime=0;
			} //end event loop
		}//file loop
		outFile->cd();
		outTree->Write();
		outFile->Close();

		return 0;
		}
