#ifndef _RUNPOSTPROCESSES_H
#define _RUNPOSTPROCESSES_H

#include <AnalysisReader.h>//This includes the classes to read the tree or copy a tree
#include <TransverseTools.h>

#include <cassert>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include <TVector3.h>

using namespace std;

string test_file = "/pnfs/minerva/persistent/users/dcoplowe/CC1P1Pi_PL13C_200317/grid/central_value/minerva/ana/v10r8p9/00/01/32/60/SIM_minerva_00013260_Subruns_0001-0002-0003-0004_CC1P1PiAnalysis_Ana_Tuple_v10r8p9-dcoplowe.root";

class RunPostProcesses {

public:
	RunPostProcesses(std::string infilemame, std::string outfilename = "some_generic_name.root", std::string rec_tree = "sel");
	~RunPostProcesses();

	void Analyse();

private:

	TFile * m_infile;
	AnalysisReader * m_reader;
	Int_t m_entries;
	TTree * m_rec;

	TFile * m_outfile;
	TTree * m_outtree;
	void CopyTree(std::string treename);

	void SetOutTree();
	void FillOutTree();

	TVector3 * GetTransVarsRec(double vtx[], double mumom[], double prmom[], double pimom[], double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT);

	TVector3 * GetTransVarsSim(double vtx[], double pdp[], double mumom[], double prmom[], double pimom[], double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT);// const;

	TransverseTools * m_TransTools;
};

#endif

RunPostProcesses::RunPostProcesses(std::string infilemame, std::string outfilename, std::string rec_tree) {

	m_infile = new TFile(infilemame.c_str(), "READ");
	assert(m_infile);

	if(!m_infile->IsOpen()){
		cout << "RunPostProcesses : Error : Could not open file named: " << infilemame << endl;
	}

	TTree * m_rec = static_cast<TTree*>(m_infile->Get( rec_tree.c_str() ) );
	assert(m_rec);

	m_outfile = new TFile(outfilename.c_str(), "RECREATE");
	m_outfile->cd();
	m_outtree = new TTree(rec_tree.c_str(),"");

	m_reader = new AnalysisReader(m_rec, m_outtree);
	m_entries = m_reader->GetEntries();

	m_TransTools = new TransverseTools();
}

RunPostProcesses::~RunPostProcesses(){
	// if(m_outtree) delete m_outtree;
	cout << "Working 1" << endl;	
	if(m_rec) delete m_rec;
	cout << "Working 2" << endl;	
	if(m_infile->IsOpen()) m_infile->Close();
	if(m_infile) delete m_infile;
	cout << "Working 3" << endl;	
	if(m_outfile->IsOpen()) m_outfile->Close();
	if(m_outfile) delete m_outfile;
	cout << "Working 4" << endl;	
	delete m_TransTools;
}

void RunPostProcesses::Analyse(){

	int percent = m_entries/20;

	m_outfile->cd();

	m_reader->SetOutTree();

	cout << "Starting to reprocess transverse variables" << endl;

	for(int ev = 0; ev < m_entries; ev++){
		if(ev % percent == 0) cout << Form("Reprocessed %.f%%", (double)100*ev/m_entries ) << endl;
		m_reader->GetEntry(ev);
		if(m_reader->accum_level[0] > 5 || m_reader->accum_level[1] > 5){
			(void)GetTransVarsRec(m_reader->sel_VTX, m_reader->sel_mu_4mom, m_reader->sel_pi_EX_4mom, m_reader->sel_pi_EX_4mom, m_reader->sel_dpTT_EX, m_reader->sel_dpT_EX, m_reader->sel_dalphaT_EX, m_reader->sel_dphiT_EX);

			// (void)GetTransVarsSim(m_reader->sel_trueVTX, m_reader->sel_truePDP, m_reader-> mumom, m_reader-> prmom, m_reader->  pimom, m_reader->  dpTT, m_reader->  dpTMag, m_reader-> dalphaT, m_reader-> dphiT);
		
					// sel_dpTT_EX;
   					// sel_dpTT_EX_tmumom;
   					// sel_dpTT_EX_tnudir;
   					// sel_dpTT_EX_tpimom;
   					// sel_dpTT_EX_tprmom;
   					// sel_dpTT_LL;
   					// sel_dpTT_LL_tmumom;
   					// sel_dpTT_LL_tnudir;
   					// sel_dpTT_LL_tpimom;
   					// sel_dpTT_LL_tprmom;
   					// sel_dpTT_pi_EX;
   					// sel_dpTT_pi_EX_tmumom;
   					// sel_dpTT_pi_EX_tnudir;
   					// sel_dpTT_pi_EX_tpimom;
   					// sel_dpTT_pi_EX_tprmom;
   					// sel_dpTT_pi_LL;
   					// sel_dpTT_pi_LL_tmumom;
   					// sel_dpTT_pi_LL_tnudir;
   					// sel_dpTT_pi_LL_tpimom;
   					// sel_dpTT_pi_LL_tprmom;
   					// sel_dpTT_pi_dir_EX;
   					// sel_dpTT_pi_dir_EX_tmumom;
   					// sel_dpTT_pi_dir_EX_tnudir;
   					// sel_dpTT_pi_dir_EX_tpidir;
   					// sel_dpTT_pi_dir_EX_tprmom;
   					// sel_dpTT_pi_dir_LL;
   					// sel_dpTT_pi_dir_LL_tmumom;
   					// sel_dpTT_pi_dir_LL_tnudir;
   					// sel_dpTT_pi_dir_LL_tpidir;
   					// sel_dpTT_pi_dir_LL_tprmom;
   					// sel_dpTT_pr_EX;
   					// sel_dpTT_pr_EX_tmumom;
   					// sel_dpTT_pr_EX_tnudir;
   					// sel_dpTT_pr_EX_tpimon;
   					// sel_dpTT_pr_EX_tprmom;
   					// sel_dpTT_pr_LL;
   					// sel_dpTT_pr_LL_tmumom;
   					// sel_dpTT_pr_LL_tnudir;
   					// sel_dpTT_pr_LL_tpimon;
   					// sel_dpTT_pr_LL_tprmom;
   					// sel_dpTT_pr_dir_EX;
   					// sel_dpTT_pr_dir_EX_tmumom;
   					// sel_dpTT_pr_dir_EX_tnudir;
   					// sel_dpTT_pr_dir_EX_tpimom;
   					// sel_dpTT_pr_dir_EX_tprdir;
   					// sel_dpTT_pr_dir_LL;
   					// sel_dpTT_pr_dir_LL_tmumom;
   					// sel_dpTT_pr_dir_LL_tnudir;
   					// sel_dpTT_pr_dir_LL_tpimom;
   					// sel_dpTT_pr_dir_LL_tprdir;

		}
		m_reader->FillOutTree();
	}
	cout << "********* Finished Processing *********" << endl;
	cout << "Copying other trees." << endl;
	m_outfile->cd();
	m_outtree->Write();
	CopyTree("Truth");
	CopyTree("Meta");
	m_outfile->Close();
	cout << "Reanalysis Complete" << endl;
}

void RunPostProcesses::CopyTree(std::string treename){
	bool good_files = true;

	if(!m_infile->IsOpen()){ 
		cout <<  "RunPostProcesses::CopyTree : Error : In file not open." << endl;
		good_files = false; 
	}
	if(!m_outfile->IsOpen()){ 
		cout <<  "RunPostProcesses::CopyTree : Error : Out file not open." << endl;
		good_files = false;
	}

	if(good_files){
		cout << "Starting to copy tree " << treename << "...";
		TTree * tree = static_cast<TTree*>(m_infile->Get( treename.c_str() ));
		assert(tree);
		m_outfile->cd();
		TTree * tree_copy = tree->CopyTree("");
		tree_copy->Write();
		cout << " Successful" << endl;
	}
}

TVector3 * RunPostProcesses::GetTransVarsRec(double vtx[], double mu[], double pr[], double pi[], double &dpTT, double &dpTMag, double &dalphaT, double &dphiT)
{
	const TVector3 * mumom = new TVector3( mu[1], mu[2], mu[3] );
	const TVector3 * prmom = new TVector3( pr[1], pr[2], pr[3] );
	const TVector3 * pimom = new TVector3( pi[1], pi[2], pi[3] );
	return m_TransTools->GetTransVarsRec(vtx, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT);
}

TVector3 * RunPostProcesses::GetTransVarsSim(double vtx[], double pdp[], double mu[], double pr[], double pi[], double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT)
{
	const TVector3 * mumom = new TVector3( mu[1], pr[2], pi[3] );
	const TVector3 * prmom = new TVector3( pr[1], pr[2], pr[3] );
	const TVector3 * pimom = new TVector3( pi[1], pi[2], pi[3] );
	return m_TransTools->GetTransVarsSim(vtx, pdp, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT);
}

int main(int argc, char *argv[])
{

	string ifile = test_file;
	string ofile = "some_generic_name.root";
	string stree = "sel";

	char cc;
	while ((cc = getopt(argc, argv, "i:o:t:")) != -1) {
		switch (cc){
			case 'i': ifile = string(optarg); break;
			case 'o': ofile = string(optarg); break;
			case 't': stree = string(optarg); break;
			default: return 1;
		}
	}

	RunPostProcesses * Run = new RunPostProcesses(ifile,ofile,stree);
	Run->Analyse();
	delete Run;

}


