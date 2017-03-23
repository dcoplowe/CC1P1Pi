#ifndef _RUNPOSTPROCESSES_H
#define _RUNPOSTPROCESSES_H

#include <AnalysisReader.h>//This includes the classes to read the tree or copy a tree
#include <TransverseTools.h>

#include <cassert>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

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
	if(m_rec) delete m_rec;
	
	if(m_infile->IsOpen()) m_infile->Close();
	if(m_infile) delete m_infile;

	if(m_outtree) delete m_outtree;

	if(m_outfile->IsOpen()) m_outfile->Close();
	if(m_outfile) delete m_outfile;

	delete m_TransTools;
}

void RunPostProcesses::Analyse(){

	m_outfile->cd();

	m_reader->SetOutTree();

	for(int ev = 0; ev < m_entries; ev++){
		m_reader->GetEntry(ev);
		cout << " Pre: dpTT = " << m_reader->sel_dpTT_pi_EX << endl;

		m_reader->sel_dpTT_pi_EX = 8008135.;

		cout << "Post: dpTT = " << m_reader->sel_dpTT_pi_EX << endl;

		m_reader->FillOutTree();

	}

	m_outfile->cd();
	m_outtree->Write();

	CopyTree("Truth");
	CopyTree("Meta");
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

int main(int argc, char *argv[])
{

	string infile = test_file;

	char cc;
	while ((cc = getopt(argc, argv, "i:o:t:")) != -1) {
		switch (cc){
			case 'i': break;
			case 'o': break;
			case 't': break;
			default: return 1;
		}
	}

	RunPostProcesses * Run = new RunPostProcesses(test_file,"some_generic_name.root","sel");
	Run->Analyse();
	delete Run;

}


