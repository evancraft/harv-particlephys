TH1F * graphmmt(const char* filename, const char* branchname, const char* branchname2, const char* branchname3, const char* branchname4, const char* title, int bins, int lower, int upper, short color, std::string type){

	// declare the histogram

	TH1F *hist = new TH1F("hist", title, bins, lower, upper);

	// declare the ttree

	TFile *input = new TFile(filename, "read");
	TTree *tree = (TTree*)input -> Get("mini");

	// set the E branch

	TBranch *branch = tree -> GetBranch(branchname);
	std::vector<float> *nrg =  0;
	tree -> SetBranchAddress(branchname, &nrg);

	// set the met branch

	TBranch *branch2 = tree -> GetBranch(branchname2);
	float mtd;
	tree -> SetBranchAddress(branchname2, &mtd);

	// set the phi branch

	TBranch *branch3 = tree -> GetBranch(branchname3);
	std::vector<float> *phd =  0;
	tree -> SetBranchAddress(branchname3, &phd);

	// set the phi2 branch

	TBranch *branch4 = tree -> GetBranch(branchname4);
	float mpi;
	tree -> SetBranchAddress(branchname4, &mpi);

	// iterate over the branch values and fill the histogram
	// each branch entry is a std:bivector <float>

	int entries = tree -> GetEntries();
	float E, met, phi, phi2, delta;

	for(int j = 0; j < entries; j++)
	{

	branch -> GetEntry(j);
	branch2 -> GetEntry(j);
	branch3 -> GetEntry(j);
	branch4 -> GetEntry(j);

	E = nrg -> at(0);
	met = mtd;
	phi = phd -> at(0);
	phi2 = mpi;

	delta = TMath::Abs(phi2-phi);
	
	//if(mass>60000)
	hist -> Fill(TMath::Sqrt(2*E*met*(1-TMath::Cos(delta)))/1000);
	//hist -> Fill(TMath::Abs(phi2-phi));
	//hist -> Fill(delta);
	//hist -> Fill(E);

	}
	
	// draw the histogram

	gStyle->SetOptTitle(0);

	hist -> SetFillColor(color);
	hist -> GetXaxis() -> SetTitle("#it{m}_{T} [Gev]");
	hist -> GetYaxis() -> SetTitle("Events");
	Double_t fact = 1.;
	hist ->Scale(fact/hist->GetEntries());
	return hist;
	
}

void project4p1(){

	TCanvas *c1 = new TCanvas();
	TH1F *lepton = 0;
	
	lepton = graphmmt("data_A.1lep.small.root", "lep_E", "met_et", "lep_phi", "met_phi", "Transverse Mass", 100, 0, 200, kBlue-9, "met");

	lepton -> Draw("HIST");

	//TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.7);
	//leg -> AddEntry(lepton, "Lepton", "f");
	//leg -> Draw();

}