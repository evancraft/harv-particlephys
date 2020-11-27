float calculateZ(float E, float Y) {

	float exp = TMath::Exp(2*Y);

	float num = exp - 1;

	float denom = exp + 1;

	float pz = (E*num)/denom;

	return pz;

}

TH1F * graphmmt(const char* filename, const char* branchname, const char* branchname2, const char* branchname3, const char* branchname4, const char* title, int bins, int lower, int upper, short color){

	// declare the histogram

	TH1F *hist = new TH1F("hist", "Mass Spectrum (Formula 1)", bins, lower, upper);

	// declare the ttree

	TFile *input = new TFile(filename, "read");
	TTree *tree = (TTree*)input -> Get("mini");

	// set the pt branch

	TBranch *branch = tree -> GetBranch(branchname);
	std::vector<float> *bivector = 0;
	tree -> SetBranchAddress(branchname, &bivector);

	// set the E branch

	TBranch *branch2 = tree -> GetBranch(branchname2);
	std::vector<float> *bivector2 = 0;
	tree -> SetBranchAddress(branchname2, &bivector2);

	// set the phi branch

	TBranch *branch3 = tree -> GetBranch(branchname3);
	std::vector<float> *bivector3 = 0;
	tree -> SetBranchAddress(branchname3, &bivector3);

	// set the eta branch

	TBranch *branch4 = tree -> GetBranch(branchname4);
	std::vector<float> *bivector4 = 0;
	tree -> SetBranchAddress(branchname4, &bivector4);

	// iterate over the branch values and fill the histogram
	// each branch entry is a std:bivector <float>

	int entries = tree -> GetEntries();
	float pt1, pt2, pt3, pt4, dotproduct, mag1, mag2;
	float E, E1, E2, E3, E4, eta1, eta2, eta3, eta4;
	float cosdelta, factor, mass;
	float px1, py1, pz1;
	float px2, py2, pz2;
	float px3, py3, pz3;
	float px4, py4, pz4;
	float phi1, phi2, phi3, phi4;
	float hx, hy, hz;
	float magsquared;

	for(int j = 0; j < entries; j++)
	{

	branch -> GetEntry(j);
	branch2 -> GetEntry(j);
	branch3 -> GetEntry(j);
	branch4 -> GetEntry(j);

	pt1 = bivector -> at(0);
	pt2 = bivector -> at(1);
	pt3 = bivector -> at(2);
	pt4 = bivector -> at(3);

	E1 = bivector2 -> at(0);
	E2 = bivector2 -> at(1);
	E3 = bivector2 -> at(2);
	E4 = bivector2 -> at(3);

	phi1 = bivector3 -> at(0);
	phi2 = bivector3 -> at(1);
	phi3 = bivector3 -> at(2);
	phi4 = bivector3 -> at(3);

	eta1 = bivector4 -> at(0);
	eta2 = bivector4 -> at(1);
	eta3 = bivector4 -> at(2);
	eta4 = bivector4 -> at(3);

	px1 = pt1*TMath::Cos(phi1);
	py1 = pt1*TMath::Sin(phi1);
	pz1 = calculateZ(E1, eta1);

	px2 = pt2*TMath::Cos(phi2);
	py2 = pt2*TMath::Sin(phi2);
	pz2 = calculateZ(E2, eta2);

	px3 = pt3*TMath::Cos(phi3);
	py3 = pt3*TMath::Sin(phi3);
	pz3 = calculateZ(E3, eta3);

	px4 = pt4*TMath::Cos(phi4);
	py4 = pt4*TMath::Sin(phi4);
	pz4 = calculateZ(E4, eta4);

	hx = px1 + px2 + px3 + px4;
	hy = py1 + py2 + py3 + py4;
	hz = pz1 + pz2 + pz3 + pz4;

	E = E1 + E2 + E3 + E4;
	
	magsquared = hx*hx + hy*hy + hz*hz;
	mass = TMath::Sqrt(E*E - magsquared);

	//std::cout << pz1 << " " << pz2 << endl;
	
	//if(mass>60000)
	hist -> Fill(mass/1000);

	}
	
	// draw the histogram

	hist -> SetFillColor(color);

	gStyle->SetOptTitle(0);

	hist -> GetXaxis() -> SetTitle("#it{m} [Gev]");
	hist -> GetYaxis() -> SetTitle("Events");

	Double_t fact = 1.;
	hist->Scale(fact/hist->GetEntries());

	return hist;
}

void project4p2(){

	TH1F *dilepton1 = 0;
	TH1F *dilepton2 = 0;
	TH1F *dilepton3 = 0;
	TH1F *dilepton4 = 0;
	TH1F *monte = 0;

	TCanvas *c1 = new TCanvas();
	
	dilepton1 = graphmmt("mc_345060.ggH125_ZZ4lep.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "Inv. Mass", 100, 0, 200, kYellow-6);

	dilepton1 -> Draw("HIST");

	TCanvas *c2 = new TCanvas();
	
	dilepton2 = graphmmt("mc_363490.llll.VV.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "Inv. Mass", 100, 0, 200, kYellow-6);

	dilepton2 -> Draw("HIST");

	TCanvas *c3 = new TCanvas();
	
	dilepton3 = graphmmt("mc_ZeeZmumu.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "Inv. Mass", 100, 0, 200, kYellow-6);

	dilepton3 -> Draw("HIST");

	TCanvas *c4 = new TCanvas();
	
	dilepton4 = graphmmt("mc_410000.ttbar_lep.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "Inv. Mass", 100, 0, 200, kYellow-6);

	dilepton4 -> Draw("HIST");


	//TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.7);
	//leg -> AddEntry(dilepton, "Higgs", "f");
	//leg -> Draw();

}