float calculateZ(float E, float Y) {

	float exp = TMath::Exp(2*Y);

	float num = exp - 1;

	float denom = exp + 1;

	float pz = (E*num)/denom;

	return pz;

}

TH1F * graphmmt(const char* filename, const char* branchname, const char* branchname2, const char* branchname3, const char* branchname4, const char* branchname5, const char* branchname6, const char* branchname7, const char* branchname8, const char* branchname9, const char* title, int bins, int lower, int upper, short color, bool weighted, int cat){

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

	// set the muon branch

	TBranch *branch5 = tree -> GetBranch(branchname5);
	float mb;
	tree -> SetBranchAddress(branchname5, &mb);

	// set the ele branch

	TBranch *branch6 = tree -> GetBranch(branchname6);
	float eleb;
	tree -> SetBranchAddress(branchname6, &eleb);

	// set the pileup branch

	TBranch *branch7 = tree -> GetBranch(branchname7);
	float plb;
	tree -> SetBranchAddress(branchname7, &plb);

	// set the trigger branch

	TBranch *branch8 = tree -> GetBranch(branchname8);
	float trb;
	tree -> SetBranchAddress(branchname8, &trb);

	// set the mc branch

	TBranch *branch9 = tree -> GetBranch(branchname9);
	float mcb;
	tree -> SetBranchAddress(branchname9, &mcb);

	// set the tight branch

	TBranch *branch10 = tree -> GetBranch("lep_isTightID");
	std::vector<bool> *tb = 0;
	tree -> SetBranchAddress("lep_isTightID", &tb);

	// set the charge branch

	TBranch *branch11 = tree -> GetBranch("lep_charge");
	std::vector<int> *cb = 0;
	tree -> SetBranchAddress("lep_charge", &cb);

	// set the type branch

	TBranch *branch12 = tree -> GetBranch("lep_type");
	std::vector<int> *tyb = 0;
	tree -> SetBranchAddress("lep_type", &tyb);

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
	float muon, ele, pileup, trigger, mc;
	float weight;
	float tight1, tight2, tight3, tight4;
	float charge1, charge2, charge3, charge4, chargesum;
	float type1, type2, type3, type4, typesum;
	std::vector<float> *weights = new vector<float>(entries);

	float sum = 0;

	for(int j = 0; j < entries; j++)
	{

	branch -> GetEntry(j);
	branch2 -> GetEntry(j);
	branch3 -> GetEntry(j);
	branch4 -> GetEntry(j);

	branch5 -> GetEntry(j);
	branch6 -> GetEntry(j);
	branch7 -> GetEntry(j);
	branch8 -> GetEntry(j);
	branch9 -> GetEntry(j);

	branch10 -> GetEntry(j);
	branch11 -> GetEntry(j);
	branch12 -> GetEntry(j);

	tight1 = tb -> at(0);
	tight2 = tb -> at(1);
	tight3 = tb -> at(2);
	tight4 = tb -> at(3);

	charge1 = cb -> at(0);
	charge2 = cb -> at(1);
	charge3 = cb -> at(2);
	charge4 = cb -> at(3);

	type1 = tyb -> at(0);
	type2 = tyb -> at(1);
	type3 = tyb -> at(2);
	type4 = tyb -> at(3);

	chargesum = charge1 + charge2 + charge3 + charge4;
	typesum = type1 + type2 + type3 + type4;

	muon = mb;
	ele = eleb;
	pileup = plb;
	trigger = trb;
	mc = mcb;

	weight = muon*ele*pileup*trigger*mc;
	weight = weight*2.16*TMath::Power(10,-6);
	weights -> at(j) = weight;
	sum = sum + weight;

	//std::cout << sum << endl;

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

	if(((mass/1000) > 100) && ((mass/1000) < 135)){
	
	if(weighted && (chargesum == 0) && (typesum != 46) && (typesum != 50)){
	hist -> Fill(mass/1000, weight);
	} else if (cat == 2) {
	//hist -> Fill(mass/10000, .000005);
	hist -> Fill(mass/1000, .0005);
	} else if ((chargesum == 0) && (typesum != 46) && (typesum != 50)) {
	//hist -> Fill(mass/10000, .000005);
	hist -> Fill(mass/1000, 1.02);
	//std::cout << chargesum << endl;
	//std::cout << typesum << endl;
	}

	}

	}
	
	// draw the histogram

	hist -> SetFillColor(color);
	gStyle->SetOptTitle(0);

	hist -> GetXaxis() -> SetTitle("#it{m} [Gev]");
	hist -> GetYaxis() -> SetTitle("Events");

	if(0==1){
	Double_t fact = 1.;
	hist->Scale(fact/hist->GetEntries());
	}

	return hist;
}

void project4weights(){

	TH1F *dilepton = 0;
	TH1F *monte = 0;
	TH1F *data = 0;
	
	dilepton = graphmmt("mc_345060.ggH125_ZZ4lep.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "scaleFactor_MUON", "scaleFactor_ELE", "scaleFactor_PILEUP", "scaleFactor_LepTRIGGER", "mcWeight", "Inv. Mass", 50, 0, 200, kRed-9, true, 1);

	monte = graphmmt("mc_345060.ggH125_ZZ4lep.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "scaleFactor_MUON", "scaleFactor_ELE", "scaleFactor_PILEUP", "scaleFactor_LepTRIGGER", "mcWeight", "Inv. Mass", 50, 0, 200, kBlue-9, false, 2);

	data = graphmmt("data.4lep.root", "lep_pt", "lep_E", "lep_phi", "lep_eta", "scaleFactor_MUON", "scaleFactor_ELE", "scaleFactor_PILEUP", "scaleFactor_LepTRIGGER", "mcWeight", "Inv. Mass", 20, 0, 200, kBlack, false, 3);

	//dilepton -> Draw("HIST");
	//monte -> Draw("SAME");

	data->SetMarkerStyle(10);
	//data->Draw("P");

	//data -> Draw("HIST");

	   THStack * hs1 = new THStack("hs1","Stacked Events");

   //signal
   hs1->Add(monte);
   hs1->Add(dilepton);
   //hs1->Add(data);

   TCanvas * c1 = new TCanvas("c1","stacked hists");
   c1->cd(1);
   hs1->Draw("hist");
	data->Draw("SAME");

	TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.7);
	leg -> AddEntry(dilepton, "Higgs (MC)", "f");
	leg -> AddEntry(monte, "Background (MC)", "f");
	leg -> AddEntry(data, "Data", "f");
	leg -> Draw();

}