#define P_ID_cxx
#include "P_ID.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <stdio.h>
#include "jsoncpp/json/json.h"
#include <fstream>
using namespace std;
double EnergyLossFitFunc2(double * x_val, double * par);
void P_ID(float sigmain, float angin)
{
	class P_ID t;
	t.Loop(sigmain, angin);
	
}
void P_ID::Loop(float sigma, float ang)
{
    vector<double> energy,data;
	Json::Reader reader;
	Json::Value value(Json::objectValue);
	Json::Value obj(Json::objectValue);
	Json::Value bin(Json::objectValue);
	Json::Value bind(Json::objectValue);
	Json::Value bint(Json::objectValue);
	Json::Value binpip(Json::objectValue);
	Json::StyledWriter styledWriter;
	int ebin,pbinend=152,dbinend=120,pbins=8,dbins=18,tbins=26;
	TCanvas *c0 = new TCanvas("Triton", "Triton",1900,1060);
	c0->Divide(5, 5,0.0001,0.0001);
	TCanvas *c1 = new TCanvas("Triton1", "Triton1",1900,1060);
	c1->Divide(6, 6,0.0001,0.0001);
	TCanvas *c2 = new TCanvas("Deuteron1", "Deuteron1",1900,1060);
	c2->Divide(7, 7,0.0001,0.0001);
	TCanvas *c3 = new TCanvas("Deuteron2", "Deuteron2",1900,1060);
	c3->Divide(8, 7,0.0001,0.0001);
	TCanvas *c4 = new TCanvas("Proton1", "Proton1",1900,1060);
	c4->Divide(9, 8,0.0001,0.0001);
	TCanvas *c5 = new TCanvas("Proton2", "Proton2",1900,1060);
	c5->Divide(9, 8,0.0001,0.0001);
	TCanvas *c6 = new TCanvas("PionPositive1", "PionPositive1",1900,1060);
	c6->Divide(6, 6,0.0001,0.0001);
	TCanvas *c7 = new TCanvas("PionPositive2", "PionPositive2",1900,1060);
	c7->Divide(6,6,0.0001,0.0001);
	ofstream myfile;
	myfile.open(TString::Format("backgroundnew%2.0f_%.1f.json", ang, sigma));
	std::ifstream background;
  	background.open(TString::Format("../fittingprogram/toffitting/Pippdttofcut_%0.1f.json", sigma));
  	reader.parse(background, value);
	ifstream myfile2;
	myfile2.open(TString::Format("../fittingprogram/toffitting/Pippdttofsystem0_%.1f.txt",sigma));
	cout<<TString::Format("../fittingprogram/toffitting/Pippdttofsystem0_%.1f.txt",sigma)<<endl;
	double mass=938.272;
	double parp[13]={ 0.999248,1.99681e+06,3.28118e+06,43.502,3838.2,15.7087,-20972.1,1.11833,mass,0.973016,-0.165029,0.01555,1.0};
	double pard[13]={ 0.999248,1.99681e+06,3.28118e+06,43.502,3838.2,15.7087,-20972.1,1.11833,2*mass,0.973016,-0.165029,0.01555,1.0};
	double part[13]={ 0.999248,1.99681e+06,3.28118e+06,43.502,3838.2,15.7087,-20972.1,1.11833,3*mass,0.973016,-0.165029,0.01555,1.0};
	double dedxp=0,dedxd=0,dedxt=0 ;
	double mom[1]={0};

	//   In a ROOT session, you can do:
	//      root> .L P_ID.C
	//      root> P_ID t
	//      root> t.GetEntry(12); // Fill t data members with entry number 12
	//      root> t.Show();       // Show values of entry 12
	//      root> t.Show(16);     // Read and show values of entry 16
	//      root> t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	//TH2F *hh =new TH2F("","",1000,0,4000,1000,0,80);
	//TH2F *hh2 =new TH2F("deutron","duetron",1000,0,4000,1000,0,80);
	//TH2F *hh3 =new TH2F("protron","protron",1000,0,4000,1000,0,80);

	TFile *fCut = TFile::Open(TString::Format("../fittingprogram/mdc_cuts_%.1f.root", sigma));
	TH2F *hh0; 
	TH2F *hh1; 
	TH2F *hh2; 
	TH2F *hh3; 
	TH2F *hh4; 
	TH2F *hh5; 
	TH2F *hh6;	
	TH2F *hh7; 
	// TH1D * mass= new TH1D("IN","IN",1000,0,4000);
	// TH1D * mass2= new TH1D("IN1","IN1",1000,0,4000);
	TH1D *mass3 = new TH1D("IN3", "IN3", 1000, 0, 4000);
	TCutG *Cutg[12];
	TGraph *g[140]; 
	int p;
	float intergral_signal = 0, intergral_Background = 0, lowlim, uplim;
	Double_t par[3000];
	Double_t par2[4];
	//TFile *fcut = TFile::Open("../Pippdtmdcnew.root");
	Cutg[0] = (TCutG *)fCut->Get("p_mdc")->Clone();
	Cutg[1] = (TCutG *)fCut->Get("d_mdc")->Clone();
	Cutg[2] = (TCutG *)fCut->Get("t_mdc")->Clone();
	Cutg[3] = (TCutG *)fCut->Get("pip_mdc")->Clone();
	if (fChain == 0)
		return;
	float mass1 = 0;
	Long64_t nentries = fChain->GetEntriesFast();
    cout<<nentries<<endl;
	Long64_t nbytes = 0, nb = 0;
	TFile file(TString::Format("mdchist%2.0f_%.1f.root", ang,sigma));
	if(file.IsZombie()){	
	hh0 = new TH2F("protron", "protron", 4000, 0, 4000, 200, 0, 15);
	hh1 = new TH2F("protron2", "protron2", 4000, 0, 4000, 200, 0, 15);
	hh2 = new TH2F("deutron", "deutron", 4000, 0, 4000, 200, 0, 15);
	hh3 = new TH2F("deutron2", "deutron2", 4000, 0, 4000, 200, 0, 15);
	hh4 = new TH2F("triton", "triton", 4000, 0, 4000, 200, 0, 15);
	hh5 = new TH2F("triton1", "triton1", 4000, 0, 4000, 200, 0, 15);
	hh6 = new TH2F("Pionp", "Pionp", 4000, 0, 4000, 200, 0, 15);
	hh7 = new TH2F("Pionp1", "Pionp1", 4000, 0, 4000, 200, 0, 15);
	TFile *file1=new TFile(TString::Format("mdchist%2.0f_%.1f.root", ang,sigma),"RECREATE");
	for (Long64_t jentry = 0; jentry <300000000; jentry++)
	{
		
		Long64_t ientry = LoadTree(jentry);
		if (jentry % 1000000 == 0)
		{
			std::cout << jentry << std::endl;
		}
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		if (Cutg[0]->IsInside(p_p, p_dedx_mdc) && p_theta < ang + 3 && p_theta > ang - 3)
		{
			if (p_system == 0)
				hh0->Fill(p_p, p_dedx_tof);
			if (p_system == 1)
				hh0->Fill(p_p, p_dedx_tof);
		}
		if (Cutg[1]->IsInside(p_p, p_dedx_mdc) && p_theta < ang + 3 && p_theta > ang - 3)
		{
			if (p_system == 0)
				hh2->Fill(p_p, p_dedx_tof);
			if (p_system == 1)
				hh2->Fill(p_p, p_dedx_tof);
		}
		if (Cutg[2]->IsInside(p_p, p_dedx_mdc) && p_theta < ang + 3 && p_theta > ang - 3)
		{
			if (p_system == 0)
				hh4->Fill(p_p, p_dedx_tof);
			if (p_system == 1)
				hh4->Fill(p_p, p_dedx_tof);
		} 
		if (Cutg[3]->IsInside(p_p, p_dedx_mdc) && p_theta < ang + 3 && p_theta > ang - 3)
		{	
			if (p_system == 0)
				hh6->Fill(p_p, p_dedx_tof);
			if (p_system == 1)
				hh6->Fill(p_p, p_dedx_tof);
		}
	}
	hh0->Write();
	hh1->Write();
	hh2->Write();
	hh3->Write();
	hh4->Write();
	hh5->Write();
	hh6->Write();
	hh7->Write();
	}
	else{
	hh0 =(TH2F *)file.Get("protron")->Clone();
	hh1 =(TH2F *)file.Get("protron2")->Clone();
	hh2 =(TH2F *)file.Get("deutron")->Clone();
	hh3 =(TH2F *)file.Get("deutron2")->Clone();
	hh4 =(TH2F *)file.Get("triton")->Clone();	
	hh5 =(TH2F *)file.Get("triton1")->Clone();
	hh6 =(TH2F *)file.Get("Pionp")->Clone();	
	hh7 =(TH2F *)file.Get("Pionp1")->Clone();}
	TFile *file1=new TFile(TString::Format("backgroundfits%2.0f_%.1f.root", ang,sigma),"RECREATE");
	if(ang>45){
		for (int i = 1; i < 434; ++i){
			myfile2.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		pbinend=120;
		dbinend=112;
		pbins=9;
		dbins=16;
		tbins=28;
	}
	float c, mean, sigl, sigr;
	TH1D *hf[161], *hf2[161], *hf3[161], *hf4[161];
	c = 700;
	mean = 9;
	sigl = 4;
	sigr = 5.9;
	float mom1, mom2;
	std::string line;
	std::getline(myfile2, line);
	cout<<line<<endl;
	std::getline(myfile2, line);
	cout<<line<<endl;
	for (int i = pbins; i < pbinend; i++)
	{   energy.clear();
		data.clear();
		myfile2 >> mom1 >> mom2 >> c >> mean >> sigl >> sigr;
		cout<<mom1<<endl;
		TF1 *f = new TF1("f", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)", mean - sigl * 2, mean + sigr * 2);
		f->SetParNames("Height", "Position", "LeftSigma", "RightSigma");
		f->SetParameters(c, mean, sigl, sigr);
		f->SetParLimits(1,mean-1.2,mean+1.2);
		f->SetParLimits(2,sigl,sigl+0.4);
		f->SetParLimits(3,sigr,sigr+0.4);
		f->SetLineColor(kGreen);
		TF1 *g1 = new TF1(TString::Format("Proton_pol1%d",i), "pol1", mean - sigl * 3.5, mean + sigr * 3.5);
		g1->SetLineColor(kRed);
		mom[0]=hh0->GetXaxis()->GetBinCenter(i*25);
		dedxp=EnergyLossFitFunc2(mom, parp);		

		hf[i] = hh0->ProjectionY(TString::Format("Proton_%d_%d", (i-1)*25,i*25), (i - 1) * 25, i * 25);
		ebin=hf[i]->FindBin(mean-sigl*3.0);
		for(int j=0;j<3;j++){
			energy.push_back(hf[i]->GetBinCenter(ebin-j));
			data.push_back(hf[i]->GetBinContent(ebin-j));
		}
		ebin=hf[i]->FindBin(mean+sigr*3.2);
		for(int j=0;j<3;j++){
			energy.push_back(hf[i]->GetBinCenter(ebin+j));
			data.push_back(hf[i]->GetBinContent(ebin+j));
		}
		hf[i]->SetTitle(TString::Format("Proton at %d -%d [MeV] ", (i - 1) * 25, i * 25));
		hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
		hf[i]->GetYaxis()->SetTitle("Counts (a.u)");
        g[i]= new TGraph(energy.size(), &energy[0], &data[0]);
		if (i-pbins < 72)
		{
			c4->cd(i-pbins+1);
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf[i]->GetYaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetTitleOffset(0.80);
			if(ang>45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf[i]->GetYaxis()->SetTitle("Counts (a.u)");
			hf[i]->GetXaxis()->SetRangeUser(0,round(hf[i]->GetBinCenter(hf[i]->FindBin(mean+sigr*4)))+1);
			hf[i]->Draw();
			TLine *l=new TLine(dedxp,0,dedxp,1000);
    		l->SetLineColor(kRed);
    		l->Draw();
			g[i]->SetMarkerStyle(24);
			g[i]->SetMarkerColor(kRed);
			g[i]->Fit(g1,"RQ");
			g[i]->Draw("p");
			hf[i]->Fit(f, "RQ");
			f->SetRange(0,14);
			f->Draw("same");
			// hf[i]->Fit(g1, "RMQ+");
			f->GetParameters(&par[0]);

			// g1->GetParameters(&par[4]);
			// int gaussup = hf[i]->GetXaxis()->FindBin(par[1] + par[3] * sigma);
			// int gausslow = hf[i]->GetXaxis()->FindBin(par[1] - par[2] * sigma);
			intergral_Background = g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
			intergral_signal = f->Integral(mean - sigl * sigma, mean + sigr * sigma) - g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
		}
		else
		{
			c5->cd(i - 71-pbins);
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf[i]->GetYaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetTitleOffset(0.80);
			if(ang>45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf[i]->GetYaxis()->SetTitle("Counts (a.u)");
			hf[i]->GetXaxis()->SetRangeUser(0,round(hf[i]->GetBinCenter(hf[i]->FindBin(mean+sigr*4)))+1);
			hf[i]->Draw();
			g[i]->SetMarkerStyle(24);
			g[i]->SetMarkerColor(kRed);
			g[i]->Fit(g1,"RQ");
			g[i]->Draw("p");
			f->SetRange(0,14);
			f->Draw("same");
			// //f->Draw("same");
			hf[i]->Fit(f, "RMQ");
			// hf[i]->Fit(g1, "RMQ+");
			f->GetParameters(&par[0]);
			// g1->SetLineColor(6);
			// g1->GetParameters(&par[4]);
			// int gaussup = hf[i]->GetXaxis()->FindBin(par[1] + par[3] * sigma);
			// int gausslow = hf[i]->GetXaxis()->FindBin(par[1] - par[2] * sigma);
			intergral_Background = g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
			intergral_signal = f->Integral(mean - sigl * sigma, mean + sigr * sigma) - g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
		}
		hf[i]->Write();
		g[i]->Write(TString::Format("Proton_%d_%d_graph", (i-1)*25,i*25));
		//myfile << (i - 1) * 25 << " " << i * 25 << " " << intergral_Background << " " << intergral_signal << " " << (intergral_Background / (intergral_Background + intergral_signal)) * 100 << endl;
		//myfile<<(i-1)*25<<" "<<i*25<<" "<<f->GetParameter(0)<<" "<<f->GetParameter(1)<<" "<<f->GetParameter(2)<<" "<<f->GetParameter(3)<<endl;
		bin["fistbin"][i-pbins]=(i-1)*25;
		bin["lastbin"][i-pbins]=(i)*25;
		if((intergral_Background / (intergral_Background + intergral_signal))<100&&(intergral_Background / (intergral_Background + intergral_signal))>-100)
		{
			bin["background"][i-pbins]=intergral_Background;
			bin["signal"][i-pbins]=intergral_signal;
			bin["backper"][i-pbins]= (intergral_Background / (intergral_Background + intergral_signal)) * 100 ;
		}
		else
		{
			bin["background"][i-pbins]=0;
			bin["signal"][i-pbins]=0;
			bin["backper"][i-pbins]= 0;
		}


	}
	c4->Write();
	c5->Write();
	obj["Proton"]=bin;
	//myfile << styledWriter.write(obj);
	c = 700;
	mean = 6.2;
	sigl = 1.6;
	sigr = 5.9;
	//myfile<<"Deutron"<<endl;
	std::getline(myfile2, line);
	cout<<line<<endl;
	std::getline(myfile2, line);
    cout<<line<<endl;
    std::getline(myfile2, line);
    cout<<line<<endl;
	//myfile<<"Momstart Momend Height Position LeftSigma RightSigma"<<endl;
	float	pmean, 	plsig ,prsig; 
	gStyle->SetOptStat(10);

	for (int i = dbins; i < dbinend; i++)
	{	if(ang<45){
			pmean=value["PSys0"]["Mean"][i-dbins+(dbins-pbins)].asFloat();
			plsig=value["PSys0"]["LeftSigma"][i-dbins+(dbins-pbins)].asFloat();
			prsig=value["PSys0"]["RightSigma"][i-dbins+(dbins-pbins)].asFloat();}
		else{
			pmean=value["PSys1"]["Mean"][i-dbins+(dbins-pbins)].asFloat();
			plsig=value["PSys1"]["LeftSigma"][i-dbins+(dbins-pbins)].asFloat();
			prsig=value["PSys1"]["RightSigma"][i-dbins+(dbins-pbins)].asFloat();
		}
		myfile2 >> mom1 >>mom2 >> c >> mean >> sigl >> sigr;
		cout<<mom1<<" "<<value["PSys0"]["Fistbin"][i-dbins+(dbins-pbins)].asFloat()<<endl;
		energy.clear();
		data.clear();
		TF1 *f0 = new TF1("f", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)", mean - sigl*1.5, mean + sigr*1.9);
		TF1 *f = new TF1("landau1", "landau", mean - sigl, mean + sigr * 4);
		TF1 *f2 = new TF1("f3", "landau", 0.0, mean - mean * 0.01);
		TF1 *f4 = new TF1("landau2", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)+[4]* TMath::Gaus(x, [5], ((x < [5]) ? [6] : [7]), 0)", 0.0, 12);
		f0->SetParNames("Height", "Position", "LeftSigma", "RightSigma");
		// f2->SetParameters(1752, 2.7, 0.8, 1.4);
		//TF1 *g1= new TF1 ("pol","pol2",0.0,12);
		f0->SetParameters(c, mean, sigl, sigr);
		f0->SetParLimits(1,mean-0.1,mean+0.1);
		f0->SetParLimits(2,sigl-0.1,sigl+0.1);
		f0->SetParLimits(3,sigr-0.1,sigr+0.1);
		f0->SetLineColor(kGreen);
		f4->SetLineColor(6);
		gStyle->SetTitleFontSize(0.1);
		gStyle->SetLabelSize(0.1);
		//gStyle->SetOptStat(11);
		hf[i] = hh2->ProjectionY(TString::Format("Deuteron_%d_%d", (i-1)*25,i*25), (i - 1) * 25, i * 25);
		hf[i]->SetTitle(TString::Format("Deuteron at %d-%d [MeV] ", (i - 1) * 25, i * 25));
		// TF1 *g1 = new TF1(TString::Format("Deuteron_pol1%d",i), "expo", mean - sigl * 3.5-1, mean + sigr * 3.5);//-1 in lower limit for angle 35
		// ebin=hf[i]->FindBin(mean-sigl*1.8);//2.5 for 35, 2.0 for 65,80,1.8 for 25 
		mom[0]=hh2->GetXaxis()->GetBinCenter(i*25+13);
		dedxp=EnergyLossFitFunc2(mom, parp);
		dedxd=EnergyLossFitFunc2(mom, pard);
		TF1 *g1 = new TF1 ("g1","[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)",0.7,14);
		g1->SetParLimits(1,pmean-0.1,pmean+0.1);
		g1->SetParLimits(2,prsig-0.1,prsig+0.1);
		g1->SetParLimits(3,prsig-0.1,prsig+0.1);
		// if(ang>45) 	f0->FixParameter(1,dedxd);
		// for(int j=0;j<3;j++){ //changed to 3 for35 and 25
		// 	energy.push_back(hf[i]->GetBinCenter(ebin-j));
		// 	data.push_back(hf[i]->GetBinContent(ebin-j));
		// }
		// ebin=hf[i]->FindBin(mean+sigr*3.0);//3.2 for 35, 3.0 for 65,25,80
		// for(int j=0;j<3;j++){
		// 	energy.push_back(hf[i]->GetBinCenter(ebin+j));
		// 	data.push_back(hf[i]->GetBinContent(ebin+j));
		// }
		
		//g[i]= new TGraph(energy.size(), &energy[0], &data[0]);
		if (i < 49+dbins)
		{	
			c2->cd(i - dbins+1);
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			// gPad->Set
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf[i]->GetYaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetTitleOffset(0.80);
			gPad->SetLogy();
			hf[i]->GetXaxis()->SetRangeUser(0,round(hf[i]->GetBinCenter(hf[i]->FindBin(mean+sigr*4))));
			hf[i]->Draw();
			// g[i]->SetMarkerStyle(24);
			// g[i]->SetMarkerColor(kRed);
			// g[i]->Fit(g1,"RQ");
			// g[i]->Draw("p");
			hf[i]->Fit(f0, "QR");
			hf[i]->Fit(g1, "QR+");
			TLine *l=new TLine(dedxp,0,dedxp,1000);
    		l->SetLineColor(kRed);
    		l->Draw();
			TLine *l2=new TLine(dedxd,0,dedxd,1000);
    		l2->SetLineColor(kRed);
    		l2->Draw();
			f0->GetParameters(&par[0]);
			if(ang>45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf[i]->GetYaxis()->SetTitle("Counts (a.u)");
			f0->GetParameters(&par[0]);
			g1->GetParameters(&par[4]);
			f4->SetParameters(&par[0]);
			f4->SetParLimits(5,pmean-0.2,pmean+0.1);
			f4->SetParLimits(1,mean,mean+2);
			f4->SetParLimits(2,sigl-0.1,sigl);
			f4->SetParLimits(3,sigr-0.2,sigr+0.2);
			f4->SetParLimits(7,par[7]-0.4,par[7]+0.1);
			hf[i]->Fit(f4, "RQ");
			f4->GetParameters(&par[0]);
			f0->SetParameters(&par[0]);
			g1->SetParameters(&par[4]);
			f0->SetRange(0,14);
			f0->Draw("same");
			g1->Draw("same");
			f4->Draw("same");
			
	// 		hf[i]->GetXaxis()->SetTitleOffset(0.90);
	// 		//f->Draw("same");
			
	// 		//hf[i]->Fit(g1,"R+");
	 		//hf[i]->Fit(f2, "QR+");
			
	// 		//g1-> GetParameters (&par[4]);
	// 		f2->GetParameters(&par[3]);
	// 		hf[i]->Fit(f0, "QR");
	// 		//int gaussup=hf[i]->GetXaxis()->FindBin(par[1]+par[3]*sigma);
	// 		//int gausslow=hf[i]->GetXaxis()->FindBin(par[1]-par[2]*sigma);
	// 		//intergral_Background=f2->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma);
	// 		//intergral_signal= f->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma)-f2->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma);
	// 		TF1 *f1 = new TF1("double_landau", "landau(0)+landau(3)", 0.0, 12); //f->Draw("same");
	// 		f1->SetParameters(par);
	// 		//f1->FixParameter(8,par[8]);
	// 		//f1->FixParameter(9,par[9]);
	// 		f1->SetLineColor(kRed);
	// 		hf[i]->Fit(f1, "+", "", 0.0, 10);
	// 		f1->GetParameters(&par[0]);
	// 		f4->SetParameters(&par[3]);
	// 		//f4->FixParameter(0,par[3]);
	// 		//f4->FixParameter(1,par[4]);
	// 		f4->Draw("same");
	// 		f0->GetParameters(&par2[0]);
	// 		//hf[i]->Fit(f4,"BR+");
			intergral_Background = g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
			intergral_signal = f0->Integral(mean - sigl * sigma, mean + sigr * sigma) - g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
			//intergral_Background = g1->Integral(par[1] - par[2]*sigma, par[1] + par[3] * sigma);
			//intergral_signal = f0->Integral(par[1] - par[2]*sigma, par[1] + par[3] * sigma) - g1->Integral(par[1] - par[2]*sigma, par[1] + par[3] * sigma);
		}
		else
		{
			c3->cd(i - (48+dbins));
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf[i]->GetYaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetLabelSize(0.11);
			hf[i]->GetXaxis()->SetTitleSize(0.11);
			hf[i]->GetYaxis()->SetTitleOffset(0.80);
			hf[i]->GetXaxis()->SetRangeUser(0,round(hf[i]->GetBinCenter(hf[i]->FindBin(mean+sigr*4)))+1);
			hf[i]->Draw();
			cout<<round(hf[i]->GetBinCenter(hf[i]->FindBin(mean+sigr*4)))<<endl;
			// g[i]->SetMarkerStyle(24);
			// g[i]->SetMarkerColor(kRed);
			// g[i]->Fit(g1,"RMQ");
			// g[i]->Draw("p");
			hf[i]->Fit(f0, "QRM");
			hf[i]->Fit(g1, "QR+");
			f0->GetParameters(&par[0]);
			if(ang>45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf[i]->GetYaxis()->SetTitle("Counts (a.u)");
			f0->GetParameters(&par[0]);
			g1->GetParameters(&par[4]);
			f4->SetParameters(&par[0]);
			f4->FixParameter(5,par[5]);
			f4->SetParLimits(1,mean,mean+2);
			f4->SetParLimits(2,sigl-0.1,sigl);
			f4->SetParLimits(3,sigr-0.2,sigr+0.2);
			hf[i]->Fit(f4, "RQ");
			f4->GetParameters(&par[0]);
			f0->SetParameters(&par[0]);
			g1->SetParameters(&par[4]);
			f0->SetRange(0,14);
			f0->Draw("same");
			g1->Draw("same");
			f4->Draw("same");
	// 		hf[i]->GetYaxis()->SetTitleSize(0.08);
	// 		hf[i]->GetYaxis()->SetLabelSize(0.08);
	// 		hf[i]->GetXaxis()->SetLabelSize(0.08);
	// 		hf[i]->GetXaxis()->SetTitleSize(0.08);
	// 		hf[i]->GetYaxis()->SetTitleOffset(0.60);
	// 		hf[i]->GetXaxis()->SetTitleOffset(0.96);
	// 		//f->Draw("same");
	// 		hf[i]->Fit(f, "QR");
	// 		//hf[i]->Fit(g1,"R+");
	 		//hf[i]->Fit(f2, "QR+");
			 
	// 		f->GetParameters(&par[0]);
	// 		//g1-> GetParameters (&par[4]);
	// 		f2->GetParameters(&par[3]);
	// 		//int gaussup=hf[i]->GetXaxis()->FindBin(par[1]+par[3]*sigma);
	// 		//int gausslow=hf[i]->GetXaxis()->FindBin(par[1]-par[2]*sigma);
	// 		hf[i]->Fit(f0, "QR");
	// 		TF1 *f1 = new TF1("double_gaus2", "landau(0)+landau(3)", 0.0, 12); //f->Draw("same");
	// 		f1->SetParameters(&par[0]);
	// 		f1->SetLineColor(kRed);
	// 		//f1->FixParameter(8,par[8]);
	// 		//f1->FixParameter(9,par[9]);
	// 		//f1->FixParameter(7,par[7]);
	// 		hf[i]->Fit(f1, "+", "", 0.0, 12);
	// 		f1->GetParameters(&par[0]);
	// 		f4->SetParameters(&par[3]);
	// 		///f4->FixParameter(0,par[3]);
	// 		//f4->FixParameter(1,par[4]);
	// 		//f4->FixParameter(2,par[5]);
	// 		f0->GetParameters(&par2[0]);
	// 		hf[i]->Fit(f4, "QBR+");
	// 		intergral_Background = f4->Integral(par2[1] - par2[2], par2[1] + par2[3] * sigma);
	// 		intergral_signal = f0->Integral(par2[1] - par2[2], par2[1] + par2[3] * sigma) - f4->Integral(par2[1] - par2[2], par2[1] + par2[3] * sigma);
	        //intergral_Background = g1->Integral(par[1] - par[2]*sigma, par[1] + par[3] * sigma);
			//intergral_signal = f0->Integral(par[1] - par[2]*sigma, par[1] + par[3] * sigma) - g1->Integral(par[1] - par[2]*sigma, par[1] + par[3] * sigma);
			intergral_Background = g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
			intergral_signal = f0->Integral(mean - sigl * sigma, mean + sigr * sigma) - g1->Integral(mean - sigl * sigma, mean + sigr * sigma);
		}
		hf[i]->Write();
		//g[i]->Write(TString::Format("Deuteron_%d_%d_graph", (i-1)*25,i*25));
		//myfile << (i - 1) * 25 << " " << i * 25 << " " << intergral_Background << " " << intergral_signal << " " << (intergral_Background / (intergral_Background + intergral_signal)) * 100 << endl;
		//myfile<<(i-1)*25<<" "<<i*25<<" "<<f->GetParameter(0)<<" "<<f->GetParameter(1)<<" "<<f->GetParameter(2)<<" "<<f->GetParameter(3)<<endl;
		bind["fistbin"][i-dbins]=(i-1)*25;
		bind["lastbin"][i-dbins]=(i)*25;
		if((intergral_Background / (intergral_Background + intergral_signal))<100&&(intergral_Background / (intergral_Background + intergral_signal))>-100)
		{
			bind["background"][i-dbins]=intergral_Background;
			bind["signal"][i-dbins]=intergral_signal;
			bind["backper"][i-dbins]= (intergral_Background / (intergral_Background + intergral_signal)) * 100 ;
		}
		else
		{
			bind["background"][i-dbins]=0;
			bind["signal"][i-dbins]=0;
			bind["backper"][i-dbins]= 0;
		}
 }
 obj["Deuteron"]=bind;
 c3->Write();
 c2->Write();
//myfile << styledWriter.write(obj);
	// //TF1 *f2 = new TF1("f", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)", 3.0, 40.0);
	c = 3314;
	mean = 12;
	sigl = 4;
	sigr = 7;
	//myfile<<"Triton"<<endl;
	std::getline(myfile2, line);
	cout<<line<<endl;
	std::getline(myfile2, line);
    cout<<line<<endl;
    std::getline(myfile2, line);
    cout<<line<<endl;
	float dmean;
	float dlsig;
	float drsig;
	float firstbin;
	//myfile<<"Momstart Momend Height Position LeftSigma RightSigma"<<endl;
	for (int i = tbins; i < 86; i++)
	{   if(ang<45){
		firstbin=value["DSys0"]["Fistbin"][i-tbins+(tbins-dbins)].asFloat();
		dmean=value["DSys0"]["Mean"][i-tbins+(tbins-dbins)].asFloat();
	    dlsig=value["DSys0"]["LeftSigma"][i-tbins+(tbins-dbins)].asFloat();
		drsig=value["DSys0"]["RightSigma"][i-tbins+(tbins-dbins)].asFloat();
		pmean=value["PSys0"]["Mean"][i-tbins+(tbins-pbins)].asFloat();
	    plsig=value["PSys0"]["LeftSigma"][i-tbins+(tbins-pbins)].asFloat();
		prsig=value["PSys0"]["RightSigma"][i-tbins+(tbins-pbins)].asFloat();}
		else{
		firstbin=value["DSys1"]["Fistbin"][i-tbins+(tbins-dbins)].asFloat();
		dmean=value["DSys1"]["Mean"][i-tbins+(tbins-dbins)].asFloat();
	    dlsig=value["DSys1"]["LeftSigma"][i-tbins+(tbins-dbins)].asFloat();
		drsig=value["DSys1"]["RightSigma"][i-tbins+(tbins-dbins)].asFloat();
		pmean=value["PSys1"]["Mean"][i-tbins+(tbins-pbins)].asFloat();
	    plsig=value["PSys1"]["LeftSigma"][i-tbins+(tbins-pbins)].asFloat();
		prsig=value["PSys1"]["RightSigma"][i-tbins+(tbins-pbins)].asFloat()	;
		}
		myfile2 >> mom1 >> mom2 >> c >> mean >> sigl >> sigr;
		cout<<mom1<<" "<<firstbin<<" "<<value["TSys1"]["Fistbin"][i-tbins].asFloat()<<endl;
		energy.clear();
		data.clear();
		TF1 *f0 = new TF1("f", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)", mean - sigl * 0.8, mean + sigr * 3);
		TF1 *f = new TF1("landau", "expo",6, 14);
		TF1 *f2 = new TF1("f3", "gaus", 0.0, mean - mean * 0.01);
		
		f0->SetParNames("Height", "Position", "LeftSigma", "RightSigma");
		// f2->SetParameters(1752, 2.7, 0.8, 1.4);
		f0->SetParameters(c, mean, sigl, sigr);
		//f0->FixParameter(1, mean);
		//f0->FixParameter(2, sigl);
		//f0->FixParameter(3, sigr);
		f0->SetLineColor(kGreen);
		f0->SetParLimits(1,mean,mean+2);
		if(ang>45)f0->SetParLimits(1,mean,mean+2);
		f0->SetParLimits(3,sigr-0.2,sigr+0.2);
		f0->SetParLimits(2,sigl-0.2,sigl+0.2);
		hf3[i] = hh4->ProjectionY(TString::Format("Triton_%d_%d", (i-1)*25,i*25), (i - 1) * 25, i * 25);
		hf3[i]->SetTitle(TString::Format("Triton at %d -%d [MeV] ", (i - 1) * 25, i * 25));
        // TF1 *g1 = new TF1(TString::Format("Triton_pol1%d",i), "expo", mean - sigl * 3.5, mean + sigr * 6);
		ebin=hf3[i]->FindBin(mean-sigl*0.8);//0.8 for lower angles
		// for(int j=0;j<5;j++){
		// 	energy.push_back(hf3[i]->GetBinCenter(ebin-j));
		// 	data.push_back(hf3[i]->GetBinContent(ebin-j));
		// }
		// ebin=hf3[i]->FindBin(mean+sigr*2.5);
		// for(int j=0;j<5;j++){
		// 	energy.push_back(hf3[i]->GetBinCenter(ebin+j));
		// 	data.push_back(hf3[i]->GetBinContent(ebin+j));
		// }
		mom[0]=hh4->GetXaxis()->GetBinCenter(i*25+13);
		dedxp=EnergyLossFitFunc2(mom, parp);
		dedxd=EnergyLossFitFunc2(mom, pard);	
		dedxt=EnergyLossFitFunc2(mom, part);	
		
		TF1 *f4;
		TF1 *g1 = new TF1 ("g1","[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)",dmean-dlsig,dmean+drsig);
		g1->SetParLimits(1,dmean-0.1,dmean+0.1);
		//f->SetParLimits(1,pmean-0.1,pmean+0.1);
		g1->SetParLimits(2,drsig-0.1,drsig+0.1);
		g1->SetParLimits(3,drsig-0.1,drsig+0.1);
		// if(ang>45){
		// f0->FixParameter(1,dedxt);
		// f0->SetParLimits(2,sigl-0.4,sigl+0.1);
		// f0->SetParLimits(3,sigr-0.4,sigr+0.2);
		// //g1->SetParameters(1752, 2.7, 50);
		// g1->FixParameter(1,dedxd);
		// f2->SetParLimits(1,dedxp-0.2,dedxp+0.2);}
		if(sigma>0.9 & ang>45){
			  f4= new TF1("landau2", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)+[4]* TMath::Gaus(x, [5], ((x < [5]) ? [6] : [7]), 0)+gaus(8)", 0.0, 12);}
		else{
			 f4= new TF1("landau2", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)+[4]* TMath::Gaus(x, [5], ((x < [5]) ? [6] : [7]), 0)", 0.0, 12);//[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)+gaus(4)", 0.0, 12);
		}
		f4->SetLineColor(6);
		// g[i]= new TGraph(energy.size(), &energy[0], &data[0]);
		if (i < 25+tbins)
			{
			c0->cd(i - tbins+1);
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf3[i]->GetYaxis()->SetTitleSize(0.11);
			hf3[i]->GetYaxis()->SetLabelSize(0.11);
			hf3[i]->GetXaxis()->SetLabelSize(0.11);
			hf3[i]->GetXaxis()->SetTitleSize(0.11);
			hf3[i]->GetYaxis()->SetTitleOffset(0.80);
			if(ang>45) hf3[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf3[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf3[i]->GetYaxis()->SetTitle("Counts (a.u)");
			hf3[i]->Draw();
		// 	if(ang>45){
		// 	for(int j=0;j<3;j++){
    	// 	hf3[i]->Fit(g1,"RQ");
		// 	hf3[i]->Fit(f2,"RQ");
    	// 	g1->GetParameters(&par[0]);
    	// 	g1->SetRange(par[1]-par[2]*0.7,par[1]+par[2]*0.7);
		// 	f2->GetParameters(&par[0]);
		// 	f2->SetRange(par[1]-par[2]*0.7,par[1]+par[2]*0.7);
    	// 	//cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<endl;
    	// 	}
		// 	g1->GetParameters(&par[0]);
		// 	g1->SetRange(par[1]-par[2]*2,par[1]+par[2]*2);
		// 	hf3[i]->Fit(g1,"RQ");
		// 	g1->SetRange(0,mean+2);
		// 	f2->SetRange(0,mean+2);
		// 	hf3[i]->Fit(f, "RQ");
		// 	f0->GetParameters(&par[0]);
		// 	g1->GetParameters(&par[4]);
		// 	f2->GetParameters(&par[7]);
		// 	f4->SetParameters(&par[0]);
		// 	f4->FixParameter(1,dedxt);
		// 	f4->FixParameter(5,dedxd);
		// 	f4->SetParLimits(8,dedxp-0.2,dedxp+0.2);
		// 	//f4->SetParLimits(2,sigl-0.4,sigl+0.2);
		// 	//f4->SetParLimits(3,sigr-0.4,sigr+0.2);
		// 	hf3[i]->Fit(f4, "RQ");
		// 	f0->SetRange(0,14);
		// 	f4->GetParameters(&par[0]);
		// 	f0->SetParameters(&par[0]);
		// 	g1->SetParameters(&par[4]);
		// 	f2->SetParameters(&par[7]);
		// 	f->Draw("same");
		// 	g1->GetParameters(&par[0]);
		// 	g1->Draw("same");
		// 	cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<endl;
		// 	cout<<g1->Integral(mean - sigl * sigma, mean + sigr * sigma)<<" "<<f->Integral(mean - sigl * sigma, mean + sigr * sigma)<<" "<<par[2]<<endl;
		// 	if(sigma>0.9)f2->Draw("same");
		// 	TLine *l=new TLine(dedxt,0,dedxt,1000);
    	// 	l->SetLineColor(kRed);
    	// 	l->Draw();
		// 	TLine *l2=new TLine(dedxd,0,dedxd,1000);
    	// 	l2->SetLineColor(kRed);
    	// 	l2->Draw();
			
		// 	}
		// else{	
				hf3[i]->Fit(g1,"RQ");
				hf3[i]->Fit(f0, "RQ+");
				f0->GetParameters(&par[0]);
				g1->GetParameters(&par[4]);
				f4->SetParameters(&par[0]);
				f4->FixParameter(5,par[5]);
				f4->SetParLimits(1,mean,mean+2);
				f4->SetParLimits(2,sigl-0.1,sigl);
				f4->SetParLimits(3,sigr-0.2,sigr+0.2);
				f4->SetParLimits(6,par[6]-0.2,par[6]+0.1);
				
				hf3[i]->Fit(f4, "RQ");
				hf3[i]->Fit(f, "RQ");
				
				f4->GetParameters(&par[0]);
				f0->SetParameters(&par[0]);
				g1->SetParameters(&par[4]);
				f0->SetRange(0,14);
				f->SetRange(0,14);
				g1->SetRange(0,14);
				f0->Draw("same");
				g1->Draw("same");
				f4->Draw("same");
				f->Draw("same");
				// g1->SetRange(0,mean+2);
			    
				
				TLine *l=new TLine(mean,0,mean,1000);
    			l->SetLineColor(kRed);
    			l->Draw();
			// }
			
			// g[i]->SetMarkerStyle(24);
			// g[i]->SetMarkerColor(kRed);
			// g[i]->Fit(g1,"RQ");
			// g[i]->Draw("p");
			// f0->GetParameters(&par2[0]);
	// 		hf3[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
	// 		hf3[i]->GetYaxis()->SetTitle("Counts (a.u)");
	// 		hf3[i]->GetYaxis()->SetTitleSize(0.1);
	// 		hf3[i]->GetYaxis()->SetLabelSize(0.1);
	// 		hf[i]->GetXaxis()->SetLabelSize(0.1);
	// 		hf3[i]->GetXaxis()->SetTitleSize(0.1);
	// 		hf3[i]->GetYaxis()->SetTitleOffset(0.70);
	// 		hf3[i]->GetXaxis()->SetTitleOffset(0.96);
	// 		//f->Draw("same");
	// 		hf3[i]->Fit(f, "R");
	// 		//hf[i]->Fit(g1,"R+");
	// 		hf3[i]->Fit(f2, "R+");
	// 		f->GetParameters(&par[0]);
	// 		//g1-> GetParameters (&par[4]);
	// 		f2->GetParameters(&par[3]);
	// 		hf3[i]->Fit(f0, "R");
	// 		//int gaussup=hf[i]->GetXaxis()->FindBin(par[1]+par[3]*sigma);
	// 		//int gausslow=hf[i]->GetXaxis()->FindBin(par[1]-par[2]*sigma);
	// 		//intergral_Background=f2->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma);
	// 		//intergral_signal= f->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma)-f2->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma);
	// 		TF1 *f1 = new TF1("double_landau", "landau(0)+landau(3)", 0.0, 12); //f->Draw("same");
	// 		f1->SetParameters(par);
	// 		//f1->FixParameter(8,par[8]);
	// 		//f1->FixParameter(9,par[9]);
	// 		f1->SetLineColor(kRed);
	// 		hf3[i]->Fit(f1, "+", "", 0.0, 10);
	// 		f1->GetParameters(&par[0]);
	// 		f4->SetParameters(&par[3]);
	// 		f0->GetParameters(&par2[0]);
	// 		hf3[i]->Fit(f4, "BR+");
			//g1->GetParameters(&par[0]);
			cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<endl;
			intergral_Background = f->Integral(mean, mean + sigr * sigma);
			intergral_signal = f0->Integral(mean - sigl * sigma, mean + sigr * sigma);
			cout<<intergral_signal<<" "<<intergral_Background <<" "<<par[2]<<endl;
		}
		else
		{
			c1->cd(i - (tbins+24));
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf3[i]->GetYaxis()->SetTitleSize(0.11);
			hf3[i]->GetYaxis()->SetLabelSize(0.11);
			hf3[i]->GetXaxis()->SetLabelSize(0.11);
			hf3[i]->GetXaxis()->SetTitleSize(0.11);
			hf3[i]->GetYaxis()->SetTitleOffset(0.80);
			if(ang>45) hf3[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf3[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf3[i]->GetYaxis()->SetTitle("Counts (a.u)");
			hf3[i]->Fit(f, "RQ");
			
	// 		gPad->SetBottomMargin(0.17);
	// 		gPad->SetTopMargin(0);
	// 		gPad->SetRightMargin(0.2);
	// 		gPad->SetLeftMargin(0.17);
			
			//hf3[i]->Fit(f0, "RQ");
			// g[i]->SetMarkerStyle(24);
			// g[i]->SetMarkerColor(kRed);
			// g[i]->Fit(g1,"RMQ");
			// g[i]->Draw("p");
			//hf3[i]->Fit(g1,"RQ");
			hf3[i]->Draw();
			// if(ang>45){
			// for(int j=0;j<3;j++){
    		// hf3[i]->Fit(g1,"RQ");
			// hf3[i]->Fit(f2,"RQ");
    		// g1->GetParameters(&par[0]);
    		// g1->SetRange(par[1]-par[2]*0.7,par[1]+par[2]*0.7);
			// f2->GetParameters(&par[0]);
			// f2->SetRange(par[1]-par[2]*0.7,par[1]+par[2]*0.7);
    		// //cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<endl;
    		// }
			// g1->GetParameters(&par[0]);
			// g1->SetRange(par[1]-par[2]*2,par[1]+par[2]*2);
			// hf3[i]->Fit(g1,"RQ");
			// g1->SetRange(0,mean+2);
			// f2->SetRange(0,mean+2);
			// hf3[i]->Fit(f, "RQ");
			// f->GetParameters(&par[0]);
			// g1->GetParameters(&par[4]);
			// f2->GetParameters(&par[7]);
			// f4->SetParameters(&par[0]);
			// f4->FixParameter(1,dedxt);
			// f4->FixParameter(5,dedxd);
			// f4->SetParLimits(8,dedxp-0.2,dedxp+0.2);
			// //f4->SetParLimits(2,sigl-0.4,sigl+0.2);
			// //f4->SetParLimits(3,sigr-0.4,sigr+0.2);
			// hf3[i]->Fit(f4, "RQ");
			// f0->SetRange(0,14);
			// f4->GetParameters(&par[0]);
			// f0->SetParameters(&par[0]);
			// g1->SetParameters(&par[4]);
			// f2->SetParameters(&par[7]);
			// f->Draw("same");
			// g1->GetParameters(&par[0]);
			// g1->Draw("same");
			// cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<endl;
			// cout<<g1->Integral(mean - sigl * sigma, mean + sigr * sigma)<<" "<<f->Integral(mean - sigl * sigma, mean + sigr * sigma)<<" "<<par[2]<<endl;
			// if(sigma>0.9)f2->Draw("same");
			// TLine *l=new TLine(dedxt,0,dedxt,1000);
    		// l->SetLineColor(kRed);
    		// l->Draw();
			// TLine *l2=new TLine(dedxd,0,dedxd,1000);
    		// l2->SetLineColor(kRed);
    		// l2->Draw();
			
			// }
			// else{	
				hf3[i]->Fit(g1,"RQ");
				hf3[i]->Fit(f0, "RQ+");
				f0->GetParameters(&par[0]);
				g1->GetParameters(&par[4]);
				f4->SetParameters(&par[0]);
				f4->FixParameter(5,par[5]);
				f4->SetParLimits(2,sigl-0.2,sigl);
				f4->SetParLimits(3,sigr-0.2,sigr+0.2);
				hf3[i]->Fit(f4, "RQ");
				f0->SetRange(0,14);
				f4->GetParameters(&par[0]);
				f0->SetParameters(&par[0]);
				f->SetRange(0,14);
				f->Draw("same");
				g1->SetParameters(&par[4]);
				f0->Draw("same");
				// g1->SetRange(0,mean+2);
			    g1->Draw("same");
				f4->Draw("same");
				TLine *l=new TLine(mean,0,mean,1000);
    			l->SetLineColor(kRed);
    			l->Draw();
			//}
	// 		hf3[i]->GetYaxis()->SetTitle("Counts (a.u)");
	// 		hf3[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
	// 		hf3[i]->GetYaxis()->SetTitleSize(0.1);
	// 		hf3[i]->GetYaxis()->SetLabelSize(0.1);
	// 		hf[i]->GetXaxis()->SetLabelSize(0.1);
	// 		hf3[i]->GetXaxis()->SetTitleSize(0.1);
	// 		hf3[i]->GetYaxis()->SetTitleOffset(0.70);
	// 		hf3[i]->GetXaxis()->SetTitleOffset(0.96);
	// 		//f->Draw("same");
	// 		hf3[i]->Fit(f, "R");
	// 		//hf[i]->Fit(g1,"R+");
	// 		hf3[i]->Fit(f2, "R+");
	// 		f->GetParameters(&par[0]);
	// 		//g1-> GetParameters (&par[4]);
	// 		f2->GetParameters(&par[3]);
	// 		//int gaussup=hf[i]->GetXaxis()->FindBin(par[1]+par[3]*sigma);
	// 		//int gausslow=hf[i]->GetXaxis()->FindBin(par[1]-par[2]*sigma);
	// 		hf3[i]->Fit(f0, "R");
	// 		TF1 *f1 = new TF1("double_gaus2", "landau(0)+landau(3)", 0.0, 12); //f->Draw("same");
	// 		f1->SetParameters(par);
	// 		f1->SetLineColor(kRed);
	// 		//f1->FixParameter(8,par[8]);
	// 		//f1->FixParameter(9,par[9]);
	// 		//f1->FixParameter(7,par[7]);
	// 		hf3[i]->Fit(f1, "+", "", 0.0, 12);
	// 		f1->GetParameters(&par[0]);
	// 		f4->SetParameters(&par[3]);
	// 		f0->GetParameters(&par2[0]);
	// 		hf3[i]->Fit(f4, "BR+");
			intergral_Background = f->Integral(mean, mean + sigr * sigma);
			intergral_signal = f0->Integral(mean - sigl * sigma, mean + sigr * sigma);
		}
		hf3[i]->Write();
		bint["fistbin"][i-tbins]=(i-1)*25;
		bint["lastbin"][i-tbins]=(i)*25;
		if((intergral_Background / intergral_signal)<1&&(intergral_Background / intergral_signal)>-1)
		{
			bint["background"][i-tbins]=intergral_Background;
			bint["signal"][i-tbins]=intergral_signal;
			bint["backper"][i-tbins]= (intergral_Background / intergral_signal) * 100 ;
		}
		else
		{
			bint["background"][i-tbins]=0;
			bint["signal"][i-tbins]=0;
			bint["backper"][i-tbins]= 0;
		}
	}
	c0->Write();
 	c1->Write();
	obj["Triton"]=bint;
	//myfile << styledWriter.write(obj);
	//myfile << "PionPositive" << endl;
	std::getline(myfile2, line);
	cout<<line<<endl;
	std::getline(myfile2, line);
    cout<<line<<endl;
    std::getline(myfile2, line);
    cout<<line<<endl;
	//myfile<<"Momstart Momend Height Position LeftSigma RightSigma"<<endl;
	for (int i = 1; i < 63; i++)
	{	if(ang<45&&pbins<=i){
			pmean=value["PSys0"]["Mean"][i-pbins].asFloat();
			plsig=value["PSys0"]["LeftSigma"][i-pbins].asFloat();
			prsig=value["PSys0"]["RightSigma"][i-pbins].asFloat();
			}
		else if(pbins<=i){
			pmean=value["PSys1"]["Mean"][i-pbins].asFloat();
			plsig=value["PSys1"]["LeftSigma"][i-pbins].asFloat();
			prsig=value["PSys1"]["RightSigma"][i-pbins].asFloat();
		}
		TF1 *f = new TF1("f", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)",mean-2.8*sigl, mean + 1.7 * sigr);
		myfile2 >> mom1 >> mom2 >> c >> mean >> sigl >> sigr;
		cout<<mom1<<endl;
		//TF1 *f = new TF1("landau", "landau", mean-1.5*sigl, mean + 1.5 * sigr);
		energy.clear();
		data.clear();
		//myfile2>>mom1>>mom2>>c>>mean>>sigl>>sigr;
		//myfile<<mom1<<mom2<<c<<mean<<sigl<<sigr<<endl;
		//f->SetParameters(c, mean, sigl, sigr);
		//f->SetLineColor(kGreen);
		//TF1 *f2 = new TF1("f2", "landau", mean + 4 * sigr, 15);
		f->SetLineColor(kGreen);
		//f->SetParNames("Height", "Position", "LeftSigma", "RightSigma");
		//f2->SetParameters(1752, 2.7, 0.8, 1.4);
		f->SetParameters(c, mean, sigl, sigr);
		f->SetParLimits(1,mean-1.2,mean);
		//f->SetParLimits(2,sigl-0.4,sigl+0.4);
		//f->SetParLimits(3,sigr,sigr);
		//f->FixParameter(1,mean);
		f->FixParameter(2,sigl);
		f->FixParameter(3,sigr);
		TF1 *f4= new TF1("twogaus", "[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)+[4]* TMath::Gaus(x, [5], ((x < [5]) ? [6] : [7]), 0)", 0.0, 12);
		hf4[i] = hh6->ProjectionY(TString::Format("PionPositive_%d_%d", (i-1)*25,i*25), (i - 1) * 25, i * 25);
		hf4[i]->SetTitle(TString::Format("Pion Positive at %d -%d [MeV] ", (i - 1) * 25, i * 25));
		// TF1 *g1 = new TF1 ("g1","landau",mean-sigl*3.2,mean-sigl*0.1);
		TF1 *g1 = new TF1("pol", "pol1", mean - sigl * 3.2, mean + sigr * 3.2);
		ebin=hf4[i]->FindBin(mean-sigl*2.8);
		for(int j=0;j<3;j++){
			energy.push_back(hf4[i]->GetBinCenter(ebin-j));
			data.push_back(hf4[i]->GetBinContent(ebin-j));
		}
		ebin=hf4[i]->FindBin(mean+sigr*2.4);
		for(int j=0;j<3;j++){
			energy.push_back(hf4[i]->GetBinCenter(ebin+j));
			data.push_back(hf4[i]->GetBinContent(ebin+j));
		}
		TF1 *g2 = new TF1 ("g2","[0] * TMath::Gaus(x, [1], ((x < [1]) ? [2] : [3]), 0)",mean,14);
		if(i-pbins>=0){
		g2->SetParLimits(1,pmean-0.1,pmean+0.1);
		g2->FixParameter(2,plsig);
		g2->SetParLimits(3,prsig-0.1,prsig+0.1);}
		g[i]= new TGraph(energy.size(), &energy[0], &data[0]);
		if (i <= 35)
		{
			c6->cd(i);
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf4[i]->GetYaxis()->SetTitleSize(0.11);
			hf4[i]->GetYaxis()->SetLabelSize(0.11);
			hf4[i]->GetXaxis()->SetLabelSize(0.11);
			hf4[i]->GetXaxis()->SetTitleSize(0.11);
			hf4[i]->GetYaxis()->SetTitleOffset(0.80);
			if(ang>45) hf4[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf4[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf4[i]->GetYaxis()->SetTitle("Counts (a.u)");
	// 		//gPad->SetBottomMargin(0.17);
	// 		//gPad->SetTopMargin(0);
	// 		//gPad->SetRightMargin(0.2);
	// 		//gPad->SetLeftMargin(0.17);
			hf4[i]->Draw();
			hf4[i]->GetXaxis()->SetRangeUser(0,hf4[i]->GetBinCenter(hf4[i]->FindBin(mean+sigr*5)));
			//hf3[i]->Fit(f0, "RMQ");
			g[i]->SetMarkerStyle(24);
			g[i]->SetMarkerColor(kRed);
			//g[i]->Fit(g1,"RMQ");
			g[i]->Draw("p");
			// hf4[i]->Fit(g1, "RMQ");
			hf4[i]->Fit(f, "RMQ");
			if(i>=pbins)hf4[i]->Fit(g2, "RMQ+");
			f->GetParameters(&par[0]);
			g2->GetParameters(&par[4]);
			f4->SetParameters(&par[0]);
			f4->FixParameter(5,par[5]);
			f4->SetParLimits(1,mean-1.2,mean);
			f4->SetParLimits(6,plsig,plsig+0.2);
			f4->SetParLimits(2,sigl-0.1,sigl+0.1);
			f4->SetParLimits(3,sigr-0.1,sigr+0.1);
			hf4[i]->Fit(f4, "RQ");
			f4->GetParameters(&par[0]);
			f->SetParameters(&par[0]);
			g2->SetParameters(&par[4]);
			f->SetRange(0,14);
			g2->SetRange(0,14);
			g2->Draw("same");
			f->Draw("same");
			f->GetParameters(&par[0]);
	// 		hf4[i]->GetYaxis()->SetTitleSize(0.05);
	// 		hf4[i]->GetYaxis()->SetLabelSize(0.05);
	// 		hf4[i]->GetXaxis()->SetLabelSize(0.05);
	// 		hf4[i]->GetXaxis()->SetTitleSize(0.05);
	// 		hf4[i]->GetYaxis()->SetTitleOffset(0.90);
	// 		hf4[i]->GetXaxis()->SetTitleOffset(0.90);
	// 		//f->Draw("same");
	// 		hf4[i]->Fit(f, "R");
	// 		//hf3[i]->Fit(g1,"R+");
	// 		//hf4[i]->Fit(f2,"R+");

	// 		f->GetParameters(&par[0]);
	// 		//	g1-> GetParameters (&par[4]);
	// 		f2->GetParameters(&par[4]);
	// 		TF1 *f1 = new TF1("double_Landau", "landau(0)+landau(4)", 0.0, 12); //f->Draw("same");
	// 		f1->SetParameters(par);
	// 		f1->SetLineColor(6);
	// 		//f1->FixParameter(8,par[8]);
	// 		//f1->FixParameter(9,par[9]);
	// 		//f1->FixParameter(7,par[7]);
	// 		hf4[i]->Fit(f1, "+", "", 0.0, 10);
	// 		//int gaussup=hf[i]->GetXaxis()->FindBin(par[1]+par[3]*sigma);
	// 		//int gausslow=hf[i]->GetXaxis()->FindBin(par[1]-par[2]*sigma);
			//intergral_Background = g1->Integral(par[1] - par[2] * sigma, par[1] + par[3] * sigma);
			intergral_Background = g2->Integral(mean - sigl * sigma, mean + sigr * sigma);
			intergral_signal = f->Integral(mean - sigl * sigma, mean + sigr * sigma) ;
			//float intergral_signal2= hf[i]->Integral(gausslow,gaussup)-f1->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma);
			//intergral_signal = f->Integral(par[1] - par[2] * sigma, par[1] + par[3] * sigma) - g1->Integral(par[1] - par[2] * sigma, par[1] + par[3] * sigma);
		}
		else
		{
			c7->cd(i - 35);
			gPad->SetLogy();;
			gPad->SetRightMargin(0.01);
			gPad->SetLeftMargin(0.2);
			gPad->SetTopMargin(0.1);
			gPad->SetBottomMargin(0.3);
			hf4[i]->GetYaxis()->SetTitleSize(0.11);
			hf4[i]->GetYaxis()->SetLabelSize(0.11);
			hf4[i]->GetXaxis()->SetLabelSize(0.11);
			hf4[i]->GetXaxis()->SetTitleSize(0.11);
			hf4[i]->GetYaxis()->SetTitleOffset(0.80);
			if(ang>45) hf4[i]->GetXaxis()->SetTitle("dE/dx signal in TOF (a.u)");
			if(ang<45) hf4[i]->GetXaxis()->SetTitle("dE/dx signal in TOFino (a.u)");
			hf4[i]->GetYaxis()->SetTitle("Counts (a.u)");
	// 		gPad->SetBottomMargin(0.17);
	// 		gPad->SetTopMargin(0);
	// 		gPad->SetRightMargin(0.2);
	// 		gPad->SetLeftMargin(0.17);
			// hf4[i]->Draw();
			g[i]->SetMarkerStyle(24);
			g[i]->SetMarkerColor(kRed);
			g[i]->Fit(g1,"RMQ");
			//hf3[i]->Fit(f0, "RMQ");
			// hf4[i]->Fit(g1, "RMQ");
			g[i]->Draw("p");
			hf4[i]->Fit(f, "RMQ");
			if(i>=pbins)hf4[i]->Fit(g2, "RMQ+");
			f->GetParameters(&par[0]);
			g2->GetParameters(&par[4]);
			f4->SetParameters(&par[0]);
			f4->FixParameter(5,par[5]);
			f4->SetParLimits(2,mean,mean+0.4);
			f4->SetParLimits(2,sigl-0.1,sigl);
			f4->SetParLimits(3,sigr-0.2,sigr+0.2);
			hf4[i]->Fit(f4, "RQ");
			f4->GetParameters(&par[0]);
			f->SetParameters(&par[0]);
			g2->SetParameters(&par[4]);
			f->SetRange(0,14);
			g2->SetRange(0,14);
			g2->Draw("same");
			f->Draw("same");
			hf4[i]->GetXaxis()->SetRangeUser(0,hf4[i]->GetBinCenter(hf4[i]->FindBin(mean+sigr*5)));
	// 		hf4[i]->GetYaxis()->SetTitleSize(0.1);
	// 		hf4[i]->GetYaxis()->SetLabelSize(0.1);
	// 		hf4[i]->GetXaxis()->SetLabelSize(0.1);
	// 		hf4[i]->GetXaxis()->SetTitleSize(0.1);
	// 		hf4[i]->GetYaxis()->SetTitleOffset(0.70);
	// 		hf4[i]->GetXaxis()->SetTitleOffset(0.96);
	// 		//f->Draw("same");
	// 		hf4[i]->Fit(f, "R");
	// 		//hf3[i]->Fit(g1,"R+");
	// 		hf4[i]->Fit(f2, "R+");

	// 		f->GetParameters(&par[0]);
	// 		//g1-> GetParameters (&par[4]);
	// 		f2->GetParameters(&par[4]);
	// 		TF1 *f1 = new TF1("double_landau", "landau(0)+landau(4)", 0.0, 12); //f->Draw("same");
	// 		f1->SetParameters(par);
	// 		f1->SetLineColor(6);
	// 		//f1->FixParameter(8,par[8]);
	// 		//f1->FixParameter(9,par[9]);
	// 		//f1->FixParameter(7,par[7]);
	// 		hf4[i]->Fit(f1, "+", "", 0.0, 10);
	// 		//int gaussup=hf[i]->GetXaxis()->FindBin(par[1]+par[3]*sigma);
	// 		//int gausslow=hf[i]->GetXaxis()->FindBin(par[1]-par[2]*sigma);
		//intergral_Background = g1->Integral(par[1] - par[2] * sigma, par[1] + par[3] * sigma);
	// 		//float intergral_signal2= hf[i]->Integral(gausslow,gaussup)-f1->Integral(par[1]-par[2]*sigma,par[1]+par[3]*sigma);
		//intergral_signal = f->Integral(par[1] - par[2] * sigma, par[1] + par[3] * sigma) - g1->Integral(par[1] - par[2] * sigma, par[1] + par[3] * sigma);
		intergral_Background = g2->Integral(mean - sigl * sigma, mean + sigr * sigma);
		intergral_signal = f->Integral(mean - sigl * sigma, mean + sigr * sigma);
		}
		hf4[i]->Write();
		g[i]->Write(TString::Format("PionPositive_%d_%d_graph", (i-1)*25,i*25));
		//myfile << (i - 1) * 25 << " " << i * 25 << " " << intergral_Background << " " << intergral_signal << " " << (intergral_Background / (intergral_Background + intergral_signal)) * 100 << endl;
		binpip["fistbin"][i-1]=(i-1)*25;
		binpip["lastbin"][i-1]=(i)*25;
		if((intergral_Background / intergral_signal)<1&&(intergral_Background / intergral_Background)>=0)
		{
			binpip["background"][i-1]=intergral_Background;
			binpip["signal"][i-1]=intergral_signal;
			binpip["backper"][i-1]= (intergral_Background / intergral_signal) * 100 ;
		}
		else
		{
			binpip["background"][i-1]=0;
			binpip["signal"][i-1]=0;
			binpip["backper"][i-1]= 0;
		}
	}
    obj["PiP"]=binpip;
	c7->Write();
 	c6->Write();
	myfile << styledWriter.write(obj);
	//c1->cd();
	// hh1->Draw("colz");
	// hh1->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh1->GetYaxis()->SetTitle("Beta (C)");
	// hh1->SetTitle("Using mdc and tof dE/dx(system 1) cut for Proton ");
	// hh1->SetLineColor(kRed);
	// c2->cd();
	// hh2->Draw("colz");
	// hh2->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh2->GetYaxis()->SetTitle("Beta (C)");
	// hh2->SetTitle("Using mdc and tof dE/dx(system 0) cut for Deuteron ");
	// hh2->SetLineColor(kRed);
	//
	// // mass->Draw();
	// // mass->GetXaxis()->SetTitle("Mass(MeV/c^{2})");
	// // mass->GetYaxis()->SetTitle("Counts(a.u)");
	// // mass->SetTitle("Mass spectra with dueteron cut for P + Nb at 3.5GeV");
	// c3->cd();
	// hh3->Draw("colz");
	// hh3->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh3->GetYaxis()->SetTitle("Beta (C)");
	// hh3->SetTitle("Using mdc and tof dE/dx(system 1) cut for deutron");
	// // mass2->Draw();
	// // mass2->GetXaxis()->SetTitle("Mass(MeV/c^{2})");
	// // mass2->GetYaxis()->SetTitle("Counts(a.u)");
	// // mass2->SetTitle("Mass spectra with triton cut for P + Nb at 3.5GeV");
	// c4->cd();
	// hh4->Draw("colz");
	// hh4->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh4->GetYaxis()->SetTitle("Beta (C)");
	// hh4->SetTitle("Using mdc and tof dE/dx(system 0) cut for Triton");
	// c5->cd();
	// hh5->Draw("colz");
	// hh5->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh5->GetYaxis()->SetTitle("Beta (C)");
	// hh5->SetTitle("Using mdc and tof dE/dx(system 1) cut for triton");
	// c6->cd();
	// hh6->Draw("colz");
	// hh6->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh6->GetYaxis()->SetTitle("Beta (C)");
	// hh6->SetTitle("Using mdc and tof dE/dx(system 0) cut for Pion Positive");
	// hh6->SetLineColor(kRed);
	//
	//c7->cd();
	// hh4->Draw("colz");
	// hh7->GetXaxis()->SetTitle("Momentum p(MeV/c)");
	// hh7->GetYaxis()->SetTitle("Beta (C)");
	// hh7->SetTitle("Using mdc and tof dE/dx(system 1) cut for Pion Positive");
	// hh7->SetLineColor(kRed);
	// /*c4->cd(1);
	// hh3_0->Draw();
	// c4->cd(2);
	// hh3_1->Draw();
	// c4->cd(3);
	// hh3_2->Draw();
	// c4->cd(4);
	// hh3_3->Draw();*/
	//deuteron->Draw();
	// mass3->Draw();
	// mass3->GetXaxis()->SetTitle("Mass(MeV/c^{2})");
	//mass3->GetYaxis()->SetTitle("Counts(a.u)");
	//mass3->SetTitle("Mass spectra with proton cut for P + Nb at 3.5GeV");
}
double EnergyLossFitFunc2(double * x_val, double * par) {
        // Calculates the dEdx (MeV cm2/g ) of an particle with GEANT ID id
        // and momentum p (MeV) for He/i-butan gas mixture with He fraction hefr
        // (He (hefr) / i-butan (1-hefr)). The fomular is taken from PDG and doesn't
        // include the density correction term. The values for the mean excitation
        // energy are taken from Sauli.

        double p,  hefr, par1, par2, par3, par4, par5, par6, par7, par8, par9, par10, par11, par12;
        p = x_val[0];
        hefr = par[0];
        par1 = par[1];
        par2 = par[2];
        par3 = par[3];
        par4 = par[4];
        par5 = par[5];
        par6 = par[6];
        par7 = par[7];
        par8 = par[8];
        par9 = par[9];
        par10 = par[10];
        par11 = par[11];
        par12 = par[12];
        if (p==0)                        return -1;
        if (hefr<0.||hefr>1.) return -1;
        double mass   = par7 * par[8];
        if (mass==0) return -1;

        // Z and A
        double Z_gas  =par3*hefr+(1.-hefr)*par1; // 14
        double A_gas  =par4*hefr+(1.-hefr)*par2; // 58

        // I_0 and I
        double I_0_gas=par5*hefr+(1.-hefr)*par6;
        double I2        =pow(I_0_gas*Z_gas*(1.e-6),2); // sauli

        double K          =0.307075; // MeV cm2 PDG, 4*pi*N_{A}*r_{e}2*m_{e}2*c2
        double mass2  =pow(mass,2);
        double m_ec2  =0.5109989461;
        double z2        =pow(par12,2);
        double p2        =pow(p,2);
        double beta2  =1./((mass2/p2)+1.);
        double gamma2 =1./(1-beta2);
        double gamma  =sqrt(gamma2);

        double Tmax   =(2.*m_ec2*beta2*gamma2)/(1.+ 2.*gamma*(m_ec2/mass)+pow((m_ec2/mass),2));
        double term1  =K*z2*(Z_gas/A_gas)*(1./beta2) + par11;
        double term2  =((2.*m_ec2*beta2*gamma2*Tmax)/I2);
        double dedx   = par9 * term1 * (0.5*log(term2)-beta2) + par10;

        return dedx;
}

