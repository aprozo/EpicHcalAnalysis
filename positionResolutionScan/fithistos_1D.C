#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>

Double_t fline(Double_t *x, Double_t *par)
{

	/*const int numHistograms = 8;
	const char *histogramNames[numHistograms] = {"diffD034w", "diffD023w", "diffD012w", "diffD01_52w", "diffD00_91_5w", "diffD001w", "diffD00_51w", "diffD000_5w" };
	double ranges[numHistograms][2] = {{1.83, 1.88}, {1.84, 1.92}, {1.84, 1.91}, {1.84, 1.91}, {1.83, 1.92}, {1.84, 1.91}, {1.84, 1.90}, {1.85, 1.91}};
	 Int_t i = 7;
	 if (x[0] >= ranges[8-i][0] && x[0] <= ranges[8-i][1])
	 */
	/////return par[0]*exp(-0.5*((x[0] - par[1])/par[2])*((x[0] - par[1])/par[2]));
	////if (x[0] >= ranges[8-i][1] && x[0] <= 2.05) return par[0] + (par[1] * x[0]);
	////if (x[0] >= 1.75 && x[0] <= ranges[8-i][0])
	return par[0] + (par[1] * x[0]);
}

void RenormalizeXY(TH2D *h, TH2D *hnorm, double a, double b)
{

	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();

	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{

			double x = h->GetXaxis()->GetBinCenter(i);
			double y = h->GetYaxis()->GetBinCenter(j);

			double counts = h->GetBinContent(i, j);

			double normx = x * a;
			double normy = y * b;

			int ixn = hnorm->GetXaxis()->FindBin(normx);
			int iyn = hnorm->GetYaxis()->FindBin(normy);

			hnorm->SetBinContent(ixn, iyn, counts);
		}
	}
}

void RemoveEntries(TH2D *h, double a, double b)
{

	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();

	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{

			double x = h->GetXaxis()->GetBinCenter(i);
			double y = h->GetYaxis()->GetBinCenter(j);

			/////////if (y > 9.5) h->SetBinContent(i, j, 0.0);
			/////////if ((x > 2.0 || x < 0.5)) h->SetBinContent(i, j, 0.0);
			if (y > a * x + b + 4)
				h->SetBinContent(i, j, 0.0);
			if (y < a * x + b - 4)
				h->SetBinContent(i, j, 0.0);
		}
	}
}

Double_t fgaus(Double_t *x, Double_t *par)
{
	return par[0] / (TMath::Sqrt(2 * TMath::Pi() * par[2] * par[2])) * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / (2 * par[2] * par[2]));
}

void fithistos_1D(TString filename)
{

	//==========Style of the plot============
	gStyle->SetPalette(1);
	gStyle->SetOptTitle(1);
	gStyle->SetTitleOffset(.85, "X");
	gStyle->SetTitleOffset(.85, "Y");
	gStyle->SetTitleSize(.04, "X");
	gStyle->SetTitleSize(.04, "Y");
	gStyle->SetLabelSize(.04, "X");
	gStyle->SetLabelSize(.04, "Y");
	gStyle->SetHistLineWidth(2);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	char *parstring = new char[50];

	TFile *file = new TFile(filename, "read");

	TH1D *h_reco_E = (TH1D *)file->Get("h_NHcal_cluster_thetaresol");
	h_reco_E->Rebin(2);
	// h_reco_E->Update();
	cout << "h_reco_E->GetMaximum()     " << h_reco_E->GetMaximum() << endl;
	Int_t bin1 = h_reco_E->FindFirstBinAbove(h_reco_E->GetMaximum() / 2);
	Int_t bin2 = h_reco_E->FindLastBinAbove(h_reco_E->GetMaximum() / 2.);
	Double_t fwhm = h_reco_E->GetBinCenter(bin2) - h_reco_E->GetBinCenter(bin1);
	TCanvas *sampfrac = new TCanvas("sampfrac");
	sampfrac->cd();
	TF1 *fitfunc = new TF1("fitfunc", fgaus, 0.02, 0.08, 3);
	fitfunc->SetParameters(1000, 5., 1.1);
	/// fitfunc->FixParameter(1,5.);
	fitfunc->SetParNames("Area", "Mean", "Sigma");

	/////h_reco_E->GetYaxis()->SetRangeUser(-10.0, 3500.0);
	h_reco_E->GetXaxis()->SetRangeUser(0.0, 15.0);
	h_reco_E->GetYaxis()->SetTitleOffset(1.25);
	h_reco_E->SetMarkerStyle(43);
	h_reco_E->SetMarkerSize(1.5);
	h_reco_E->SetMarkerColor(kBlue);
	h_reco_E->Fit("fitfunc", "MR", "", ((h_reco_E->GetBinCenter(bin1)) - 0.), ((h_reco_E->GetBinCenter(bin2)) + 0.));
	h_reco_E->Draw("p");

	TLatex tl;
	tl.SetTextSize(0.1);
	/// tl.SetTextStyle(4);
	tl.SetTextColor(1);
	tl.DrawLatex(0.05, 5000, parstring);
	sprintf(parstring, "Mean ~ %.2f", fitfunc->GetParameter(1));
	tl.DrawLatex(0.05, 0.55, parstring);
	sprintf(parstring, "FWHM ~ %.2f", fwhm);
	tl.DrawLatex(0.05, 3000, parstring);
	////sprintf(parstring,"#chi^{2}/ndf = %.2f/%i",fitfunc->GetChisquare(),fitfunc->GetNDF());

	// TString outputfile = filename;
	// outputfile.ReplaceAll(".root", "");
	// outputfile.ReplaceAll("outputFill", "");

	sampfrac->SaveAs("resolution.pdf");
}