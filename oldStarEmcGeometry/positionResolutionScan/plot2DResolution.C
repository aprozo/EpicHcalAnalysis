#include <iostream>
#include <fstream>
#include <set>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TGraph2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

const Double_t globalZ = 383.95; // first HCAL layer
const Double_t pi = 3.14159265358979323846;

Double_t getEta(const Double_t &theta) { return -log(tan((pi * theta / 180) / 2)); }
Double_t getTheta(const Double_t &x, const Double_t &y) { return atan2(sqrt(x * x + y * y), globalZ); }
Double_t getPhi(const Double_t &x, const Double_t &y) { return atan2(y, x); }

Double_t getX(const Double_t &theta, const Double_t &phi) { return globalZ * tan(pi * theta / 180) * cos(pi * phi / 180); }
Double_t getY(const Double_t &theta, const Double_t &phi) { return globalZ * tan(pi * theta / 180) * sin(pi * phi / 180); }

const double etaFullRange[26] =
    {3.48733, 3.33696, 3.19285, 3.05481, 2.92264, 2.79614, 2.67511, 2.55935, 2.44868, 2.34289,
     2.24167, 2.14506, 2.05272, 1.96465, 1.88032, 1.79973, 1.72258, 1.64865,
     1.57776, 1.50981, 1.44473, 1.38221, 1.3219, 1.26371};

struct TileHists
{
    TileHists()
    {
    }

    TileHists(TString name)
    {
        tilePhiTheta = new TGraph2D();
        tilePhiEta = new TGraph2D();
        tileXY = new TGraph2D();
        phiScan = new TGraph();
        thetaScan = new TGraph();

        TString titleZ;
        if (name.Contains("phi"))
            titleZ = "#Delta#phi/#phi";
        else
            titleZ = "#Delta#theta/#theta";

        tilePhiTheta->SetTitle(name + " ;#phi, deg; #theta, deg ;" + titleZ);
        tilePhiEta->SetTitle(name + " ;#phi, deg; #eta;" + titleZ);
        tileXY->SetTitle(name + "  ;x,mm;y,mm;" + titleZ);
        phiScan->SetTitle(name + " ;#phi, deg;" + titleZ);
        thetaScan->SetTitle(name + " ;#theta, deg;" + titleZ);
        tilePhiTheta->SetName(name + "_PhiTheta");
        tilePhiEta->SetName(name + "_PhiEta");
        tileXY->SetName(name + "_XY");
        phiScan->SetName(name + "PhiScan");
        thetaScan->SetName(name + "ThetaScan");

        tilePhiThetaHist = new TH2D(name + "PhiTheta", name + "PhiTheta", 21, 58, 90, 60, 0, 60);
    }

    void Fill(const Double_t phi, const Double_t theta, const Double_t fwhm)
    {
        Double_t eta = getEta(theta);
        Double_t x = globalZ * tan((180 - theta) * pi / 180) * cos(phi * pi / 180);
        Double_t y = globalZ * tan((180 - theta) * pi / 180) * sin(phi * pi / 180);

        tilePhiEta->SetPoint(tilePhiEta->GetN(), phi, eta, fwhm);
        tilePhiTheta->SetPoint(tilePhiTheta->GetN(), phi, theta, fwhm);
        tileXY->SetPoint(tileXY->GetN(), x, y, fwhm);
        if (phi < 60)
            phiScan->SetPoint(phiScan->GetN(), phi, fwhm);
    }

    void Draw(TString outPdf, TCanvas *can)
    {
        // tilePhiTheta->GetYaxis()->SetNdivisions(505);
        // tilePhiTheta->GetZaxis()->SetNdivisions(505);
        tilePhiTheta->Draw("colz");
        can->SaveAs(outPdf);

        tilePhiEta->Draw("colz");
        can->SaveAs(outPdf);

        tileXY->GetYaxis()->SetNdivisions(305);
        tileXY->Draw("colz");
        can->SaveAs(outPdf);
    }

    void Write()
    {
        tilePhiTheta->Write();
        tilePhiEta->Write();
        tileXY->Write();
    }

    TGraph2D *tilePhiTheta;
    TGraph2D *tilePhiEta;
    TGraph2D *tileXY;

    TH2D *tilePhiThetaHist;
    TH2D *tilePhiEtaHist;

    TGraph *phiScan;
    TGraph *thetaScan;
};

void plot2DResolution(TString outFile = "sectorResolution")
{

    // gStyle->SetTitleOffset(1.2, "y");
    set<pair<TString, TString>> tileMap;

    fstream fin;
    fin.open("tileMap.txt", ios::in);
    if (fin.is_open())
    {
        cout << "File Opened successfully!!!. Reading data from file into array" << endl;
        while (!fin.eof())
        {
            TString phi, theta;
            fin >> phi >> theta;
            tileMap.insert({phi, theta});
        }
        fin.close();
    }
    else
    {
        cout << "File Not Found" << endl;
    }

    cout << "Number of points : " << tileMap.size() << endl;

    TFile *inputFile = new TFile("output.root", "READ");

    const Int_t nBins = 12;
    TileHists resolution[nBins];
    TString labels[nBins] = {"hcal+ecal #theta", "hcal+ecal #phi",
                             "hcal #theta", "hcal #phi",
                             "ecal #theta", "ecal #phi",
                             "hcal+ecal #theta Truth", "hcal+ecal #phi Truth",
                             "hcal #theta Truth", "hcal #phi Truth",
                             "ecal #theta Truth", "ecal #phi Truth"};

    for (int i = 0; i < nBins; i++)
    {
        resolution[i] = TileHists(labels[i]);
    }

    const double PI = TMath::Pi();

    for (auto &angle : tileMap)
    {
        cout << angle.first << " " << angle.second << endl;

        TString angleStr = "Phi" + angle.first + "Theta" + angle.second;

        TH1D *hist = (TH1D *)inputFile->Get("hResolution" + angleStr);
        if (!hist)
        {
            cout << "Histogram not found" << ("hResolution" + angleStr) << endl;
            continue;
        }

        Double_t phi = angle.first.Atof();
        Double_t theta = angle.second.Atof();

        // Double_t x = getX(theta, phi);
        // Double_t y = getY(theta, phi);
        // x -= 1; // remove disksGap
        // phi = 180 * getPhi(x, y) / PI;
        // theta = 180 * getTheta(x, y) / PI;
        // if (phi < 0)
        //     phi += 180;

        cout << phi << " " << theta << endl;

        for (int i = 0; i < nBins; i++)
        {
            Float_t fwhm = hist->GetBinContent(i + 1);
            Double_t upperCut = 2.;
            if (i == 0 || i == 6) // theta hcal+ecal Resolution
                upperCut = 0.1;

            if (fwhm > 0. && fwhm < upperCut)
                resolution[i].Fill(phi, theta, hist->GetBinContent(i + 1));
        }
    }
    TCanvas *can = new TCanvas("can", "can", 1200, 1000);
    can->cd();
    TString outPdf = outFile + ".pdf";
    can->SaveAs(outPdf + "[");
    gPad->SetRightMargin(0.2);

    for (int i = 0; i < nBins; i++)
    {
        resolution[i].Draw(outPdf, can);
    }

    can->SaveAs(outPdf + "]");

    TFile *output = new TFile("sectorResolution.root", "RECREATE");
    output->cd();
    for (int i = 0; i < nBins; i++)
    {
        resolution[i].Write();
    }

    output->Close();
}