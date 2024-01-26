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

const double globalZ = -395; // first HCAL layer

const double rmin = 14;  // first HCAL layer
const double rmax = 270; // first HCAL layer

const Double_t pi = 3.14159265358979323846;
double getTheta(const double &eta) { return 2 * atan(exp(-eta)); }
Double_t getEta(const Double_t &theta) { return -log(tan((pi * theta / 180) / 2)); }
double getEta(const double &z, const double &r) { return asinh(z / r); }
Double_t getTheta(const Double_t &x, const Double_t &y) { return atan2(sqrt(x * x + y * y), globalZ); }
Double_t getPhi(const Double_t &x, const Double_t &y) { return atan2(y, x); }

Double_t getX(const Double_t &theta, const Double_t &phi) { return globalZ * tan(pi * theta / 180) * cos(pi * phi / 180); }
Double_t getY(const Double_t &theta, const Double_t &phi) { return globalZ * tan(pi * theta / 180) * sin(pi * phi / 180); }

const int nPhiBins = 30;
const int nEtaBins = 20;

const double phiMin = 0;
const double phiMax = 90;
const double phiStep = (phiMax - phiMin) / (double)nPhiBins;

const double etaMin = getEta(globalZ, rmax);
const double etaMax = getEta(globalZ, rmin);
const double etaStep = (etaMax - etaMin) / (double)nEtaBins;

vector<double> thetaBins;
vector<double> phiBins;

struct ClusterHists
{

    ClusterHists()
    {
    }
    ClusterHists(TString _name, TFile *inputFile, TString folder)
    {
        hTheta = (TH1D *)inputFile->Get(folder + "/" + _name + "/hTheta_" + _name);
        if (!hTheta)
            cout << "hTheta not found" << endl;
        hPhi = (TH1D *)inputFile->Get(folder + "/" + _name + "/hPhi_" + _name);
        if (!hPhi)
            cout << "hPhi not found" << endl;

        hThetaResol = (TH1D *)inputFile->Get(folder + "/" + _name + "/hThetaResolution_" + _name);
        if (!hThetaResol)
            cout << "hThetaResol not found" << endl;
        hPhiResol = (TH1D *)inputFile->Get(folder + "/" + _name + "/hPhiResolution_" + _name);
        if (!hPhiResol)
            cout << "hPhiResol not found" << endl;
        hEnergy = (TH1D *)inputFile->Get(folder + "/" + _name + "/hEnergy_" + _name);
        if (!hEnergy)
            cout << "hEnergy not found" << endl;

        name = _name;
    }
    TString name;
    TH1D *hTheta;
    TH1D *hPhi;
    TH1D *hThetaResol;
    TH1D *hPhiResol;
    TH1D *hEnergy;
};

struct AnglePoint
{
    AnglePoint(const pair<TString, TString> _angle, TFile *_inputFile)
    {
        phi = _angle.first.Atof();
        theta = _angle.second.Atof();
        angleStr = "Phi" + _angle.first + "Theta" + _angle.second;
        inputFile = _inputFile;

        hcalRecoHist = {"HCal_Reco", inputFile, angleStr};
        ecalRecoHist = {"ECal_Reco", inputFile, angleStr};
        //  ebarellRecHist={"EBarell_Reco", inputFile, angleStr};
        hcalAndEcalSumHist = {"HcalAndEcalSum_Reco", inputFile, angleStr};
        scatteredHcal = {"scattered_Hcal", inputFile, angleStr};

        hcalTruthHist = {"HCal_Truth", inputFile, angleStr};
        ecalTruthHist = {"ECal_Truth", inputFile, angleStr};
        //  ebarellTruthHist={"EBarell_Truth", inputFile, angleStr};
        hcalAndEcalSumTruthHist = {"HcalAndEcalSum_Truth", inputFile, angleStr};
        scatteredHcalTruth = {"scattered_Hcal_Truth", inputFile, angleStr};
    }

    Double_t getHcalOnlyToScatteredHcal() const
    {
        if (scatteredHcal.hTheta->GetEntries() == 0)
            return 0;
        Double_t ratio = hcalRecoHist.hTheta->GetEntries() / scatteredHcal.hTheta->GetEntries();
        return ratio;
    }

    TFile *inputFile;
    Double_t phi;
    Double_t theta;
    TString angleStr;
    TString folder;

    ClusterHists hcalRecoHist;
    ClusterHists ecalRecoHist;
    // ClusterHists ebarellRecHist;
    ClusterHists hcalAndEcalSumHist;
    ClusterHists scatteredHcal;

    ClusterHists hcalTruthHist;
    ClusterHists ecalTruthHist;
    // ClusterHists ebarellTruthHist;
    ClusterHists hcalAndEcalSumTruthHist;
    ClusterHists scatteredHcalTruth;
};

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

        tileRatioHcalOnlyToScatteredHcal = new TGraph2D();

        // TH2D *tileMap = new TH2D("tileMap", "tileMap; phi ;theta", phiBins.size() - 1, &phiBins[0], thetaBins.size() - 1, &thetaBins[0]);
        TString titleZ;
        if (name.Contains("phi"))
            titleZ = "#Delta#phi, deg  ";
        else
            titleZ = "#Delta#theta, deg ";

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

        tileRatioHcalOnlyToScatteredHcal->SetName(name + "RatioHcalOnlyToScatteredHcal");
        tileRatioHcalOnlyToScatteredHcal->SetTitle(name + "HcalOnly/Scattered; #phi, deg;#theta, deg;ratio");
    }

    void Fill(const AnglePoint &anglePoint)
    {
        Double_t ratio = anglePoint.getHcalOnlyToScatteredHcal();
        cout << "ratio : " << ratio << endl;
        // Int_t binx = tileRatioHcalOnlyToScatteredHcal->GetXaxis()->FindBin(anglePoint.phi);
        // Int_t biny = tileRatioHcalOnlyToScatteredHcal->GetYaxis()->FindBin(anglePoint.theta);
        tileRatioHcalOnlyToScatteredHcal->SetPoint(tileRatioHcalOnlyToScatteredHcal->GetN(), anglePoint.phi, anglePoint.theta, ratio);
        // tileRatioHcalOnlyToScatteredHcal->SetBinContent(binx, biny, ratio);
    }

    void Fill(const Double_t phi, const Double_t theta, const Double_t sigma)
    {
        Double_t eta = getEta(theta);
        Double_t x = globalZ * tan((180 - theta) * pi / 180) * cos(phi * pi / 180);
        Double_t y = globalZ * tan((180 - theta) * pi / 180) * sin(phi * pi / 180);

        tilePhiEta->SetPoint(tilePhiEta->GetN(), phi, eta, sigma);
        tilePhiTheta->SetPoint(tilePhiTheta->GetN(), phi, theta, sigma);
        tileXY->SetPoint(tileXY->GetN(), x, y, sigma);
        // if (phi < 60)
        //     phiScan->SetPoint(phiScan->GetN(), phi, sigma);
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

    TGraph2D *tileRatioHcalOnlyToScatteredHcal;

    // TH2D *tilePhiThetaHist;
    // TH2D *tilePhiEtaHist;
    // TH2D *tileRatioHcalOnlyToScatteredHcal;

    TGraph *phiScan;
    TGraph *thetaScan;
};

void plot2DResolution(TString outFile = "sectorResolution")
{

    for (double iEta = 0; iEta <= nEtaBins; iEta++) // 3 degree scan
    {
        double eta = etaMin + iEta * etaStep;
        double theta = 180 * getTheta(eta) / TMath::Pi();
        thetaBins.push_back(theta);
    }

    

    for (double iPhi = 0; iPhi <= nPhiBins; iPhi++) // 3 degree scan
    {
        double phi = phiMin + iPhi * phiStep;
        phiBins.push_back(phi);
    }

    // gStyle->SetTitleOffset(1.2, "y");
    vector<pair<TString, TString>> tileMap;

    fstream fin;
    fin.open("tileMap.txt", ios::in);
    if (fin.is_open())
    {
        cout << "File Opened successfully! Reading data from file into array" << endl;
        while (!fin.eof())
        {
            TString phi, theta;
            fin >> phi >> theta;
            tileMap.push_back({phi, theta});
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

        TString angleStr = "Phi" + angle.first + "Theta" + angle.second;
        cout << angleStr << endl;

        TH1D *hist = (TH1D *)inputFile->Get("hResolution" + angleStr);
        if (!hist)
        {
            cout << "Histogram not found" << ("hResolution" + angleStr) << endl;
            continue;
        }

        Double_t phi = angle.first.Atof();
        Double_t theta = angle.second.Atof();

        AnglePoint point(angle, inputFile);

        resolution[0].Fill(point);

        for (int i = 0; i < nBins; i++)
        {
            Float_t sigma = 180 * hist->GetBinContent(i + 1) / 3.14;
            if (sigma < 0.001)
                continue;
            resolution[i].Fill(phi, theta, sigma);
        }
    }
    TCanvas *can = new TCanvas("can", "can", 1200, 1000);
    can->cd();
    TString outPdf = outFile + ".pdf";
    can->SaveAs(outPdf + "[");
    gPad->SetRightMargin(0.2);
    can->SetLogz();
    resolution[0].tileRatioHcalOnlyToScatteredHcal->Draw("colz");
    can->SaveAs(outPdf);

    can->SetLogz(0);
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