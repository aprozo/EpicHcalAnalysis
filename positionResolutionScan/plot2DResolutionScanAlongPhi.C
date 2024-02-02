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
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"

const double alongPhiAngle = 48;

// const double globalZ = -395; // first HCAL layer
const double globalZ = -395 + 55 / 2; // middle HCAL layer
//const double globalZ = -395 + 55 / 2; // 1 interaction length HCAL layer

const double rmin = 14;  // first HCAL layer
const double rmax = 270; // first HCAL layer
const double pi = TMath::Pi();
double getTheta(const double &eta) { return 2 * atan(exp(-eta)); }
double getEta(const double &theta) { return -log(tan((pi * theta / 180) / 2)); }
double getEta(const double &z, const double &r) { return asinh(z / r); }
// double getTheta(const double &x, const double &y) { return atan2(sqrt(x * x + y * y), globalZ); }
// double getPhi(const double &x, const double &y) { return atan2(y, x); }
// double getX(const double &theta, const double &phi) { return globalZ * tan(pi * theta / 180) * cos(pi * phi / 180); }
// double getY(const double &theta, const double &phi) { return globalZ * tan(pi * theta / 180) * sin(pi * phi / 180); }

const int nPhiBins = 30;
const int nEtaBins = 20;

const double maxEcalTheta = 177.758;
const double minEcalTheta = 160.81035;

const double maxEcalEta = getEta(maxEcalTheta);
const double minEcalEta = getEta(minEcalTheta);

const double phiMin = 0;
const double phiMax = 90;
const double phiStep = (phiMax - phiMin) / (double)nPhiBins;

const double etaMin = getEta(globalZ, rmax);
const double etaMax = getEta(globalZ, rmin);
const double etaStep = (etaMax - etaMin) / (double)nEtaBins;

vector<double> thetaBins;
vector<double> phiBins;
vector<double> etaBins;

ostream &operator<<(ostream &os, const vector<double> &v)
{
    for (auto &i : v)
        os << i << " ";
    return os;
}

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
        hcalRecoHist.hThetaResol->GetXaxis()->SetRangeUser(-1, 1);
        scatteredHcal.hThetaResol->GetXaxis()->SetRangeUser(-1, 1);
        Double_t ratio = hcalRecoHist.hThetaResol->GetEntries() / scatteredHcal.hThetaResol->GetEntries();
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

    TileHists(TString _name)
    {
        name = _name;
        TString titleY;
        if (name.Contains("phi"))
            titleY = "#Delta#phi, deg  ";
        else
            titleY = "#Delta#theta, deg ";

        ratioHcalOnlyToScatteredHcalVsEta = new TGraph();
        ratioHcalOnlyToScatteredHcalVsEta->SetName(name + " ratioEta");
        ratioHcalOnlyToScatteredHcalVsEta->SetTitle("HcalOnly / Scattered; #eta ; ratio");

        ratioHcalOnlyToScatteredHcalVsTheta = new TGraph();
        ratioHcalOnlyToScatteredHcalVsTheta->SetName(name + " ratioTheta");
        ratioHcalOnlyToScatteredHcalVsTheta->SetTitle("HcalOnly / Scattered; #theta, deg; ratio");

        resolutionVsEta = new TGraph();
        resolutionVsEta->SetName(name + "resolutionVsEta");
        resolutionVsEta->SetTitle(name + "; #eta ;" + titleY);

        resolutionVsTheta = new TGraph();
        resolutionVsTheta->SetName(name + "resolutionVsTheta");
        resolutionVsTheta->SetTitle(name + "; #theta, deg;" + titleY);
    }

    void setStyle(const Int_t &color = 2001, const Int_t &marker = 20)
    {
        ratioHcalOnlyToScatteredHcalVsEta->SetMarkerColor(color);
        ratioHcalOnlyToScatteredHcalVsEta->SetMarkerSize(2);
        ratioHcalOnlyToScatteredHcalVsEta->SetMarkerStyle(marker);
        ratioHcalOnlyToScatteredHcalVsEta->SetLineColor(color);
        ratioHcalOnlyToScatteredHcalVsEta->SetLineWidth(2);
        ratioHcalOnlyToScatteredHcalVsTheta->SetMarkerColor(color);
        ratioHcalOnlyToScatteredHcalVsTheta->SetMarkerSize(2);
        ratioHcalOnlyToScatteredHcalVsTheta->SetMarkerStyle(marker);
        ratioHcalOnlyToScatteredHcalVsTheta->SetLineColor(color);
        ratioHcalOnlyToScatteredHcalVsTheta->SetLineWidth(2);
        resolutionVsEta->SetMarkerColor(color);
        resolutionVsEta->SetMarkerStyle(marker);
        resolutionVsEta->SetMarkerSize(2);
        resolutionVsEta->SetLineColor(color);
        resolutionVsEta->SetLineWidth(2);

        resolutionVsTheta->SetMarkerColor(color);
        resolutionVsTheta->SetMarkerSize(2);
        resolutionVsTheta->SetMarkerStyle(marker);
        resolutionVsTheta->SetLineColor(color);
        resolutionVsTheta->SetLineWidth(2);
    }

    void Fill(const AnglePoint &anglePoint)
    {
        Double_t ratio = anglePoint.getHcalOnlyToScatteredHcal();
        if (ratio == 0)
            return;
        Double_t eta = getEta(anglePoint.theta);
        ratioHcalOnlyToScatteredHcalVsEta->SetPoint(ratioHcalOnlyToScatteredHcalVsEta->GetN(), eta, ratio);
        ratioHcalOnlyToScatteredHcalVsTheta->SetPoint(ratioHcalOnlyToScatteredHcalVsTheta->GetN(), anglePoint.theta, ratio);
    }

    void Fill(const Double_t theta, const Double_t sigma)
    {
        Double_t eta = getEta(theta);

        resolutionVsEta->SetPoint(resolutionVsEta->GetN(), eta, sigma);
        resolutionVsTheta->SetPoint(resolutionVsTheta->GetN(), theta, sigma);
    }

    void Draw(TString outPdf, TCanvas *can)
    {

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);

        TLine *line = new TLine();
        line->SetLineColor(kRed);
        line->SetLineStyle(2);

        can->cd();
        resolutionVsEta->Draw("APL");
        Double_t min = TMath::MinElement(resolutionVsEta->GetN(), resolutionVsEta->GetY());
        Double_t max = TMath::MaxElement(resolutionVsEta->GetN(), resolutionVsEta->GetY());

        line->DrawLine(minEcalEta, min, minEcalEta, max);
        line->DrawLine(maxEcalEta, min, maxEcalEta, max);

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));

        can->SaveAs(outPdf);
        resolutionVsTheta->Draw("APL");
        min = TMath::MinElement(resolutionVsTheta->GetN(), resolutionVsTheta->GetY());
        max = TMath::MaxElement(resolutionVsTheta->GetN(), resolutionVsTheta->GetY());
        line->DrawLine(minEcalTheta, min, minEcalTheta, max);
        line->DrawLine(maxEcalTheta, min, maxEcalTheta, max);

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        can->SaveAs(outPdf);

        line->Draw("same");
    }

    void DrawComparisonTruth(TString outPdf, TCanvas *can, TileHists &other)
    {

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);

        TLine *line = new TLine();
        line->SetLineColor(kRed);
        line->SetLineStyle(2);

        TLegend *leg = new TLegend(0.2, 0.15, 0.4, 0.35);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        leg->AddEntry(resolutionVsEta, "Reco", "lp");
        leg->AddEntry(other.resolutionVsEta, "Truth", "lp");

        leg->AddEntry(line, "ecal acceptance", "l");

        TString title;
        if (name.Contains("phi"))
            title = "#phi resolution";
        else
            title = "#theta resolution";

        can->cd();

        resolutionVsEta->Draw("APL");
        Double_t min = TMath::MinElement(resolutionVsEta->GetN(), resolutionVsEta->GetY());
        Double_t max = TMath::MaxElement(resolutionVsEta->GetN(), resolutionVsEta->GetY());

        line->DrawLine(minEcalEta, min, minEcalEta, max);
        line->DrawLine(maxEcalEta, min, maxEcalEta, max);

        other.resolutionVsEta->Draw("PL same");

        leg->Draw();

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        tex->DrawLatex(0.2, 0.6, title);

        can->SaveAs(outPdf);

        resolutionVsTheta->Draw("APL");
        min = TMath::MinElement(resolutionVsTheta->GetN(), resolutionVsTheta->GetY());
        max = TMath::MaxElement(resolutionVsTheta->GetN(), resolutionVsTheta->GetY());
        line->DrawLine(minEcalTheta, min, minEcalTheta, max);
        line->DrawLine(maxEcalTheta, min, maxEcalTheta, max);

        other.resolutionVsTheta->Draw("PL");

        leg->Draw();

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        tex->DrawLatex(0.2, 0.6, title);
        can->SaveAs(outPdf);
    }

    void DrawComparison(TString outPdf, TCanvas *can, TileHists &other, TileHists &other2)
    {

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);

        TLine *line = new TLine();
        line->SetLineColor(kRed);
        line->SetLineStyle(2);

        TLegend *leg = new TLegend(0.2, 0.15, 0.4, 0.35);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        leg->AddEntry(resolutionVsEta, "Hcal only", "lp");
        leg->AddEntry(other.resolutionVsEta, "Ecal only", "lp");
        leg->AddEntry(other2.resolutionVsEta, "Hcal+Ecal both", "lp");

        leg->AddEntry(line, "ecal acceptance", "l");

        TString title;
        if (name.Contains("phi"))
            title = "#phi resolution";
        else
            title = "#theta resolution";

        if (name.Contains("Truth"))
            title += " Truth";

        can->cd();
        resolutionVsEta->SetTitle("");
        resolutionVsEta->Draw("APL");
        Double_t min = TMath::MinElement(resolutionVsEta->GetN(), resolutionVsEta->GetY());
        Double_t max = TMath::MaxElement(resolutionVsEta->GetN(), resolutionVsEta->GetY());

        line->DrawLine(minEcalEta, min, minEcalEta, max);
        line->DrawLine(maxEcalEta, min, maxEcalEta, max);

        other.resolutionVsEta->Draw("PL same");
        other2.resolutionVsEta->Draw("PL same");

        leg->Draw();

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        tex->DrawLatex(0.2, 0.6, title);

        can->SaveAs(outPdf);
        resolutionVsTheta->SetTitle("");
        resolutionVsTheta->Draw("APL");
        min = TMath::MinElement(resolutionVsTheta->GetN(), resolutionVsTheta->GetY());
        max = TMath::MaxElement(resolutionVsTheta->GetN(), resolutionVsTheta->GetY());
        line->DrawLine(minEcalTheta, min, minEcalTheta, max);
        line->DrawLine(maxEcalTheta, min, maxEcalTheta, max);

        other.resolutionVsTheta->Draw("PL");
        other2.resolutionVsTheta->Draw("PL");

        leg->Draw();

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        tex->DrawLatex(0.2, 0.6, title);
        can->SaveAs(outPdf);

        line->Draw("same");
    }

    void DrawScatteringRatio(TString outPdf, TCanvas *can)
    {

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        TLine *line = new TLine();
        line->SetLineColor(kRed);
        line->SetLineStyle(2);

        can->cd();
        can->SetLogy();
        can->SetGridy();

        ratioHcalOnlyToScatteredHcalVsEta->Draw("APL");
        Double_t min = TMath::MinElement(ratioHcalOnlyToScatteredHcalVsEta->GetN(), ratioHcalOnlyToScatteredHcalVsEta->GetY());
        Double_t max = TMath::MaxElement(ratioHcalOnlyToScatteredHcalVsEta->GetN(), ratioHcalOnlyToScatteredHcalVsEta->GetY());

        line->DrawLine(minEcalEta, min, minEcalEta, max);
        line->DrawLine(maxEcalEta, min, maxEcalEta, max);
        can->SaveAs(outPdf);

        can->cd();
        ratioHcalOnlyToScatteredHcalVsTheta->Draw("APL");
        min = TMath::MinElement(ratioHcalOnlyToScatteredHcalVsTheta->GetN(), ratioHcalOnlyToScatteredHcalVsTheta->GetY());
        max = TMath::MaxElement(ratioHcalOnlyToScatteredHcalVsTheta->GetN(), ratioHcalOnlyToScatteredHcalVsTheta->GetY());

        line->DrawLine(minEcalTheta, min, minEcalTheta, max);
        line->DrawLine(maxEcalTheta, min, maxEcalTheta, max);

        tex->DrawLatex(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        can->SaveAs(outPdf);
        can->SetLogy(0);
    }

    void Write()
    {
    }

    TGraph *ratioHcalOnlyToScatteredHcalVsEta;
    TGraph *ratioHcalOnlyToScatteredHcalVsTheta;
    TString name;

    TGraph *resolutionVsEta;
    TGraph *resolutionVsTheta;
};

void plot2DResolutionScanAlongPhi(TString outFile = "scanAlongPhi")
{

    for (double iEta = 0; iEta <= nEtaBins; iEta++) // 3 degree scan
    {
        double eta = etaMin + iEta * etaStep;
        double theta = 180 * getTheta(eta) / TMath::Pi();
        thetaBins.push_back(theta);
    }

    for (double iEta = 0; iEta <= nEtaBins; iEta++) // 3 degree scan
    {
        double eta = etaMax - iEta * etaStep;
        etaBins.push_back(eta);
    }

    cout << "Theta bins : " << thetaBins << endl;
    cout << "Eta bins : " << etaBins << endl;

    for (double iPhi = 0; iPhi <= nPhiBins; iPhi++) // 3 degree scan
    {
        double phi = phiMin + iPhi * phiStep;
        phiBins.push_back(phi);
    }

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
    vector<TString> labels = {"hcal+ecal #theta", "hcal+ecal #phi",
                              "hcal only #theta", "hcal only #phi",
                              "ecal only #theta", "ecal only #phi",
                              "hcal+ecal #theta Truth", "hcal+ecal #phi Truth",
                              "hcal only #theta Truth", "hcal only #phi Truth",
                              "ecal only #theta Truth", "ecal only #phi Truth"};

    for (int i = 0; i < nBins; i++)
    {
        resolution[i] = TileHists(labels[i]);
        resolution[i].setStyle(2001 + (int)(i / 2), 20 + (int)(i / 2));
    }

    for (auto &angle : tileMap)
    {

        TString angleStr = "Phi" + angle.first + "Theta" + angle.second;

        Double_t phi = angle.first.Atof();
        Double_t theta = angle.second.Atof();
        if (phi != alongPhiAngle)
            continue;

        TH1D *hist = (TH1D *)inputFile->Get("hResolution" + angleStr);
        if (!hist)
        {
            cout << "Histogram not found" << ("hResolution" + angleStr) << endl;
            continue;
        }

        AnglePoint point(angle, inputFile);
        resolution[0].Fill(point);

        for (int i = 0; i < nBins; i++)
        {
            Float_t sigma = 180 * hist->GetBinContent(i + 1) / pi;
            if (sigma == 0. || sigma > 40)
                continue;
            resolution[i].Fill(theta, sigma);
        }
    }
    TCanvas *can = new TCanvas("can", "can", 1200, 1000);
    can->cd();
    TString outPdf = outFile + ".pdf";
    can->SaveAs(outPdf + "[");

    resolution[0].DrawScatteringRatio(outPdf, can);

    for (int i = 0; i < nBins; i++)
    {
        resolution[i].Draw(outPdf, can);
    }

    // truth vs reco
    for (int i = 0; i < nBins / 2; i++)
    {
        resolution[i].DrawComparisonTruth(outPdf, can, resolution[i + 6]);
    }

    // theta
    resolution[2].DrawComparison(outPdf, can, resolution[4], resolution[0]);
    // phi
    resolution[3].DrawComparison(outPdf, can, resolution[5], resolution[1]);
    // theta truth
    resolution[8].DrawComparison(outPdf, can, resolution[10], resolution[6]);
    // phi truth
    resolution[9].DrawComparison(outPdf, can, resolution[11], resolution[7]);

    can->SaveAs(outPdf + "]");

    // TFile *output = new TFile("sectorResolution.root", "RECREATE");
    // output->cd();
    // for (int i = 0; i < nBins; i++)
    // {
    //     resolution[i].Write();
    // }

    // output->Close();
}


//TO do:   delta x , delta y, delta r 

// event display ? 

