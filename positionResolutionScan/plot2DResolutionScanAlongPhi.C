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
#include "TGaxis.h"

const double alongPhiAngle = 48;

const double globalZ = -395; // first HCAL layer
// const double globalZ = -395 + 55 / 2; // middle HCAL layer
// const double globalZ = -395 + 55 / 2; // 1 interaction length HCAL layer

const double rmin = 14;  // first HCAL layer
const double rmax = 270; // first HCAL layer
const double pi = TMath::Pi();
double getTheta(const double &eta) { return 2 * atan(exp(-eta)); }
double getEta(const double &theta) { return -log(tan((pi * theta / 180) / 2)); }
double getEta(const double &z, const double &r) { return asinh(z / r); }
// double getTheta(const double &x, const double &y) { return atan2(sqrt(x * x + y * y), globalZ); }
// double getPhi(const double &x, const double &y) { return atan2(y, x); }
double getX(const double &theta, const double &phi) { return globalZ * tan(pi * theta / 180) * cos(pi * phi / 180); }
double getY(const double &theta, const double &phi) { return globalZ * tan(pi * theta / 180) * sin(pi * phi / 180); }

const int nPhiBins = 30;
const int nEtaBins = 20;

const double maxEcalTheta = 177.758;
const double minEcalTheta = 160.81035;

const double maxEcalEta = getEta(maxEcalTheta);
const double minEcalEta = getEta(minEcalTheta);

const double minEcalR = globalZ * tan(maxEcalTheta * pi / 180);
const double maxEcalR = globalZ * tan(minEcalTheta * pi / 180);

const double phiMin = 0;
const double phiMax = 90;
const double phiStep = (phiMax - phiMin) / (double)nPhiBins;

const double etaMin = getEta(globalZ, rmax);
const double etaMax = getEta(globalZ, rmin);
const double etaStep = (etaMax - etaMin) / (double)nEtaBins;

vector<double> thetaBins;
vector<double> phiBins;
vector<double> etaBins;

Double_t getMin(TGraph *g1, TGraph *g2)
{
    Double_t min = TMath::MinElement(g1->GetN(), g1->GetY());
    min = TMath::Min(min, TMath::MinElement(g2->GetN(), g2->GetY()));
    return min;
}

Double_t getMax(TGraph *g1, TGraph *g2)
{
    Double_t max = TMath::MaxElement(g1->GetN(), g1->GetY());
    max = TMath::Max(max, TMath::MaxElement(g2->GetN(), g2->GetY()));
    return max;
}

Double_t getMin(TGraph *g1, TGraph *g2, TGraph *g3)
{
    Double_t min = TMath::MinElement(g1->GetN(), g1->GetY());
    min = TMath::Min(min, TMath::MinElement(g2->GetN(), g2->GetY()));
    min = TMath::Min(min, TMath::MinElement(g3->GetN(), g3->GetY()));
    return min;
}

Double_t getMax(TGraph *g1, TGraph *g2, TGraph *g3)
{
    Double_t max = TMath::MaxElement(g1->GetN(), g1->GetY());
    max = TMath::Max(max, TMath::MaxElement(g2->GetN(), g2->GetY()));
    max = TMath::Max(max, TMath::MaxElement(g3->GetN(), g3->GetY()));
    return max;
}

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

        hRxyResol = (TH1D *)inputFile->Get(folder + "/" + _name + "/hRxyResolution_" + _name);
        if (!hRxyResol)
            cout << "hRxyResol not found" << endl;

        name = _name;
    }
    TString name;
    TH1D *hTheta;
    TH1D *hPhi;
    TH1D *hThetaResol;
    TH1D *hPhiResol;
    TH1D *hRxyResol;
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
        ecalRecoHist = {"EMCal_Reco", inputFile, angleStr};
        //  ebarellRecHist={"EBarell_Reco", inputFile, angleStr};
        hcalAndEcalSumHist = {"HcalAndEcalSum_Reco", inputFile, angleStr};
        scatteredHcal = {"scattered_Hcal", inputFile, angleStr};

        hcalTruthHist = {"HCal_Truth", inputFile, angleStr};
        ecalTruthHist = {"EMCal_Truth", inputFile, angleStr};
        //  ebarellTruthHist={"EBarell_Truth", inputFile, angleStr};
        hcalAndEcalSumTruthHist = {"HcalAndEcalSum_Truth", inputFile, angleStr};
        scatteredHcalTruth = {"scattered_Hcal_Truth", inputFile, angleStr};

        //   vector<TString> labels = {"hcal+emcal #theta", "hcal+emcal #phi", "hcal+emcal R_{xy}",
        //                             "hcal only #theta", "hcal only #phi", "hcal only R_{xy}",
        //                             "emcal only #theta", "emcal only #phi", "emcal only R_{xy}",
        //                             "hcal+emcal #theta Truth", "hcal+emcal #phi Truth", "hcal+emcal R_{xy} Truth",
        //                             "hcal only #theta Truth", "hcal only #phi Truth", "hcal only R_{xy} Truth",
        //                             "emcal only #theta Truth", "emcal only #phi Truth", "emcal only R_{xy} Truth"};

        histMap.push_back(hcalAndEcalSumHist.hThetaResol);
        histMap.push_back(hcalAndEcalSumHist.hPhiResol);
        histMap.push_back(hcalAndEcalSumHist.hRxyResol);
        histMap.push_back(hcalRecoHist.hThetaResol);
        histMap.push_back(hcalRecoHist.hPhiResol);
        histMap.push_back(hcalRecoHist.hRxyResol);
        histMap.push_back(ecalRecoHist.hThetaResol);
        histMap.push_back(ecalRecoHist.hPhiResol);
        histMap.push_back(ecalRecoHist.hRxyResol);
        histMap.push_back(hcalAndEcalSumTruthHist.hThetaResol);
        histMap.push_back(hcalAndEcalSumTruthHist.hPhiResol);
        histMap.push_back(hcalAndEcalSumTruthHist.hRxyResol);
        histMap.push_back(hcalTruthHist.hThetaResol);
        histMap.push_back(hcalTruthHist.hPhiResol);
        histMap.push_back(hcalTruthHist.hRxyResol);
        histMap.push_back(ecalTruthHist.hThetaResol);
        histMap.push_back(ecalTruthHist.hPhiResol);
        histMap.push_back(ecalTruthHist.hRxyResol);
    }
    Double_t getEntries(const Int_t &i)
    {
        return histMap[i]->GetEntries();
    }

    Double_t getHcalOnlyToScatteredHcal() const
    {
        // if (scatteredHcal.hTheta->GetEntries() == 0)
        //     return 0;
        // hcalRecoHist.hThetaResol->GetXaxis()->SetRangeUser(-1, 1);
        // hcalAndEcalSumHist.hThetaResol->GetXaxis()->SetRangeUser(-1, 1);
        Double_t hcalOnlyEntries = hcalRecoHist.hThetaResol->GetEntries();
        Double_t ecalOnlyEntries = ecalRecoHist.hThetaResol->GetEntries();
        Double_t hcalAndEcalEntries = hcalAndEcalSumHist.hThetaResol->GetEntries();

        Double_t ratio = hcalOnlyEntries / (hcalAndEcalEntries + ecalOnlyEntries);
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

    vector<TH1D *> histMap;
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
        else if (name.Contains("theta"))
            titleY = "#Delta#theta, deg ";
        else
            titleY = "#Delta R_{xy}, cm ";

        line = new TLine();
        tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);

        ratioHcalOnlyToScatteredHcalVsEta = new TGraph();
        ratioHcalOnlyToScatteredHcalVsEta->SetName(name + " ratioEta");
        ratioHcalOnlyToScatteredHcalVsEta->SetTitle("'Hcal Only' /('Hcal+Ecal both' + 'Ecal only'); #eta ; ratio");

        ratioHcalOnlyToScatteredHcalVsTheta = new TGraph();
        ratioHcalOnlyToScatteredHcalVsTheta->SetName(name + " ratioTheta");
        ratioHcalOnlyToScatteredHcalVsTheta->SetTitle("'Hcal Only' /('Hcal+Ecal both' + 'Ecal only'); #theta, deg; ratio");

        resolutionVsEta = new TGraph();
        resolutionVsEta->SetName(name + "resolutionVsEta");
        resolutionVsEta->SetTitle(name + "; #eta ;" + titleY);

        resolutionVsR = new TGraph();
        resolutionVsR->SetName(name + "resolutionVsR");
        resolutionVsR->SetTitle(name + "; r=#sqrt{x^{2}+y^{2}}, cm ;" + titleY);

        resolutionVsTheta = new TGraph();
        resolutionVsTheta->SetName(name + "resolutionVsTheta");
        resolutionVsTheta->SetTitle(name + "; #theta, deg;" + titleY);

        numberOfEntriesVsEta = new TGraph();
        numberOfEntriesVsEta->SetName(name + "numberOfEntriesVsEta");
        numberOfEntriesVsEta->SetTitle(name + "; #eta ; Number of Entries (out of 50k)");
    }

    void setStyle(TGraph *graph, const Int_t &color = 2001, const Int_t &marker = 20)
    {
        graph->SetMarkerColor(color);
        graph->SetMarkerSize(2);
        graph->SetMarkerStyle(marker);
        graph->SetLineColor(color);
        graph->SetLineWidth(2);
    }

    void setStyle(const Int_t &color = 2001, const Int_t &marker = 20)
    {
        setStyle(ratioHcalOnlyToScatteredHcalVsEta, color, marker);
        setStyle(ratioHcalOnlyToScatteredHcalVsTheta, color, marker);
        setStyle(resolutionVsEta, color, marker);
        setStyle(resolutionVsTheta, color, marker);
        setStyle(resolutionVsR, color, marker);
        setStyle(numberOfEntriesVsEta, color, marker);
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

    void Fill(const Double_t &theta, const Double_t &sigma)
    {
        Double_t eta = getEta(theta);

        resolutionVsEta->SetPoint(resolutionVsEta->GetN(), eta, sigma);
        resolutionVsTheta->SetPoint(resolutionVsTheta->GetN(), theta, sigma);
        Double_t x = getX(theta, alongPhiAngle);
        Double_t y = getY(theta, alongPhiAngle);
        Double_t r = TMath::Sqrt(x * x + y * y);
        resolutionVsR->SetPoint(resolutionVsR->GetN(), r, sigma);
    }

    void FillEntries(const Double_t theta, const Double_t entries)
    {
        Double_t eta = getEta(theta);
        numberOfEntriesVsEta->SetPoint(numberOfEntriesVsEta->GetN(), eta, entries);
    }

    void Draw(TGraph *graph, TString outPdf, TCanvas *can, Double_t minEcal, Double_t maxEcal, TGraph *graph2 = NULL, TGraph *graph3 = NULL)
    {
        can->cd();
        if (graph->GetN() == 0)
            return;
        graph->Draw("APL");
        Double_t min = 0;
        Double_t max = TMath::MaxElement(graph->GetN(), graph->GetY());

        TLegend *leg = new TLegend(0.2, 0.15, 0.4, 0.35);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        if (graph2 && !graph3)
        {
            leg->AddEntry(graph, "Reco", "lp");
            leg->AddEntry(graph2, "Truth", "lp");
            leg->AddEntry(line, "ecal acceptance", "l");
            if (graph2->GetN() != 0)
            {
                graph2->Draw("PL same");
                min = TMath::Min(min, TMath::MinElement(graph2->GetN(), graph2->GetY()));
                max = TMath::Max(max, TMath::MaxElement(graph2->GetN(), graph2->GetY()));
            }
        }

        if (graph3)
        {
            leg->AddEntry(graph, "Hcal only", "lp");
            leg->AddEntry(graph2, "Ecal only", "lp");
            leg->AddEntry(graph3, "Hcal+Ecal both", "lp");
            leg->AddEntry(line, "ecal acceptance", "l");

            if (graph2->GetN() != 0)
            {
                graph2->Draw("PL same");
                min = TMath::Min(min, TMath::MinElement(graph2->GetN(), graph2->GetY()));
                max = TMath::Max(max, TMath::MaxElement(graph2->GetN(), graph2->GetY()));
            }

            if (graph3->GetN() != 0)
            {
                graph3->Draw("PL same");
                min = TMath::Min(min, TMath::MinElement(graph3->GetN(), graph3->GetY()));
                max = TMath::Max(max, TMath::MaxElement(graph3->GetN(), graph3->GetY()));
            }
        }

        TString nameS = graph->GetName();
        if (nameS.Contains("numberOfEntriesVsEta"))
        {
            min = 0;
            max = 45000;
        }

        graph->GetYaxis()->SetRangeUser(min, max);
        line->DrawLine(minEcal, min, minEcal, max);
        line->DrawLine(maxEcal, min, maxEcal, max);
        TString title;
        if (nameS.Contains("phi"))
            title = "#phi resolution";
        else if (nameS.Contains("Truth"))
            title = "#theta resolution";
        else
            title = "R_{xy} resolution";
        if (nameS.Contains("Truth"))
            title += " Truth";

        tex->DrawLatexNDC(0.2, 0.7, Form("Scan Along #phi=%.0f #circ", alongPhiAngle));
        tex->DrawLatexNDC(0.2, 0.6, title);

        if (graph2 || graph3)
            leg->Draw();
        can->SaveAs(outPdf);
    }

    void Draw(TString outPdf, TCanvas *can)
    {
        Draw(ratioHcalOnlyToScatteredHcalVsEta, outPdf, can, minEcalEta, maxEcalEta);
        Draw(ratioHcalOnlyToScatteredHcalVsTheta, outPdf, can, minEcalTheta, maxEcalTheta);
        Draw(resolutionVsEta, outPdf, can, minEcalEta, maxEcalEta);
        Draw(numberOfEntriesVsEta, outPdf, can, minEcalEta, maxEcalEta);
    }

    void DrawComparisonTruth(TString outPdf, TCanvas *can, TileHists &other)
    {
        if (resolutionVsEta->GetN() == 0)
            return;
        Draw(resolutionVsEta, outPdf, can, minEcalEta, maxEcalEta, other.resolutionVsEta);
        Draw(resolutionVsTheta, outPdf, can, minEcalTheta, maxEcalTheta, other.resolutionVsTheta);
        Draw(resolutionVsR, outPdf, can, minEcalR, maxEcalR, other.resolutionVsR);
        Draw(numberOfEntriesVsEta, outPdf, can, minEcalEta, maxEcalEta, other.numberOfEntriesVsEta);
    }

    void
    DrawComparison(TString outPdf, TCanvas *can, TileHists &other, TileHists &other2)
    {
        if (resolutionVsEta->GetN() == 0)
            return;

        can->cd();
        gPad->SetRightMargin(0.2);

        Draw(resolutionVsEta, outPdf, can, minEcalEta, maxEcalEta, other.resolutionVsEta, other2.resolutionVsEta);
        Draw(resolutionVsTheta, outPdf, can, minEcalTheta, maxEcalTheta, other.resolutionVsTheta, other2.resolutionVsTheta);
        Draw(resolutionVsR, outPdf, can, minEcalR, maxEcalR, other.resolutionVsR, other2.resolutionVsR);
        Draw(numberOfEntriesVsEta, outPdf, can, minEcalEta, maxEcalEta, other.numberOfEntriesVsEta, other2.numberOfEntriesVsEta);
    }

    void DrawScatteringRatio(TString outPdf, TCanvas *can)
    {
        if (ratioHcalOnlyToScatteredHcalVsEta->GetN() == 0)
            return;

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

    TGraph *resolutionVsR;
    TGraph *resolutionVsEta;
    TGraph *resolutionVsTheta;

    TGraph *numberOfEntriesVsEta;
    TLatex *tex;
    TLine *line;
};

void plot2DResolutionScanAlongPhi(TString outFile = "scanAlongPhi")
{

    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(0.9, "y");
    TGaxis::SetMaxDigits(3);

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

    const Int_t nBins = 18;
    TileHists resolution[nBins];
    vector<TString> labels = {"hcal+emcal #theta", "hcal+emcal #phi", "hcal+emcal #Rxy",
                              "hcal only #theta", "hcal only #phi", "hcal only #Rxy",
                              "emcal only #theta", "emcal only #phi", "emcal only #Rxy",
                              "hcal+emcal #theta Truth", "hcal+emcal #phi Truth", "hcal+emcal #Rxy Truth",
                              "hcal only #theta Truth", "hcal only #phi Truth", "hcal only #Rxy Truth",
                              "emcal only #theta Truth", "emcal only #phi Truth", "emcal only #Rxy Truth"};

    for (int i = 0; i < nBins; i++)
    {
        resolution[i] = TileHists(labels[i]);
        resolution[i].setStyle(2001 + (int)(i / 3), 20 + (int)(i / 3));
    }

    for (auto &angle : tileMap)
    {

        TString angleStr = "Phi" + angle.first + "Theta" + angle.second;

        Double_t phi = angle.first.Atof();
        Double_t theta = angle.second.Atof();
        if (phi != alongPhiAngle)
            continue;
        cout << "Angle : " << angleStr << endl;

        TH1D *hist = (TH1D *)inputFile->Get("hResolution" + angleStr);
        if (!hist)
        {
            cout << "Histogram not found" << ("hResolution" + angleStr) << endl;
            continue;
        }

        AnglePoint point(angle, inputFile);
        resolution[0].Fill(point);
        for (int i = 0; i < nBins; i++)
        { // 0-8 reco, 9-17 truth
            if (point.getEntries(i) != 0)
                resolution[i].FillEntries(theta, point.getEntries(i));

            if (point.getEntries(i) < 2000)
                continue;
            Float_t sigma = 180 * hist->GetBinContent(i + 1) / pi;
            if ((i + 1) % 3 == 0) // not an angle but Rxy
                sigma = hist->GetBinContent(i + 1);
            resolution[i].Fill(theta, sigma);
        }
    }
    TCanvas *can = new TCanvas("can", "can", 1200, 1000);
    can->cd();
    TString outPdf = outFile + ".pdf";
    can->SaveAs(outPdf + "[");
    cout << "DrawScatteringRatio" << endl;
    resolution[0].DrawScatteringRatio(outPdf, can);

    cout << "Draw resolutions" << endl;
    for (int i = 0; i < nBins; i++)
    {
        resolution[i].Draw(outPdf, can);
    }

    cout << "DrawComparisonTruth" << endl;
    // truth vs reco
    for (int i = 0; i < nBins / 3; i++)
    {
        resolution[i].DrawComparisonTruth(outPdf, can, resolution[i + 9]);
    }

    cout << "DrawComparison Reco resolutions" << endl;
    for (int i = 0; i < 3; i++)
    {
        resolution[i + 3].DrawComparison(outPdf, can, resolution[i + 6], resolution[i]);
    }

    cout << "DrawComparison Truth resolutions" << endl;
    for (int i = 9; i < 12; i++)
    {
        cout << resolution[i].name << endl;
        cout << resolution[i + 3].name << endl;
        cout << resolution[i + 6].name << endl;
        resolution[i + 3].DrawComparison(outPdf, can, resolution[i + 6], resolution[i]);
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

// TO do:   delta x , delta y, delta r

// event display ?
