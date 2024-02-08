#include <cmath>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"

using namespace std;

const double globalZ = -395; // first HCAL layer

const double rmin = 14;  // first HCAL layer
const double rmax = 270; // first HCAL layer

// double getR(const double &z, const double &theta) { return z / (sinh(eta)); }
// double getTheta(const double &r) { return atan(r / globalZ); }

double getTheta(const double &eta) { return 2 * atan(exp(-eta)); }

// // double getR(const double &z, const double &eta) { return z / (sinh(eta)); }
double getEta(const double &z, const double &r) { return asinh(z / r); }
// double getTheta(const double &x, const double &y) { return atan2(sqrt(x * x + y * y), globalZ); }
// double getPhi(const double &x, const double &y) { return atan2(y, x); }

struct PointHit
{
    PointHit(const Double_t &phi, const Double_t &theta) : phi(phi), theta(theta){};
    Double_t phi;
    Double_t theta;
    bool operator<(const PointHit &rhs) const
    {
        if (phi < rhs.phi)
            return true;
        if (rhs.phi < phi)
            return false;
        return theta < rhs.theta;
    }
};

ostream &operator<<(ostream &os, const set<PointHit> &set)
{
    for (auto &element : set)
    {
        os << element.phi << " " << element.theta << " " << endl;
    }
    return os;
}

ostream &operator<<(ostream &os, const set<Double_t> &set)
{
    for (auto &element : set)
    {
        os << element << endl;
    }
    return os;
}

void calculateAngles()
{
    TGraph *tilePhiTheta = new TGraph();
    TGraph *tileXY = new TGraph();
    TGraph *tilePhiEta = new TGraph();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetTitleOffset(1.0, "y");

    set<PointHit> tileValues;

    const int nPhiBins = 30;
    const int nEtaBins = 20;

    const double phiMin = 0;
    const double phiMax = 90;
    const double phiStep = (phiMax - phiMin) / (double)nPhiBins;

    const double etaMin = getEta(globalZ, rmax);
    const double etaMax = getEta(globalZ, rmin);
    const double etaStep = (etaMax - etaMin) / (double)nEtaBins;

    vector<double> thetaBins;
    for (double iEta = 0; iEta <= nEtaBins; iEta++) // 3 degree scan
    {
        double eta = etaMin + iEta * etaStep;
        double theta = 180 * getTheta(eta) / TMath::Pi();
        thetaBins.push_back(theta);
    }

    vector<double> phiBins;
    for (double iPhi = 0; iPhi <= nPhiBins; iPhi++) // 3 degree scan
    {
        double phi = phiMin + iPhi * phiStep;
        phiBins.push_back(phi);
    }

    TH2D *tileMap = new TH2D("tileMap", "tileMap; phi ;theta", phiBins.size() - 1, &phiBins[0], thetaBins.size() - 1, &thetaBins[0]);

    cout << "etaMin : " << etaMin << endl;
    cout << "etaMax : " << etaMax << endl;
    cout << "etaStep : " << etaStep << endl;
    Double_t counter = 1;
    for (double iEta = 0; iEta < nEtaBins; iEta++) // 3 degree scan
    {
        for (double iPhi = 0; iPhi < nPhiBins; iPhi++) // 3 degree scan
        {

            double eta = etaMin + iEta * etaStep;
            double phi = phiMin + iPhi * phiStep;
            if (phi == 0)
                phi = 0.001;

            double theta = 180 * getTheta(eta) / TMath::Pi();
            double x = globalZ * tan(getTheta(eta)) * cos(phi * TMath::Pi() / 180);
            double y = globalZ * tan(getTheta(eta)) * sin(phi * TMath::Pi() / 180);

            tileXY->SetPoint(tileXY->GetN(), x, y);
            tileMap->Fill(phi, theta, counter);
            tileValues.insert(PointHit(phi, theta));
            tilePhiTheta->SetPoint(tilePhiTheta->GetN(), phi, theta);
            tilePhiEta->SetPoint(tilePhiEta->GetN(), phi, eta);

            counter += 1;
        }
    }

    cout << tileValues;
    cout << "Number of points : " << tileValues.size() << endl;

    fstream fout;
    fout.open("tileMap.txt", ios::out);
    fout << tileValues;
    fout.close();

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
    c1->SaveAs("tilePhiTheta.pdf[");
    tilePhiTheta->SetMarkerStyle(20);

    TAxis *axis = tilePhiTheta->GetXaxis();

    tilePhiTheta->SetTitle("Sector Hit Points  ;#phi, deg;#theta,deg");

    gPad->SetMargin(0.15, 0.01, 0.1, 0.01);
    tilePhiTheta->Draw("AP");
    tilePhiTheta->SetMarkerStyle(20);
    tilePhiTheta->SetLineColor(kRed);

    c1->SaveAs("tilePhiTheta.pdf");

    tileMap->Draw("colz");
    c1->SaveAs("tilePhiTheta.pdf");

    tileXY->SetTitle("Sector Hit Points ;x, cm;y , cm");
    tileXY->SetMarkerStyle(20);
    tileXY->Draw("AP");
    c1->SaveAs("tilePhiTheta.pdf");

    tilePhiEta->SetTitle("Sector Hit Points ;#phi, deg;#eta");
    tilePhiEta->SetMarkerStyle(20);
    tilePhiEta->Draw("AP");
    c1->SaveAs("tilePhiTheta.pdf");

    c1->SaveAs("tilePhiTheta.pdf]");
}