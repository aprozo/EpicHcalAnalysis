#include <cmath>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

const double globalZ = 383.95; // first HCAL layer
double getR(const double &z, const double &eta) { return z / (sinh(eta)); }
double getEta(const double &z, const double &r) { return asinh(z / r); }
double getTheta(const double &x, const double &y) { return atan2(sqrt(x * x + y * y), globalZ); }
double getPhi(const double &x, const double &y) { return atan2(y, x); }

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

void calculateMap()
{
    const int nPointsPerPhi = 2; // actually it is 3 points per tile but not to overlap with the next tile and avoid double counting
    const int nPointsPerEta = 2; // actually it is 3points per tile but not to overlap with the next tile and avoid double counting

    set<PointHit> tileValues;
    set<Double_t> etaValues;
    TGraph *tilePhiTheta = new TGraph();
    TGraph *tileSizes = new TGraph();
    const double PI = TMath::Pi();
    const double phiStart = PI / 2;
    const int nHcalSectors = 12;
    const double dphiSector = 2 * PI / nHcalSectors; // 30  deg

    const int nOuterRows = 2;
    const int nStarEmcRows = 12;
    const int nInnerRowsGroup1 = 4; // innermost // 4 rows
    const int nInnerRowsGroup2 = 2; // central inner //2 rows
    const int nInnerRowsGroup3 = 3; // star emc size continuation //3 rows

    vector<double> outerHcalEtaBins = {1.38221, 1.3219, 1.26371};                                                                                          //[nOuterRows + 1]
    vector<double> starEmcEtaBins = {2.34289, 2.24167, 2.14506, 2.05272, 1.96465, 1.88032, 1.79973, 1.72258, 1.64865, 1.57776, 1.50981, 1.44473, 1.38221}; //[nStarEmcRows + 1]
    vector<double> innerHcalEtaBins3 = {2.67511, 2.55935, 2.44868, 2.34289};                                                                               //[nInnerRowsGroup3 + 1]
    vector<double> innerHcalEtaBins2 = {2.92264, 2.79614, 2.67511};                                                                                        //[nInnerRowsGroup2 + 1]
    vector<double> innerHcalEtaBins1 = {3.48733, 3.33696, 3.19285, 3.05481, 2.92264};                                                                      //[nInnerRowsGroup1 + 1]

    const double etaRangeWithHalves[47]=
    {1.26371, 1.29218, 1.3219, 1.35141, 1.38221, 1.4128, 1.44473, 1.47657, 1.50981, 1.54305, 1.57776, 1.61243, 1.64865, 1.68479, 1.72258, 1.76028, 1.79973, 1.83909, 1.88032,
    1.92148, 1.96465, 2.00761, 2.05272, 2.09773, 2.14506, 2.19211, 2.24167, 2.29092, 2.34289, 2.39432, 2.44868, 2.50242, 2.55935, 2.6155, 2.67511, 2.73375, 2.79614, 2.85735,
    2.92264, 2.98651, 3.05481, 3.12142, 3.19285, 3.26229, 3.33696, 3.4093, 3.48733};

    const double etaFullRange[nInnerRowsGroup1 + nInnerRowsGroup2 + nInnerRowsGroup3 + nStarEmcRows + nOuterRows + 1] =
        {3.48733, 3.33696, 3.19285, 3.05481, 2.92264, 2.79614, 2.67511, 2.55935, 2.44868, 2.34289,
         2.24167, 2.14506, 2.05272, 1.96465, 1.88032, 1.79973, 1.72258, 1.64865,
         1.57776, 1.50981, 1.44473, 1.38221, 1.3219, 1.26371};

    vector<vector<double>>
        allEtaBins = {innerHcalEtaBins1, innerHcalEtaBins2, innerHcalEtaBins3, starEmcEtaBins, outerHcalEtaBins};

    vector<pair<int, int>> sectorTileMap = {
        {nInnerRowsGroup1, 2}, // 2 tile per sector - 15 degree tiles
        {nInnerRowsGroup2, 3}, // 3 tile per sector - 10 degree tiles
        {nInnerRowsGroup3, 5}, // 5 tile per sector - 6 degree tiles
        {nStarEmcRows, 5},     // 5 tile per sector - 6 degree tiles
        {nOuterRows, 10}};     // 10 tile per sector - 3 degree tiles};

    Int_t rowCounter = -1;

    for (int iGroup = 0; iGroup < sectorTileMap.size(); iGroup++)
    {
        int nRows = sectorTileMap[iGroup].first;
        int nTiles = sectorTileMap[iGroup].second;
        vector<double> etaBins = allEtaBins[iGroup];
        double oldEmcGap = 0;
        if (iGroup < 3)
            oldEmcGap = 2.18;
        else if (iGroup == 3)
            oldEmcGap = 0;
        else
            oldEmcGap = -2.18;

        for (int iRow = 0; iRow < nRows; iRow++)
        {
            rowCounter++;
            double dPhiTile = dphiSector / nTiles;
            double dphi = dPhiTile / nPointsPerPhi;
            double rTop = getR(globalZ, etaBins[iRow + 1]) - oldEmcGap; // closest tile distance to the beam pipe

            // last row in 2 most inner groups - a  solution for avoiding overlaps
            if (iRow == (nRows - 1) && iGroup == 1)
            {
                rTop = rTop * cos(dPhiTile / 2) * cos(dPhiTile / 3);
            }
            if (iRow == (nRows - 1) && iGroup == 0)
            {
                rTop = rTop * 1.002 * cos(dPhiTile / 2) * cos(dPhiTile / 3);
            }

            double rBottom = getR(globalZ, etaBins[iRow]) - oldEmcGap;
            double deltaR = (rTop - rBottom) / nPointsPerEta;
            for (int iTile = 0; iTile < nTiles; iTile++)
            {
                for (int iEta = 0; iEta < nPointsPerEta + 1; iEta++)
                {
                    Int_t iGlobalEtaBin;

                    iGlobalEtaBin = rowCounter * nPointsPerEta + iEta;
                    cout << "iGlobalEtaBin " << iGlobalEtaBin << endl;

                    for (int iPhi = 0; iPhi < nPointsPerPhi + 1; iPhi++)
                    {
                        double phi = phiStart - iTile * dPhiTile - dphi * iPhi;
                        double r = rBottom + deltaR * iEta;
                        // r *= cos(dPhiTile / 2) / cos(dphi * iPhi - dPhiTile / 2); // using sine theorem to reduce radius - a trapezoid, not a tube segment along constant radius

                        double x = r * cos(phi);
                        double y = r * sin(phi);

                        Double_t eta = getEta(globalZ, r);

                        etaValues.insert(eta);
                        x += 1; // disksGap
                        if (iEta == 0 || iEta == nPointsPerEta || iPhi == 0 || iPhi == nPointsPerPhi)
                            tileSizes->SetPoint(tileSizes->GetN(), x, y);
                        phi = 180 * getPhi(x, y) / PI;
                        double theta = 180 * getTheta(x, y) / PI;
                        tilePhiTheta->SetPoint(tilePhiTheta->GetN(), x, y);

                        tileValues.insert({phi, theta});
                    }
                }
            }
        }
    }
    cout << tileValues;
    cout << "Number of points : " << tileValues.size() << endl;

    cout << "Number of eta bins : " << etaValues.size() << endl;
    cout << etaValues;

    fstream fout;
    fout.open("tileMap.txt", ios::out);
    fout << tileValues;
    fout.close();

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 1000);
    tilePhiTheta->SetMarkerStyle(20);
    TAxis *axis = tileSizes->GetXaxis();
    axis->SetLimits(0., 125.);                   // along X
    tileSizes->GetHistogram()->SetMaximum(250.); // along
    tileSizes->GetHistogram()->SetMinimum(0.);   //   Y
    tileSizes->SetTitle("Sector Hit Points ;x,mm;y,mm");
    tileSizes->GetYaxis()->SetTitleOffset(1.1);
    tileSizes->GetXaxis()->SetTitleOffset(0.6);
    gPad->SetMargin(0.15, 0.01, 0.1, 0.01);
    tileSizes->Draw("APL");
    tileSizes->SetMarkerStyle(20);
    tileSizes->SetLineColor(kRed);
    tilePhiTheta->Draw("P SAME");
    c1->SaveAs("tilePhiTheta.pdf");
}