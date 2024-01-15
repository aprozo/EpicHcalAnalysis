
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLine.h"
#include "TF1.h"
#include "TLatex.h"
#include <iostream>

using namespace std;
const int nEtaBins = 12;
const int nOldLayers = 24;

const double oldEtaBin[nEtaBins + 1] = {
    2.0, 1.9008, 1.8065, 1.7168, 1.6317, 1.5507, 1.4738,
    1.4007, 1.3312, 1.2651, 1.2023, 1.1427, 1.086};

const double oldZPositions[nOldLayers] = {
    270.19, 271.695, 273.15, 274.555, 275.96, 277.365,
    282.363, 283.768, 285.173, 286.578, 287.983, 289.388,
    290.793, 292.198, 293.603, 295.008, 296.413, 297.818,
    299.223, 300.628, 302.033, 303.438, 304.843, 306.158};

// const double hCalLength = 329.6;
const double hcalStart = 383.95; // ~!!!!! need to subtract 3 mm
const double steelWidth = 4;
const double polystereneWidth = 0.4;
const double airGapWidth = 0.15;
const double hcalTotalLength = 10 * (steelWidth + airGapWidth / 2 + polystereneWidth + airGapWidth / 2) + (airGapWidth / 2 + polystereneWidth + airGapWidth / 2);
// const double hcalEnd = 433.00;
//  const double hcalStart = hcalEnd - 10 * (steelWidth + airGapWidth + polystereneWidth) - 1 * (airGapWidth + polystereneWidth);

const double newZstart = hcalStart + airGapWidth / 2 + polystereneWidth / 2;
const double zShift = newZstart - oldZPositions[0];
const double newLayerWidth = steelWidth + airGapWidth + polystereneWidth;
const int firstNewLayer = 0;
const int nNewLayers = 11;
// const int tileMap[10] = {0, 2, 5, 6, 8, 10, 13, 16, 18, 21}; //start at 330
// const int tileMap[10] = {0, 2, 4, 5, 6, 8, 10, 12, 15, 17};
int tileMap[nNewLayers];
// eta = - log(tan(theta/2.0));
double getEta(const double &z, const double &r) { return asinh(z / r); }
double getR(const double &z, const double &eta) { return z / (sinh(eta)); }
double getNewEta(const double &zOld, const double &etaOld, double zNew) { return getEta(zNew, getR(zOld, etaOld)); }
double getNewEta(const int &layer, const int &etaBin, double zNew) { return getNewEta(oldZPositions[layer], oldEtaBin[etaBin], zNew); }

int zToEtaEpicHcal(void)
{
    // find the best match for each layer
    tileMap[0] = firstNewLayer;
    for (int iLayer = 0; iLayer < nNewLayers; iLayer++)
    {
        double newZPosition = newZstart + iLayer * newLayerWidth;
        cout << "layer " << iLayer << "  newZPosition = " << newZPosition << endl;
        vector<double> chi2;

        for (int iOldLayer = 0; iOldLayer < nOldLayers; iOldLayer++)
        {
            double sumSquare = 0;
            for (int iEta = 0; iEta < nEtaBins + 1; iEta++)
            {
                double newEtaLow = getNewEta(firstNewLayer, iEta, newZstart);
                double newRLow = getR(newZPosition, newEtaLow);
                double newEtaUp = getNewEta(firstNewLayer, iEta + 1, newZstart);
                double newRUp = getR(newZPosition, newEtaUp);
                double newTileHeight = (newRUp - newRLow);

                double oldRLow = getR(oldZPositions[iOldLayer], oldEtaBin[iEta]);
                double oldRUp = getR(oldZPositions[iOldLayer], oldEtaBin[iEta + 1]);
                double oldTileHeight = (oldRUp - oldRLow);
                sumSquare += pow(newTileHeight - oldTileHeight, 2);
            }
            chi2.push_back(sumSquare);
        }
        tileMap[iLayer] = std::distance(chi2.begin(), min_element(chi2.begin(), chi2.end()));
        std::cout << "hcal layer " << iLayer << " is " << tileMap[iLayer] << " emc layer ,  sqrt(chi2) =" << sqrt(chi2.at(tileMap[iLayer])) << endl;
    }
    for (int iLayer = 0; iLayer < nNewLayers; iLayer++)
    {
        cout << tileMap[iLayer] << ", ";
    }
    cout << endl;

    TGraph *oldEtaZ = new TGraph();
    oldEtaZ->SetMarkerColor(kOrange + 1);
    oldEtaZ->SetMarkerSize(1.5);
    oldEtaZ->SetMarkerStyle(20);
    oldEtaZ->GetXaxis()->SetRangeUser(260, 360);
    oldEtaZ->GetXaxis()->SetTitle("Z, cm (shifted Hcal)");
    oldEtaZ->GetYaxis()->SetTitle("Eta");

    TGraph *newEtaZ = new TGraph();
    newEtaZ->SetMarkerColor(kBlue);
    newEtaZ->SetMarkerSize(1.5);
    newEtaZ->SetMarkerStyle(22);
    newEtaZ->GetXaxis()->SetTitle("Z, cm ");
    newEtaZ->GetYaxis()->SetTitle("Eta");

    TGraph *oldRZ = (TGraph *)oldEtaZ->Clone();
    oldRZ->GetYaxis()->SetTitle("R, cm");
    TGraph *newRZ = (TGraph *)newEtaZ->Clone();
    TGraph *tileHeightVsOldZ = (TGraph *)oldRZ->Clone();
    tileHeightVsOldZ->GetYaxis()->SetTitle("Tile Height, cm");
    tileHeightVsOldZ->GetXaxis()->SetTitle("Z, cm (shifted Hcal)");

    TGraph *tileHeightVsNewZ = (TGraph *)newEtaZ->Clone();
    TGraph *selectedTileHeightVsOldZ = (TGraph *)newEtaZ->Clone();
    TGraph *newEtaOfOldLayersVsNewZ = (TGraph *)newEtaZ->Clone();

    for (int iLayer = 0; iLayer < nOldLayers; iLayer++)
    {
        for (int iEta = 0; iEta < nEtaBins + 1; iEta++)
        {
            double newZPosition = newZstart + iLayer * newLayerWidth;

            double newEta = getNewEta(firstNewLayer, iEta, newZstart);
            double oldR = getR(oldZPositions[iLayer], oldEtaBin[iEta]);
            double newR = getR(newZPosition, newEta);

            oldEtaZ->SetPoint(oldEtaZ->GetN(), oldZPositions[iLayer], oldEtaBin[iEta]);
            oldRZ->SetPoint(oldRZ->GetN(), oldZPositions[iLayer], oldR);
            if (iLayer < nNewLayers)
            { // only 10 layers in new hcal

                newEtaZ->SetPoint(newEtaZ->GetN(), newZPosition - zShift, newEta);
                newRZ->SetPoint(newRZ->GetN(), oldZPositions[tileMap[iLayer]], newR);
            }

            if (iEta < nEtaBins) // tile heights needs 2 eta values - low and up
            {

                double oldRLow = getR(oldZPositions[iLayer], oldEtaBin[iEta]);
                double oldRUp = getR(oldZPositions[iLayer], oldEtaBin[iEta + 1]);
                double oldTileHeight = (oldRUp - oldRLow);

                double newEtaLow = getNewEta(firstNewLayer, iEta, newZstart);
                double newRLow = getR(newZPosition, newEtaLow);
                double newEtaUp = getNewEta(firstNewLayer, iEta + 1, newZstart);
                double newRUp = getR(newZPosition, newEtaUp);

                double newTileHeight = (newRUp - newRLow);

                tileHeightVsOldZ->SetPoint(tileHeightVsOldZ->GetN(), oldZPositions[iLayer], oldTileHeight);

                if (iLayer < nNewLayers)
                {
                    tileHeightVsNewZ->SetPoint(tileHeightVsNewZ->GetN(), newZPosition - zShift, newTileHeight);
                    selectedTileHeightVsOldZ->SetPoint(selectedTileHeightVsOldZ->GetN(), oldZPositions[tileMap[iLayer]], newTileHeight);
                    int newLayerNumber = tileMap[iLayer];
                    double newEtaOfOldLayer = getEta(newZPosition, getR(oldZPositions[newLayerNumber], oldEtaBin[iEta]));
                    newEtaOfOldLayersVsNewZ->SetPoint(newEtaOfOldLayersVsNewZ->GetN(), newZPosition, newEtaOfOldLayer);
                }
            }
        }
    }

    TH1D *etaRanges = new TH1D("etaRanges", "etaRanges; etaBin; #Delta #eta", nEtaBins, 0, nEtaBins);

    TH1D *sizeRanges = new TH1D("sizeRanges", "tile heigth; etaBin; #Delta R", nEtaBins, 0, nEtaBins);

    TGraph *tileSizeVsEta = new TGraph();
    sizeRanges->SetMarkerColor(kBlue);
    sizeRanges->SetMarkerSize(2);
    sizeRanges->SetMarkerStyle(20);
    etaRanges->SetMarkerColor(kBlue);
    etaRanges->SetMarkerSize(2);
    etaRanges->SetMarkerStyle(20);

    double tileSizes[20];

    for (int iEta = 0; iEta < nEtaBins; iEta++)
    {
        double newEta = getNewEta(firstNewLayer, iEta, newZstart);
        double newEta2 = getNewEta(firstNewLayer, iEta + 1, newZstart);
        double deltaEta = newEta2 - newEta;
        double newR = getR(newZstart, newEta);
        double newR2 = getR(newZstart, newEta2);
        double deltaR = newR2 - newR;
        tileSizes[iEta] = deltaR;
        tileSizeVsEta->SetPoint(tileSizeVsEta->GetN(), iEta, deltaEta);
        etaRanges->SetBinContent(iEta + 1, deltaEta);
        sizeRanges->SetBinContent(iEta + 1, deltaR);
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->cd();
    tileSizeVsEta->SetTitle("tile size vs eta; eta bin; Tile size , cm");
    tileSizeVsEta->SetMarkerColor(kBlue);
    tileSizeVsEta->SetMarkerSize(2);
    tileSizeVsEta->SetMarkerStyle(20);
    // tileSizeVsEta->GetXaxis()->SetLimits(0., 290.);
    // tileSizeVsEta->GetHistogram()->SetMaximum(20);
    // tileSizeVsEta->GetHistogram()->SetMinimum(0);

    // TF1 *f1 = new TF1("f1", "pol1", 12.4, 275);
    // tileSizeVsEta->Fit(f1, "REM");
    // f1->Draw("same");

    TF1 *f2 = new TF1("f2", "pol2", -10, 15);
    tileSizeVsEta->GetXaxis()->SetLimits(-10, 15);
    tileSizeVsEta->GetHistogram()->SetMaximum(-0.05);
    tileSizeVsEta->GetHistogram()->SetMinimum(-0.2);
    tileSizeVsEta->Fit(f2, "REM");
    tileSizeVsEta->Draw("AP");

    TGraph *newEtaLayers = new TGraph();

    newEtaLayers->SetMarkerColor(kBlue);
    newEtaLayers->SetMarkerSize(2);
    newEtaLayers->SetMarkerStyle(20);
    TGraph *innerLayersPlot = new TGraph();
    TGraph *outerLayersPlot = new TGraph();
    innerLayersPlot->SetMarkerColor(kAzure + 1);
    innerLayersPlot->SetMarkerSize(2);
    innerLayersPlot->SetMarkerStyle(21);
    outerLayersPlot->SetMarkerColor(kViolet);
    outerLayersPlot->SetMarkerSize(2);
    outerLayersPlot->SetMarkerStyle(22);

    const int nInnerLayers = 10;
    const int nOuterLayers = 3;
    const int nEtaLayers = 13;
    const int nTotalLayers = nInnerLayers + nOuterLayers + nEtaLayers;

    double newEmcLayers[nEtaLayers] = {0};
    double innerLayers[nInnerLayers] = {0};
    double outerLayers[nOuterLayers] = {0};

    for (int iEta = 0; iEta < nEtaLayers; iEta++)
    {
        double newEta = getNewEta(firstNewLayer, iEta, newZstart);
        newEmcLayers[iEta] = newEta;
    }

    innerLayers[nInnerLayers - 1] = newEmcLayers[0];
    for (int iEta = nInnerLayers - 2; iEta >= 0; iEta--)
    {
        innerLayers[iEta] = innerLayers[iEta + 1] - f2->Eval(iEta - nInnerLayers + 1);
    }

    outerLayers[0] = newEmcLayers[nEtaLayers - 1];
    for (int iEta = 1; iEta < nOuterLayers; iEta++)
    {
        outerLayers[iEta] = outerLayers[iEta - 1] + f2->Eval(iEta + nEtaLayers - 2);
    }

    // printing and drawing

    cout << "new emc eta [0 layer]: " << endl;
    for (int iEta = 0; iEta < nEtaLayers; iEta++)
    {
        cout << newEmcLayers[iEta] << ", ";
        if (iEta < nEtaLayers - 1)
            newEtaLayers->SetPoint(newEtaLayers->GetN(), iEta, newEmcLayers[iEta + 1] - newEmcLayers[iEta]);
    }
    cout << endl;

    cout << "innerLayers eta: " << endl;
    for (int iEta = 0; iEta < nInnerLayers; iEta++)
    {
        cout << innerLayers[iEta] << ", ";
        if (iEta < nInnerLayers - 1)
            innerLayersPlot->SetPoint(innerLayersPlot->GetN(), -nInnerLayers + iEta + 1, innerLayers[iEta + 1] - innerLayers[iEta]);
    }
    cout << endl;
    cout << "outerLayers eta: " << endl;
    for (int iEta = 0; iEta < nOuterLayers; iEta++)
    {
        cout << outerLayers[iEta] << ", ";
        if (iEta < nOuterLayers - 1)
            outerLayersPlot->SetPoint(outerLayersPlot->GetN(), nEtaLayers + iEta - 1, outerLayers[iEta + 1] - outerLayers[iEta]);
    }
    cout << endl;
    newEtaLayers->GetXaxis()->SetLimits(-7., 14.);
    newEtaLayers->GetHistogram()->SetMaximum(-0.04);
    newEtaLayers->GetHistogram()->SetMinimum(-0.15);
    newEtaLayers->SetTitle("; eta bin (row); #Delta#eta");
    newEtaLayers->Draw("AP");
    f2->Draw("same");
    innerLayersPlot->Draw("P same");
    outerLayersPlot->Draw("P same");
    TLine *lin = new TLine();
    lin->SetLineColor(kOrange + 1);
    lin->SetLineWidth(1);
    lin->SetLineStyle(9);
    lin->SetLineColorAlpha(kRed, 0.3);
    lin->DrawLine(-0.5, -0.15, -0.5, -0.04);
    lin->DrawLine(11.5, -0.15, 11.5, -0.04);

    TLegend *legEta = new TLegend(0.4, 0.17, 0.74, 0.45);

    legEta->AddEntry(newEtaLayers, "old STAR EMC", "p");
    legEta->AddEntry(f2, "pol2 fit", "l");
    legEta->AddEntry(innerLayersPlot, "new Inner rows", "p");
    legEta->AddEntry(outerLayersPlot, "new Outer rows", "p");
    legEta->Draw("same");
    c2->SaveAs("newEtaRowsCheck.pdf");

    TLegend *leg = new TLegend(0.1, 0.8, 0.48, 0.9);
    leg->AddEntry(oldEtaZ, "old STAR EMC", "p");
    leg->AddEntry(newEtaZ, "new  HCAL( HCAL layer 0 = EMC layer 0 at z'_{0})", "p");

    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
    c1->cd();
    oldEtaZ->Draw("AP");
    oldEtaZ->GetHistogram()->SetMaximum(2.7); // along
    // oldEtaZ->GetHistogram()->SetMinimum(0.9); //   Y
    newEtaZ->Draw("P same");
    leg->DrawClone();
    c1->SaveAs("etaZ.pdf");

    oldRZ->Draw("AP");
    newRZ->Draw("P same");
    leg->DrawClone();
    c1->SaveAs("radiiZ.pdf");

    // tileHeightVsOldZ->GetXaxis()->SetLimits(260., 450.);
    tileHeightVsOldZ->Draw("AP");
    tileHeightVsNewZ->Draw("P same");
    leg->DrawClone();
    c1->SaveAs("tilesVsZ.pdf");

    tileHeightVsOldZ->Draw("AP");
    selectedTileHeightVsOldZ->Draw("P same");
    leg->DrawClone();
    c1->SaveAs("selectedTilesVsZ.pdf");

    TCanvas *can2 = new TCanvas("can2", "can2", 800, 800);
    can2->cd();
    newEtaOfOldLayersVsNewZ->GetYaxis()->SetTitle("recalculated #eta for Star Emc tiles");
    newEtaOfOldLayersVsNewZ->GetXaxis()->SetTitle("new Z for Star Emc, cm");
    newEtaOfOldLayersVsNewZ->Draw("AP");
    newEtaOfOldLayersVsNewZ->GetHistogram()->SetMaximum(2.7); // along
    // leg->DrawClone();
    TLine *l = new TLine();
    l->SetLineColor(newEtaOfOldLayersVsNewZ->GetMarkerColor());
    l->SetLineWidth(1);
    l->SetLineColorAlpha(newEtaOfOldLayersVsNewZ->GetMarkerColor(), 0.3);
    l->SetLineStyle(10);

    TLatex *tex = new TLatex();
    tex->SetTextSize(0.02);
    tex->SetTextAngle(90);
    tex->SetTextColor(tileHeightVsOldZ->GetMarkerColor());

    for (size_t iEta = 0; iEta < nEtaBins; iEta++)
    {
        double y = getEta(newZstart, getR(oldZPositions[0], oldEtaBin[iEta]));
        l->DrawLine(newZstart - 1.5, y, newZstart + 46, y);
    }

    for (size_t iLayer = 0; iLayer < nNewLayers; iLayer++)
    {
        tex->DrawLatex(newZstart + iLayer * newLayerWidth, 2.45, Form("emc layer %d", tileMap[iLayer]));
    }
    can2->SaveAs("checkForConstantEtaOfChosenEmcLayers.pdf");

    return 0;
}
