#include <iostream>
#include <fstream>
#include <set>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TGraph2D.h"
#include "TH1D.h"
#include "TCanvas.h"

const Double_t globalZ = 383.95; // first HCAL layer
const Double_t pi = 3.14159265358979323846;

Double_t getEta(const Double_t &theta) { return -log(tan((3.14 * theta / 180) / 2)); }
Double_t getTheta(const Double_t &x, const Double_t &y) { return atan2(sqrt(x * x + y * y), globalZ); }
Double_t getPhi(const Double_t &x, const Double_t &y) { return atan2(y, x); }

struct TileHists
{

    TileHists(TString name)
    {
        tilePhiTheta = new TGraph2D();
        tilePhiEta = new TGraph2D();
        tileXY = new TGraph2D();

        tilePhiTheta->SetTitle(name + " ;#phi, deg; #theta, deg ; fwhm");
        tilePhiEta->SetTitle(name + " ;#phi, deg; #eta; fwhm");
        tileXY->SetTitle(name + "  ;x,mm;y,mm; fwhm");
    }

    void Fill(const Double_t phi, const Double_t theta, const Double_t fwhm)
    {
        Double_t eta = getEta(theta);
        Double_t x = globalZ * tan((180 - theta) * pi / 180) * cos(phi * pi / 180);
        Double_t y = globalZ * tan((180 - theta) * pi / 180) * sin(phi * pi / 180);

        tilePhiEta->SetPoint(tilePhiEta->GetN(), phi, eta, fwhm);
        tilePhiTheta->SetPoint(tilePhiTheta->GetN(), phi, theta, fwhm);
        tileXY->SetPoint(tileXY->GetN(), x, y, fwhm);
    }

    void Draw(TString outPdf, TCanvas *can)
    {

        tilePhiTheta->Draw("TRI1");
        can->SaveAs(outPdf);

        tilePhiEta->Draw("TRI1");
        can->SaveAs(outPdf);

        tileXY->Draw("TRI1");
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
};

void plot2DResolution(TString outFile = "sectorResolution")
{
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

    TFile *inputFile = new TFile("output_hists.root", "READ");

    TileHists thetaHcal("Theta Hcal");
    TileHists phiHcal("Phi Hcal");

    // TileHists thetaEcal("Theta Ecal");
    // TileHists phiEcal("Phi Ecal");

    TileHists phiSum("Phi Hcal+Ecal");
    TileHists thetaSum("Theta Hcal+Ecal");

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

        Double_t fwhmThetaSum = hist->GetBinContent(1);
        Double_t fwhmPhiSum = hist->GetBinContent(2);
        Double_t fwhmThetaHcal = hist->GetBinContent(3);
        Double_t fwhmPhiHcal = hist->GetBinContent(4);

        if (fwhmThetaSum < 0.2 && fwhmThetaSum >= 0.0)
        {
            thetaSum.Fill(phi, theta, fwhmThetaSum);
        }

        if (fwhmPhiSum < 10 && fwhmPhiSum >= 0.0)
        {
            phiSum.Fill(phi, theta, fwhmPhiSum);
        }

        if (fwhmThetaHcal < 0.2 && fwhmThetaHcal >= 0.0)
        {
            thetaHcal.Fill(phi, theta, fwhmThetaHcal);
        }

        if (fwhmPhiHcal < 10 && fwhmPhiHcal >= 0.0)
        {
            phiHcal.Fill(phi, theta, fwhmPhiHcal);
        }
    }
    TCanvas *can = new TCanvas("can", "can", 2000, 1000);
    can->cd();
    TString outPdf = outFile + ".pdf";
    can->SaveAs(outPdf + "[");
    //can->SetLogz();
    thetaHcal.Draw(outPdf, can);
    phiHcal.Draw(outPdf, can);
    phiSum.Draw(outPdf, can);
    thetaSum.Draw(outPdf, can);

    can->SaveAs(outPdf + "]");

    TFile *output = new TFile("sectorResolution.root", "RECREATE");
    output->cd();
    thetaHcal.Write();
    phiHcal.Write();
    phiSum.Write();
    thetaSum.Write();

    output->Close();
}