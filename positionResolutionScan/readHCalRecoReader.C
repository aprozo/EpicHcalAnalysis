
#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"

Int_t getMaximumEnergyIndex(TTreeReaderArray<Float_t> &energyVector)
{
  Int_t maxIndex = -1;
  Float_t maxEnergy = 0;
  for (int i = 0; i < energyVector.GetSize(); ++i)
  {
    if (energyVector[i] > maxEnergy)
    {
      maxEnergy = energyVector[i];
      maxIndex = i;
    }
  }
  return maxIndex;
}

struct Cluster
{
  Cluster() : theta(0), phi(0), energy(0), x(0), y(0), phiResolution(0), thetaResolution(0){};
  Cluster addEcalWithSampleCoefficient(const Float_t &sampleCoefficient, Cluster &ecal)
  {
    Cluster result;
    result.energy = (sampleCoefficient * ecal.energy) + energy;
    Double_t HcalContribution = energy / result.energy;
    Double_t EcalContribution = 1 - HcalContribution;
    result.theta = (EcalContribution * ecal.theta) + (HcalContribution * theta);
    result.phi = (EcalContribution * ecal.phi) + (HcalContribution * phi);
    result.x = (EcalContribution * ecal.x) + (HcalContribution * x);
    result.y = (EcalContribution * ecal.y) + (HcalContribution * y);
    return result;
  }

  Cluster &operator-(const Cluster &mc_position)
  {
    thetaResolution = (theta - mc_position.theta);
    phiResolution = (phi - mc_position.phi);
    return *this;
  }

  Double_t theta;
  Double_t phi;
  Double_t energy;
  Double_t x;
  Double_t y;
  Double_t phiResolution;
  Double_t thetaResolution;
};
Double_t getMaximum(TH1D *h1, TH1D *h2, TH1D *h3)
{
  Double_t max = h1->GetMaximum();
  if (h2->GetMaximum() > max)
    max = h2->GetMaximum();
  if (h3->GetMaximum() > max)
    max = h3->GetMaximum();
  return max;
}
struct ClusterHists
{
  ClusterHists(TString _name)
  {
    hTheta = new TH1D("hTheta_" + _name, _name + " Cluster #theta_{cluster}; #theta_{cluster} [rad]; Entries", 500, 2, 4);
    hPhi = new TH1D("hPhi_" + _name, _name + " Cluster #phi_{cluster}; #phi_{cluster} [rad]; Entries", 500, -TMath::Pi(), TMath::Pi());
    hThetaResol = new TH1D("hThetaResolution_" + _name, _name + " Cluster #theta Resolution ; #Delta#theta; Entries", 500, 1, 1);
    hPhiResol = new TH1D("hPhiResolution_" + _name, _name + " Cluster #phi Resolution; #Delta#phi; Entries", 500, 1, 1);
    hEnergy = new TH1D("hEnergy_" + _name, _name + " Cluster energy; E [GeV]; Entries", 10000, 0, 10);
    hPos = new TH2D("hPosition_" + _name, _name + " Cluster position x,y; x [mm]; y [mm]; Entries", 300, -1500, 1500, 300, -1500, 1500);

    hPhiEnergy = new TH2D("hPhiEnergy_" + _name, _name + " Cluster #phi vs Energy; #phi; Energy", 100, -TMath::Pi(), TMath::Pi(), 100, 0, 10);
    name = _name;
  }

  void Fill(const Cluster &cluster)
  {
    hTheta->Fill(cluster.theta);
    hPhi->Fill(cluster.phi);
    hEnergy->Fill(cluster.energy);
    hPos->Fill(cluster.x, cluster.y);
    hThetaResol->Fill(cluster.thetaResolution);
    hPhiResol->Fill(cluster.phiResolution);
  }
  void Draw(TCanvas *can, TString outPdf)
  {
    can->cd();
    hTheta->Draw("hist");
    can->SaveAs(outPdf);
    hPhi->Draw("hist");
    can->SaveAs(outPdf);
    hEnergy->Draw("hist");
    can->SaveAs(outPdf);
    hPos->Draw("colz");
    can->SaveAs(outPdf);
  }

  void DrawTogether(TCanvas *can, TString outPdf, ClusterHists &ecalHists, ClusterHists &sumHists)
  {
    can->cd();
    sumHists.hPhiResol->SetTitle("ECal, HCal, Sum #phi Resolution; #Delta#phi; Entries");
    sumHists.hPhiResol->SetMarkerColor(kBlack);
    sumHists.hPhiResol->SetLineColor(kBlack);
    sumHists.hPhiResol->GetYaxis()->SetRangeUser(0, 1.2 * getMaximum(sumHists.hPhiResol, ecalHists.hPhiResol, hPhiResol));
    sumHists.hPhiResol->Draw("hist");
    ecalHists.hPhiResol->Draw("same");
    hPhiResol->SetMarkerColor(kRed);
    hPhiResol->SetLineColor(kRed);
    hPhiResol->Draw("same");

    TLegend leg(0.4, 0.7, 0.55, 0.9);
    leg.AddEntry(ecalHists.hPhiResol, "ECal", "l");
    leg.AddEntry(hPhiResol, "HCal", "l");
    leg.AddEntry(sumHists.hPhiResol, "Sum", "l");
    leg.Draw();

    can->SaveAs(outPdf);
    sumHists.hThetaResol->SetTitle("ECal, HCal, Sum #theta Resolution; #Delta#theta; Entries");
    sumHists.hThetaResol->SetMarkerColor(kBlack);
    sumHists.hThetaResol->SetLineColor(kBlack);
    sumHists.hThetaResol->GetYaxis()->SetRangeUser(0, 1.2 * getMaximum(sumHists.hThetaResol, ecalHists.hThetaResol, hThetaResol));
    sumHists.hThetaResol->Draw("hist");
    ecalHists.hThetaResol->Draw("same");
    hThetaResol->SetMarkerColor(kRed);
    hThetaResol->SetLineColor(kRed);
    hThetaResol->Draw("same");

    leg.Draw();
    // TLatex *tl = new TLatex();
    // tl->SetTextSize(0.05);
    // tl->DrawLatexNDC(0.5, 0.5, "No ECal in");
    can->SaveAs(outPdf);

    sumHists.hTheta->SetTitle("ECal, HCal, Sum #theta; #theta; Entries");
    sumHists.hTheta->SetMarkerColor(kBlack);
    sumHists.hTheta->SetLineColor(kBlack);
    sumHists.hTheta->GetYaxis()->SetRangeUser(0, 1.2 * getMaximum(sumHists.hTheta, ecalHists.hTheta, hTheta));
    sumHists.hTheta->Draw("hist");
    ecalHists.hTheta->Draw("same");
    hTheta->SetMarkerColor(kRed);
    hTheta->SetLineColor(kRed);
    hTheta->Draw("same");

    leg.Draw();
    can->SaveAs(outPdf);

    sumHists.hPhi->SetTitle("ECal, HCal, Sum  #phi; #phi; Entries");
    sumHists.hPhi->SetMarkerColor(kBlack);
    sumHists.hPhi->SetLineColor(kBlack);
    sumHists.hPhi->GetYaxis()->SetRangeUser(0, 1.2 * getMaximum(sumHists.hPhi, ecalHists.hPhi, hPhi));
    sumHists.hPhi->Draw("hist");
    ecalHists.hPhi->Draw("same");
    hPhi->SetMarkerColor(kRed);
    hPhi->SetLineColor(kRed);
    hPhi->Draw("same");

    leg.Draw();
    can->SaveAs(outPdf);
  }

  Double_t calculateFWHM(TString varname)
  {
    TH1D *hist;
    TF1 *fitfunc = new TF1("fitfunc", "gaus", -4, 4);
    fitfunc->SetParNames("Area", "Mean", "Sigma");
    fitfunc->SetLineColor(kViolet);

    if (varname.Contains("theta"))
    {
      hist = hThetaResol;
      // fitfunc->SetRange(-0.02, 0.02);
      // fitfunc->SetParLimits(1, -0.1, 0.1);
      // fitfunc->SetParLimits(2, 0, 0.1);
    }
    else if (varname.Contains("phi"))
    {
      hist = hPhiResol;

      // fitfunc->SetRange(-0.5, 0.5);
      // fitfunc->SetParLimits(1, -0.5, 0.5);
      // fitfunc->SetParLimits(2, 0, 1);
    }
    else
    {
      cout << "Wrong variable name" << endl;
      return -1;
    }

    // Double_t x_min = hist->GetMean() - hist->GetStdDev();
    // Double_t x_max = hist->GetMean() + hist->GetStdDev();
    // hist->GetXaxis()->SetRangeUser(x_min, x_max);

    fitfunc->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetStdDev());
    // fitfunc->SetParLimits(0, 0, 1.2 * hist->GetMaximum());

    hist->GetYaxis()->SetTitleOffset(1.25);
    hist->SetMarkerStyle(43);
    hist->SetMarkerSize(1.5);
    hist->SetLineColor(kBlue);
    hist->SetMarkerColor(kBlack);
    hist->Fit("fitfunc", "MR", "");
    hist->Draw("p");
    Double_t sigma = fitfunc->GetParameter(2);

    Double_t fwhm = hist->GetStdDev(); // sigma * 2.355;
    // Double_t fwhm = sigma * 2.355;
    TLatex *tl = new TLatex();
    tl->SetTextSize(0.05);
    tl->DrawLatexNDC(0.12, 0.5, "#" + varname + " " + name);
    tl->DrawLatexNDC(0.12, 0.4, Form("Mean ~ %.2f", hist->GetMean()));
    tl->DrawLatexNDC(0.12, 0.3, Form("Sigma ~ %.2f", fwhm));

    return fwhm;
  }

  void Write(TDirectory *output)
  {
    cout << "Writing histograms" << name << endl;
    output->cd();
    TDirectory *dir = output->mkdir(name);
    dir->cd();
    hTheta->Write();
    hPhi->Write();
    hThetaResol->Write();
    hPhiResol->Write();
    hEnergy->Write();
    hPos->Write();
  }

  TString name;
  TH1D *hTheta;
  TH1D *hPhi;
  TH1D *hThetaResol;
  TH1D *hPhiResol;
  TH1D *hEnergy;
  TH2D *hPos;
  TH2D *hPhiEnergy;
};

// void readHCalRecoReader(TString inFileName = "eicrecon_neutron_5GeV.edm4eic.root", TString outFileName = "test.root")
void readHCalRecoReader(TString inFileName = "../output_eicrecon.edm4eic.root", TString outFileName = "test.root")
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

  //=======Reading the root file DD4HEP===========
  TFile *file = new TFile(inFileName);  // Tree with tracks and hits
                                        // Create the tree reader and its data containers
  TTreeReader myReader("events", file); // name of tree and file

  TTreeReaderArray<Float_t> charge(myReader, "MCParticles.charge");
  TTreeReaderArray<Double_t> vx_mc(myReader, "MCParticles.vertex.x");
  TTreeReaderArray<Double_t> vy_mc(myReader, "MCParticles.vertex.y");
  TTreeReaderArray<Double_t> vz_mc(myReader, "MCParticles.vertex.z");
  TTreeReaderArray<Float_t> px_mc(myReader, "MCParticles.momentum.x");
  TTreeReaderArray<Float_t> py_mc(myReader, "MCParticles.momentum.y");
  TTreeReaderArray<Float_t> pz_mc(myReader, "MCParticles.momentum.z");
  TTreeReaderArray<Int_t> status(myReader, "MCParticles.generatorStatus");
  TTreeReaderArray<Int_t> pdg(myReader, "MCParticles.PDG");

  TTreeReaderArray<Float_t> hcal_truth_E(myReader, "HcalEndcapNTruthClusters.energy");
  TTreeReaderArray<Float_t> hcal_truth_theta(myReader, "HcalEndcapNTruthClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> hcal_truth_phi(myReader, "HcalEndcapNTruthClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> hcal_truth_x(myReader, "HcalEndcapNTruthClusters.position.x");
  TTreeReaderArray<Float_t> hcal_truth_y(myReader, "HcalEndcapNTruthClusters.position.y");

  TTreeReaderArray<Float_t> ecal_truth_E(myReader, "EcalEndcapNTruthClusters.energy");
  TTreeReaderArray<Float_t> ecal_truth_theta(myReader, "EcalEndcapNTruthClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> ecal_truth_phi(myReader, "EcalEndcapNTruthClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> ecal_truth_x(myReader, "EcalEndcapNTruthClusters.position.x");
  TTreeReaderArray<Float_t> ecal_truth_y(myReader, "EcalEndcapNTruthClusters.position.y");

  TTreeReaderArray<Float_t> hcal_reco_E(myReader, "HcalEndcapNClusters.energy");
  TTreeReaderArray<Float_t> hcal_reco_theta(myReader, "HcalEndcapNClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> hcal_reco_phi(myReader, "HcalEndcapNClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> hcal_reco_x(myReader, "HcalEndcapNClusters.position.x");
  TTreeReaderArray<Float_t> hcal_reco_y(myReader, "HcalEndcapNClusters.position.y");

  TTreeReaderArray<Float_t> ecal_reco_E(myReader, "EcalEndcapNClusters.energy");
  TTreeReaderArray<Float_t> ecal_reco_theta(myReader, "EcalEndcapNClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> ecal_reco_phi(myReader, "EcalEndcapNClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> ecal_reco_x(myReader, "EcalEndcapNClusters.position.x");
  TTreeReaderArray<Float_t> ecal_reco_y(myReader, "EcalEndcapNClusters.position.y");

  TTreeReaderArray<Float_t> ebarell_truth_E(myReader, "EcalBarrelTruthClusters.energy");
  TTreeReaderArray<Float_t> ebarell_truth_theta(myReader, "EcalBarrelTruthClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> ebarell_truth_phi(myReader, "EcalBarrelTruthClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> ebarell_truth_x(myReader, "EcalBarrelTruthClusters.position.x");
  TTreeReaderArray<Float_t> ebarell_truth_y(myReader, "EcalBarrelTruthClusters.position.y");

  TTreeReaderArray<Float_t> ebarell_reco_E(myReader, "EcalBarrelClusters.energy");
  TTreeReaderArray<Float_t> ebarell_reco_theta(myReader, "EcalBarrelClusters.intrinsicTheta");
  TTreeReaderArray<Float_t> ebarell_reco_phi(myReader, "EcalBarrelClusters.intrinsicPhi");
  TTreeReaderArray<Float_t> ebarell_reco_x(myReader, "EcalBarrelClusters.position.x");
  TTreeReaderArray<Float_t> ebarell_reco_y(myReader, "EcalBarrelClusters.position.y");

  TCanvas *can = new TCanvas("can", "can", 1200, 1000);
  can->SetMargin(0.09, 0.1, 0.1, 0.06);

  TFile *output = new TFile(outFileName, "recreate");
  TH2D *hPtEtaMc = new TH2D("hPtEtaMc", "hPtEtaMc;p_{gen} (GeV/c);#eta_{gen} ", 1500, 0., 10.0, 200, -5, 5);

  ClusterHists hcalRecoHist("HCal_Reco");
  ClusterHists ecalRecoHist("ECal_Reco");
  ClusterHists ebarellRecHist("EBarell_Reco");
  ClusterHists sumRecoHist("Sum_Reco");

  ClusterHists hcalTruthHist("HCal_Truth");
  ClusterHists ecalTruthHist("ECal_Truth");
  ClusterHists ebarellTruthHist("EBarell_Truth");
  ClusterHists sumTruthHist("Sum_Truth");

  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////

  Int_t nEvents = myReader.GetEntries() / 10;
  cout << "Total Events: " << nEvents << endl;

  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    myReader.SetEntry(iEvent);
    if (iEvent % 1000 == 0)
      cout << "Event No: " << iEvent << endl;
    // MC Particle
    Cluster mc;
    for (int iParticle = 0; iParticle < charge.GetSize(); ++iParticle)
    {
      if (status[iParticle] == 1)
      {
        Double_t p_mc = sqrt(px_mc[iParticle] * px_mc[iParticle] + py_mc[iParticle] * py_mc[iParticle] + pz_mc[iParticle] * pz_mc[iParticle]);
        Double_t pt_mc = sqrt(px_mc[iParticle] * px_mc[iParticle] + py_mc[iParticle] * py_mc[iParticle]);
        Double_t eta_mc = -1.0 * TMath::Log(TMath::Tan((TMath::ACos(pz_mc[iParticle] / p_mc)) / 2));
        Double_t theta_mc = TMath::ACos(pz_mc[iParticle] / p_mc);
        Double_t phi_mc = TMath::ATan(py_mc[iParticle] / px_mc[iParticle]);
        mc.theta = theta_mc;
        mc.phi = phi_mc;
        mc.energy = p_mc;
        mc.x = vx_mc[iParticle];
        mc.y = vy_mc[iParticle];
        hPtEtaMc->Fill(pt_mc, eta_mc);
      }
    }

    Int_t indexHClusterTruth = getMaximumEnergyIndex(hcal_truth_E);
    Int_t indexHClusterReco = getMaximumEnergyIndex(hcal_reco_E);
    Int_t indexEClusterReco = getMaximumEnergyIndex(ecal_reco_E);
    Int_t indexEClusterTruth = getMaximumEnergyIndex(ecal_truth_E);
    Int_t indexEBarellClusterReco = getMaximumEnergyIndex(ebarell_reco_E);
    Int_t indexEBarellClusterTruth = getMaximumEnergyIndex(ebarell_truth_E);

    // if (indexHClusterReco < 0 && indexEClusterReco < 0)
    //   continue;

    Cluster hcalReco;
    Cluster ecalReco;
    Cluster ebarellReco;

    Cluster hcalTruth;
    Cluster ecalTruth;
    Cluster ebarellTruth;

    if (indexHClusterReco >= 0) // only hcal
    {
      hcalReco.energy = hcal_reco_E[indexHClusterReco];
      hcalReco.theta = hcal_reco_theta[indexHClusterReco];
      hcalReco.phi = hcal_reco_phi[indexHClusterReco];
      hcalReco.x = hcal_reco_x[indexHClusterReco];
      hcalReco.y = hcal_reco_y[indexHClusterReco];
      hcalReco = hcalReco - mc;

      if (indexEClusterReco < 0 && indexEBarellClusterReco < 0)
        hcalRecoHist.Fill(hcalReco);
    }

    if (indexEBarellClusterReco >= 0)
    {
      ebarellReco.energy = ebarell_reco_E[indexEBarellClusterReco];
      ebarellReco.theta = ebarell_reco_theta[indexEBarellClusterReco];
      ebarellReco.phi = ebarell_reco_phi[indexEBarellClusterReco];
      ebarellReco.x = ebarell_reco_x[indexEBarellClusterReco];
      ebarellReco.y = ebarell_reco_y[indexEBarellClusterReco];
      ebarellReco = ebarellReco - mc;

      if (indexHClusterReco < 0)
        ebarellRecHist.Fill(ebarellReco);
    }

    if (indexEClusterReco >= 0) // only hcal
    {
      ecalReco.energy = ecal_reco_E[indexEClusterReco];
      ecalReco.theta = ecal_reco_theta[indexEClusterReco];
      ecalReco.phi = ecal_reco_phi[indexEClusterReco];
      ecalReco.x = ecal_reco_x[indexEClusterReco];
      ecalReco.y = ecal_reco_y[indexEClusterReco];
      ecalReco = ecalReco - mc;

      if (indexHClusterReco < 0)
        ecalRecoHist.Fill(ecalReco);
    }

    if (indexHClusterReco >= 0 && indexEClusterReco >= 0)
    {
      Cluster sumReco = hcalReco;
      sumReco = hcalReco.addEcalWithSampleCoefficient(1.58, ecalReco);
      sumReco = sumReco - mc;
      sumRecoHist.Fill(sumReco);
    }

    if (indexHClusterTruth >= 0) // only hcal
    {
      hcalTruth.energy = hcal_truth_E[indexHClusterTruth];
      hcalTruth.theta = hcal_truth_theta[indexHClusterTruth];
      hcalTruth.phi = hcal_truth_phi[indexHClusterTruth];
      hcalTruth.x = hcal_truth_x[indexHClusterTruth];
      hcalTruth.y = hcal_truth_y[indexHClusterTruth];
      hcalTruth = hcalTruth - mc;

      if (indexEClusterTruth < 0 && indexEBarellClusterTruth < 0)
        hcalTruthHist.Fill(hcalTruth);
    }

    if (indexEBarellClusterTruth >= 0)
    {
      ebarellTruth.energy = ebarell_truth_E[indexEBarellClusterTruth];
      ebarellTruth.theta = ebarell_truth_theta[indexEBarellClusterTruth];
      ebarellTruth.phi = ebarell_truth_phi[indexEBarellClusterTruth];
      ebarellTruth.x = ebarell_truth_x[indexEBarellClusterTruth];
      ebarellTruth.y = ebarell_truth_y[indexEBarellClusterTruth];
      ebarellTruth = ebarellTruth - mc;

      if (indexHClusterTruth < 0)
        ebarellTruthHist.Fill(ebarellTruth);
    }

    if (indexEClusterTruth >= 0) // only ecal
    {
      ecalTruth.energy = ecal_truth_E[indexEClusterTruth];
      ecalTruth.theta = ecal_truth_theta[indexEClusterTruth];
      ecalTruth.phi = ecal_truth_phi[indexEClusterTruth];
      ecalTruth.x = ecal_truth_x[indexEClusterTruth];
      ecalTruth.y = ecal_truth_y[indexEClusterTruth];
      ecalTruth = ecalTruth - mc;

      if (indexHClusterTruth < 0)
        ecalTruthHist.Fill(ecalTruth);
    }

    if (indexHClusterTruth >= 0 && indexEClusterTruth >= 0)
    {
      Cluster sumTruth = hcalTruth;
      sumTruth = hcalTruth.addEcalWithSampleCoefficient(1.58, ecalTruth);
      sumTruth = sumTruth - mc;
      sumTruthHist.Fill(sumTruth);
    }

  } // Event For loop

  outFileName.ReplaceAll(".root", "");
  TString outPdf = outFileName + "QaCheck.pdf";

  can->cd();
  can->SaveAs(outPdf + "[");
  hPtEtaMc->Draw("colz");
  can->SaveAs(outPdf);

  Double_t fwhmThetaSum = sumRecoHist.calculateFWHM("theta");
  can->SaveAs(outPdf);
  Double_t fwhmPhiSum = sumRecoHist.calculateFWHM("phi");
  can->SaveAs(outPdf);
  Double_t fwhmThetaHcal = hcalRecoHist.calculateFWHM("theta");
  can->SaveAs(outPdf);
  Double_t fwhmPhiHcal = hcalRecoHist.calculateFWHM("phi");
  can->SaveAs(outPdf);

  Double_t fwhmThetaEcal = ecalRecoHist.calculateFWHM("theta");
  can->SaveAs(outPdf);
  Double_t fwhmPhiEcal = ecalRecoHist.calculateFWHM("phi");
  can->SaveAs(outPdf);

  Double_t fwhmThetaSumTruth = sumTruthHist.calculateFWHM("theta");
  can->SaveAs(outPdf);
  Double_t fwhmPhiSumTruth = sumTruthHist.calculateFWHM("phi");
  can->SaveAs(outPdf);
  Double_t fwhmThetaHcalTruth = hcalTruthHist.calculateFWHM("theta");
  can->SaveAs(outPdf);
  Double_t fwhmPhiHcalTruth = hcalTruthHist.calculateFWHM("phi");
  can->SaveAs(outPdf);

  Double_t fwhmThetaEcalTruth = ecalTruthHist.calculateFWHM("theta");
  can->SaveAs(outPdf);
  Double_t fwhmPhiEcalTruth = ecalTruthHist.calculateFWHM("phi");
  can->SaveAs(outPdf);

  hcalRecoHist.Draw(can, outPdf);
  ebarellRecHist.Draw(can, outPdf);
  ecalRecoHist.Draw(can, outPdf);
  sumRecoHist.Draw(can, outPdf);
  hcalRecoHist.DrawTogether(can, outPdf, ecalRecoHist, sumRecoHist);

  hcalTruthHist.Draw(can, outPdf);
  ebarellTruthHist.Draw(can, outPdf);
  ecalTruthHist.Draw(can, outPdf);
  sumTruthHist.Draw(can, outPdf);
  hcalTruthHist.DrawTogether(can, outPdf, ecalTruthHist, sumTruthHist);

  output->cd();

  string phiString = "";
  string thetaString = "";

  // string filename = "eicrecon_neutron_50000events_p5gev_phi78_theta155.872.edm4eic.root"
  std::string filename = inFileName.Data();
  std::regex phiRegex("phi(\\d+\\.?\\d*)");
  std::regex thetaRegex("theta(\\d+\\.?\\d*)");
  std::smatch phiMatch, thetaMatch;

  if (std::regex_search(filename, phiMatch, phiRegex) && std::regex_search(filename, thetaMatch, thetaRegex))
  {
    phiString = phiMatch.str(1);
    thetaString = thetaMatch.str(1);
    cout << "Phi =" << phiString << endl;
    cout << "Theta =" << thetaString << endl;
  }

  TH1D *hResolution = new TH1D(Form("hResolutionPhi%sTheta%s", phiString.data(), thetaString.data()), "Sigma values", 12, 0, 12);
  hResolution->GetXaxis()->SetBinLabel(1, "hcal+ecal #theta");
  hResolution->SetBinContent(1, fwhmThetaSum);
  hResolution->GetXaxis()->SetBinLabel(2, "hcal+ecal #phi");
  hResolution->SetBinContent(2, fwhmPhiSum);
  hResolution->GetXaxis()->SetBinLabel(3, "hcal #theta");
  hResolution->SetBinContent(3, fwhmThetaHcal);
  hResolution->GetXaxis()->SetBinLabel(4, "hcal #phi");
  hResolution->SetBinContent(4, fwhmPhiHcal);
  hResolution->GetXaxis()->SetBinLabel(5, "ecal #theta");
  hResolution->SetBinContent(5, fwhmThetaEcal);
  hResolution->GetXaxis()->SetBinLabel(6, "ecal #phi");
  hResolution->SetBinContent(6, fwhmPhiEcal);

  hResolution->GetXaxis()->SetBinLabel(7, "hcal+ecal #theta Truth");
  hResolution->SetBinContent(7, fwhmThetaSumTruth);
  hResolution->GetXaxis()->SetBinLabel(8, "hcal+ecal #phi Truth");
  hResolution->SetBinContent(8, fwhmPhiSumTruth);
  hResolution->GetXaxis()->SetBinLabel(9, "hcal #theta Truth");
  hResolution->SetBinContent(9, fwhmThetaHcalTruth);
  hResolution->GetXaxis()->SetBinLabel(10, "hcal #phi Truth");
  hResolution->SetBinContent(10, fwhmPhiHcalTruth);
  hResolution->GetXaxis()->SetBinLabel(11, "ecal #theta Truth");
  hResolution->SetBinContent(11, fwhmThetaEcalTruth);
  hResolution->GetXaxis()->SetBinLabel(12, "ecal #phi Truth");
  hResolution->SetBinContent(12, fwhmPhiEcalTruth);

  can->cd();
  hResolution->Draw("hist");
  can->SaveAs(outPdf);
  can->SaveAs(outPdf + "]");

  output->cd();
  hResolution->Write();
  TDirectory *dir = output->mkdir(Form("Phi%sTheta%s", phiString.data(), thetaString.data()));
  dir->cd();

  hcalRecoHist.Write(dir);
  ecalRecoHist.Write(dir);
  sumRecoHist.Write(dir);
  hcalTruthHist.Write(dir);
  ecalTruthHist.Write(dir);
  sumTruthHist.Write(dir);

  output->Save();
  output->Close();
}
