/*
 * calculateLimits.C
 *
 *  Created on: 21 wrz 2022
 *      Author: Khaless
 */


#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>


//#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TVector3.h>

using namespace std;
using namespace ROOT;
using namespace TMath;

void printArray(double *arr, int size);

int calculateLimits()
{
	const int nEtaBins = 12;
	const int nEtaLimits = nEtaBins+1;

	double z = 270; // [cm]

	static const float defaultEtaBin[nEtaLimits] = {
		2.0,
		1.9008, 1.8065, 1.7168, 1.6317, 1.5507, 1.4738,
		1.4007, 1.3312, 1.2651, 1.2023, 1.1427, 1.086
	};

	double rBins[nEtaLimits];

	for (int i = 0; i < nEtaLimits; ++i) {

		double eta = defaultEtaBin[i];

		double r = z/(sinh(eta));

		rBins[i] = r;

	}

	printArray(rBins, nEtaLimits);


	return 1;
}
////////////////////////////

void printArray(double *arr, int size)
{


	for (int i = 0; i < size; ++i) {

		cout<<arr[i]<<", ";

	}

	cout<<endl;

}
