#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/Plotting.hh"

using namespace std;
using namespace Garfield;

int main(int argc, char *argv[])
{
	TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();

	// Load a gas file which includes excitation and ionisation rates.
	MediumMagboltz *gas = new MediumMagboltz();
	gas->LoadGasFile(string(app.Argv(1)));
	gas->PrintGas();

	const std::string path = getenv("GARFIELD_INSTALL");
	gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
	gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ne+_Ne.txt");
	gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_He+_He.txt");

	ViewMedium *view = new ViewMedium();
	view->SetMedium(gas);
	// view->SetMagneticField(2);

	TCanvas *c1 = new TCanvas("c1", "", 1000, 1000);
	c1->Divide(2, 2);
	view->SetCanvas(dynamic_cast<TPad*>(c1->cd(1)));
	view->PlotElectronVelocity('e');
	view->SetCanvas(dynamic_cast<TPad*>(c1->cd(2)));
	view->PlotElectronDiffusion('e');
	view->SetCanvas(dynamic_cast<TPad*>(c1->cd(3)));
	view->PlotElectronTownsend('e');
	view->SetCanvas(dynamic_cast<TPad*>(c1->cd(4)));
	view->PlotElectronAttachment('e');

	TCanvas *c2 = new TCanvas("c2", "",1000, 500);
	c2->Divide(2, 1);
	// Plot the Townsend coefficient (without Penning transfer).
	view->SetCanvas(dynamic_cast<TPad*>(c2->cd(1)));
	view->PlotElectronTownsend('e');
	// Switch on Penning transfer for all excitation levels in the mixture.
	gas->EnablePenningTransfer(0.426, 0.);
	// Plot the Townsend coefficient with Penning transfer.
	view->SetCanvas(dynamic_cast<TPad*>(c2->cd(2)));
	view->PlotElectronTownsend('e');

	app.Run(kTRUE);
}
