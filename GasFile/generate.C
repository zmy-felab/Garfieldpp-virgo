#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace std;
using namespace Garfield;

namespace
{
	void PrintUsage()
	{
		cerr << " Usage: " << endl;
		cerr << " ./generate [-g1 Gas1] [-f1 Fraction1] [-g2 Gas2] [-f2 Fraction2] [-g3 Gas3] [-f3 Fraction3] [-p Pressure] [-k Temperature]" << endl;
	}
} // namespace

int main(int argc, char *argv[])
{
	TApplication app("app", &argc, argv);

	if (argc < 3)
	{
		PrintUsage();
		return 1;
	}

	// Default parameters
	string gas1 = "", gas2 = "", gas3 = ""; //
	double fra1 = 0., fra2 = 0., fra3 = 0.; //
	double pre = 1.;						// atm
	double tem = 20.;						// C

	double fraction = 0.;

	string filename = "./result/";

	for (int i = 1; i < argc; i = i + 2)
	{
		if (string(argv[i]) == "-g1")
		{
			gas1 = string(argv[i + 1]);
			filename += gas1;
		}
		else if (string(argv[i]) == "-f1")
		{
			fra1 = atof(argv[i + 1]);
			fraction += fra1;
			char temp[10];
			sprintf(temp, "%.1f", fra1);
			filename += "_" + string(temp);
		}
		else if (string(argv[i]) == "-g2")
		{
			gas2 = string(argv[i + 1]);
			filename += "_" + gas2;
		}
		else if (string(argv[i]) == "-f2")
		{
			fra2 = atof(argv[i + 1]);
			fraction += fra2;
			char temp[10];
			sprintf(temp, "%.1f", fra2);
			filename += "_" + string(temp);
		}
		else if (string(argv[i]) == "-g3")
		{
			gas3 = string(argv[i + 1]);
			filename += "_" + gas3;
		}
		else if (string(argv[i]) == "-f3")
		{
			fra3 = atof(argv[i + 1]);
			fraction += fra3;
			char temp[10];
			sprintf(temp, "%.1f", fra3);
			filename += "_" + string(temp);
		}
		else if (string(argv[i]) == "-p")
		{
			pre = atof(argv[i + 1]);
			char temp[10];
			sprintf(temp, "%.1fatm", pre);
			filename += "_" + string(temp);
		}
		// "-t" can not work.
		else if (string(argv[i]) == "-k")
		{
			tem = atof(argv[i + 1]);
			char temp[10];
			sprintf(temp, "%.1fC", tem);
			filename += "_" + string(temp);
		}
		else
		{
			PrintUsage();
			return 1;
		}
	}
	if (fraction != 100)
	{
		printf("Please set the correct gas fractions.\n");
		return 2;
	}
	filename += ".gas";

	// Setup the gas.
	MediumMagboltz *gas = new MediumMagboltz();

	gas->SetComposition(gas1, fra1, gas2, fra2, gas3, fra3);
	gas->SetTemperature(273.15 + tem);
	gas->SetPressure(pre * AtmosphericPressure);
	// gas->SetMaxElectronEnergy(500);

	// Set the field range to be covered by the gas table.
	const int nFields = 50;
	const double emin = 100.;
	const double emax = 1000000.;

	// Flag to request logarithmic spacing.
	const bool useLog = true;
	gas->SetFieldGrid(emin, emax, nFields, useLog);

	// Set the angle between E- and B-field to 0 instead of 90 degrree.
	vector<double> efields, bfields, angles;
	gas->GetFieldGrid(efields, bfields, angles);
	angles = {0.};
	gas->SetFieldGrid(efields, bfields, angles);

	const int ncoll = 10;
	// Run Magboltz to generate the gas table.
	gas->GenerateGasTable(ncoll);
	// Save the table.
	gas->WriteGasFile(filename);

	// app.Run(kTRUE);
}