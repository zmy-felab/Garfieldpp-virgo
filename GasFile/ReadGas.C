#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <numeric>
#include <cmath>
#include <dirent.h>

#include "TFile.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TColor.h"
#include "TStyle.h"

using namespace std;

void ReadRecord(ifstream &gasfile, vector<double> &ve, vector<double> &dl, vector<double> &dt,
				vector<double> &alpha, vector<double> &alpha0, vector<double> &eta, int &nexc, int &nion);
void ReadHeader(ifstream &gasfile, vector<double> &mixture, vector<double> &efields, int &nexc, int &nion);
void ReadFooter(ifstream &gasfile, double &pgas, double &tgas);

void GetMixture(const vector<double> &mixture, vector<string> &gasnames, vector<double> &percentages);
string GetGasName(const int gasnumber);

bool GetPath(string path, vector<string> &name);

void ReadGas(const string dataDir)
{
	gStyle->SetTextFont(22);
	gStyle->SetLabelFont(22);
	gStyle->SetTitleFont(22);
	gStyle->SetLegendFont(22);
	gStyle->SetPadGridX(true);
	gStyle->SetPadGridY(true);
	gStyle->SetPadTickX(true);
	gStyle->SetPadTickY(true);
	gStyle->SetOptLogx();
	gStyle->SetFillStyle(0);

	TFile *f = new TFile((dataDir+".root").c_str(), "recreate");

	vector<string> filename; // gas file
	if(GetPath(dataDir, filename)) 
	{
		const int phy = 5;
		const int num = filename.size();
		TGraph *gr[phy][num];

		for(Int_t i = 0; i < num; i++)
		{			
			ifstream gasfile;
			gasfile.open(filename[i].c_str());
			if (!gasfile.is_open())
			{
				cerr << "Cannot open file " << filename[i] << endl;
				exit(-1);
			}
			cout << "\nReading " << filename[i] << endl;

			const int nMagboltzGases = 60;
			vector<double> mixture(nMagboltzGases, 0);
			vector<double> EF;

			int nexc; // excitations
			int nion; // ionisations

			// Check header
			ReadHeader(gasfile, mixture, EF, nexc, nion);

			const int nE = EF.size();
			cout << "LoadGasFile: " << nE << " electric field(s)\n";
			cout << "Found " << nexc << " excitations, " << nion << " ionisations.\n";

			// Get mixture gas information (optional)
			// vector<string> gasnames;
			// vector<double> percentages;
			// GetMixture(mixture, gasnames, percentages);

			// Physics parameters 
			vector<double> DV(nE, 0); // Drift velocity
			vector<double> TD(nE, 0); // Transverse diffusion
			vector<double> LD(nE, 0); // Longitudinal diffusion
			vector<double> TC(nE, 0); // Townsend coefficient
			vector<double> AC(nE, 0); // Attachment coefficient
			// Table of Townsend coefficients without Penning transfer
			vector<double> TC0(nE, 0);

			// Get data
			ReadRecord(gasfile, DV, LD, TD, TC, TC0, AC, nexc, nion);

			// Get pressure
			double pressure;
			double temperature;
			ReadFooter(gasfile, pressure, temperature);

			gasfile.close();

			for(int j = 0; j < phy; j++)
				gr[j][i] = new TGraph();

			const double sqrp = sqrt(pressure);
			const double logp = log(pressure);
			for (int j = 0; j < nE; ++j)
			{
				EF[j] *= pressure; //
				DV[j] *= 1.e-3;	   // from cm/us to cm/ns
				LD[j] /= sqrp;	   //
				TD[j] /= sqrp;	   //

				TC[j] = exp(TC[j] + logp);	 //
				TC0[j] = exp(TC0[j] + logp); //
				AC[j] = exp(AC[j] + logp);	 //
				// cout << EF[j] << setw(15) << DV[j] << setw(15) << TD[j] << setw(15) << LD[j] << setw(15) << TC[j] << setw(15) << AC[j] << endl;

				gr[0][i]->SetPoint(j, EF[j], DV[j]);
				gr[1][i]->SetPoint(j, EF[j], LD[j]);
				gr[2][i]->SetPoint(j, EF[j], TD[j]);
				gr[3][i]->SetPoint(j, EF[j], TC[j]);
				gr[4][i]->SetPoint(j, EF[j], AC[j]);
			}
		}

		char ytitle[phy][100] = {"Drift Velocity [cm/ns]", "Longitudinal Diffusion Coefficient [#kern[-0.1]{#sqrt{cm}}#kern[0.1]{]}", "Transverse Diffusion Coefficient [#kern[-0.1]{#sqrt{cm}}#kern[0.1]{]}", "Townsend coefficient #it{#alpha} [1/cm]", "Attachment coefficient #it{#eta} [1/cm]"};
		TMultiGraph *mg[phy];
		TLegend *lg[phy];
		TCanvas *c[phy];
		for(int i = 0; i < phy; i++)
		{
			mg[i] = new TMultiGraph();
			lg[i] = new TLegend();
			int color[12] = {1,2,3,4,kOrange+7,6,7,8,9,28,38,50};
			for(int j = 0; j < num; j++)
			{
				gr[i][j]->SetLineWidth(2);
				gr[i][j]->SetLineColor(color[j]);
				// gr[i][j]->SetLineStyle(j+1);
				gr[i][j]->SetMarkerColor(color[j]);
				gr[i][j]->SetMarkerStyle(20+j);

				mg[i]->Add(gr[i][j]);
				lg[i]->AddEntry(gr[i][j], filename[j].substr(filename[j].size() - 10, 6).c_str(), "lpf");
			}
			mg[i]->GetXaxis()->SetTitle("Electric Field [V/cm]");
			mg[i]->GetXaxis()->CenterTitle();
			mg[i]->GetXaxis()->SetTitleOffset(1.2);
			mg[i]->GetXaxis()->SetLimits(90, 1.1e6);
			mg[i]->GetYaxis()->SetTitle(ytitle[i]);
			mg[i]->GetYaxis()->CenterTitle();
			mg[i]->GetYaxis()->SetLabelFont(22);
			mg[i]->GetYaxis()->SetTitleFont(22);
			mg[i]->GetYaxis()->SetTitleOffset(1.5);

			c[i] = new TCanvas();
			mg[i]->Draw("ACP");
			lg[i]->Draw();
			c[i]->Write();
			c[i]->Update();
		}
	}
	f->Close();
}
void ReadHeader(ifstream &gasfile, vector<double> &mixture, vector<double> &efields, int &nexc, int &nion)
{
	while (1)
	{
		char line[256];
		gasfile.getline(line, 256);
		const bool quotes = (strstr(line, "\"") != NULL);
		// Head End
		if (strncmp(line, " The gas tables follow:", 8) == 0 || strncmp(line, "The gas tables follow:", 7) == 0)
			break;

		char *token = strtok(line, " :,%");
		while (token)
		{
			if (strcmp(token, "Identifier") == 0)
			{
				// Get the identification string.
				string identifier = "";
				token = strtok(NULL, "\n");
				if (token != NULL)
					identifier += token;
				cout << "Condition : " << identifier << "\n";
			}
			else if (strcmp(token, "Dimension") == 0)
			{
				token = strtok(NULL, " :,%\t");
				token = strtok(NULL, " :,%\t");

				const int nE = atoi(token);
				// Check the number of E points.
				if (nE <= 0)
					cerr << "Number of E fields out of range.\n";
				efields.resize(nE);

				token = strtok(NULL, " :,%\t");
				token = strtok(NULL, " :,%\t");

				// Excitation
				token = strtok(NULL, " :,%\t");
				nexc = atoi(token);

				// Ionization
				token = strtok(NULL, " :,%\t");
				nion = atoi(token);
			}
			else if (strcmp(token, "E") == 0)
			{
				token = strtok(NULL, " :,%");
				if (strncmp(token, "fields", 6) == 0)
				{
					const int nE = efields.size();
					// ELECTRIC FIELD
					for (int i = 0; i < nE; ++i)
						gasfile >> efields[i];
				}
			}
			else if (strcmp(token, "Mixture") == 0)
			{
				const unsigned int nMagboltzGases = mixture.size();
				// Mixture Gas
				for (unsigned int i = 0; i < nMagboltzGases; ++i)
					gasfile >> mixture[i];
			}
			token = strtok(NULL, " :,%");
		}
	}
}
void ReadRecord(ifstream &gasfile, vector<double> &ve, vector<double> &dl, vector<double> &dt,
				vector<double> &alpha, vector<double> &alpha0, vector<double> &eta, int &nexc, int &nion)
{
	double waste = 0.;
	const int nE = ve.size();

	for (int i = 0; i < nE; i++)
	{
		gasfile >> ve[i] >> waste >> waste >> waste >> waste >> waste;
		gasfile >> dl[i] >> waste >> dt[i] >> waste >> alpha[i] >> waste >> alpha0[i] >> eta[i] >> waste;

		// ion mobility, lorentz angle, dissociation, excitations, ionisations
		int res = 3 * 2 + 6 * 2 + nexc * 2 + nion * 2;
		for (int i = 0; i < res; ++i)
			gasfile >> waste;
	}
}
void ReadFooter(ifstream &gasfile, double &pgas, double &tgas)
{
	while (1)
	{
		char line[256];
		gasfile.getline(line, 256);
		if (gasfile.fail())
			break;
		char *token = strtok(line, " :,%=\t");
		while (token)
		{
			if (strcmp(token, "PGAS") == 0)
			{
				// Pressure [Torr]
				token = strtok(NULL, " :,%=\t");
				if (token != NULL)
					pgas = atof(token);
			}
			else if (strcmp(token, "TGAS") == 0)
			{
				// Temperature [K]
				token = strtok(NULL, " :,%=\t");
				if (token != NULL)
					tgas = atof(token);
				break;
			}
			token = strtok(NULL, " :,%=\t");
		}
	}
}
void GetMixture(const vector<double> &mixture, vector<string> &gasnames, vector<double> &percentages)
{
	gasnames.clear();
	percentages.clear();
	const double Small = 1.e-20;
	const unsigned int m_nMaxGases = 6;
	const unsigned int nMagboltzGases = mixture.size();
	for (unsigned int i = 0; i < nMagboltzGases; ++i)
	{
		if (mixture[i] < Small)
			continue;
		const string gasname = GetGasName(i + 1);
		if (gasname.empty())
		{
			cerr << "GetMixture: Unknown gas (gas number " << i + 1 << ").\n";
			exit(1);
		}
		gasnames.push_back(gasname);
		percentages.push_back(mixture[i]);
	}
	if (gasnames.size() > m_nMaxGases)
	{
		cerr << "GetMixture: Gas mixture has " << gasnames.size() << " components.\n"
			 << "    Number of gases is limited to " << m_nMaxGases << ".\n";
		exit(1);
	}
	else if (gasnames.empty())
	{
		cerr << "GetMixture: Gas mixture is not defined (zero components).\n";
		exit(1);
	}
	double sum = accumulate(percentages.begin(), percentages.end(), 0.);
	if (sum != 100.)
	{
		cout << "GetMixture: Renormalizing the percentages.\n";
		for (auto &percentage : percentages)
			percentage *= 100. / sum;
		exit(1);
	}
}
string GetGasName(const int gasnumber)
{
	switch (gasnumber)
	{
	case 1:
		return "CF4";
		break;
	case 2:
		return "Ar";
		break;
	case 3:
		return "He";
		break;
	case 4:
		return "He-3";
		break;
	case 5:
		return "Ne";
		break;
	case 6:
		return "Kr";
		break;
	case 7:
		return "Xe";
		break;
	case 8:
		return "CH4";
		break;
	case 9:
		return "C2H6";
		break;
	case 10:
		return "C3H8";
		break;
	case 11:
		return "iC4H10";
		break;
	case 12:
		return "CO2";
		break;
	case 13:
		return "neoC5H12";
		break;
	case 14:
		return "H2O";
		break;
	case 15:
		return "O2";
		break;
	case 16:
		return "N2";
		break;
	case 17:
		return "NO";
		break;
	case 18:
		return "N2O";
		break;
	case 19:
		return "C2H4";
		break;
	case 20:
		return "C2H2";
		break;
	case 21:
		return "H2";
		break;
	case 22:
		return "D2";
		break;
	case 23:
		return "CO";
		break;
	case 24:
		return "Methylal";
		break;
	case 25:
		return "DME";
		break;
	case 26:
		return "Reid-Step";
		break;
	case 27:
		return "Maxwell-Model";
		break;
	case 28:
		return "Reid-Ramp";
		break;
	case 29:
		return "C2F6";
		break;
	case 30:
		return "SF6";
		break;
	case 31:
		return "NH3";
		break;
	case 32:
		return "C3H6";
		break;
	case 33:
		return "cC3H6";
		break;
	case 34:
		return "CH3OH";
		break;
	case 35:
		return "C2H5OH";
		break;
	case 36:
		return "C3H7OH";
		break;
	case 37:
		return "Cs";
		break;
	case 38:
		return "F2";
		break;
	case 39:
		return "CS2";
		break;
	case 40:
		return "COS";
		break;
	case 41:
		return "CD4";
		break;
	case 42:
		return "BF3";
		break;
	case 43:
		return "C2H2F4";
		break;
	case 44:
		return "TMA";
		break;
	case 45:
		return "paraH2";
		break;
	case 46:
		return "nC3H7OH";
		break;
	case 47:
		return "Ar";
		break;
	case 48:
		return "Kr";
		break;
	case 49:
		return "Xe";
		break;
	case 50:
		return "CHF3";
		break;
	case 51:
		return "CF3Br";
		break;
	case 52:
		return "C3F8";
		break;
	case 53:
		return "O3";
		break;
	case 54:
		return "Hg";
		break;
	case 55:
		return "H2S";
		break;
	case 56:
		return "nC4H10";
		break;
	case 57:
		return "nC5H12";
		break;
	case 58:
		return "N2";
		break;
	case 59:
		return "GeH4";
		break;
	case 60:
		return "SiH4";
		break;
	default:
		return "";
		break;
	}
	return "";
}
bool GetPath(string path, vector<string> &name)
{
	DIR *dir = opendir(path.c_str());
	if(!dir)
	{
		cout << "opendir " << path << " error." << endl;
		return false;
	}
	dirent *p = NULL;
	while((p = readdir(dir)) != NULL)
		if(p->d_name[0] != '.')
			name.push_back(path + "/" + string(p->d_name));

	closedir(dir);
	return true;
}