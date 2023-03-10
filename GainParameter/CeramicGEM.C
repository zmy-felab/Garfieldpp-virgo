#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"

using namespace Garfield;
using namespace std;
namespace
{
    void PrintUsage()
    {
        cerr << " Usage: " << endl;
        cerr << " ./CeramicGEM [-n nEvents] [-p Pressure] [-t Temperature] [-v Voltage] [-d Drift] [-i Induction] [-r Rim]" << endl;
    }
} // namespace
double GetPenning(double p, string gas1, string gas2, double c2);
int main(int argc, char *argv[])
{
    if (argc > 15)
    {
        PrintUsage();
        return 1;
    }

    // Default parameters
    int nEvents = 10;         // num
    double pressure = 1.;     // atm
    double temperature = 20.; // C
    int voltage = 700;        // Voltage
    double driftE = 1.;       // kV/cm
    double inductionE = 3.;   // kV/cm
    int rim = 80;             // um

    string rootname = "./result/";

    for (int i = 1; i < argc; i = i + 2)
    {
        if (string(argv[i]) == "-n")
        {
            nEvents = atoi(argv[i + 1]);
            rootname += "Ele_" + to_string(nEvents);
        }
        else if (string(argv[i]) == "-p")
        {
            pressure = atof(argv[i + 1]);
            char temp[10];
            sprintf(temp, "%.1fatm", pressure);
            rootname += "_P_" + string(temp);
        }
        else if (string(argv[i]) == "-t")
        {
            temperature = atof(argv[i + 1]);
            char temp[10];
            sprintf(temp, "%.1fC", temperature);
            rootname += "_T_" + string(temp);
        }
        else if (string(argv[i]) == "-v")
        {
            voltage = atoi(argv[i + 1]);
            rootname += "_VGEM_" + to_string(voltage);
        }
        else if (string(argv[i]) == "-d")
        {
            driftE = atof(argv[i + 1]);
            char temp[10];
            sprintf(temp, "%.1fkV", driftE);
            rootname += "_D_" + string(temp);
        }
        else if (string(argv[i]) == "-i")
        {
            inductionE = atof(argv[i + 1]);
            char temp[10];
            sprintf(temp, "%.1fkV", inductionE);
            rootname += "_I_" + string(temp);
        }
        else if (string(argv[i]) == "-r")
        {
            rim = atoi(argv[i + 1]);
            rootname += "_R_" + to_string(rim);
        }
        else
        {
            PrintUsage();
            return 1;
        }
    }

    char ansysPath[50];
    sprintf(ansysPath, "./ansys/D_%.1fkV_I_%.1fkV/%dV/", driftE, inductionE, voltage);

    // Start time
    time_t t;
    time(&t);
    struct tm *lt = localtime(&t);

    char timeRecond[255];
    sprintf(timeRecond, "Start   time: %d/%02d/%02d %02d:%02d:%02d", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    const bool saveData = true;
    const bool plotDrift = false;
    const bool plotField = false;
    const bool plotLine = false;
    const bool plotMesh = false;
    const bool driftIon = false;

    // Information of detector [cm]
    const double pitch = 0.06;
    // const double dia = 0.02;
    const double ceramic = 168.e-4;
    const double metal = 18.e-4;
    const double drift = 0.2;
    const double induct = 0.2;

    // Load the field map.
    ComponentAnsys123 *thgem = new ComponentAnsys123();
    thgem->Initialise(string(ansysPath) + "ELIST.lis", string(ansysPath) + "NLIST.lis", string(ansysPath) + "MPLIST.lis", string(ansysPath) + "PRNSOL.lis", "mm");
    thgem->EnableMirrorPeriodicityX();
    thgem->EnableMirrorPeriodicityY();
    thgem->PrintRange();

    // Setup the gas.
    string gas1 = "ne", gas2 = "co2";
    double f1 = 90., f2 = 10.;
    string mixgas = gas1 + "_" + to_string(int(f1)) + "_" + gas2 + "_" + to_string(int(f2));
    string gasfilePath = "./GasFile/" + mixgas + ".gas";
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition(gas1, f1, gas2, f2);
    gas->SetTemperature(temperature + 273.15);
    gas->SetPressure(pressure * AtmosphericPressure);
    gas->LoadGasFile(gasfilePath);
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    // gas->PrintGas();
    // Set the Penning transfer efficiency.
    const double rPenning = GetPenning(pressure, gas1, gas2, f2 / 100.);
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, gas1);
    // Load the ion mobilities.
    if (driftIon)
    {
        if (gas1 == "he")
            gas->LoadIonMobility(string(getenv("GARFIELD_INSTALL")) + "/share/Garfield/Data/IonMobility_He+_He.txt");
        else if (gas1 == "ne")
            gas->LoadIonMobility(string(getenv("GARFIELD_INSTALL")) + "/share/Garfield/Data/IonMobility_Ne+_Ne.txt");
        else if (gas1 == "ar")
            gas->LoadIonMobility(string(getenv("GARFIELD_INSTALL")) + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
        else
            cout << "Please set correct nobe gas." << endl;

        if (gas2 == "co2")
            gas->LoadIonMobility(string(getenv("GARFIELD_INSTALL")) + "/share/Garfield/Data/IonMobility_CO2+_CO2.txt");
    }

    // Associate the gas with the corresponding field map material.
    const unsigned int nMaterials = thgem->GetNumberOfMaterials();
    for (unsigned int i = 0; i < nMaterials; ++i)
    {
        const double eps = thgem->GetPermittivity(i);
        if (eps == 1.)
            thgem->SetMedium(i, gas);
    }
    thgem->PrintMaterials();

    // Create the sensor.
    Sensor *sensor = new Sensor();
    sensor->AddComponent(thgem);
    sensor->SetArea(-5 * pitch, -5 * pitch, -induct - metal - ceramic / 2., 5 * pitch, 5 * pitch, drift + metal + ceramic / 2.);

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);

    AvalancheMC *aval_mc = nullptr;
    if (driftIon)
    {
        aval_mc = new AvalancheMC();
        aval_mc->SetSensor(sensor);
        aval_mc->SetDistanceSteps(2.e-4);
    }

    ViewDrift *driftView = new ViewDrift();
    if (plotDrift)
    {
        aval->EnablePlotting(driftView);
        if (driftIon)
            aval_mc->EnablePlotting(driftView);
    }

    int ne = 0, ni = 0, np = 0, npp = 0, ntotal = 0, ntotaleff = 0;
    double xe0 = 0., ye0 = 0., ze0 = 0.21, te0 = 0., ee0 = 0.1;
    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.;
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.;
    double xi1 = 0., yi1 = 0., zi1 = 0., ti1 = 0.;
    double xi2 = 0., yi2 = 0., zi2 = 0., ti2 = 0.;
    int statuse, statusi;

    TFile *ff;
    TTree *tt_pri, *tt_ele, *tt_ion;
    if (saveData)
    {
        rootname = rootname + "_" + mixgas + ".root";

        ff = new TFile(rootname.c_str(), "RECREATE");
        tt_pri = new TTree("pri", "Primary electrons");
        tt_pri->Branch("xe0", &xe0, "xe0/D");
        tt_pri->Branch("ye0", &ye0, "ye0/D");
        tt_pri->Branch("ze0", &ze0, "ze0/D");
        tt_pri->Branch("ne", &ne, "ne/I");
        tt_pri->Branch("ni", &ni, "ni/I");
        tt_pri->Branch("np", &np, "np/I");
        tt_pri->Branch("npp", &npp, "npp/I");
        // tt_gain->AutoSave();
        tt_ele = new TTree("ele", "Electrons information");
        tt_ele->Branch("xe1", &xe1, "xe1/D");
        tt_ele->Branch("ye1", &ye1, "ye1/D");
        tt_ele->Branch("ze1", &ze1, "ze1/D");
        tt_ele->Branch("te1", &te1, "te1/D");
        tt_ele->Branch("ee1", &ee1, "ee1/D");
        tt_ele->Branch("xe2", &xe2, "xe2/D");
        tt_ele->Branch("ye2", &ye2, "ye2/D");
        tt_ele->Branch("ze2", &ze2, "ze2/D");
        tt_ele->Branch("ee2", &ee2, "ee2/D");
        tt_ele->Branch("te2", &te2, "te2/D");
        tt_ele->Branch("statuse", &statuse, "statuse/I");
        // tt_ele->AutoSave();

        if (driftIon)
        {
            tt_ion = new TTree("ion", "Ions information");
            tt_ion->Branch("xi1", &xi1, "xe1/D");
            tt_ion->Branch("yi1", &yi1, "yi1/D");
            tt_ion->Branch("zi1", &zi1, "zi1/D");
            tt_ion->Branch("ti1", &ti1, "ti1/D");
            tt_ion->Branch("xi2", &xi2, "xi2/D");
            tt_ion->Branch("yi2", &yi2, "yi2/D");
            tt_ion->Branch("zi2", &zi2, "zi2/D");
            tt_ion->Branch("ti2", &ti2, "ti2/D");
            tt_ion->Branch("statusi", &statusi, "statusi/I");
            // tt_ion->AutoSave();
        }
    }
    for (int i = 0; i < nEvents; i++)
    {
        // Randomize the initial position. RndmUniform->[0,1) RndmUniformPos->(0,1)
        xe0 = -pitch / 2. + RndmUniform() * pitch;
        ye0 = -sqrt(3) * pitch / 2. + RndmUniform() * sqrt(3) * pitch;
        // ze0 = RndmUniformPos() * 0.2 + ceramic / 2.;

        aval->AvalancheElectron(xe0, ye0, ze0, te0, ee0, 0., 0., 0.);
        aval->GetAvalancheSize(ne, ni);
        np = aval->GetNumberOfElectronEndpoints();

        for (int j = 0; j < np; j++)
        {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, ee1, xe2, ye2, ze2, te2, ee2, statuse);
            if (saveData)
                tt_ele->Fill();

            // arrive to the readout plane
            if (ze2 <= -induct - metal - ceramic / 2.)
                npp++;

            if (driftIon)
            {
                aval_mc->DriftIon(xe1, ye1, ze1, te1);
                aval_mc->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, statusi);
                if (saveData)
                    tt_ion->Fill();
            }
        }
        if (saveData)
            tt_pri->Fill();

        ntotal += np;
        ntotaleff += npp;

        printf("%d/%d: %10.1lfum %10.1lfum  %10.1fum %10d %10d %10d %10d\n", i, nEvents, xe0 * 10000, ye0 * 10000, ze0 * 10000, ni, ne, np, npp);

        npp = 0;
    }

    if (saveData)
    {
        ff->Write();
        ff->Close();
    }

    printf("Average   Gain: %d / %d = %.2lf\n", ntotal, nEvents, (double)ntotal / nEvents);
    printf("Effective Gain: %d / %d = %.2lf\n", ntotaleff, nEvents, (double)ntotaleff / nEvents);

    if (plotField)
    {
        ViewField *fieldView = new ViewField();
        fieldView->SetComponent(thgem);
        fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
        double xmin = -3 * pitch / 2.;
        double xmax = 3 * pitch / 2.;
        double zmin = -0.05;
        double zmax = 0.07;
        fieldView->SetArea(xmin, zmin, xmax, zmax);
        fieldView->EnableAutoRange();
        TCanvas *cf = new TCanvas();
        fieldView->SetCanvas(cf);
        fieldView->PlotContour(); // e v p
        // fieldView->Plot("v", "CONT1"); // e v p; SCAT Box ARR COLZ TEXT CONT4Z CONT1 CONT2 CONT3
        // fieldView->PlotProfile(0., 0., -0.21, 0., 0., 0.41, "e"); // e v p

        if (plotLine)
        {
            vector<double> xf, yf, zf;
            fieldView->EqualFluxIntervals(xmin, 0, 0.99 * zmax, xmax, 0, 0.99 * zmax, xf, yf, zf, 60);
            fieldView->PlotFieldLines(xf, yf, zf, true, false);
        }
    }
    if (plotDrift)
    {
        TCanvas *cd = new TCanvas();
        if (plotMesh)
        {
            ViewFEMesh *meshView = new ViewFEMesh();
            meshView->SetComponent(thgem);
            meshView->SetViewDrift(driftView);
            meshView->SetPlane(0, -1, 0, 0, 0, 0);
            meshView->SetPlaneXZ(); // master
            meshView->SetArea(-3 * pitch, -3 * pitch, -induct - metal - ceramic / 2., 3 * pitch, 3 * pitch, drift + metal + ceramic / 2.);
            meshView->SetFillMesh(true);
            meshView->SetColor(0, kBlack);
            // set the color of ceramic
            meshView->SetColor(2, kYellow + 2);
            meshView->EnableAxes();
            meshView->SetYaxisTitle("z[cm]");
            meshView->SetCanvas(cd);
            meshView->Plot();
        }
        else
        {
            driftView->SetArea(-3 * pitch, -3 * pitch, -induct - metal - ceramic / 2., 3 * pitch, 3 * pitch, drift + metal + ceramic / 2.);
            driftView->SetCanvas(cd);
            driftView->Plot();
        }
    }

    // Print start and end time
    printf("%s\n", timeRecond);
    time(&t);
    lt = localtime(&t);
    printf("End     time: %d/%02d/%02d %02d:%02d:%02d\n", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    if (plotDrift || plotField)
        app.Run(kTRUE);
}
// calculate penning coefficient of the mixture gas.
// p was the pressure (atm), gas1 (he, ne, ar), gas2 (co2, n2, ch4, ic4h10...), c was the fraction of gas2 (co2 or n2).
double GetPenning(double p, string gas1, string gas2, double c2)
{
    double penning = 0.5;
    if (gas1 == "ar")
    {
        if (gas2 == "co2")
        {
            double a1 = 0.627898;
            double a2 = 0.041394;
            double a3 = 0.004716;
            double a4 = 0.001562;
            double a5 = 0.002422;
            double a6 = 0.027115;

            penning = (a5 * p * p * (1 - c2) * (1 - c2) + a1 * p * c2 + a4 * c2 + a3) / (a6 * p * p * (1 - c2) * (1 - c2) + p * c2 + a2);
        }
    }
    else if (gas1 == "ne")
    {
        if (gas2 == "co2")
        {
            double a1 = 0.71104;
            double a2 = 0.06323;
            double a3 = 0.03085;
            double a4 = 4.20089;
            double a5 = 0.07831;
            double a6 = 0.13235;
            double a7 = 1.47470;

            penning = (a5 * p * p * (1 - c2) * (1 - c2) + a7 * c2 * c2 + a1 * p * c2 + a3) / (a6 * p * p * (1 - c2) * (1 - c2) + a4 * c2 * c2 + p * c2 + a2);
        }
        else if (gas2 == "n2")
        {
            double a1 = 0.55802;
            double a2 = 0.00514;
            double a3 = 0.00206;
            double a4 = 0.55385;
            double a5 = 0.01153;
            double a6 = 0.02073;
            double a7 = 0.01;

            penning = (a5 * p * p * (1 - c2) * (1 - c2) + a7 * c2 * c2 + a1 * p * c2 + a3) / (a6 * p * p * (1 - c2) * (1 - c2) + a4 * c2 * c2 + p * c2 + a2);
        }
        else
        {
            cout << "Default Penning Coefficient" << endl;
        }
    }
    else
    {
        cout << "Default Penning Coefficient" << endl;
    }
    return penning;
}
