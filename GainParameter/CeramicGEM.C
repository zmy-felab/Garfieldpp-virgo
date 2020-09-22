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
#include "Garfield/ViewSignal.hh"
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
int main(int argc, char *argv[])
{
    if (argc > 15)
    {
        PrintUsage();
        return 1;
    }

    // Default parameters
    int nEvents = 10;            // num
    double pressure = 760.;      // Torr
    double temperature = 293.15; // k
    int voltage = 900;           // Voltage
    double driftE = 1.;          // kV/cm
    double inductionE = 3.;      // kV/cm
    int rim = 80;                // um

    string ansysPath = "./ansys/" + to_string(voltage) + "V/";
    string rootname = "./result/EleInformation";

    for (int i = 1; i < argc; i = i + 2)
    {
        if (string(argv[i]) == "-n")
        {
            nEvents = atoi(argv[i + 1]);
            rootname += "_" + to_string(nEvents) + "Ele";
        }
        else if (string(argv[i]) == "-p")
        {
            pressure = atof(argv[i + 1]);
            rootname += "_" + to_string(pressure) + "Torr";
        }
        else if (string(argv[i]) == "-t")
        {
            temperature = atof(argv[i + 1]);
            rootname += "_" + to_string(temperature) + "K";
        }
        else if (string(argv[i]) == "-v")
        {
            voltage = atoi(argv[i + 1]);
            rootname += "_" + to_string(voltage) + "V";
        }
        else if (string(argv[i]) == "-d")
        {
            driftE = atof(argv[i + 1]);
            rootname += "_" + to_string(driftE) + "kV_D";
            ansysPath = "./ansys/" + to_string(driftE) + "kV/";
        }
        else if (string(argv[i]) == "-i")
        {
            inductionE = atof(argv[i + 1]);
            rootname += "_" + to_string(inductionE) + "kV_I";
            ansysPath = "./ansys/" + to_string(inductionE) + "kV/";
        }
        else if (string(argv[i]) == "-r")
        {
            rim = atoi(argv[i + 1]);
            rootname += "_" + to_string(rim) + "Rim";
            ansysPath = "./ansys/" + to_string(rim) + "um/";
        }
        else
        {
            PrintUsage();
            return 1;
        }
    }
    rootname += ".root";

    // Start time
    time_t t;
    time(&t);
    struct tm *lt = localtime(&t);

    char timeRecond[255];
    sprintf(timeRecond, "Start   time: %d/%02d/%02d %02d:%02d:%02d", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    const bool plotDrift = false;
    const bool plotField = false;
    const bool plotMesh = false;
    const bool driftIon = false;
    // const bool calculateSignal = false;

    // Information of detector [cm]
    const double pitch = 0.06;
    // const double dia = 0.02;
    const double ceramic = 168.e-4;
    const double metal = 18.e-4;
    const double drift = 0.2;
    const double induct = 0.2;

    // Load the field map.
    ComponentAnsys123 *thgem = new ComponentAnsys123();
    thgem->Initialise(ansysPath + "ELIST.lis", ansysPath + "NLIST.lis", ansysPath + "MPLIST.lis", ansysPath + "PRNSOL.lis", "mm");
    thgem->EnableMirrorPeriodicityX();
    thgem->EnableMirrorPeriodicityY();
    thgem->PrintRange();

    // Setup the gas.
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(temperature);
    gas->SetPressure(pressure);
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    // gas->PrintGas();
    // Set the Penning transfer efficiency.
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    // Load the ion mobilities.
    if (driftIon)
    {
        const string path = getenv("GARFIELD_HOME");
        gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
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
    sensor->SetArea(-42 * pitch, -42 * pitch, -induct - metal - ceramic / 2., 42 * pitch, 42 * pitch, drift + metal + ceramic / 2.);

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);

    AvalancheMC *aval_mc = new AvalancheMC();
    aval_mc->SetSensor(sensor);
    aval_mc->SetDistanceSteps(2.e-4);

    ViewDrift *driftView = new ViewDrift();
    if (plotDrift)
    {
        driftView->SetArea(-3 * pitch, -3 * pitch, -induct - metal - ceramic / 2., 3 * pitch, 3 * pitch, drift + metal + ceramic / 2.);
        aval->EnablePlotting(driftView);

        if (driftIon)
            aval_mc->EnablePlotting(driftView);
    }

    int ne = 0, ni = 0, np = 0, npp = 0, ntotal = 0;
    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., e1 = 0.;
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., e2 = 0.;
    double xi1 = 0., yi1 = 0., zi1 = 0., ti1 = 0.;
    double xi2 = 0., yi2 = 0., zi2 = 0., ti2 = 0.;
    int status;

    TFile *ff = new TFile(rootname.c_str(), "RECREATE");
    TTree *tt_ele = new TTree("ele", "Electrons information");
    tt_ele->Branch("xe1", &xe1, "xe1/D");
    tt_ele->Branch("ye1", &ye1, "ye1/D");
    tt_ele->Branch("ze1", &ze1, "ze1/D");
    tt_ele->Branch("te1", &te1, "te1/D");
    tt_ele->Branch("xe2", &xe2, "xe2/D");
    tt_ele->Branch("ye2", &ye2, "ye2/D");
    tt_ele->Branch("ze2", &ze2, "ze2/D");
    tt_ele->Branch("te2", &te2, "te2/D");
    tt_ele->Branch("status", &status, "status/I");
    TTree *tt_gain = new TTree("gain", "Electrons avalanche");
    tt_gain->Branch("ne", &ne, "ne/I");
    tt_gain->Branch("ni", &ni, "ni/I");
    tt_gain->Branch("np", &np, "np/I");
    tt_gain->Branch("npp", &npp, "npp/I");

    for (int i = 0; i < nEvents; i++)
    {
        // Randomize the initial position.
        double x0 = -pitch / 2. + RndmUniform() * pitch;
        double y0 = -sqrt(3) * pitch / 2. + RndmUniform() * sqrt(3) * pitch;
        // double x0 = 0.;
        // double y0 = 0.;
        double z0 = 0.21;
        double t0 = 0.;
        double e0 = 0.1;

        aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
        aval->GetAvalancheSize(ne, ni);
        np = aval->GetNumberOfElectronEndpoints();

        for (int j = 0; j < np; j++)
        {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            tt_ele->Fill();

            if (ze2 <= -induct - metal - ceramic / 2.)
                npp++;

            if (driftIon)
            {
                aval_mc->DriftIon(xe1, ye1, ze1, te1);
                aval_mc->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
            }
        }
        tt_gain->Fill();

        ntotal += np;

        printf("%d/%d: %10.1lfum %10.1lfum %10d %10d %10d %10d\n", i, nEvents, x0 * 10000, y0 * 10000, ni, ne, np, npp);

        npp = 0;
    }
    ff->Write();
    ff->Close();

    printf("Average Gain: %d / %d = %.2lf\n", ntotal, nEvents, (double)ntotal / nEvents);

    if (plotDrift)
    {
        TCanvas *cd = new TCanvas();
        driftView->SetCanvas(cd);
        driftView->Plot();
    }
    if (plotField)
    {
        ViewField *fieldView = new ViewField();
        fieldView->SetComponent(thgem);
        fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
        fieldView->SetArea(-3 * pitch / 2., -0.05, 3 * pitch / 2., 0.07);
        double vmin = 0., vmax = 0.;
        thgem->GetVoltageRange(vmin, vmax);
        fieldView->SetVoltageRange(vmin, vmax);
        // fieldView->SetElectricFieldRange(0., 10000.);
        TCanvas *cf = new TCanvas();
        fieldView->SetCanvas(cf);
        fieldView->PlotContour();                                 // e v p
        fieldView->PlotSurface("e");                              // e v p
        fieldView->Plot("v", "CONT1");                            // e v p; SCAT Box ARR COLZ TEXT CONT4Z CONT1 CONT2 CONT3
        fieldView->PlotProfile(0., 0., -0.21, 0., 0., 0.41, "e"); // e v p
    }
    if (plotMesh)
    {
        ViewFEMesh *meshView = new ViewFEMesh();
        meshView->SetComponent(thgem);
        meshView->SetPlane(0, -1, 0, 0, 0, 0);
        meshView->SetViewDrift(driftView);
        meshView->SetArea(-3 * pitch, -induct - metal - ceramic / 2., -3 * pitch, 3 * pitch, drift + metal + ceramic / 2., 3 * pitch);
        meshView->SetFillMesh(false);
        meshView->EnableAxes();
        meshView->SetYaxisTitle("z");
        TCanvas *cm = new TCanvas();
        meshView->SetCanvas(cm);
        meshView->Plot();
    }

    // Print start and end time
    printf("%s\n", timeRecond);
    time(&t);
    lt = localtime(&t);
    printf("End     time: %d/%02d/%02d %02d:%02d:%02d\n", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    if (plotDrift || plotField || plotMesh)
        app.Run(kTRUE);
}
