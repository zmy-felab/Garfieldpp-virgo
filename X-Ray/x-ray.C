#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
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
        cerr << " ./x-ray [-n nEvents] [-v Voltage]" << endl;
    }
} // namespace
int main(int argc, char *argv[])
{
    if (argc > 5)
    {
        PrintUsage();
        return 1;
    }

    // Default parameters
    int nEvents = 10;  // num
    int voltage = 900; // Voltage

    string rootname = "./result/Information";

    for (int i = 1; i < argc; i = i + 2)
    {
        if (string(argv[i]) == "-n")
        {
            nEvents = atoi(argv[i + 1]);
            rootname += "_" + to_string(nEvents) + "Ele";
        }
        else if (string(argv[i]) == "-v")
        {
            voltage = atoi(argv[i + 1]);
            rootname += "_" + to_string(voltage) + "V";
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
    // const double rim = 0.008;

    // Load the field map.
    const string ansysPath = "./ansys/" + to_string(voltage) + "V/";
    ComponentAnsys123 *thgem = new ComponentAnsys123();
    thgem->Initialise(ansysPath + "ELIST.lis", ansysPath + "NLIST.lis", ansysPath + "MPLIST.lis", ansysPath + "PRNSOL.lis", "mm");
    thgem->EnableMirrorPeriodicityX();
    thgem->EnableMirrorPeriodicityY();
    thgem->PrintRange();

    // Setup the gas.
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.0);
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    // gas->PrintGas();
    // Set the Penning transfer efficiency.
    const double rPenning = 0.426;
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
    sensor->SetArea(-10 * pitch, -10 * pitch, -induct - metal - ceramic / 2., 10 * pitch, 42 * pitch, drift + metal + ceramic / 2.);

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);

    AvalancheMC *aval_mc = new AvalancheMC();
    aval_mc->SetSensor(sensor);
    aval_mc->SetDistanceSteps(2.e-4);

    //
    TrackHeed *track = new TrackHeed();
    track->SetSensor(sensor);

    ViewDrift *driftView = new ViewDrift();
    if (plotDrift)
    {
        driftView->SetArea(-3 * pitch, -3 * pitch, -induct - metal - ceramic / 2., 3 * pitch, 3 * pitch, drift + metal + ceramic / 2.);
        aval->EnablePlotting(driftView);

        if (driftIon)
            aval_mc->EnablePlotting(driftView);
    }

    int nex = 0., nix = 0., netotal = 0., netotaleff = 0.;
    int ne = 0, ni = 0, np = 0, npp = 0;
    double xe0 = 0., ye0 = 0., ze0 = 0., te0 = 0., ee0 = 0., dx0 = 0., dy0 = 0., dz0 = 0.;
    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.;
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.;
    int statuse;
    double xi0 = 0., yi0 = 0., zi0 = 0., ti0 = 0.;
    double xi1 = 0., yi1 = 0., zi1 = 0., ti1 = 0.;
    double xi2 = 0., yi2 = 0., zi2 = 0., ti2 = 0.;
    int statusi;

    TFile *ff = new TFile(rootname.c_str(), "RECREATE");
    TTree *tt_x = new TTree("x_ray", "number of electrons and ions");
    tt_x->Branch("nex", &nex, "nex/I");
    tt_x->Branch("nix", &nix, "nix/I");
    tt_x->Branch("netotal", &netotal, "netotal/I");
    tt_x->Branch("netotaleff", &netotaleff, "netotaleff/I");
    TTree *tt_pri_e = new TTree("pri_e", "Primary electrons");
    tt_pri_e->Branch("xe0", &xe0, "xe0/D");
    tt_pri_e->Branch("ye0", &ye0, "ye0/D");
    tt_pri_e->Branch("ze0", &ze0, "ze0/D");
    tt_pri_e->Branch("te0", &te0, "te0/D");
    tt_pri_e->Branch("ee0", &ee0, "ee0/D");
    TTree *tt_ele = new TTree("ele", "Electrons information");
    tt_ele->Branch("xe1", &xe1, "xe1/D");
    tt_ele->Branch("ye1", &ye1, "ye1/D");
    tt_ele->Branch("ze1", &ze1, "ze1/D");
    tt_ele->Branch("te1", &te1, "te1/D");
    tt_ele->Branch("ee1", &ee1, "ee1/D");
    tt_ele->Branch("xe2", &xe2, "xe2/D");
    tt_ele->Branch("ye2", &ye2, "ye2/D");
    tt_ele->Branch("ze2", &ze2, "ze2/D");
    tt_ele->Branch("te2", &te2, "te2/D");
    tt_ele->Branch("ee2", &ee2, "ee2/D");
    tt_ele->Branch("statuse", &statuse, "statuse/I");
    TTree *tt_pri_i = new TTree("pri_i", "Primary electrons");
    tt_pri_i->Branch("xi0", &xi0, "xi0/D");
    tt_pri_i->Branch("yi0", &yi0, "yi0/D");
    tt_pri_i->Branch("zi0", &zi0, "zi0/D");
    tt_pri_i->Branch("ti0", &ti0, "ti0/D");
    TTree *tt_ion = new TTree("ion", "Ions information");
    tt_ion->Branch("xi1", &xi1, "xe1/D");
    tt_ion->Branch("yi1", &yi1, "yi1/D");
    tt_ion->Branch("zi1", &zi1, "zi1/D");
    tt_ion->Branch("ti1", &ti1, "ti1/D");
    tt_ion->Branch("xi2", &xi2, "xi2/D");
    tt_ion->Branch("yi2", &yi2, "yi2/D");
    tt_ion->Branch("zi2", &zi2, "zi2/D");
    tt_ion->Branch("ti2", &ti2, "ti2/D");
    tt_ion->Branch("statusi", &statusi, "statusi/I");
    TTree *tt_gain = new TTree("gain", "Electrons avalanche");
    tt_gain->Branch("ne", &ne, "ne/I");
    tt_gain->Branch("ni", &ni, "ni/I");
    tt_gain->Branch("np", &np, "np/I");
    tt_gain->Branch("npp", &npp, "npp/I");

    for (int i = 0; i < nEvents; i++)
    {
        printf("----> Event %d/%d Start:\n", i, nEvents);

        // Randomize the initial position.
        const double x0 = 0.;
        const double y0 = 0.;
        const double z0 = drift + metal + ceramic / 2.;
        const double t0 = 0.;
        // Sample the x ray energy, using the relative intensities according to XDB.
        const double e0 = 5900.; // eV
        // const double r = 167. * RndmUniform();
        // const double e0 = r < 100. ? 5898.8 : r < 150. ? 5887.6 : 6490.4;

        while (1)
        {
            track->TransportPhoton(x0, y0, z0, t0, e0, 0., 0., -1, nex, nix);
            if (nex != 0)
            {
                track->GetElectron(0, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
                if (ze0 > ceramic / 2.)
                    break;
            }
        }
        if (driftIon)
        {
            // primary ions
            for (int j = 0; j < nix; j++)
            {
                track->GetIon(i, xi0, yi0, zi0, ti0);
                aval_mc->DriftIon(xi0, yi0, zi0, ti0);
                tt_pri_i->Fill();
            }
        }
        for (int j = 0; j < nex; j++)
        {
            track->GetElectron(j, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
            tt_pri_e->Fill();

            // aval->DriftElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
            aval->AvalancheElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);

            np = aval->GetNumberOfElectronEndpoints();
            aval->GetAvalancheSize(ne, ni);

            for (int k = 0; k < np; k++)
            {
                aval->GetElectronEndpoint(k, xe1, ye1, ze1, te1, ee1, xe2, ye2, ze2, te2, ee2, statuse);
                tt_ele->Fill();

                if (ze2 <= -induct - metal - ceramic / 2.)
                    npp++;

                if (driftIon)
                {
                    // ions by electrons avalanche
                    aval_mc->DriftIon(xe1, ye1, ze1, te1);
                    aval_mc->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, statusi);
                    tt_ion->Fill();
                }
            }
            tt_gain->Fill();

            netotal += np;
            netotaleff += npp;

            // print information of the primary electrons avalanche
            printf("Ele: %d/%d: %10.1lfum %10.1lfum %10.1fum %6d %6d %6d %6d %10d %10d\n", j, nex, xe0 * 10000, ye0 * 10000, ze0 * 10000, ni, ne, np, npp, netotal, netotaleff);
            npp = 0;
        }

        tt_x->Fill();

        printf("Event %d Average    Gain: %d / %d = %.2lf\n", i, netotal, nex, (double)netotal / nex);
        printf("Event %d Efficiency Gain: %d / %d = %.2lf\n", i, netotaleff, nex, (double)netotaleff / nex);

        // Reset
        netotal = 0;
        netotaleff = 0;
    }
    ff->Write();
    ff->Close();

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
        fieldView->SetArea(-3 * pitch / 2., -0.1, 3 * pitch / 2., 0.1);
        double vmin = 0., vmax = 0.;
        thgem->GetVoltageRange(vmin, vmax);
        fieldView->SetVoltageRange(vmin, vmax);
        // fieldView->SetElectricFieldRange(0., 10000.);
        TCanvas *cf = new TCanvas();
        fieldView->SetCanvas(cf);
        fieldView->PlotContour(); // e v p
        // fieldView->Plot("v", "CONT1");                            // e v p; SCAT Box ARR COLZ TEXT CONT4Z CONT1 CONT2 CONT3
        // fieldView->PlotProfile(0., 0., -0.21, 0., 0., 0.41, "e"); // e v p
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
