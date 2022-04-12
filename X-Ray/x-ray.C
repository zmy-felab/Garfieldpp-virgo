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
double transfer(double t)
{
    constexpr double tau = 25.;
    return (t / tau) * exp(1 - t / tau);
}
int main(int argc, char *argv[])
{
    constexpr bool plotField = true;
    constexpr bool plotDrift = true;
    constexpr bool driftIon = false;
    constexpr bool calSignal = true;
    constexpr bool plotSignal = true;

    // Start time
    time_t t;
    time(&t);
    struct tm *lt = localtime(&t);

    char timeRecond[255];
    sprintf(timeRecond, "Start   time: %d/%02d/%02d %02d:%02d:%02d", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    int nEvents = atoi(argv[1]); // event num

    string rootname = "./result/x-ray.root";

    // information of detector [cm]
    const double pitch = 0.06;
    // const double dia = 0.02;
    const double ceramic = 168.e-4;
    const double metal = 18.e-4;
    const double drift = 0.2;
    const double induct = 0.2;
    // const double rim = 0.008;

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    // load the field map
    const string ansysPath = "./ansys/";
    ComponentAnsys123 *thgem = new ComponentAnsys123();
    thgem->Initialise(ansysPath + "ELIST.lis", ansysPath + "NLIST.lis", ansysPath + "MPLIST.lis", ansysPath + "PRNSOL.lis", "mm");
    thgem->EnableMirrorPeriodicityX();
    thgem->EnableMirrorPeriodicityY();
    thgem->PrintRange();
    if (calSignal)
    {
        thgem->SetWeightingField(ansysPath + "ANODE.lis", "anode");
        thgem->SetWeightingField(ansysPath + "GEMDOWN.lis", "gemdown");
        thgem->SetWeightingField(ansysPath + "GEMUP.lis", "gemup");
        thgem->SetWeightingField(ansysPath + "CATHODE.lis", "cathode");
    }

    // setup the gas
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.0);
    gas->LoadGasFile("./GasFile/ar_90_co2_10.gas");
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
        const string path = getenv("GARFIELD_INSTALL");
        gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
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

    if (calSignal)
    {
        sensor->AddElectrode(thgem, "anode");
        sensor->AddElectrode(thgem, "gemdown");
        sensor->AddElectrode(thgem, "gemup");
        sensor->AddElectrode(thgem, "cathode");

        const double tMin = -1.;
        const double tMax = 300.;
        const double tStep = 0.1;
        const int nTimeBins = int((tMax - tMin) / tStep);
        sensor->SetTimeWindow(0, tStep, nTimeBins);
        sensor->SetTransferFunction(transfer);
    }

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    if (calSignal)
        aval->EnableSignalCalculation();

    AvalancheMC *aval_mc;
    if (driftIon)
    {
        aval_mc = new AvalancheMC();
        aval_mc->SetSensor(sensor);
        aval_mc->SetDistanceSteps(2.e-4);
        if (calSignal)
            aval_mc->EnableSignalCalculation();
    }
    // use Heed for simulating the photon absorption
    TrackHeed *track = new TrackHeed();
    track->SetSensor(sensor);

    ViewDrift *driftView;
    ViewFEMesh *meshView = nullptr;
    ViewField *fieldView;
    ViewSignal *signalView;
    TCanvas *cd, *cf, *cs;

    if (plotField)
    {
        fieldView = new ViewField();
        fieldView->SetComponent(thgem);
        // Set the normal vector of the viewing plane (xz plane).
        fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
        // Set the plot limits in the current viewing plane.
        double xmin = -0.5 * pitch, xmax = 0.5 * pitch;
        double zmin = -0.07, zmax = 0.07;
        fieldView->SetArea(xmin, zmin, xmax, zmax);
        fieldView->EnableAutoRange();
        cf = new TCanvas("Field", "", 500, 500);
        fieldView->SetCanvas(cf);
        // fieldView->PlotContour();
        fieldView->Plot("v", "CONT4Z");

        if (!meshView)
            meshView = new ViewFEMesh();

        meshView->SetComponent(thgem);
        meshView->SetFillMesh(true);
        meshView->SetColor(0, kYellow + 2);
        meshView->SetColor(2, kGray);
        meshView->SetPlane(0, -1, 0, 0, 0, 0);
        meshView->SetArea(xmin, zmin, xmax, zmax);
        meshView->SetCanvas(cf);
        meshView->Plot(true);
    }
    if (plotDrift)
    {
        driftView = new ViewDrift();
        aval->EnablePlotting(driftView);
        track->EnablePlotting(driftView);
        if (driftIon)
            aval_mc->EnablePlotting(driftView);

        cd = new TCanvas("DriftLine", "DriftLine", 500, 500);

        if (!meshView)
            meshView = new ViewFEMesh();

        meshView->SetComponent(thgem);
        meshView->SetPlane(0, -1, 0, 0, 0, 0);
        meshView->SetArea(-3 * pitch, -induct - metal - ceramic / 2., 3 * pitch, drift + metal + ceramic / 2.);
        meshView->SetFillMesh(true);
        meshView->SetColor(0, kGray);
        meshView->SetColor(2, kYellow + 2);
        meshView->EnableAxes();
        meshView->SetViewDrift(driftView);
    }
    if (calSignal)
    {
        signalView = new ViewSignal();
        signalView->SetSensor(sensor);
        if (plotSignal)
        {
            cs = new TCanvas("Signal", "Signal", 1200, 600);
            cs->Divide(4, 2);
        }
    }

    double x0 = 0, y0 = 0, z0 = 0.21, t0 = 0., e0 = 5900., dx = 0, dy = 0, dz = -1;        // x ray information
    int nex = 0., nix = 0., netotal = 0., netotaleff = 0.;                                 // x ray
    double xe0 = 0., ye0 = 0., ze0 = 0., te0 = 0., ee0 = 0., dx0 = 0., dy0 = 0., dz0 = 0.; // primary electron
    int ne = 0, ni = 0, np = 0, npp = 0;
    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.; // avalanche electron
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.;
    int statuse;
    double xi0 = 0., yi0 = 0., zi0 = 0., ti0 = 0.; // primary ion
    double xi1 = 0., yi1 = 0., zi1 = 0., ti1 = 0.; // avalanche electron ion
    double xi2 = 0., yi2 = 0., zi2 = 0., ti2 = 0.;
    int statusi;
    double raw_sa = 0., raw_sd = 0., raw_su = 0., raw_sc = 0.; // raw signal
    double con_sa = 0., con_sd = 0., con_su = 0., con_sc = 0.; // convolute signal

    TFile *ff = new TFile(rootname.c_str(), "RECREATE");
    TTree *tt_x = new TTree("x_ray", "number of electrons and ions");
    tt_x->Branch("x0", &x0, "x0/D");
    tt_x->Branch("y0", &y0, "y0/D");
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
    tt_pri_e->Branch("ne", &ne, "ne/I");
    tt_pri_e->Branch("ni", &ni, "ni/I");
    tt_pri_e->Branch("np", &np, "np/I");
    tt_pri_e->Branch("npp", &npp, "npp/I");
    TTree *tt_ele = new TTree("ele", "Avalanche electrons information");
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

    TTree *tt_pri_i, *tt_ion, *tt_s;
    if (driftIon)
    {
        tt_pri_i = new TTree("pri_i", "Primary ions");
        tt_pri_i->Branch("xi0", &xi0, "xi0/D");
        tt_pri_i->Branch("yi0", &yi0, "yi0/D");
        tt_pri_i->Branch("zi0", &zi0, "zi0/D");
        tt_pri_i->Branch("ti0", &ti0, "ti0/D");

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
    }
    if (calSignal)
    {
        tt_s = new TTree("signal", "Signal information");
        tt_s->Branch("AnodeRaw", &raw_sa, "raw_sa/D");
        tt_s->Branch("GEMDownRaw", &raw_sd, "raw_sd/D");
        tt_s->Branch("GEMUpRaw", &raw_su, "raw_su/D");
        tt_s->Branch("CathodeRaw", &raw_sc, "raw_sc/D");
        tt_s->Branch("AnodeCon", &con_sa, "con_sa/D");
        tt_s->Branch("GEMDownCon", &con_sd, "con_sd/D");
        tt_s->Branch("GEMUpCon", &con_su, "con_su/D");
        tt_s->Branch("CathodeCon", &con_sc, "con_sc/D");
    }

    for (int i = 0; i < nEvents; i++)
    {
        printf("----> Event %d/%d Start:\n", i, nEvents);
        if (calSignal)
            sensor->ClearSignal();
        if (plotDrift)
            driftView->Clear();

        // Randomize the initial position.
        x0 = -pitch / 2. + RndmUniform() * pitch;
        y0 = -sqrt(3) * pitch / 2. + RndmUniform() * sqrt(3) * pitch;
        while (1)
        {
            track->TransportPhoton(x0, y0, z0, t0, e0, dx, dy, dz, nex, nix);
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
                track->GetIon(j, xi0, yi0, zi0, ti0);
                aval_mc->DriftIon(xi0, yi0, zi0, ti0);
                // aval_mc->GetIonEndpoint(0,);
                tt_pri_i->Fill();
            }
        }
        for (int j = 0; j < nex; j++)
        {
            track->GetElectron(j, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);

            // aval->DriftElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
            aval->AvalancheElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);

            np = aval->GetNumberOfElectronEndpoints();
            aval->GetAvalancheSize(ne, ni);

            for (int k = 0; k < np; k++)
            {
                aval->GetElectronEndpoint(k, xe1, ye1, ze1, te1, ee1, xe2, ye2, ze2, te2, ee2, statuse);
                tt_ele->Fill();

                if (ze2 <= -induct)
                    npp++;

                if (driftIon)
                {
                    // ions by electrons avalanche
                    aval_mc->DriftIon(xe1, ye1, ze1, te1);
                    aval_mc->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, statusi);
                    tt_ion->Fill();
                }
            }
            tt_pri_e->Fill();

            netotal += np;
            netotaleff += npp;

            // print information of the primary electrons avalanche
            printf("Ele: %d/%d: %10.1lfum %10.1lfum %10.1fum %6d %6d %6d %6d\n", j, nex, xe0 * 10000, ye0 * 10000, ze0 * 10000, ni, ne, np, npp);
            npp = 0;
        }

        tt_x->Fill();

        printf("Event %d Primary position: %10.1lfum %10.1lfum\n", i, x0 * 10000, y0 * 10000);
        printf("Event %d Average    Gain: %d / %d = %.2lf\n", i, netotal, nex, (double)netotal / nex);
        printf("Event %d Efficiency Gain: %d / %d = %.2lf\n", i, netotaleff, nex, (double)netotaleff / nex);

        // Reset
        netotal = 0;
        netotaleff = 0;

        if (calSignal)
        {
            double tstart, tstep;
            unsigned int nsteps;
            sensor->GetTimeWindow(tstart, tstep, nsteps);
            for (unsigned int j = 0; j < nsteps; j++)
            {
                raw_sa = sensor->GetSignal("anode", j);
                raw_sd = sensor->GetSignal("gemdown", j);
                raw_su = sensor->GetSignal("gemup", j);
                raw_sc = sensor->GetSignal("cathode", j);
            }
            // Exporting induced signal to a csv file
            char name[50];
            sprintf(name, "./result/signal_anode_raw_%d", i);
            sensor->ExportSignal("anode", name);

            if (plotSignal)
            {
                cs->Clear("D");
                signalView->SetCanvas((TPad *)cs->cd(1));
                signalView->PlotSignal("anode", true, false, false);
                signalView->SetCanvas((TPad *)cs->cd(2));
                signalView->PlotSignal("gemdown", true, false, false);
                signalView->SetCanvas((TPad *)cs->cd(3));
                signalView->PlotSignal("gemup", true, false, false);
                signalView->SetCanvas((TPad *)cs->cd(4));
                signalView->PlotSignal("cathode", true, false, false);
            }
            sensor->ConvoluteSignals();
            for (unsigned int j = 0; j < nsteps; j++)
            {
                con_sa = sensor->GetSignal("anode", j);
                con_sd = sensor->GetSignal("gemdown", j);
                con_su = sensor->GetSignal("gemup", j);
                con_sc = sensor->GetSignal("cathode", j);
            }
            // Exporting induced signal to a csv file
            sprintf(name, "./result/signal_anode_con_%d", i);
            sensor->ExportSignal("anode", name);

            if (plotSignal)
            {
                signalView->SetCanvas((TPad *)cs->cd(5));
                signalView->PlotSignal("anode", true, false, false);
                signalView->SetCanvas((TPad *)cs->cd(6));
                signalView->PlotSignal("gemdown", true, false, false);
                signalView->SetCanvas((TPad *)cs->cd(7));
                signalView->PlotSignal("gemup", true, false, false);
                signalView->SetCanvas((TPad *)cs->cd(8));
                signalView->PlotSignal("cathode", true, false, false);

                sprintf(name, "./result/signal_%d.pdf", i);
                cs->SaveAs(name);
            }
            tt_s->Fill();
        }
        if (plotDrift)
        {
            cd->Clear();
            meshView->SetCanvas(cd);
            meshView->Plot();
            char name[50];
            sprintf(name, "./result/driftline_%d.pdf", i);
            cd->SaveAs(name);
        }
    }
    ff->Write();
    ff->Close();

    // Print start and end time
    printf("%s\n", timeRecond);
    time(&t);
    lt = localtime(&t);
    printf("End     time: %d/%02d/%02d %02d:%02d:%02d\n", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    if (plotDrift || plotField || plotSignal)
        app.Run(kTRUE);
}
