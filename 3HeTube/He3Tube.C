#include <iostream>
#include <fstream>
#include <sstream>

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ComponentConstant.hh"
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
    constexpr bool plotDrift = true;
    constexpr bool plotSignal = true;
    constexpr bool driftIon = true;

    const double atm = atof(argv[1]);   // gas pressure [atm]
    const double rWire = (15.e-4) / 2;  // Wire radius [cm]
    const double rTube = 2.54 / 2;      // Outer radius of the tube [cm]
    const double vWire = atof(argv[2]); // Voltages [V]
    const double vTube = 0.;
    const int nEvents = atoi(argv[3]); // event num

    // set x ray position, direction and energy.
    const double x0 = 0.5, y0 = 0.5, z0 = 0, t0 = 0; // cm
    const double dx = 0, dy = 0, dz = 1;
    const double e0 = 5900.; // eV

    // set the time window [ns] for the signal calculation.
    const double tMin = 0.;
    const double tMax = 2000. * atm;
    const double tStep = 1 * atm;
    const int nTimeBins = int((tMax - tMin) / tStep);

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    // make a gas medium.
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90, "co2", 10);
    gas->SetTemperature(293.15);
    gas->SetPressure(atm * AtmosphericPressure);
    gas->LoadGasFile("ar_90.0_co2_10.0_1.0atm.gas");
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    // gas->PrintGas();
    // set the Penning transfer efficiency.
    const double rPenning = 0.45;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    const std::string path = std::getenv("GARFIELD_INSTALL");
    gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

    // make a component with analytic electric field.
    ComponentAnalyticField *cmp = new ComponentAnalyticField();
    cmp->SetMedium(gas);
    cmp->AddWire(0, 0, 2 * rWire, vWire, "s"); // Add the wire in the centre.
    cmp->AddTube(rTube, vTube, 0, "t");        // Add the tube.
    cmp->AddReadout("s");                      // Request calculation of the weighting field.

    // make a sensor.
    Sensor *sensor = new Sensor();
    sensor->AddComponent(cmp);
    sensor->AddElectrode(cmp, "s");
    sensor->SetTimeWindow(tMin, tStep, nTimeBins);
    sensor->SetTransferFunction(transfer);

    // set up Heed.
    TrackHeed *track = new TrackHeed();
    track->SetSensor(sensor);

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    aval->EnableSignalCalculation();

    AvalancheMC *aval_mc = new AvalancheMC();
    aval_mc->SetSensor(sensor);
    aval_mc->SetDistanceSteps(2.e-4);
    aval_mc->EnableSignalCalculation();

    TCanvas *cD = nullptr;
    ViewCell *cellView = new ViewCell();
    ViewDrift *driftView = new ViewDrift();
    if (plotDrift)
    {
        cD = new TCanvas("cD", "", 500, 500);
        cellView->SetCanvas(cD);
        cellView->SetComponent(cmp);
        driftView->SetCanvas(cD);

        aval->EnablePlotting(driftView);
        aval_mc->EnablePlotting(driftView);
        track->EnablePlotting(driftView);
    }

    TCanvas *cS = nullptr;
    ViewSignal *signalView = new ViewSignal();
    if (plotSignal)
    {
        cS = new TCanvas("cS", "", 1200, 800);
        cS->Divide(3, 2);
        signalView->SetSensor(sensor);
    }

    int nex = 0., nix = 0., netotal = 0.;
    double rawSignal[nTimeBins], rawEleSignal[nTimeBins], rawIonSignal[nTimeBins];
    double conSignal[nTimeBins], conEleSignal[nTimeBins], conIonSignal[nTimeBins];
    int ne = 0, ni = 0, np = 0;
    double xe0 = 0., ye0 = 0., ze0 = 0., te0 = 0., ee0 = 0., dx0 = 0., dy0 = 0., dz0 = 0.;
    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.;
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.;
    int statuse;
    double xi0 = 0., yi0 = 0., zi0 = 0., ti0 = 0.;
    double xi1 = 0., yi1 = 0., zi1 = 0., ti1 = 0.;
    double xi2 = 0., yi2 = 0., zi2 = 0., ti2 = 0.;
    int statusi;

    string rootname = "EleIonInfo_" + string(argv[1]) + "atm_" + string(argv[2]) + "V.root";
    TFile *ff = new TFile(rootname.c_str(), "RECREATE");
    TTree *tt_x = new TTree("x_ray", "X ray information");
    tt_x->Branch("nex", &nex, "nex/I");
    tt_x->Branch("nix", &nix, "nix/I");
    tt_x->Branch("netotal", &netotal, "netotal/I");
    tt_x->Branch("rawSignal", rawSignal, "rawSignal[2000]/D");
    tt_x->Branch("rawEleSignal", rawEleSignal, "rawEleSignal[2000]/D");
    tt_x->Branch("rawIonSignal", rawIonSignal, "rawIonSignal[2000]/D");
    tt_x->Branch("conSignal", conSignal, "conSignal[2000]/D");
    tt_x->Branch("conEleSignal", conEleSignal, "conEleSignal[2000]/D");
    tt_x->Branch("conIonSignal", conIonSignal, "conIonSignal[2000]/D");
    TTree *tt_pri_e = new TTree("pri_ele", "Primary electrons infomation");
    tt_pri_e->Branch("xe0", &xe0, "xe0/D");
    tt_pri_e->Branch("ye0", &ye0, "ye0/D");
    tt_pri_e->Branch("ze0", &ze0, "ze0/D");
    tt_pri_e->Branch("te0", &te0, "te0/D");
    tt_pri_e->Branch("ee0", &ee0, "ee0/D");
    tt_pri_e->Branch("ne", &ne, "ne/I");
    tt_pri_e->Branch("ni", &ni, "ni/I");
    tt_pri_e->Branch("np", &np, "np/I");
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

    TTree *tt_pri_i, *tt_ion;
    if (driftIon)
    {
        tt_pri_i = new TTree("pri_ion", "Primary ions information");
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

    for (int i = 0; i < nEvents; i++)
    {
        printf("----> Event %d/%d Start:\n", i, nEvents);
        sensor->ClearSignal();

        while (1)
        {
            track->TransportPhoton(x0, y0, z0, t0, e0, dx, dy, dz, nex, nix);
            if (nex != 0)
                break;
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
            aval->AvalancheElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
            np = aval->GetNumberOfElectronEndpoints();
            aval->GetAvalancheSize(ne, ni);

            for (int k = 0; k < np; k++)
            {
                aval->GetElectronEndpoint(k, xe1, ye1, ze1, te1, ee1, xe2, ye2, ze2, te2, ee2, statuse);
                tt_ele->Fill();

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

            // print information of the primary electrons avalanche
            printf("Ele: %d/%d: %6.2lfmm %6.3lfmm %6.3lfmm %6.3lfns %6d %6d %6d %10d\n", j, nex, xe0 * 10, ye0 * 10, ze0 * 10, te0, ni, ne, np, netotal);
        }

        printf("Event %d/%d: AveGain %d / %d = %.2lf PriIon: %d\n", i, nEvents, netotal, nex, (double)netotal / nex, nix);

        // Reset
        netotal = 0;

        // get electron, ion, total signal in every time bin
        for (int j = 0; j < nTimeBins; j++)
        {
            rawSignal[j] = sensor->GetSignal("s", j);
            rawEleSignal[j] = sensor->GetElectronSignal("s", j);
            rawIonSignal[j] = sensor->GetIonSignal("s", j);
        }
        //
        if (plotDrift)
        {
            cD->Clear();
            cellView->Plot2d();
            constexpr bool twod = true;
            constexpr bool drawaxis = false;
            driftView->Plot(twod, drawaxis);
        }
        // save currnt signal of every event to csv file
        char name[50];
        sprintf(name, "./signals/raw_%.1lfatm_%d_%d", atm, int(vWire), i);
        sensor->ExportSignal("s", name);
        // plot raw current signal
        if (plotSignal)
        {
            cS->Clear("D");
            signalView->SetCanvas((TPad *)cS->cd(1));
            signalView->PlotSignal("s", true, false, false);
            signalView->SetCanvas((TPad *)cS->cd(2));
            signalView->PlotSignal("s", false, true, false);
            signalView->SetCanvas((TPad *)cS->cd(3));
            signalView->PlotSignal("s", false, false, true);
        }
        //
        sensor->ConvoluteSignals();
        // get evergy convolute signal (electron, ion, total)
        for (int j = 0; j < nTimeBins; j++)
        {
            conSignal[j] = sensor->GetSignal("s", j);
            conEleSignal[j] = sensor->GetElectronSignal("s", j);
            conIonSignal[j] = sensor->GetIonSignal("s", j);
        }
        // save voltage signal of every event to csv file
        sprintf(name, "./signals/con_%.1lfatm_%d_%d", atm, int(vWire), i);
        sensor->ExportSignal("s", name);
        // plot voltage signal
        if (plotSignal)
        {
            signalView->SetCanvas((TPad *)cS->cd(4));
            signalView->PlotSignal("s", true, false, false);
            signalView->SetCanvas((TPad *)cS->cd(5));
            signalView->PlotSignal("s", false, true, false);
            signalView->SetCanvas((TPad *)cS->cd(6));
            signalView->PlotSignal("s", false, false, true);
        }
        tt_x->Fill();
    }
    ff->Write();
    ff->Close();

    if (plotSignal || plotDrift)
        app.Run(kTRUE);
}