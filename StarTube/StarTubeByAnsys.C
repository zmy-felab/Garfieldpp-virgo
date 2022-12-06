#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TLine.h>

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
void DrawDetector(double tubeR, double wireR, double starX, int starN)
{
    TEllipse *cir = new TEllipse(0, 0, tubeR);
    TEllipse *wire = new TEllipse(0, 0, wireR);
    cir->SetFillStyle(0);
    cir->SetLineWidth(2);
    cir->Draw();
    wire->SetFillColor(kBlack);
    wire->Draw();
    TLine *line[starN];
    for (int i = 0; i < starN; i++)
    {
        double x1 = (tubeR - starX) * cos(2 * Pi / starN * i);
        double y1 = (tubeR - starX) * sin(2 * Pi / starN * i);
        double x2 = tubeR * cos(2 * Pi / starN * i);
        double y2 = tubeR * sin(2 * Pi / starN * i);
        line[i] = new TLine(x1, y1, x2, y2);
        line[i]->SetLineColor(kBlack);
        line[i]->SetLineWidth(2);
        line[i]->Draw();
    }
}
int main(int argc, char *argv[])
{
    constexpr bool plotField = true;
    constexpr bool plotDrift = true;
    constexpr bool plotSignal = true;
    constexpr bool driftIon = true;

    // Start time
    time_t t;
    time(&t);
    struct tm *lt = localtime(&t);

    char timeRecond[255];
    sprintf(timeRecond, "Start   time: %d/%02d/%02d %02d:%02d:%02d", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    const int nEvents = atoi(argv[1]); // event num

    // information of detector [cm]
    const double zLength = 0.01;
    const int zPeriod = 1000;
    const double tubeR = 2.54 / 2;
    // const double tubeWallR = tubeR + 0.03;
    const double wireR = 0.01 / 2;
    const double starX = tubeR / 2.;
    // const double starY = 0.2;
    const int starN = 6;

    TApplication app("app", &argc, argv);

    plottingEngine.SetDefaultStyle();

    // load the field map
    const string ansysPath = "./ansys/";
    ComponentAnsys123 *tube = new ComponentAnsys123();
    tube->Initialise(ansysPath + "ELIST.lis", ansysPath + "NLIST.lis", ansysPath + "MPLIST.lis", ansysPath + "PRNSOL.lis", "mm");
    tube->EnableMirrorPeriodicityZ();
    tube->PrintRange();

    // setup the gas
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.0);
    gas->LoadGasFile("./gasFile/ar_90.0_co2_10.0_1.0atm.gas");
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
    const unsigned int nMaterials = tube->GetNumberOfMaterials();
    for (unsigned int i = 0; i < nMaterials; ++i)
    {
        const double eps = tube->GetPermittivity(i);
        if (eps == 1.)
            tube->SetMedium(i, gas);
    }
    tube->PrintMaterials();

    // Create the sensor.
    Sensor *sensor = new Sensor();
    sensor->AddComponent(tube);
    sensor->SetArea(-tubeR, -tubeR, -zLength * zPeriod / 2., tubeR, tubeR, zLength * zPeriod / 2.);

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);

    AvalancheMC *aval_mc;
    if (driftIon)
    {
        aval_mc = new AvalancheMC();
        aval_mc->SetSensor(sensor);
        aval_mc->SetDistanceSteps(2.e-4);
    }

    // Use Heed for simulating the photon absorption
    TrackHeed *track = new TrackHeed();
    track->SetSensor(sensor);

    ViewDrift *driftView;
    ViewField *fieldView;
    ViewSignal *signalView;
    TCanvas *cd, *cf, *cs;

    if (plotField)
    {
        fieldView = new ViewField();
        fieldView->SetComponent(tube);
        fieldView->SetPlane(0., 0., 1., 0., 0., 0.);
        fieldView->SetArea(-tubeR - 0.1, -tubeR - 0.1, tubeR + 0.1, tubeR + 0.1);
        fieldView->EnableAutoRange();
        cf = new TCanvas("Field", "Field", 550, 500);
        fieldView->SetCanvas(cf);
        // fieldView->PlotContour();
        fieldView->Plot("v", "CONT0Z");

        DrawDetector(tubeR, wireR, starX, starN);
    }
    if (plotDrift)
    {
        driftView = new ViewDrift();
        aval->EnablePlotting(driftView);
        track->EnablePlotting(driftView);
        if (driftIon)
            aval_mc->EnablePlotting(driftView);

        cd = new TCanvas("DriftLine", "DriftLine", 500, 500);
        driftView->SetCanvas(cd);
        driftView->SetArea(-tubeR - 0.1, -tubeR - 0.1, -zLength * zPeriod / 2, tubeR + 0.1, tubeR + 0.1, zLength * zPeriod / 2);
    }
    if (plotSignal)
    {
        tube->SetWeightingField(ansysPath + "ANODE.lis", "anode");

        const double tMin = -1.;
        const double tMax = 2000.;
        const double tStep = 1;
        const int nTimeBins = int((tMax - tMin) / tStep);
        sensor->AddElectrode(tube, "anode");
        sensor->SetTimeWindow(0, tStep, nTimeBins);
        sensor->SetTransferFunction(transfer);

        aval->EnableSignalCalculation();
        if (driftIon)
            aval_mc->EnableSignalCalculation();

        signalView = new ViewSignal();
        signalView->SetSensor(sensor);
        signalView->SetColourTotal(kMagenta);
        signalView->SetColourElectrons(kViolet);
        signalView->SetColourIons(kRed);

        cs = new TCanvas("Signal", "Signal", 1000, 500);
        cs->Divide(2, 1);
    }

    double x0 = 0, y0 = 0, z0 = 1., t0 = 0., e0 = 5900., dx = 0, dy = 0, dz = -1; // x ray information
    int nex = 0, nix = 0, netotal = 0;
    double xe0 = 0., ye0 = 0., ze0 = 0., te0 = 0., ee0 = 0., dx0 = 0., dy0 = 0., dz0 = 0.; // primary electron
    int ne = 0, ni = 0, np = 0;
    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.; // avalanche electron
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.;
    int statuse;
    double xi0 = 0., yi0 = 0., zi0 = 0., ti0 = 0.; // primary ion
    double xi1 = 0., yi1 = 0., zi1 = 0., ti1 = 0.; // avalanche electron ion
    double xi2 = 0., yi2 = 0., zi2 = 0., ti2 = 0.;
    int statusi;

    TFile *ff = new TFile("./result/x-ray.root", "RECREATE");
    TTree *tt_x = new TTree("x_ray", "number of electrons and ions");
    tt_x->Branch("x0", &x0, "x0/D");
    tt_x->Branch("y0", &y0, "y0/D");
    tt_x->Branch("nex", &nex, "nex/I");
    tt_x->Branch("nix", &nix, "nix/I");
    tt_x->Branch("netotal", &netotal, "netotal/I");
    TTree *tt_pri_e = new TTree("pri_e", "Primary electrons");
    tt_pri_e->Branch("xe0", &xe0, "xe0/D");
    tt_pri_e->Branch("ye0", &ye0, "ye0/D");
    tt_pri_e->Branch("ze0", &ze0, "ze0/D");
    tt_pri_e->Branch("te0", &te0, "te0/D");
    tt_pri_e->Branch("ee0", &ee0, "ee0/D");
    tt_pri_e->Branch("ne", &ne, "ne/I");
    tt_pri_e->Branch("ni", &ni, "ni/I");
    tt_pri_e->Branch("np", &np, "np/I");
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

    TTree *tt_ion;
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
    }

    for (int i = 0; i < nEvents; i++)
    {
        printf("-----> Event %d/%d Start:\n", i, nEvents);

        // Reset
        ne = 0, ni = 0, np = 0, netotal = 0;

        if (plotSignal)
            sensor->ClearSignal();
        if (plotDrift)
            driftView->Clear();

        // Randomize the initial position.
        while (1)
        {
            x0 = -tubeR / 2. + RndmUniform() * tubeR;
            y0 = -tubeR / 2. + RndmUniform() * tubeR;

            double r = sqrt(x0 * x0 + y0 * y0);
            if (r < tubeR && r > (tubeR - starX))
                if (tube->GetMedium(x0, y0, z0)->GetName() == "Ar/CO2")
                    break;
        }

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
                track->GetIon(j, xi0, yi0, zi0, ti0);
                aval_mc->DriftIon(xi0, yi0, zi0, ti0);
                aval_mc->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, statusi);
                tt_ion->Fill();
            }
        }
        for (int j = 0; j < nex; j++)
        {
            track->GetElectron(j, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);

            if (xe0 > 0 && xe0 < wireR)
                xe0 += 2 * wireR;
            else if (xe0 < 0 && xe0 > -wireR)
                xe0 += -2 * wireR;
            if (ye0 > 0 && ye0 < wireR)
                ye0 += 2 * wireR;
            else if (ye0 < 0 && ye0 > -wireR)
                ye0 += -2 * wireR;

            // aval->DriftElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
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
            printf("Ele: %d/%d: %8.2lfmm %8.2lfmm %8.2fmm %5d %5d %5d\n", j, nex, xe0 * 10, ye0 * 10, ze0 * 10, ni, ne, np);
        }

        tt_x->Fill();

        printf("Event %d Primary position: %8.2lfmm %8.2lfmm\n", i, x0 * 10, y0 * 10);
        printf("Event %d Average     Gain: %d / %d = %.2lf\n", i, netotal, nex, (double)netotal / nex);

        if (plotSignal)
        {

            cs->Clear("D");
            signalView->SetCanvas((TPad *)cs->cd(1));
            signalView->PlotSignal("anode", "t");
            signalView->PlotSignal("anode", "e");
            if (driftIon)
                signalView->PlotSignal("anode", "i");

            sensor->ConvoluteSignals();

            signalView->SetCanvas((TPad *)cs->cd(2));
            signalView->PlotSignal("anode", "t");
            signalView->PlotSignal("anode", "e");
            if (driftIon)
                signalView->PlotSignal("anode", "i");

            char name[50];
            sprintf(name, "./result/signal_%d.pdf", i);
            cs->SaveAs(name);
        }
        if (plotDrift)
        {
            cd->Clear();
            constexpr bool twod = true;
            constexpr bool drawaxis = true;
            driftView->Plot(twod, drawaxis);

            DrawDetector(tubeR, wireR, starX, starN);

            char name[50];
            sprintf(name, "./result/driftline_%d.pdf", i);
            cd->SaveAs(name);
        }
    }

    tt_x->Write();
    tt_pri_e->Write();
    tt_ele->Write();
    if (driftIon)
        tt_ion->Write();

    ff->Close();

    // Print start and end time
    printf("%s\n", timeRecond);
    time(&t);
    lt = localtime(&t);
    printf("End     time: %d/%02d/%02d %02d:%02d:%02d\n", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    if (plotDrift || plotField || plotSignal)
        app.Run(kTRUE);
}
