#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char *argv[])
{
    // drift line
    const bool plotDrift = true;
    const bool plotMesh = true;
    // electric field
    const bool plotField = true;
    const bool plotLine = true;
    // calculation method
    const bool useRKF = false;
    const bool useMc = false;
    const bool useMicro = true;

    // information of detector [cm]
    const double pitch = 0.06;
    // const double dia = 0.02;
    const double ceramic = 168.e-4;
    const double metal = 18.e-4;
    const double drift = 0.2;
    const double induct = 0.2;

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    // load the field map.
    string ansysPath = "./ansys/900V/";
    ComponentAnsys123 *thgem = new ComponentAnsys123();
    thgem->Initialise(ansysPath + "ELIST.lis", ansysPath + "NLIST.lis", ansysPath + "MPLIST.lis", ansysPath + "PRNSOL.lis", "mm");
    thgem->EnableMirrorPeriodicityX();
    thgem->EnableMirrorPeriodicityY();
    thgem->PrintRange();

    // setup the gas.
    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.0);
    gas->LoadGasFile("./ar_90.0_co2_10.0_1.0atm.gas");
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    // gas->PrintGas();
    // Set the Penning transfer efficiency.
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");

    // associate the gas with the corresponding field map material.
    const unsigned int nMaterials = thgem->GetNumberOfMaterials();
    for (unsigned int i = 0; i < nMaterials; ++i)
    {
        const double eps = thgem->GetPermittivity(i);
        if (eps == 1.)
            thgem->SetMedium(i, gas);
    }
    thgem->PrintMaterials();

    // create the sensor.
    Sensor *sensor = new Sensor();
    sensor->AddComponent(thgem);
    sensor->SetArea(-5 * pitch, -5 * pitch, -induct - metal - ceramic / 2., 5 * pitch, 5 * pitch, drift + metal + ceramic / 2.);

    DriftLineRKF *driftRKF;
    if(useRKF)
    {
        driftRKF = new DriftLineRKF();
        driftRKF->SetSensor(sensor);
        // driftRKF->SetGainFluctuationsPolya(0.4, 10);
    }
    AvalancheMC *avalMC;
    if(useMc)
    {
        avalMC = new AvalancheMC();
        avalMC->SetSensor(sensor);
    }
    AvalancheMicroscopic *avalMicro;
    if(useMicro)
    {
        avalMicro = new AvalancheMicroscopic();
        avalMicro->SetSensor(sensor);
    }

    ViewDrift *driftView = new ViewDrift();
    if (plotDrift)
    {
        driftView->SetArea(-3 * pitch, -3 * pitch, -induct - metal - ceramic / 2., 3 * pitch, 3 * pitch, drift + metal + ceramic / 2.);

        if(useRKF)
            driftRKF->EnablePlotting(driftView);
        if(useMc)
            avalMC->EnablePlotting(driftView);
        if(useMicro)
            avalMicro->EnablePlotting(driftView);
    }

    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0.;
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0.;
    int status;

    TFile *ff = new TFile("result.root", "RECREATE");
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

    const int nEvents = 300;
    for (int i = 0; i < nEvents; i++)
    {
        // Randomize the initial position.
        double x0 = -pitch / 2. + RndmUniform() * pitch;
        double y0 = -sqrt(3) * pitch / 2. + RndmUniform() * sqrt(3) * pitch;
        double z0 = 0.2;
        double t0 = 0.;

        int p = 0;
        if(useRKF)
        {
            xe1 = x0; ye1 = y0; ze1 = z0; te1 = t0;
            driftRKF->DriftElectron(x0, y0, z0, t0);
            driftRKF->GetEndPoint(xe2, ye2, ze2, te2, status);
            p = driftRKF->GetNumberOfDriftLinePoints(); // 可以对drifline上每一点的进一步进行追踪
        }
        if(useMc)
        {
            avalMC->DriftElectron(x0,y0, z0,t0);
            p = avalMC->GetNumberOfDriftLinePoints();
            avalMC->GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2, status);
        }
        if(useMicro)
        {
            double e0 = 0.1, ee1, ee2; // eV
            avalMicro->DriftElectron(x0, y0, z0, t0, e0);
            p = avalMicro->GetNumberOfElectronDriftLinePoints();
            avalMicro->GetElectronEndpoint(0, xe1, ye1, ze1, te1, ee1, xe2, ye2, ze2, te2, ee2, status);

        }
        printf("%d/%d: %10.1lfum %10.1lfum %10d\n", i, nEvents, x0 * 10000, y0 * 10000, p);

        tt_ele->Fill();
    }
    ff->Write();
    ff->Close();

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
    if (plotDrift || plotField || plotMesh)
        app.Run(kTRUE);
}
