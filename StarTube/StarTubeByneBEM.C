#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TEllipse.h>
#include <TLine.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/SolidWire.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

using namespace std;
using namespace Garfield;

void DrawDetector(double tubeR, double wireR, double starX, int starN);

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Please set the star number, for example: ./StarTube 6" << endl;
        return 1;
    }
    constexpr bool photo = true;

    const double wireV = 300.;

    const double tubeR = 2.54 / 8; // [cm]
    const double wireR = 0.01 / 2;
    const double zHalfL = 1 / 2.;

    const int starN = atoi(argv[1]);
    const double starX = tubeR / 2.;
    const double starY = 0.01;
    const double starXPos = tubeR - 0.5 * starX;

    TApplication app("app", &argc, argv);

    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.);
    gas->LoadGasFile("./gasFile/ar_90.0_co2_10.0_1.0atm.gas");
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();

    const double rPenning = 0.426;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");

    const string path = getenv("GARFIELD_INSTALL");
    gas->LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

    //
    MediumConductor *Cu = new MediumConductor();

    GeometrySimple *geo = new GeometrySimple();

    SolidWire *wire = new SolidWire(0, 0, 0, wireR, zHalfL);
    wire->SetBoundaryPotential(wireV);
    wire->SetLabel("anode");
    wire->SetColour(kRed);
    geo->AddSolid(wire, Cu);

    SolidTube *tube = new SolidTube(0, 0, 0, tubeR, zHalfL);
    tube->SetBoundaryPotential(0);
    tube->SetSectors(10);
    tube->SetColour(kBlack);
    geo->AddSolid(tube, gas);

    vector<SolidBox *> box(starN, nullptr);
    for (int i = 0; i < starN; i++)
    {
        double angle = 2 * Pi / starN * i;
        box[i] = new SolidBox(starXPos * cos(angle), starXPos * sin(angle), 0, zHalfL, starY / 2., starX / 2., cos(angle), sin(angle), 0);
        box[i]->SetBoundaryPotential(0);
        box[i]->SetColour(kBlue);
        geo->AddSolid(box[i], Cu);
    }

    ComponentNeBem3d *nebem = new ComponentNeBem3d();
    nebem->SetGeometry(geo);
    nebem->SetTargetElementSize(0.1);
    nebem->SetMinMaxNumberOfElements(3, 15);
    nebem->EnableDebugging();
    nebem->Initialise();
    nebem->DisableDebugging();

    ViewField *fieldView = new ViewField();
    fieldView->SetComponent(nebem);
    fieldView->SetPlane(0, 0, 1, 0, 0, 0);
    fieldView->SetArea(-tubeR - 0.1, -tubeR - 0.1, -zHalfL - 0.1, tubeR + 0.1, tubeR + 0.1, zHalfL + 0.1);
    fieldView->EnableAutoRange();
    TCanvas *cf = new TCanvas("Field", "Field", 550, 500);
    fieldView->SetCanvas(cf);
    // fieldView->PlotContour();
    fieldView->Plot("v", "CONT0Z");
    DrawDetector(tubeR, wireR, starX, starN);

    ViewGeometry *geoView = new ViewGeometry();
    geoView->SetGeometry(geo);
    geoView->SetArea(-tubeR, -tubeR, -zHalfL, tubeR, tubeR, zHalfL);
    geoView->SetPlane(0, 1, 0, 0, 0, 0);
    TCanvas *cg = new TCanvas("Geometry", "Geometry", 500, 500);
    geoView->SetCanvas(cg);
    geoView->Plot();

    //
    Sensor *sensor = new Sensor();
    sensor->AddComponent(nebem);
    sensor->SetArea(-tubeR, -tubeR, -zHalfL, tubeR, tubeR, zHalfL);

    DriftLineRKF *drift = new DriftLineRKF();
    drift->SetSensor(sensor);
    // drift->SetGainFluctuationsPolya(0.4, 10);

    AvalancheMicroscopic *aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    aval->EnableSignalCalculation();

    AvalancheMC *aval_mc;
    aval_mc = new AvalancheMC();
    aval_mc->SetSensor(sensor);
    aval_mc->SetDistanceSteps(2.e-4);
    aval_mc->EnableSignalCalculation();

    // use Heed for simulating the photon absorption
    TrackHeed *track = new TrackHeed();
    track->SetSensor(sensor);

    ViewDrift *driftView = new ViewDrift();
    drift->EnablePlotting(driftView);
    aval->EnablePlotting(driftView);
    aval_mc->EnablePlotting(driftView);
    track->EnablePlotting(driftView);
    driftView->SetArea(-tubeR - 0.1, -tubeR - 0.1, -zHalfL - 0.1, tubeR + 0.1, tubeR + 0.1, zHalfL + 0.1);
    TCanvas *cd = new TCanvas("Drift", "Drift", 500, 500);
    driftView->SetCanvas(cd);

    int nex = 0, nix = 0;
    double xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0;
    while (photo)
    {
        track->TransportPhoton(tubeR / 2., tubeR / 2., 0, 0, 5900, 0, 0, 1, nex, nix);
        if (nex > 0)
            break;
    }
    for (int i = 0; i < nex; i++)
    {
        track->GetElectron(i, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);

        cout << i << "/" << nex << ": " << xe0 << " " << ye0 << " " << ze0 << endl;

        // aval->DriftElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
        aval_mc->DriftElectron(xe0, ye0, ze0, te0);
        // drift->DriftElectron(xe0, ye0, ze0, te0);
    }

    driftView->Plot(true, true);
    DrawDetector(tubeR, wireR, starX, starN);

    app.Run(true);

    return 0;
}
void DrawDetector(double tubeR, double wireR, double starX, int starN)
{
    TEllipse *cir = new TEllipse(0, 0, tubeR);
    TEllipse *wire = new TEllipse(0, 0, wireR);
    cir->SetFillStyle(0);
    cir->SetLineWidth(3);
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
        line[i]->SetLineWidth(3);
        line[i]->Draw();
    }
}