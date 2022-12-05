#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>


#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/SolidWire.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

using namespace std;
using namespace Garfield;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout << "Please set the star number, for example: ./StarTube 6" << endl;
        return 1;
    }
    
    const double wireV = 500.; 

    const double tubeR = 2.54/2; // [cm]
    const double wireR = 50e-4;
    const double zHalfL = 3.;

    const int starN = atoi(argv[1]);
    const double starX = 20e-3;
    const double starY = tubeR /4.;


    MediumMagboltz *gas = new MediumMagboltz();
    gas->SetComposition("ar", 90., "co2", 10.);
    gas->SetTemperature(293.15);
    gas->SetPressure(760.);
    gas->LoadGasFile("../ar_90.0_co2_10.0_1.0atm.gas");
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

    // vector<SolidBox*> box(starN, nullptr);
    // for(int i = 0; i < starN; i++)
    // {
    //     double angle = 2*Pi/starN * i;
    //     box[i] = new SolidBox(0.75*tubeR*sin(angle), 0.75*tubeR*cos(angle), 0, zHalfL, starX, starY, sin(angle), cos(angle), 0);
    //     box[i]->SetBoundaryPotential(0);
    //     box[i]->SetColour(kBlue);
    //     geo->AddSolid(box[i], Cu);
    // }

    SolidTube *tube = new SolidTube(0, 0, 0, tubeR, zHalfL);
    tube->SetBoundaryPotential(0);
    tube->SetColour(kViolet);
    geo->AddSolid(tube, gas);

    ComponentNeBem3d *nebem = new ComponentNeBem3d();
    nebem->SetGeometry(geo);
    nebem->SetTargetElementSize(0.1);
    nebem->SetMinMaxNumberOfElements(3, 15);
    nebem->EnableDebugging();
    nebem->Initialise();
    nebem->DisableDebugging();

    ViewField *fieldView = new ViewField();
    fieldView->SetComponent(nebem);
    fieldView->SetArea(-tubeR/2, -tubeR/2, -zHalfL, tubeR/2, tubeR/2, zHalfL);
    fieldView->SetPlane(0, 0, 1, 0, 0, 0);
    fieldView->Plot("e", "CONT4Z");

    //
    Sensor * sensor = new Sensor();
    sensor->AddComponent(nebem);
    sensor->SetArea(-tubeR, -tubeR, -zHalfL, tubeR, tubeR, zHalfL);

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

    TCanvas *cd = new TCanvas("Drift", "Drift", 600, 600);

    ViewGeometry *geoView = new ViewGeometry();
    geoView->SetGeometry(geo);
    geoView->SetArea(-tubeR, -tubeR, -zHalfL, tubeR, tubeR, zHalfL);
    geoView->SetPlane(0, 0, 1, 0, 0, 0);
    // geoView->Rotate(Pi/2.);
    geoView->Plot();

    ViewDrift *driftView = new ViewDrift();
    aval->EnablePlotting(driftView);
    aval_mc->EnablePlotting(driftView);
    track->EnablePlotting(driftView);
    driftView->SetCanvas(cd);
    driftView->SetArea(-tubeR, -tubeR, -zHalfL, tubeR, tubeR, zHalfL);
    driftView->SetPlaneXY();

    int nex = 0, nix = 0;
    double xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0;
    while (1)
    {
        track->TransportPhoton(0.5, 0.5, 0, 0, 5900, 0, 0, 1, nex, nix);
        if(nex >0) break;
    }
    cout << "zhou\t" << nex << "\t" << nix << endl;
    for(int i = 0; i < nex; i++)
    {
        track->GetElectron(i, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);

        cout << i << "\t" << xe0 << "\t" << ye0 << "\t" << ze0 << endl; 

        aval->DriftElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0);
    }
    
    // aval->AvalancheElectron(1, 1, 1, 0, 0.1, 0, 0 ,0);

    driftView->Plot(true, false);

    TApplication app("app", &argc, argv);
    app.Run(true);

    return 0;
}

