import ROOT
import os
import sys
import ctypes
import math

path = os.getenv('GARFIELD_INSTALL')
if sys.platform == 'darwin':
    ROOT.gSystem.Load(path + '/lib/libmagboltz.dylib')
    ROOT.gSystem.Load(path + '/lib/libGarfield.dylib')
else:
    ROOT.gSystem.Load(path + '/lib/libmagboltz.so')
    ROOT.gSystem.Load(path + '/lib/libGarfield.so')


def transfer(t):
    tau = 25.
    return (t/tau)*math.exp(1-t/tau)


#
plotField = True
plotDrift = True
driftIon = False
calSignal = True
plotSignal = True
nEvents = int(sys.argv[1])

# information of detector [cm]
pitch = 0.06
dia = 0.02
ceramic = 168.e-4
metal = 18.e-4
rim = 0.008
driftR = 0.2
inductR = 0.2

# Load the field map.
ansysPath = "./ansys/"
thgem = ROOT.Garfield.ComponentAnsys123()
thgem.Initialise(ansysPath+"ELIST.lis", ansysPath+"NLIST.lis",
                 ansysPath+"MPLIST.lis", ansysPath+"PRNSOL.lis", "mm")
thgem.EnableMirrorPeriodicityX()
thgem.EnableMirrorPeriodicityY()
thgem.PrintRange()

if calSignal:
    thgem.SetWeightingField(ansysPath+"ANODE.lis", "anode")

# Setup the gas
gas = ROOT.Garfield.MediumMagboltz("ar", 90, "co2", 10)
gas.SetTemperature(293.15)
gas.SetPressure(760.)
gas.LoadGasFile("./GasFile/ar_90_co2_10.gas")
gas.Initialise(True)

# Set Penning efficiency
rPenning = 0.426
lambdaPenning = 0.
gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar")

# Load the ion mobilities file
if driftIon:
    gas.LoadIonMobility(path + '/share/Garfield/Data/IonMobility_Ar+_Ar.txt')

thgem.SetGas(gas)
thgem.PrintMaterials()

#
sensor = ROOT.Garfield.Sensor()
sensor.AddComponent(thgem)
sensor.SetArea(-5 * pitch, -5 * pitch, -inductR - metal - ceramic /
               2., 5 * pitch, 5 * pitch, driftR + metal + ceramic / 2.)

#
aval = ROOT.Garfield.AvalancheMicroscopic()
aval.SetSensor(sensor)

#
track = ROOT.Garfield.TrackHeed()
track.SetSensor(sensor)

if driftIon:
    aval_mc = ROOT.Garfield.AvalancheMC()
    aval_mc.SetSensor(sensor)
    aval_mc.SetDistanceSteps(2.e-4)
    if calSignal:
        aval_mc.EnableSignalCalculation()

if calSignal:
    tMin = -1
    tMax = 300
    tStep = 0.1
    nTimeBins = int((tMax-tMin)/tStep)

    sensor.AddElectrode(thgem, "anode")
    sensor.SetTimeWindow(0, tStep, nTimeBins)
    sensor.SetTransferFunction(transfer)

    aval.EnableSignalCalculation()

    if plotSignal:
        signal = ROOT.Garfield.ViewSignal()
        signal.SetSensor(sensor)

if plotDrift:
    drift = ROOT.Garfield.ViewDrift()
    track.EnablePlotting(drift)
    aval.EnablePlotting(drift)
    if driftIon:
        aval_mc.EnablePlotting(drift)

    mesh = ROOT.Garfield.ViewFEMesh()
    mesh.SetComponent(thgem)
    mesh.SetViewDrift(drift)
    mesh.SetPlane(0, -1, 0, 0, 0, 0)
    mesh.SetArea(-3 * pitch, -3 * pitch, -inductR - metal - ceramic /
                 2., 3 * pitch, 3 * pitch, driftR + metal + ceramic / 2.)
    mesh.SetFillMesh(True)
    mesh.SetColor(0, ROOT.kYellow+2)
    mesh.SetColor(2, ROOT.kGray)
    mesh.EnableAxes()

if plotField:
    field = ROOT.Garfield.ViewField()
    field.SetComponent(thgem)
    field.SetPlane(0, -1, 0, 0, 0, 0)
    field.SetArea(-0.5*pitch, -0.07, 0.5*pitch, 0.07)
    field.EnableAutoRange()
    cf = ROOT.TCanvas("ElectricField", "ElectricField", 500, 500)
    field.SetCanvas(cf)
    # field.PlotContour()
    field.Plot("v", "CONT4Z")

# x ray information
x0 = ctypes.c_double(0.)
y0 = ctypes.c_double(0.)
z0 = ctypes.c_double(0.2)
e0 = ctypes.c_double(5900.)
#
nex = ctypes.c_int(0)
nix = ctypes.c_int(0)
netot = ctypes.c_int(0)
netoteff = ctypes.c_int(0)

# primary electron information
xe0 = ctypes.c_double(0.)
ye0 = ctypes.c_double(0.)
ze0 = ctypes.c_double(0.)
te0 = ctypes.c_double(0.)
ee0 = ctypes.c_double(0.)
dx0 = ctypes.c_double(0.)
dy0 = ctypes.c_double(0.)
dz0 = ctypes.c_double(0.)

# primary electron avalanche size
ne = ctypes.c_int(0)
ni = ctypes.c_int(0)
np = ctypes.c_int(0)
npp = ctypes.c_int(0)

# secondary(avalanche) electron informations
xe1 = ctypes.c_double(0.)
ye1 = ctypes.c_double(0.)
ze1 = ctypes.c_double(0.)
te1 = ctypes.c_double(0.)
ee1 = ctypes.c_double(0.)
xe2 = ctypes.c_double(0.)
ye2 = ctypes.c_double(0.)
ze2 = ctypes.c_double(0.)
te2 = ctypes.c_double(0.)
ee2 = ctypes.c_double(0.)
statuse = ctypes.c_int(0)

# primary ion
xi0 = ctypes.c_double(0.)
yi0 = ctypes.c_double(0.)
zi0 = ctypes.c_double(0.)
ti0 = ctypes.c_double(0.)

# secondary ion
xi1 = ctypes.c_double(0.)
yi1 = ctypes.c_double(0.)
zi1 = ctypes.c_double(0.)
ti1 = ctypes.c_double(0.)
xi2 = ctypes.c_double(0.)
yi2 = ctypes.c_double(0.)
zi2 = ctypes.c_double(0.)
ti2 = ctypes.c_double(0.)
statusi = ctypes.c_int(0)

for i in range(nEvents):
    print("-----> Event %d/%d start:" % (i, nEvents))

    # Reset
    ne.value = 0
    ni.value = 0
    if calSignal:
        sensor.ClearSignal()
    if plotDrift:
        drift.Clear()

    # Randomize the initial position
    x0.value = -0.5*pitch+ROOT.Garfield.RndmUniform()*pitch
    y0.value = -0.5*math.sqrt(3)*pitch + \
        ROOT.Garfield.RndmUniform()*math.sqrt(3)*pitch
    while 1:
        track.TransportPhoton(x0, y0, z0, 0, e0, 0, 0, -1, nex, nix)
        # x ray interact with gas
        if nex.value > 0:
            track.GetElectron(0, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0)
            # Determine if the electron is in the drift region.
            if ze0.value > ceramic/2.:
                break
    if driftIon:
        for j in range(nix.value):
            track.GetIon(j, xi0, yi0, zi0, ti0)
            aval_mc.DriftIon(xi0, yi0, zi0, ti0)
            aval_mc.GetIonEndpoint(0, xi1, yi1, zi1, ti1,
                                   xi2, yi2, zi2, ti2, statusi)

    for j in range(nex.value):
        npp.value = 0

        track.GetElectron(j, xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0)

        # aval.DriftElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0)
        aval.AvalancheElectron(xe0, ye0, ze0, te0, ee0, dx0, dy0, dz0)
        aval.GetAvalancheSize(ne, ni)
        np.value = aval.GetNumberOfElectronEndpoints()

        for k in range(np.value):
            aval.GetElectronEndpoint(
                k, xe1, ye1, ze1, te1, ee1, xe2, ye2, ze2, te2, ee2, statuse)
            if ze2.value <= -inductR:
                npp.value += 1
            if driftIon:
                aval_mc.DriftIon(xe1, ye1, ze1, te1)
                aval_mc.GetIonEndpoint(
                    0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, statusi)

        netot.value += np.value
        netoteff.value += npp.value

        print("Ele: %d/%d: %6.1fum %6.1fum %6.1fum %5d%5d%5d%5d" % (j, nex.value, xe0.value *
              1000, ye0.value*1000, ze0.value*1000, ne.value, ni.value, np.value, npp.value))

    print("Event %d Primary position: %6.1fum %6.1fum" %
          (i, x0.value*1000, y0.value*1000))
    print("Event %d Average     gain: %d / %d = %.2f" %
          (i, netot.value, nex.value, netot.value/nex.value))
    print("Event %d Efficiency  gain: %d / %d = %.2f" %
          (i, netoteff.value, nex.value, netoteff.value/nex.value))

    if calSignal:
        name = "./result/signal_anode_raw_{}".format(i)
        sensor.ExportSignal("anode", name)

        if plotSignal:
            cs = ROOT.TCanvas("Signal", "Signal", 1000, 500)
            cs.Divide(2, 1)
            signal.SetCanvas(cs.cd(1))
            signal.PlotSignal("anode", "t")

        sensor.ConvoluteSignals()
        name = "./result/signal_anode_con_{}".format(i)
        sensor.ExportSignal("anode", name)

        if plotSignal:
            signal.SetCanvas(cs.cd(2))
            signal.PlotSignal("anode", "t")
            name = "./result/signal_anode_{}.pdf".format(i)
            cs.SaveAs(name)

    if plotDrift:
        cd = ROOT.TCanvas("DriftLine", "DriftLine", 500, 500)
        mesh.SetCanvas(cd)
        mesh.Plot(True)
        name = "./result/driftLine_{}.pdf".format(i)
        cd.SaveAs(name)
