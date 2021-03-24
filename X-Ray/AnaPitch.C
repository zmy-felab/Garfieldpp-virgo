#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

void OutputTime();
void AnaPitch()
{
    // Start time
    time_t t;
    time(&t);
    struct tm *lt = localtime(&t);

    const int nFiles = 5;                                            // root文件数
    const int nPitch = 8;                                            // 200 400 600 800 1000 1200 1400 1600
    const int nBins = 40 - 1;                                        //
    const double zPlane = -0.2 - (18.e-4) - (168.e-4) / 2.;          // cm,-induct - metal - ceramic / 2.
    double PitchWidth[nPitch] = {0.};                                // cm，每个直方图的bin宽
    double xStart[nPitch] = {0.}, xEnd[nPitch] = {0.};               // cm
    double axis[nPitch][nBins];                                      // cm，不同条宽下，各直方图每个bin的中心位置
    TH1F *hx[nPitch], *hy[nPitch];                                   // 不同条宽的直方图，用来寻找单个事件的位置
    TH1F *hx_m[nPitch], *hy_m[nPitch], *hx_c[nPitch], *hy_c[nPitch]; // 不同周期重心法与中心法直方图

    for (int i = 0; i < nPitch; i++)
    {
        PitchWidth[i] = (i + 1) * 200 / 10000.; // cm

        xStart[i] = -nBins / 2. * PitchWidth[i];
        xEnd[i] = nBins / 2. * PitchWidth[i];

        hx[i] = new TH1F("", "; x[cm]; Counts", nBins, xStart[i], xEnd[i]);
        hy[i] = new TH1F("", "; x[cm]; Counts", nBins, xStart[i], xEnd[i]);

        for (int j = 0; j < nBins; j++)
            axis[i][j] = xStart[i] + (j + 0.5) * PitchWidth[i];

        char namex_m[20], namey_m[20], namex_c[20], namey_c[20];
        sprintf(namex_m, "X_Mean_%dum", (i + 1) * 200);
        sprintf(namey_m, "Y_Mean_%dum", (i + 1) * 200);
        sprintf(namex_c, "X_Center_%dum", (i + 1) * 200);
        sprintf(namey_c, "Y_Center_%dum", (i + 1) * 200);

        hx_m[i] = new TH1F(namex_m, "; x[mm]; Counts", nBins, xStart[i] * 10, xEnd[i] * 10);
        hy_m[i] = new TH1F(namey_m, "; y[mm]; Counts", nBins, xStart[i] * 10, xEnd[i] * 10);
        hx_c[i] = new TH1F(namex_c, "; x[mm]; Counts", nBins, xStart[i] * 10, xEnd[i] * 10);
        hy_c[i] = new TH1F(namey_c, "; y[mm]; Counts", nBins, xStart[i] * 10, xEnd[i] * 10);
    }

    int nex = 0., nix = 0., netotal = 0., netotaleff = 0.; // electrons that x-ray ion gas
    double XRayCenterX = 0., XRayCenterY = 0.;             // x-ray reconstruction position
    double xMeanPitch[nPitch], yMeanPitch[nPitch];         // reconstruction position by gravity method
    double xCenterPitch[nPitch], yCenterPitch[nPitch];     // reconstruction position by center method
    TTree *tt_x = new TTree("x_ray", "X Ray information");
    tt_x->Branch("nex", &nex, "nex/I");
    tt_x->Branch("nix", &nix, "nix/I");
    tt_x->Branch("netotal", &netotal, "netotal/I");
    tt_x->Branch("netotaleff", &netotaleff, "netotaleff/I");
    tt_x->Branch("XRayCenterX", &XRayCenterX, "XRayCenterX/D");
    tt_x->Branch("XRayCenterY", &XRayCenterY, "XRayCenterY/D");
    tt_x->Branch("xMeanPitch", xMeanPitch, "xMeanPitch[8]/D");
    tt_x->Branch("yMeanPitch", yMeanPitch, "yMeanPitch[8]/D");
    tt_x->Branch("xCenterPitch", xCenterPitch, "xCenterPitch[8]/D");
    tt_x->Branch("yCenterPitch", yCenterPitch, "yCenterPitch[8]/D");

    double xe0 = 0., ye0 = 0., ze0 = 0., te0 = 0., ee0 = 0.; // electron information
    int ne = 0, ni = 0, np = 0, npp = 0;                     // electron gain
    double eleMeanX = 0., eleMeanY = 0.;                 // electron reconstruction position
    TTree *tt_pri_e = new TTree("pri_e", "Primary Electrons information");
    tt_pri_e->Branch("xe0", &xe0, "xe0/D");
    tt_pri_e->Branch("ye0", &ye0, "ye0/D");
    tt_pri_e->Branch("ze0", &ze0, "ze0/D");
    tt_pri_e->Branch("te0", &te0, "te0/D");
    tt_pri_e->Branch("ee0", &ee0, "ee0/D");
    tt_pri_e->Branch("ne", &ne, "ne/I");
    tt_pri_e->Branch("ni", &ni, "ni/I");
    tt_pri_e->Branch("np", &np, "np/I");
    tt_pri_e->Branch("npp", &npp, "npp/I");
    tt_pri_e->Branch("eleMeanX", &eleMeanX, "eleCenterX/D");
    tt_pri_e->Branch("eleMeanY", &eleMeanY, "eleCenterY/D");

    double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.; // electron information begin
    double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.; // electron information end

    TChain *ch_x_ray = new TChain("x_ray");
    TChain *ch_pri_e = new TChain("pri_e");
    TChain *ch_pri_g = new TChain("gain");
    TChain *ch_ele = new TChain("ele");
    for (int i = 0; i < nFiles; i++)
    {
        char temp[100];
        sprintf(temp, "./result/Information_%dEvents_900V.root", 185 + i);
        ch_x_ray->Add(temp);
        ch_pri_e->Add(temp);
        ch_pri_g->Add(temp);
        ch_ele->Add(temp);
    }
    ch_pri_e->AddFriend(ch_pri_g);

    ch_x_ray->SetBranchAddress("nex", &nex);
    ch_x_ray->SetBranchAddress("nix", &nix);
    ch_x_ray->SetBranchAddress("netotal", &netotal);
    ch_x_ray->SetBranchAddress("netotaleff", &netotaleff);

    ch_pri_e->SetBranchAddress("xe0", &xe0);
    ch_pri_e->SetBranchAddress("ye0", &ye0);
    ch_pri_e->SetBranchAddress("ze0", &ze0);
    ch_pri_e->SetBranchAddress("te0", &te0);
    ch_pri_e->SetBranchAddress("ee0", &ee0);
    ch_pri_e->SetBranchAddress("ne", &ne);
    ch_pri_e->SetBranchAddress("ni", &ni);
    ch_pri_e->SetBranchAddress("np", &np);
    ch_pri_e->SetBranchAddress("npp", &npp);

    ch_ele->SetBranchAddress("xe1", &xe1);
    ch_ele->SetBranchAddress("ye1", &ye1);
    ch_ele->SetBranchAddress("ze1", &ze1);
    ch_ele->SetBranchAddress("te1", &te1);
    ch_ele->SetBranchAddress("ee1", &ee1);
    ch_ele->SetBranchAddress("xe2", &xe2);
    ch_ele->SetBranchAddress("ye2", &ye2);
    ch_ele->SetBranchAddress("ze2", &ze2);
    ch_ele->SetBranchAddress("te2", &te2);
    ch_ele->SetBranchAddress("ee2", &ee2);

    Long64_t PriEleID = 0, AvaEleID = 0;
    Long64_t nXRayEntries = ch_x_ray->GetEntries();
    for (Long64_t i = 0; i < nXRayEntries; i++)
    {
        ch_x_ray->GetEntry(i);
        double xRaySum = 0, yRaySum = 0;
        for (Long64_t j = 0; j < nex; j++)
        {
            ch_pri_e->GetEntry(PriEleID);
            double xEleSum = 0., yEleSum = 0.;
            for (Long64_t k = 0; k < np; k++)
            {
                ch_ele->GetEntry(AvaEleID);
                if (ze2 <= zPlane)
                {
                    xEleSum += xe2;
                    yEleSum += ye2;
                    for (int l = 0; l < nPitch; l++)
                    {
                        hx[l]->Fill(xe2);
                        hy[l]->Fill(ye2);
                    }
                }
                AvaEleID++;
            }
            eleMeanX = xEleSum / npp;
            eleMeanY = yEleSum / npp;

            tt_pri_e->Fill();

            xRaySum += xEleSum;
            yRaySum += yEleSum;

            PriEleID++;
        }
        XRayCenterX = xRaySum / netotaleff;
        XRayCenterY = yRaySum / netotaleff;

        // position caculation
        for (int k = 0; k < nPitch; k++)
        {
            double xSumPitch = 0., ySumPitch = 0.;
            vector<int> xCenter, yCenter;
            for (int m = 0; m < nBins; m++)
            {
                int temp_x = hx[k]->GetBinContent(m + 1);
                if (temp_x != 0)
                {
                    xSumPitch += temp_x * axis[k][m];
                    xCenter.push_back(m);
                }

                int temp_y = hy[k]->GetBinContent(m + 1);
                if (temp_y != 0)
                {
                    ySumPitch += temp_y * axis[k][m];
                    yCenter.push_back(m);
                }
            }

            // Gravity method
            xMeanPitch[k] = xSumPitch / hx[k]->GetEntries();
            yMeanPitch[k] = ySumPitch / hy[k]->GetEntries();

            // Center method
            if (xCenter.size() % 2 == 0)
                xCenterPitch[k] = axis[k][xCenter[xCenter.size() / 2.]];
            else
                xCenterPitch[k] = axis[k][xCenter[xCenter.size() / 2.]];

            if (yCenter.size() % 2 == 0)
                yCenterPitch[k] = axis[k][yCenter[yCenter.size() / 2.]];
            else
                yCenterPitch[k] = axis[k][yCenter[yCenter.size() / 2.]];

            hx_m[k]->Fill(xMeanPitch[k] * 10);
            hy_m[k]->Fill(yMeanPitch[k] * 10);
            hx_c[k]->Fill(xCenterPitch[k] * 10);
            hy_c[k]->Fill(yCenterPitch[k] * 10);
        }
        tt_x->Fill();

        cout << "--> Event " << i << "/" << nXRayEntries << ": pri_ele = " << nex << setw(15) << " ave_gain = " << double(netotaleff) / nex << setw(12) << "pos_x = " << XRayCenterX * 10e4 << setw(12) << "pos_y = " << XRayCenterY * 10e4 << endl;

        // reset hist
        for (int k = 0; k < nPitch; k++)
        {
            hx[k]->Reset();
            hy[k]->Reset();
        }
    }

    TFile *fWrite = new TFile("./result_x_ray.root", "RECREATE");
    fWrite->mkdir("Analyze");
    fWrite->mkdir("Draw/X");
    fWrite->mkdir("Draw/Y");
    fWrite->cd("Analyze");
    tt_x->Write();
    tt_pri_e->Write();
    for (int i = 0; i < nPitch; i++)
    {
        fWrite->cd("Draw/X");
        hx_m[i]->Write();
        hx_c[i]->Write();

        fWrite->cd("Draw/Y");
        hy_m[i]->Write();
        hy_c[i]->Write();
    }
    fWrite->Close();

    // Print start and end time
    printf("Start   time: %d/%02d/%02d %02d:%02d:%02d\n", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);

    time(&t);
    lt = localtime(&t);
    printf("End     time: %d/%02d/%02d %02d:%02d:%02d\n", lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);
}
