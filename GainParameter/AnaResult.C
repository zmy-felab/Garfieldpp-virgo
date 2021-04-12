#include <iostream>
#include <fstream>
#include <dirent.h>

#include "TFile.h"
#include "TTree.h"

using namespace std;

bool GetPath(string path, vector<string> &name, string cut);

int AnaResult(string path)
{
    vector<string> filename;
    if (GetPath(path, filename, ".root"))
    {
        TFile *ff = new TFile("./data/pos_res.root", "recreate");
        TTree *tt;

        TFile *f;
        TTree *t_pri, *t_gain, *t_ele;
        ofstream data;

        for (int i = 0; i < filename.size(); i++)
        {
            cout << filename[i] << endl;

            string subdir = filename[i].substr(filename[i].find_last_of("/") + 1, filename[i].size() - filename[i].find_last_of("/") - 5);

            ff->mkdir(subdir.c_str());
            tt = new TTree("FWHM", "Position Resolution");
            double pos_x = 0, pos_y = 0;
            tt->Branch("X", &pos_x, "pos_x/D");
            tt->Branch("Y", &pos_y, "pos_y/D");

            f = new TFile(filename[i].c_str());

            if(f->IsOpen())
            {
                double xe0 = 0., ye0 = 0.;
                t_pri = (TTree *)f->Get("pri");
                t_pri->SetBranchAddress("xe0", &xe0);
                t_pri->SetBranchAddress("ye0", &ye0);

                int ne = 0, ni = 0, np = 0, npp = 0; // electron gain
                t_gain = (TTree *)f->Get("gain");
                // t_gain->SetBranchAddress("ne", &ne);
                // t_gain->SetBranchAddress("ni", &ni);
                t_gain->SetBranchAddress("np", &np);
                // t_gain->SetBranchAddress("npp", &npp);

                double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.; // electron information begin
                double xe2 = 0., ye2 = 0., ze2 = 0., te2 = 0., ee2 = 0.; // electron information end
                t_ele = (TTree *)f->Get("ele");
                t_ele->SetBranchAddress("xe1", &xe1);
                t_ele->SetBranchAddress("ye1", &ye1);
                t_ele->SetBranchAddress("ze1", &ze1);
                // t_ele->SetBranchAddress("te1", &te1);
                // t_ele->SetBranchAddress("ee1", &ee1);
                t_ele->SetBranchAddress("xe2", &xe2);
                t_ele->SetBranchAddress("ye2", &ye2);
                t_ele->SetBranchAddress("ze2", &ze2);
                // t_ele->SetBranchAddress("te2", &te2);
                // t_ele->SetBranchAddress("ee2", &ee2);

                Long64_t ele_pri = t_gain->GetEntries();
                Long64_t ele_ava = t_ele->GetEntries();

                Long64_t ele_pri_id = 0;
                Long64_t ele_ava_id = 0;

                Long64_t ele_start_drift = 0;
                Long64_t ele_start_copper_up = 0;
                Long64_t ele_start_ceramic = 0;
                Long64_t ele_start_copper_down = 0;
                Long64_t ele_start_induction = 0;

                Long64_t ele_end_drift = 0;
                Long64_t ele_end_copper_up = 0;
                Long64_t ele_end_ceramic = 0;
                Long64_t ele_end_copper_down = 0;
                Long64_t ele_end_induction = 0;
                Long64_t ele_end_readplane = 0;

                Long64_t gain_total = 0;
                Long64_t gain_total_drift = 0;
                Long64_t gain_total_gem = 0;
                Long64_t gain_total_induction = 0;
                Long64_t gain_eff = 0;
                Long64_t gain_eff_drift = 0;
                Long64_t gain_eff_gem = 0;
                Long64_t gain_eff_induction = 0;

                Long64_t ele_total = 0;
                Long64_t ele_tran = 0;

                for (int j = 0; j < ele_pri; j++)
                {
                    t_pri->GetEntry(j);
                    t_gain->GetEntry(j);

                    double sum_x = 0.;
                    double sum_y = 0.;
                    int ele_num = 0;

                    for (int k = 0; k < np; k++)
                    {
                        t_ele->GetEntry(ele_ava_id);

                        gain_total++;

                        if (ze1 == 0.21)
                        {
                            ele_total++;
                            if (ze2 < 0.0102 && np != 1)
                                ele_tran++;
                        }
                        else if (ze1 < 0.21 && ze1 > 0.0102)
                        {
                            ele_start_drift++;
                            gain_total_drift++;
                        }
                        else if (ze1 <= 0.0102 && ze1 > 0.0084)
                        {
                            ele_start_copper_up++;
                            gain_total_gem++;
                        }
                        else if (ze1 <= 0.0084 && ze1 > -0.0084)
                        {
                            ele_start_ceramic++;
                            gain_total_gem++;
                        }
                        else if (ze1 <= -0.0084 && ze1 > -0.0102)
                        {
                            ele_start_copper_down++;
                            gain_total_gem++;
                        }
                        else if (ze1 <= -0.0102 && ze1 > -0.2102)
                        {
                            ele_start_induction++;
                            gain_total_induction++;
                        }
                        else
                            ;

                        if (ze2 <= 0.21 && ze2 > 0.0102)
                            ele_end_drift++;
                        else if (ze2 <= 0.0102 && ze2 > 0.0084)
                            ele_end_copper_up++;
                        else if (ze2 <= 0.0084 && ze2 > -0.0084)
                            ele_end_ceramic++;
                        else if (ze2 <= -0.0084 && ze2 > -0.0102)
                            ele_end_copper_down++;
                        else if (ze2 <= -0.0102 && ze2 > -0.2102)
                            ele_end_induction++;
                        else
                        {
                            ele_end_readplane++;

                            gain_eff++;
                            if (ze1 < 0.21 && ze1 > 0.0102)
                                gain_eff_drift++;
                            else if (ze1 <= 0.0102 && ze1 > -0.0102)
                                gain_eff_gem++;
                            else if (ze1 <= -0.0102 && ze1 > -0.2102)
                                gain_eff_induction++;
                            else
                                ;

                            sum_x += xe2;
                            sum_y += ye2;
                            ele_num++;
                        }
                        ele_ava_id++;
                    }
                    if (ele_num != 0)
                    {
                        pos_x = (sum_x / ele_num - xe0) * 10;
                        pos_y = (sum_y / ele_num - ye0) * 10;
                        
                        tt->Fill();
                    }
                }
                ff->cd(subdir.c_str());
                tt->Write();

                data.open("./data/gain.txt", ios::app | ios::out);
                data << filename[i] << "\t" << gain_total << "\t" << gain_total_drift << "\t" << gain_total_gem << "\t" << gain_total_induction << "\t" << gain_eff << "\t" << gain_eff_drift << "\t" << gain_eff_gem << "\t" << gain_eff_induction << endl;
                data.close();
                data.open("./data/distribute.txt", ios::app | ios::out);
                data << filename[i] << "\t" << ele_start_drift << "\t" << ele_start_copper_up << "\t" << ele_start_ceramic << "\t" << ele_start_copper_down << "\t" << ele_start_induction << "\t" << ele_end_drift << "\t" << ele_end_copper_up << "\t" << ele_end_ceramic << "\t" << ele_end_copper_down << "\t" << ele_end_induction << "\t" << ele_end_readplane << endl;
                data.close();
                data.open("./data/transparency.txt", ios::app | ios::out);
                data << filename[i] << "\t" << ele_total << "\t" << ele_tran << endl;
                data.close();
            }            
            f->Close();
        }
        ff->Close();
    }
    return 0;
}
bool GetPath(string path, vector<string> &name, string cut)
{
    DIR *dir = opendir(path.c_str());
    if (!dir)
    {
        cout << "opendir " << path << " error." << endl;
        return false;
    }
    dirent *p = NULL;
    while ((p = readdir(dir)) != NULL)
    {
        if (p->d_name[0] != '.')
        {
            char *filetype = strrchr(p->d_name, '.');
            if (strcmp(filetype, cut.c_str()) == 0)
                name.push_back(path + "/" + string(p->d_name));
        }
    }
    closedir(dir);
    if (name.size() == 0)
    {
        cout << cut << " file tyle can not be found." << endl;
        return false;
    }
    return true;
}