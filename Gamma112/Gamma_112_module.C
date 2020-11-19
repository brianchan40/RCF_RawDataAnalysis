using namespace std;

/// C++ headers

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include <string>
#include <cstring>
#include <TLorentzVector.h>

/// PicoDst headers
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoDstReader.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoDst.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoEvent.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoTrack.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTofHit.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTowHit.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoHelix.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoPhysicalHelix.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTofPidTraits.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/StRefMultCorr.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/CentralityMaker.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/particle.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/event.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoEmcTrigger.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoTrackCovMatrix.h"

/// Analysis headers
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StV0Type.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StDcaService.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/structs.h"
#include "./namespaces/gv_gamma.h"
#include "./Gamma_112_module.h"

//void __attribute__((constructor)) LoadLib();
/////////////////////////////main program starts here/////////////////////////////////////
void Gamma_112_module(int cen = 1, int opt_weight = 1, const Char_t *inFile = "./st_physics_12160020_raw_1010001.picoDst.root", const TString JobIDName = "1234")      //main_function
{
    set_debug();

    delete gRandom;
    gRandom = new TRandom3(0);

    char fname_old[200];
    TString fname_new;
    sprintf(fname_old, "cen%d.weight_112_module.root", cen);
    Weight_Read = ReadWeight(fname_old);

    TString Name = "sched";
    Name.Append(JobIDName);
    Name.Append(TString::Format("cen%d.weight_112_module_new.root", cen));
    std::cout << "Name = " << Name << endl;
    fname_new = Name;
    std::cout << "fname_new = " << fname_new << endl;

    TChain *corr_tree = new TChain ("corr_tree");
    std::ifstream fin(inFile);
    string line1;
    char line[200];
    while(getline(fin, line1)) //loop wiill run till end of file
    {
        strcpy(line, line1.c_str());
        corr_tree->AddFile(line);
    }
    fin.close();


    corr_tree->SetBranchAddress("p_lambda", &p_lambda);
    corr_tree->SetBranchAddress("p_proton", &p_proton);
    corr_tree->SetBranchAddress("p_all", &p_all);
    corr_tree->SetBranchAddress("p_dau", &p_dau);
    corr_tree->SetBranchAddress("details", &details);

    char fname[200];
    char title[200];
    for(int j = 0; j < 9; j++)
    {
        sprintf(fname, "lam_pt_%d", j);
        sprintf(title, "Lambda vs Pt for Cen %d", j);
        lam_pt[j] = new TH1D(fname, title, 300, 0, 15);
        sprintf(fname, "p_pt_%d", j);
        sprintf(title, "Proton vs Pt for Cen %d", j);
        p_pt[j] = new TH1D(fname, title, 300, 0, 15);
    }

    for(int k = 0; k < 5; k++)
    {
        sprintf(fname, "Hist_v2_pt_obs2_caysm_%d", k);
        sprintf(title, "Hist_v2_pt_obs2_caysm %d", k);
        Hist_v2_pt_obs2_caysm[k] = new TProfile(fname, title, 300, 0, 15, -100, 100, "");
        sprintf(fname, "Hist_v2_pt_obs2_p_caysm_%d", k);
        sprintf(title, "Hist_v2_pt_obs2_p_caysm %d", k);
        Hist_v2_pt_obs2_p_caysm[k] = new TProfile(fname, title, 300, 0, 15, -100, 100, "");
    }

    for(int sp_pt = 0; sp_pt < 15; sp_pt++){
        sprintf(fname, "Parity_int_obs3_splitpt_%d", sp_pt);
        sprintf(title, "Parity_int_obs3_splitpt %d", sp_pt);
        pParity_int_obs3_splitpt[sp_pt] = new TProfile(fname, title, 4, 0.5, 4.5, -100, 100, "");
        sprintf(fname, "Parity_int_ss_obs3_splitpt_%d", sp_pt);
        sprintf(title, "Parity_int_ss_obs3_splitpt %d", sp_pt);
        pParity_int_ss_obs3_splitpt[sp_pt] = new TProfile(fname, title, 4, 0.5, 4.5, -100, 100, "");
        sprintf(fname, "Delta_int_ss_obs3_splitpt_%d", sp_pt);
        sprintf(title, "Delta_int_ss_obs3_splitpt %d", sp_pt);
        pDelta_int_ss_obs3_splitpt[sp_pt] = new TProfile(fname, title, 4, 0.5, 4.5, -100, 100, "");
    }

    Int_t nentries = corr_tree->GetEntries();

    //loop through events
    for(int i = 0; i < nentries; i++)
    {
        if((i + 1) % 1000 == 0) cout << "Processing entry == " << i + 1 << " == out of " << nentries << ".\n";

        corr_tree->GetEntry(i);

        Run     = details->Run;
        //pV      = event->primaryVertex();
        pVz     = details->PVtxz;
        pVx     = details->PVtxx;
        pVx     = details->PVtxy;
        VPDvz   = details->VPDvz;
        RefMult = details->RefMult;
        TOFMult = details->TOFMult;
        NPTracks = details->num_trk;
        Day     = (int)((Run % 1000000) / 1000);
        Day2    = (int)((Run % 1000000) / 10);
        Day3    = (int)((Run % 1000000) / 1);
        Centrality = details->cent;
        Eweight = details->Eweight;

        determine_particlelists();

        n_particle1 = particle1->size();
        n_particle2 = particle2->size();
        num_lam_tree->Fill(n_particle1);
        num_proton_tree->Fill(n_particle2);

        //if(debug_1) std::cout << "Centrality = " << Centrality << endl;

        if(!IsGoodEvent(cen)) continue;

        if(debug_1) std::cout << "debug_1 1" << endl;

        if(!CountCharge()) continue;

        if(debug_1) std::cout << "debug_1 2" << endl;

        //shuffle tracks for random EPs
        gv_gamma::iTrack.clear();
        Scount = Fcount / 2 - 1;
        for(int q = 0; q < Fcount; q++) gv_gamma::iTrack.push_back(q);
        random_shuffle(gv_gamma::iTrack.begin(), gv_gamma::iTrack.end());

        if(debug_1) std::cout << "debug_1 3" << endl;

        //TPC EP reconstruction
        MakeTPC_EP();

        //flatten EPs
        ShiftPsi();
        FillEP_resolution();

        //Store the flattened phi for POI (Lambdas)
        if(debug_1 && (n_particle1 != 0)) std::cout << "n_particle1 1 = " << n_particle1 << endl;

        Phi_new.resize(n_particle1);
        for(int trki = 0; trki < n_particle1; trki++)
        {
            //if((particle_option1 == 0) && ((particle1->at(trki).mass < 1.113) || (particle1->at(trki).mass > 1.119))) continue;
            gv_gamma::trk_mom->SetXYZ(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz);
            Eta   = gv_gamma::trk_mom->Eta();
            Pt    = gv_gamma::trk_mom->Pt();
            Charge = particle1->at(trki).Charge; //1 for lambda, -1 for antilambdas
            //DCAGlobal = p_all[trk].DCAGlobal();
            Phi   = gv_gamma::trk_mom->Phi();

            // if((Pt < 0.5) || (Pt > 2.0)) continue;

            TLorentzVector trk_lorentz(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz, sqrt(pow(gv_gamma::trk_mom->Mag(), 2) + pow(particle1->at(trki).mass, 2)));
            rapidity = trk_lorentz.Rapidity();
            // if(fabs(rapidity) > 1) continue;
            trk_lorentz.Clear();

            if(debug_1) std::cout << "Eta = " << Eta << endl;

            //if(!IsGoodPOI(Pt, Eta, DCAGlobal)) continue;
            // if((particle_option1 == 1) && (Eta > EtaCut || Eta < -EtaCut)) continue;
            if((Eta > EtaCut || Eta < -EtaCut)) continue;

            if(debug_1) std::cout << "debug_1 5" << endl;

            FillPhiPOI();       //Charge is needed here

            if(debug_1) std::cout << "debug_1 6" << endl;

            ShiftPhiPOI(trki);
        }

        //Store the flattened phi for POI (Protons)
        if(debug_1 && (n_particle2 != 0)) std::cout << "n_particle2 1 = " << n_particle2 << endl;

        Phi_pnew.resize(n_particle2);
        for(int trki = 0; trki < n_particle2; trki++)
        {
            gv_gamma::trk_mom->SetXYZ(particle2->at(trki).px, particle2->at(trki).py, particle2->at(trki).pz);
            Eta2   = gv_gamma::trk_mom->Eta();
            Pt2   = gv_gamma::trk_mom->Pt();
            Charge2 = particle2->at(trki).Charge;
            //DCAGlobal = p_all[trk].DCAGlobal();
            Phi2   = gv_gamma::trk_mom->Phi();

            if(debug_1) std::cout << "Eta2 = " << Eta2 << endl;

            // if((Pt2 < 0.4) || (Pt2 > 2.0)) continue;

            //if(!IsGoodPOI(Pt, Eta, DCAGlobal)) continue;
            if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;
            TLorentzVector trk_lorentz(particle2->at(trki).px, particle2->at(trki).py, particle2->at(trki).pz, sqrt(pow(gv_gamma::trk_mom->Mag(), 2) + pow(particle2->at(trki).mass, 2)));
            rapidity2 = trk_lorentz.Rapidity();
            // if(fabs(rapidity2) > 1) continue;
            trk_lorentz.Clear();

            if(debug_1) std::cout << "debug_1 5" << endl;

            FillPhiPOI_p();       //Charge is needed here

            if(debug_1) std::cout << "debug_1 6" << endl;

            ShiftPhiPOI_p(trki);
        }
        //////////Real analysis begins here//////////////////////////////////////////////////////////////////////////////////////
        Fcount = 0;

        int n_lam_used = 0;
        int n_proton_used = 0;
        int n_gamma = 0;
        //loop for the real analysis
        for(int trki = 0; trki < n_particle1; trki++)
        {
            //if((particle_option1 == 0) && ((particle1->at(trki).mass < 1.113) || (particle1->at(trki).mass > 1.119))) continue;
            gv_gamma::trk_mom->SetXYZ(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz);
            Eta   = gv_gamma::trk_mom->Eta();
            Pt    = gv_gamma::trk_mom->Pt();
            Charge = particle1->at(trki).Charge; //1 for lambda, -1 for antilambdas
            cout << "Lambda charge = " << Charge << endl;
            //DCAGlobal = p_all[trk].DCAGlobal();
            Phi   = gv_gamma::trk_mom->Phi();
            Theta     = 2.*atan(exp(-Eta));
            //ndEdx         = picoTrack->nHitsDedx();
            //eff = (cen > 0) ? PP0[cen - 1] * exp(-pow(PP1[cen - 1] / Pt, PP2[cen - 1])) : PP0[0] * exp(-pow(PP1[0] / Pt, PP2[0]));
            eff = 1; //Efficiency information later
            float eff_tof = 1;
            //                      if(TOF_eff && TOF_eff->GetEntries()) eff_tof = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt));
            eff *= eff_tof;

            // if((Pt < 0.5) || (Pt > 2.0)) continue;

            TLorentzVector trk_lorentz(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz, sqrt(pow(gv_gamma::trk_mom->Mag(), 2) + pow(particle1->at(trki).mass, 2)));
            rapidity = trk_lorentz.Rapidity();
        
            //if(DCAGlobal < DcaCut)
            //{
            hEtaPtDist->Fill(rapidity, Pt, Eweight);
            hEtaPhiDist->Fill(Eta, Phi, Eweight);
            hPhiPtDist->Fill(Phi, Pt, Eweight);
            //}
            //if(!IsGoodPOI(Pt, Eta, DCAGlobal)) continue;
            
            if(Eta > EtaCut || Eta < -EtaCut) continue;
            // if(fabs(rapidity) > 1) continue;
            trk_lorentz.Clear();

            //if((particle_option1 == 1) && (Eta > EtaCut || Eta < -EtaCut)) continue;
            //if(DCAGlobal < 1 && fabs(Eta) < 0.9)        // && Pt>0.2 && Pt*cosh(Eta)<2) {
            //{
            Hist_Pt->Fill(Pt, Eweight);
            /*if(picoTrack->isTofTrack() && dst->btofPidTraits(picoTrack->bTofPidTraitsIndex())->btof() > 0)*/ Hist_Pt_TOF->Fill(Pt, Eweight);
            //}

            //                  if(!IsGoodPion(dst,picoTrack,1)) continue;
            v2_sub = (Eta > 0) ? cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_new) * 100 : cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_new) * 100;
            FillCMW();

            float mQx_i = mQx, mQy_i = mQy;

            if(particle_option1 == 0)
            {
                //Take out protons/pions if they are included in the event plane reconstruction
                for(int ind = 0; ind < NPTracks; ind++)
                {
                    if(p_all->at(ind).trk_id == p_dau->at(2 * trki).trk_id)
                    {
                        TVector3 trk_mom_proton(p_dau->at(2 * trki).px, p_dau->at(2 * trki).py, p_dau->at(2 * trki).pz);
                        if(IsGoodAsso(trk_mom_proton.Pt(), trk_mom_proton.Eta(), p_dau->at(2 * trki).dcaglobal))
                        {
                            float Pt_proton = trk_mom_proton.Pt();
                            mQx_i -= Pt_proton * cos(PhiAsso_new[ind] * nHar);
                            mQy_i -= Pt_proton * sin(PhiAsso_new[ind] * nHar);
                        }
                    }
                    if(p_all->at(ind).trk_id == p_dau->at(2 * trki + 1).trk_id)
                    {
                        TVector3 trk_mom_pion(p_dau->at(2 * trki + 1).px, p_dau->at(2 * trki + 1).py, p_dau->at(2 * trki + 1).pz);
                        if(IsGoodAsso(trk_mom_pion.Pt(), trk_mom_pion.Eta(), p_dau->at(2 * trki + 1).dcaglobal))
                        {
                            float Pt_pion = trk_mom_pion.Pt();
                            mQx_i -= Pt_pion * cos(PhiAsso_new[ind] * nHar);
                            mQy_i -= Pt_pion * sin(PhiAsso_new[ind] * nHar);
                        }
                    }
                }
            }
            else if(particle_option1 == 1)
            {
                for(int ind2 = 0; ind2 < NPTracks; ind2++)
                {
                    if(p_all->at(ind2).trk_id == particle1->at(trki).trk_id)
                    {
                        cout << "Hit an EP proton!" << endl;
                        if(IsGoodAsso(Pt, Eta, DCAGlobal))
                        {
                            // mQx_i -= Pt * cos(Phi_new[ind2] * nHar);
                            // mQy_i -= Pt * sin(Phi_new[ind2] * nHar);
                            mQx_i -= Pt * cos(Phi_new[trki] * nHar);
                            mQy_i -= Pt * sin(Phi_new[trki] * nHar);
                        }
                    }
                }
            }

            TVector2 mQ_i(mQx_i, mQy_i);
            float psi_F = mQ_i.Phi() / nHar;
            Hist_TPC_EP_full_m1->Fill(psi_F);
            float psi_F_new = psi_F;
            for(int jj = 0; jj < order; jj++)
                psi_F_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F) / nHar / (jj + 1);
            Hist_TPC_EP_full_m1_flat->Fill(psi_F_new);

            v2  = cos(nHar * Phi_new[trki] - nHar * psi_F_new) * 100;
            v2e = cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_new) * 100;
            v2w = cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_new) * 100;
            /*Hist_v2_pt_obs1->Fill(Pt, v2);
            if(net_Nch_Asym > -99) Hist_v2_Ach->Fill(net_Nch_Asym, v2);
            //Hist_v2_pt_obs2->Fill(Pt, v2, Eweight / eff);
            Hist_v2_pt_obs2->Fill(Pt, v2, Eweight);
            //Hist_v2_eta_obs1->Fill(Eta, v2, 1. / eff);
            //Hist_v2_eta_obs2->Fill(Eta, v2, Eweight / eff);
            Hist_v2_eta_obs1->Fill(Eta, v2, 1.);
            Hist_v2_eta_obs2->Fill(Eta, v2, Eweight);
            lam_pt[Centrality]->Fill(Pt);*/

            float v2_final = 0;

            if(Eta > 0) v2_final = v2w;
            else if(Eta < 0) v2_final = v2e;

            if(particle_option1 == 0)
            {
                Hist_v2_pt_obs1->Fill(Pt, v2_final);
                if(net_Nch_Asym > -99) Hist_v2_Ach->Fill(net_Nch_Asym, v2_final);
                //Hist_v2_pt_obs2->Fill(Pt, v2, Eweight / eff);
                Hist_v2_pt_obs2->Fill(Pt, v2_final, Eweight);
            }
            else if(particle_option1 == 1)
            {
                if(Charge == 1)
                {
                    Hist_v2_pt_obs1->Fill(Pt, v2_final);
                    if(net_Nch_Asym > -99) Hist_v2_Ach->Fill(net_Nch_Asym, v2_final);
                    //Hist_v2_pt_obs2->Fill(Pt, v2, Eweight / eff);
                    Hist_v2_pt_obs2->Fill(Pt, v2_final, Eweight);
                }
            }
            //Hist_v2_eta_obs1->Fill(Eta, v2, 1. / eff);
            //Hist_v2_eta_obs2->Fill(Eta, v2, Eweight / eff);
            Hist_v2_eta_obs1->Fill(Eta, v2_final, 1.);
            Hist_v2_eta_obs2->Fill(Eta, v2_final, Eweight);
            Hist_v2_pt_obs2_caysm[net_charge_asym_bin]->Fill(Pt, v2_final, Eweight);
            //cout << "net_charge_asym_bin = " << net_charge_asym_bin << endl;

            /*if(fabs(Eta) < 0.5)
            {
                pTemp_v2->Fill(1, v2e, 1. / eff);
                pTemp_v2->Fill(2, v2w, 1. / eff);
                if(Charge > 0)
                {
                    pTemp_v2->Fill(3, v2e, 1. / eff);
                    pTemp_v2->Fill(4, v2w, 1. / eff);
                }
                if(Charge < 0)
                {
                    pTemp_v2->Fill(5, v2e, 1. / eff);
                    pTemp_v2->Fill(6, v2w, 1. / eff);
                }
            }*/

            Fcount++;
            n_lam_used++;
            if(opt_weight == 1) continue;
            for(int trkj = 0; trkj < n_particle2; trkj++)
            {
                Charge2 = particle2->at(trkj).Charge; //1 for protons, -1 for antiprotons
                cout << "Proton charge = " << Charge2 << endl;
                gv_gamma::trk_mom2->SetXYZ(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz);
                Eta2   = gv_gamma::trk_mom2->Eta();
                Pt2    = gv_gamma::trk_mom2->Pt();
                DCAGlobal2 = particle2->at(trkj).dcaglobal;
                Phi2   = gv_gamma::trk_mom2->Phi();

                // if((Pt2 < 0.4) || (Pt2 > 2.0)) continue;

                //if((fabs(p_proton->at(trkj).px - p_dau->at(2 * trki).px) < 1e-5) && (fabs(p_proton->at(trkj).py - p_dau->at(2 * trki).py) < 1e-5) && (fabs(p_proton->at(trkj).pz - p_dau->at(2 * trki).pz) < 1e-5 ) && (Charge2 == Charge)) {
                if((particle_option1 == 0) && (p_dau->at(2 * trki).trk_id == particle2->at(trkj).trk_id))
                {
                    cout << "Hit the Daughter Proton!" << endl;
                    continue;
                }

                //eff2 = (cen > 0) ? PP0[cen - 1] * exp(-pow(PP1[cen - 1] / Pt2, PP2[cen - 1])) : PP0[0] * exp(-pow(PP1[0] / Pt2, PP2[0]));
                eff2 = 1; //efficiency information acquired later
                float eff_tof2 = 1;
                //                                if(TOF_eff && TOF_eff->GetEntries()) eff_tof2 = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt2));
                eff2 *= eff_tof2;

                //if(!IsGoodPOI(Pt2, Eta2, DCAGlobal2)) continue;
                if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;
                TLorentzVector trk_lorentz2(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz, sqrt(pow(gv_gamma::trk_mom2->Mag(), 2) + pow(particle2->at(trkj).mass, 2)));
                rapidity2 = trk_lorentz2.Rapidity();
                // if(fabs(rapidity2) > 1) continue;
                trk_lorentz2.Clear();
                //                              if(!IsGoodPion(dst,picoTrack2,1)) continue;
                hDpt->Fill(fabs(Pt - Pt2), Eweight);

                float mQx_j = mQx_i, mQy_j = mQy_i;
                for(int ind2 = 0; ind2 < NPTracks; ind2++)
                {
                    if(p_all->at(ind2).trk_id == particle2->at(trkj).trk_id)
                    {
                        cout << "Hit an EP proton!" << endl;
                        if(IsGoodAsso(Pt2, Eta2, DCAGlobal2))
                        {
                            mQx_j -= Pt2 * cos(Phi_pnew[ind2] * nHar);
                            mQy_j -= Pt2 * sin(Phi_pnew[ind2] * nHar);
                            // mQx_j -= Pt2 * cos(Phi_pnew[trkj] * nHar);
                            // mQy_j -= Pt2 * sin(Phi_pnew[trkj] * nHar);
                        }
                    }
                }

                TVector2 mQ_j(mQx_j, mQy_j);
                float psi_F2 = mQ_j.Phi() / nHar;
                Hist_TPC_EP_full_m2->Fill(psi_F2);
                float psi_F2_new = psi_F2;
                for(int jj = 0; jj < order; jj++)
                    psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
                if(psi_F2_new < 0) psi_F2_new += PI;
                if(psi_F2_new > PI) psi_F2_new -= PI;
                Hist_TPC_EP_full_m2_flat->Fill(psi_F2_new);

                // if(trki == 0)
                // {
                //     n_proton_used++;

                //     if(Charge2 == Charge)
                //     {
                //         v2p  = cos(nHar * Phi_pnew[trkj] - nHar * psi_F2_new) * 100;
                //         /*Hist_v2_pt_obs1_p->Fill(Pt2, v2p);
                //         if(net_Nch_Asym > -99) Hist_v2_Ach_p->Fill(net_Nch_Asym, v2p);
                //         Hist_v2_pt_obs2_p->Fill(Pt2, v2p, Eweight);
                //         p_pt[Centrality]->Fill(Pt2);*/
                //         v2pe = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_for_new) * 100;
                //         v2pw = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_bac_new) * 100;

                //         float v2_final = 0;

                //         if(Eta2 > 0) v2_final = v2pw;
                //         else if(Eta2 < 0) v2_final = v2pe;

                //         Hist_v2_pt_obs1_p->Fill(Pt2, v2_final);
                //         if(net_Nch_Asym > -99) Hist_v2_Ach_p->Fill(net_Nch_Asym, v2_final);
                //         Hist_v2_pt_obs2_p->Fill(Pt2, v2_final, Eweight);

                //         Hist_v2_pt_obs2_p_caysm[net_charge_asym_bin]->Fill(Pt2, v2_final, Eweight);
                //     }
                // }

                //correlations
                //correlator0 = cos(Phi + Phi2 - 2 * psi_F2_new) * 100;
                correlator3 = cos(Phi_new[trki] - Phi_pnew[trkj]) * 100;
                correlator4e = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * TPC_EP_for_new) * 100;
                correlator4w = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * TPC_EP_bac_new) * 100;
                correlator4 = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * psi_F2_new) * 100;
                correlator0 = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + 2 * psi_F2_new) * 100;
                n_gamma++;

                cout << "Gamma112 = (correlator4) " << correlator4 << endl;
                cout << "Gamma112 = (correlator4e) " << correlator4e << endl;
                cout << "Gamma112 = (correlator4w) " << correlator4w << endl;

                if(Charge > 0 && Charge2 > 0) FillGamma(1);
                if(Charge < 0 && Charge2 < 0) FillGamma(2);
                if(Charge * Charge2 > 0)
                {
                    FillGamma(3);
                    cout << "filled 3, same sign" << endl;
                }
                if(Charge * Charge2 < 0)
                {
                    FillGamma(4);
                    cout << "filled 4, opposite sign" << endl;
                }

                if(fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
                {
                    Hist_v2_eta_obs3->Fill(Eta, v2, Eweight / eff);
                    if(fabs(Eta) < 0.5 && fabs(Eta2) < 0.5)
                    {
                        pTemp_v2_noHBT->Fill(1, v2e, 1. / eff);
                        pTemp_v2_noHBT->Fill(2, v2w, 1. / eff);
                        if(Charge > 0)
                        {
                            pTemp_v2_noHBT->Fill(3, v2e, 1. / eff);
                            pTemp_v2_noHBT->Fill(4, v2w, 1. / eff);
                        }
                        if(Charge < 0)
                        {
                            pTemp_v2_noHBT->Fill(5, v2e, 1. / eff);
                            pTemp_v2_noHBT->Fill(6, v2w, 1. / eff);
                        }
                    }
                }
            } // 2nd track

        }  //1st Track

        if(opt_weight != 1)
        {
            for(int trkj = 0; trkj < n_particle2; trkj++)
            {
                Charge2 = particle2->at(trkj).Charge; //1 for protons, -1 for antiprotons
                cout << "Proton charge = " << Charge2 << endl;
                gv_gamma::trk_mom2->SetXYZ(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz);
                Eta2   = gv_gamma::trk_mom2->Eta();
                Pt2    = gv_gamma::trk_mom2->Pt();
                DCAGlobal2 = particle2->at(trkj).dcaglobal;
                Phi2   = gv_gamma::trk_mom2->Phi();

                //eff2 = (cen > 0) ? PP0[cen - 1] * exp(-pow(PP1[cen - 1] / Pt2, PP2[cen - 1])) : PP0[0] * exp(-pow(PP1[0] / Pt2, PP2[0]));
                eff2 = 1; //efficiency information acquired later
                float eff_tof2 = 1;
                //                                if(TOF_eff && TOF_eff->GetEntries()) eff_tof2 = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt2));
                eff2 *= eff_tof2;

                if((Pt2 < 0.4) || (Pt2 > 2.0)) continue;
                TLorentzVector trk_lorentz2(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz, sqrt(pow(gv_gamma::trk_mom2->Mag(), 2) + pow(particle2->at(trkj).mass, 2)));
                rapidity2 = trk_lorentz2.Rapidity();
                if(fabs(rapidity2) > 1) continue;
                trk_lorentz2.Clear();

                //if(!IsGoodPOI(Pt2, Eta2, DCAGlobal2)) continue;
                //if(Eta2 > EtaCut || Eta2 < -EtaCut) continue;
                //                              if(!IsGoodPion(dst,picoTrack2,1)) continue;

                float mQx_j = mQx, mQy_j = mQy;
                for(int ind2 = 0; ind2 < NPTracks; ind2++)
                {
                    if(p_all->at(ind2).trk_id == particle2->at(trkj).trk_id)
                    {
                        cout << "Hit an EP proton!" << endl;
                        if(IsGoodAsso(Pt2, Eta2, DCAGlobal2))
                        {
                            // mQx_j -= Pt2 * cos(Phi_pnew[ind2] * nHar);
                            // mQy_j -= Pt2 * sin(Phi_pnew[ind2] * nHar);
                            mQx_j -= Pt2 * cos(Phi_pnew[trkj] * nHar);
                            mQy_j -= Pt2 * sin(Phi_pnew[trkj] * nHar);
                        }
                    }
                }

                TVector2 mQ_j(mQx_j, mQy_j);
                float psi_F2 = mQ_j.Phi() / nHar;
                Hist_TPC_EP_full_m2->Fill(psi_F2);
                float psi_F2_new = psi_F2;
                for(int jj = 0; jj < order; jj++)
                    psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
                if(psi_F2_new < 0) psi_F2_new += PI;
                if(psi_F2_new > PI) psi_F2_new -= PI;
                Hist_TPC_EP_full_m2_flat->Fill(psi_F2_new);

                n_proton_used++;

                if(((Charge2 == Charge) && (particle_option1 == 0)) || ((Charge2 == -1) && (particle_option1 == 1)))
                {
                    v2p  = cos(nHar * Phi_pnew[trkj] - nHar * psi_F2_new) * 100;
                    /*Hist_v2_pt_obs1_p->Fill(Pt2, v2p);
                    if(net_Nch_Asym > -99) Hist_v2_Ach_p->Fill(net_Nch_Asym, v2p);
                    Hist_v2_pt_obs2_p->Fill(Pt2, v2p, Eweight);
                    p_pt[Centrality]->Fill(Pt2);*/
                    v2pe = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_for_new) * 100;
                    v2pw = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_bac_new) * 100;

                    float v2_final = 0;

                    if(Eta2 > 0) v2_final = v2pw;
                    else if(Eta2 < 0) v2_final = v2pe;

                    Hist_v2_pt_obs1_p->Fill(Pt2, v2_final);
                    if(net_Nch_Asym > -99) Hist_v2_Ach_p->Fill(net_Nch_Asym, v2_final);
                    Hist_v2_pt_obs2_p->Fill(Pt2, v2_final, Eweight);

                    Hist_v2_pt_obs2_p_caysm[net_charge_asym_bin]->Fill(Pt2, v2_final, Eweight);
                }

            }
        }

        num_lam_final->Fill(n_lam_used);
        num_proton_final->Fill(n_proton_used);
        num_gamma_final->Fill(n_gamma);

        PhiAsso_new.clear();
        Phi_new.clear();
        Phi_pnew.clear();
        pTemp_v2->Reset();
        pTemp_v2_noHBT->Reset();
        pTemp_parity_e->Reset();
        pTemp_parity_e_noHBT->Reset();
        pTemp_parity_w->Reset();
        pTemp_parity_w_noHBT->Reset();
        pTemp_delta->Reset();
        pTemp_delta_noHBT->Reset();

    } // Event

    rc = (TH1D *)Hist_Pt_TOF->Clone();
    rc->SetName("rc");
    rc->Divide(Hist_Pt);
    WriteHistogram(cen, opt_weight, JobIDName);
    if(opt_weight == 1)
    {
        WriteWeight(fname_new);
    }
    return;
}

//////////////////////////////////
void set_debug()
{
    debug_1 = false;
    debug_2 = false;

    particle_option1 = 0; // 0 - Lambda, 1 - proton
    particle_option2 = 0; // 0 - proton
}
//////////////////////////////////
bool IsGoodEvent(int c)
{
    hTally->Fill(1);
    //hZDCcoin->Fill(ZDCcoin);
    //  if(ZDCcoin<4000) return false;
    //hTally->Fill(2);
    //hVertexXY->Fill(pVx, pVy);

    //if((pVx - 0.06) * (pVx - 0.06) + (pVy + 0.13) * (pVy + 0.13) > 1) return false;
    hTally->Fill(3);
    hVzDiff->Fill(pVz - VPDvz);
    if((pVz - VPDvz) > 4 || (pVz - VPDvz) < -4) return false;


    //hMult_Vz->Fill(RefMult, pVz);
    hVertexZ->Fill(pVz);
    hTally->Fill(4);
    //if(TMath::Abs(pVz) > Vz_cut) return false;
    hMult_Vz_new->Fill(RefMult, pVz, Eweight);
    Ref_TOF->Fill(RefMult, TOFMult);
    //if(TOFMult > (4.52*RefMult + 47.54)) return false;
    //if(TOFMult < (3.4375*RefMult - 31.09375)) return false;
    if(RefMult > (11.070446 + 1.384789 * TOFMult - 0.002402 * pow(TOFMult, 2) + 7.514934 * pow(10, -6)*pow(TOFMult, 3) - 9.148756 * pow(10, -9)*pow(TOFMult, 4))) return false;
    if(RefMult < (-12.764095 + 0.594566 * TOFMult + 0.001040 * pow(TOFMult, 2) - 2.006525 * pow(10, -6)*pow(TOFMult, 3) + 1.127352 * pow(10, -9)*pow(TOFMult, 4))) return false;
    Ref_Day3->Fill(Day3, RefMult);
    TOF_Day3->Fill(Day3, TOFMult);
    NPT_Day3->Fill(Day3, NPTracks);

    //hTally->Fill(5);
    /*int Bad = 0;
    for(int jj = 0; jj < Nrun_MB; jj++) if(Day3 == bad_day3_MB[jj])
        {
            Bad = 1;
            break;
        }*/
    //hTally->Fill(6);
    //Centrality = details->cent;
    //for(int j = 0; j < 9; j++) if(RefMult > cenDef[j]) Centrality = j + 1;
    //if(RefMult < cenDef[0]) return false;
    hTally->Fill(7);

    hCentrality->Fill(Centrality);//, Eweight);

    hMult_Vz->Fill(RefMult, TOFMult);
    //if(c && Centrality != c && Centrality != (c + 1) && Centrality != (c - 1) && Centrality != (c + 2) && Centrality != (c - 2) ) return false;
    if((Centrality != c)) return false;

    hTally->Fill(8);

    return true;
}
//////////////////////////////////////////////
bool CountCharge()
{
    int Ntof = 0, Npos = 0, Nneg = 0;
    Fcount = 0;
    for(int trk = 0; trk < NPTracks; trk++)
    {
        //StPicoTrack *picoTrack = d->track(trk);
        //if(p_all->at(trk).prim != 1) continue;
        gv_gamma::trk_mom_temp->SetXYZ(p_all->at(trk).px, p_all->at(trk).py, p_all->at(trk).pz);
        EtaAsso   = gv_gamma::trk_mom_temp->Eta();
        PtAsso    = gv_gamma::trk_mom_temp->Pt();
        ChargeAsso = p_all->at(trk).Charge;
        DCAGlobalAsso = p_all->at(trk).dcaglobal;
        nSigma_p  = p_all->at(trk).nSigmaProton;

        if(ChargeAsso > 0 && fabs(EtaAsso) < 1 && PtAsso > 0.15 && DCAGlobalAsso < 1 && !(fabs(nSigma_p) < 3 && PtAsso < 0.4)) Npos++;
        if(ChargeAsso < 0 && fabs(EtaAsso) < 1 && PtAsso > 0.15 && DCAGlobalAsso < 1 && !(fabs(nSigma_p) < 3 && PtAsso < 0.4)) Nneg++;
        if(p_all->at(trk).isTofTrack == 1) Ntof++;
        Hist_DCA->Fill(DCAGlobalAsso);
        if(!IsGoodAsso(PtAsso, EtaAsso, DCAGlobalAsso)) continue;
        Fcount++;
    }
    if(Ntof < 2) return false; //at least 2 tracks match TOF
    NPA_Day3->Fill(Day3, Fcount);
    hTally->Fill(10);
    int net_Nch = Npos - Nneg;
    net_Nch_Asym = -99;
    if((Npos + Nneg) > 0) net_Nch_Asym = (Npos - Nneg) / float(Npos + Nneg);
    Hist_positive->Fill(Npos);
    Hist_negative->Fill(Nneg);
    Hist_Ch->Fill(Npos + Nneg);
    Hist_netCh->Fill(net_Nch);
    if (net_Nch_Asym > -99)
    {
        Hist_netChAsym->Fill(net_Nch_Asym);
        p_netChAsym_RefMult->Fill(net_Nch_Asym, RefMult);
        if(net_Nch_Asym < (MeanNetChargeAsym - 0.8 * StdDevNetChargeAsym))
        {
            Hist_netChAsym_bin->Fill(1, net_Nch_Asym);
            net_charge_asym_bin = 0;
        }
        else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3 * StdDevNetChargeAsym)))
        {
            Hist_netChAsym_bin->Fill(2, net_Nch_Asym);
            net_charge_asym_bin = 1;
        }
        else if(net_Nch_Asym < (MeanNetChargeAsym + (0.2 * StdDevNetChargeAsym)))
        {
            Hist_netChAsym_bin->Fill(3, net_Nch_Asym);
            net_charge_asym_bin = 2;
        }
        else if(net_Nch_Asym < (MeanNetChargeAsym + 0.8 * StdDevNetChargeAsym))
        {
            Hist_netChAsym_bin->Fill(4, net_Nch_Asym);
            net_charge_asym_bin = 3;
        }
        else
        {
            Hist_netChAsym_bin->Fill(5, net_Nch_Asym);
            net_charge_asym_bin = 4;
        }

    }
    return true;
}
/////////////////////////////////////////////
void MakeTPC_EP()
{
    TVector2 mQ, mQ1, mQ2, mQ3, mQ4;
    mQx = 0., mQy = 0.;
    float mQQx = 0., mQQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0., mQx3 = 0., mQy3 = 0., mQx4 = 0., mQy4 = 0.;
    Fcount = 0;
    int Qcount = 0;
    PhiAsso_new.resize(NPTracks);
    for(int trk = 0; trk < NPTracks; trk++)
    {
        if(p_all->at(trk).prim != 1) continue;
        gv_gamma::trk_mom_temp->SetXYZ(p_all->at(trk).px, p_all->at(trk).py, p_all->at(trk).pz);
        EtaAsso   = gv_gamma::trk_mom_temp->Eta();
        PtAsso    = gv_gamma::trk_mom_temp->Pt();
        ChargeAsso = p_all->at(trk).Charge;
        DCAGlobalAsso = p_all->at(trk).dcaglobal;
        PhiAsso   = gv_gamma::trk_mom_temp->Phi();

        if(!IsGoodAsso(PtAsso, EtaAsso, DCAGlobalAsso)) continue;

        //if((PtAsso >= 0.4) && (p_all[trk].TOFflag >= 1) && (DCAGlobalAsso < 1) && (fabs(p_all[trk].nSigmaProton) <= 2) && (fabs(ChargeAsso) == 1)) continue;

        FillPhiAsso();  //ChargeAsso is needed here
        ShiftPhiAsso(trk);

        mQx += PtAsso * cos(PhiAsso_new[trk] * nHar);
        mQy += PtAsso * sin(PhiAsso_new[trk] * nHar);
        if(fabs(EtaAsso) < 0.5)
        {
            mQQx += cos(PhiAsso_new[trk] * nHar);
            mQQy += sin(PhiAsso_new[trk] * nHar);
            Qcount++;
        }
        if(gv_gamma::iTrack[Fcount] > Scount)
        {
            mQx1 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy1 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        else
        {
            mQx2 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy2 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        if(EtaAsso > Eta_EP_Cut)
        {
            mQx3 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy3 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        if(EtaAsso < -Eta_EP_Cut)
        {
            mQx4 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy4 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        Fcount++;
    }
    mQ.Set(mQx, mQy);
    mQ1.Set(mQx1, mQy1);
    mQ2.Set(mQx2, mQy2);
    mQ3.Set(mQx3, mQy3);
    mQ4.Set(mQx4, mQy4);
    TPC_EP_full = mQ.Phi() / nHar;
    TPC_EP_east = mQ1.Phi() / nHar;
    TPC_EP_west = mQ2.Phi() / nHar;
    TPC_EP_for  = mQ3.Phi() / nHar;
    TPC_EP_bac  = mQ4.Phi() / nHar;
    Q2 = (mQQx * mQQx + mQQy * mQQy) / float(Qcount);
    TPC_Day3_cos2->Fill(Day3, cos(2 * TPC_EP_full));
    TPC_Day3_sin2->Fill(Day3, sin(2 * TPC_EP_full));
}

//////////////////////////////////////////////
void FillEP_resolution()
{
    float cos_ew = cos(nHar * TPC_EP_east_new - nHar * TPC_EP_west_new);
    float cos_fb = cos(nHar * TPC_EP_for_new - nHar * TPC_EP_bac_new);
    Hist_cos->Fill(1, cos_fb, Eweight);
    Hist_cos->Fill(2, cos_ew, Eweight);
    Hist_cos->Fill(3, cos(nHar * TPC_EP_for_new - nHar * TPC_EP_full_new));
    Hist_cos->Fill(4, cos(nHar * TPC_EP_bac_new - nHar * TPC_EP_full_new));
    Hist_Q2->Fill(Q2);
    p_cos_Q2->Fill(Q2, cos_fb, Eweight);
    p_RefMult_Q2->Fill(Q2, Fcount);
    Hist_RefMult_Q2->Fill(Q2, Fcount);

    if(net_Nch_Asym > -99)
    {
        p_netChAsym_cos->Fill(net_Nch_Asym, cos_fb);
        if(net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym)) Hist_cos_Ach->Fill(1, cos_fb, Eweight);
        else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym))) Hist_cos_Ach->Fill(2, cos_fb, Eweight);
        else if(net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym))) Hist_cos_Ach->Fill(3, cos_fb, Eweight);
        else if(net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym)) Hist_cos_Ach->Fill(4, cos_fb, Eweight);
        else Hist_cos_Ach->Fill(5, cos_fb, Eweight);
    }
}
/////////////////////////////////////////////////
bool IsGoodAsso(float p, float e, float d)
{
    if(p > pt_asso_up || p < pt_asso_lo) return false;
    if(d > DcaCut) return false;
    if(e > EtaCut || e < -EtaCut) return false;
    return true;
}
/////////////////////////////////////////////////////
bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt)
{
    //hTall->Fill(1);
    if(p->gDCA(pV).Mag() > 1)   return false;
    //hTall->Fill(2);
    float eta_i = p->pMom().Eta();
    if(fabs(eta_i) > 0.9) return false;
    //hTall->Fill(3);
    float pt_i =  p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if(pt_i < 0.2 || p_i > 1.6) return false;
    //hTall->Fill(4);
    float nSig_i = p->nSigmaPion();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2) return false;
    //hTall->Fill(5);
    if(opt == 0) return true;
    if(!(p->isTofTrack())) return false;
    //hTall->Fill(6);
    StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
    if(!trait) return false;
    //hTall->Fill(7);
    if(trait->btof() <= 0) return false;
    //hTall->Fill(8);
    if(fabs(trait->btofYLocal()) > 1.8) return false;
    //hTall->Fill(9);
    float beta_i = trait->btofBeta();
    if(beta_i == 0) return false;
    //hTall->Fill(10);
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if(mass2_i < -0.01 || mass2_i > 0.1) return false;
    //hTall->Fill(11);
    return true;
}
///////////////////////////////////////////////////
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt)
{
    if(p->gDCA(pV).Mag() > 1) return false;
    float eta_i = p->pMom().Eta();
    if(fabs(eta_i) > 0.9) return false;
    float pt_i =  p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if(pt_i < 0.2 || p_i > 1.6) return false;

    float nSig_i = p->nSigmaKaon();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2) return false;

    if(opt == 0) return true;
    if(!(p->isTofTrack())) return false;
    StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
    if(!trait) return false;
    if(trait->btof() <= 0) return false;
    if(fabs(trait->btofYLocal()) > 1.8) return false;
    float beta_i = trait->btofBeta();
    if(beta_i == 0) return false;
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if(mass2_i < 0.2 || mass2_i > 0.35) return false;
    return true;
}
////////////////////////////////////////////////////
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt)
{
    if(p->gDCA(pV).Mag() > 1) return false;
    float eta_i = p->pMom().Eta();
    if(fabs(eta_i) > 0.9) return false;
    float pt_i =  p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if(pt_i < 0.4 || p_i > 2) return false;

    float nSig_i = p->nSigmaProton();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2) return false;

    if(opt == 0) return true;
    if(!(p->isTofTrack())) return false;
    StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
    if(!trait) return false;
    if(trait->btof() <= 0) return false;
    if(fabs(trait->btofYLocal()) > 1.8) return false;
    float beta_i = trait->btofBeta();
    if(beta_i == 0) return false;
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if(mass2_i < 0.8 || mass2_i > 1) return false;
    return true;
}
////////////////////////////////////
void ShiftPsi()
{
    Hist_TPC_EP_full->Fill(TPC_EP_full, Day);
    Hist_TPC_EP_east->Fill(TPC_EP_east, Day);
    Hist_TPC_EP_west->Fill(TPC_EP_west, Day);
    Hist_TPC_EP_for->Fill(TPC_EP_for, Day);
    Hist_TPC_EP_bac->Fill(TPC_EP_bac, Day);

    for(int kk = 0; kk < order; kk++)
    {
        pTPC_EP_full->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_full), Eweight);
        pTPC_EP_full->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_full), Eweight);
        pTPC_EP_east->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_east), Eweight);
        pTPC_EP_east->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_east), Eweight);
        pTPC_EP_west->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_west), Eweight);
        pTPC_EP_west->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_west), Eweight);
        pTPC_EP_for->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_for), Eweight);
        pTPC_EP_for->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_for), Eweight);
        pTPC_EP_bac->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_bac), Eweight);
        pTPC_EP_bac->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_bac), Eweight);
    }

    if(Weight_Read && Read_TPC_EP_full->GetEntries())
    {
        for(int k = 0; k < 2 * order; k++)
        {
            PsiMean_F[k] = Read_TPC_EP_full->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_E[k] = Read_TPC_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_W[k] = Read_TPC_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_f[k] = Read_TPC_EP_for->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_b[k] = Read_TPC_EP_bac->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
        }
    }

    TPC_EP_full_new = TPC_EP_full, TPC_EP_east_new = TPC_EP_east, TPC_EP_west_new = TPC_EP_west;
    TPC_EP_for_new = TPC_EP_for, TPC_EP_bac_new = TPC_EP_bac;

    for(int jj = 0; jj < order; jj++)
    {
        TPC_EP_full_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_full) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_full) / nHar / (jj + 1);
        TPC_EP_east_new += -2 * PsiMean_E[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_east) / nHar / (jj + 1) + 2 * PsiMean_E[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_east) / nHar / (jj + 1);
        TPC_EP_west_new += -2 * PsiMean_W[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_west) / nHar / (jj + 1) + 2 * PsiMean_W[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_west) / nHar / (jj + 1);
        TPC_EP_for_new  += -2 * PsiMean_f[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_for) / nHar / (jj + 1)  + 2 * PsiMean_f[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_for) / nHar / (jj + 1);
        TPC_EP_bac_new  += -2 * PsiMean_b[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_bac) / nHar / (jj + 1)  + 2 * PsiMean_b[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_bac) / nHar / (jj + 1);
    }
    if(TPC_EP_full_new > PI) TPC_EP_full_new -= PI;
    if(TPC_EP_full_new < 0) TPC_EP_full_new += PI;
    if(TPC_EP_east_new > PI) TPC_EP_east_new -= PI;
    if(TPC_EP_east_new < 0) TPC_EP_east_new += PI;
    if(TPC_EP_west_new > PI) TPC_EP_west_new -= PI;
    if(TPC_EP_west_new < 0) TPC_EP_west_new += PI;
    if(TPC_EP_for_new > PI)  TPC_EP_for_new  -= PI;
    if(TPC_EP_for_new < 0)  TPC_EP_for_new  += PI;
    if(TPC_EP_bac_new > PI)  TPC_EP_bac_new  -= PI;
    if(TPC_EP_bac_new < 0)  TPC_EP_bac_new  += PI;

    Hist_TPC_EP_east_flat->Fill(TPC_EP_east_new, Day);
    Hist_TPC_EP_west_flat->Fill(TPC_EP_west_new, Day);
    Hist_TPC_EP_for_flat->Fill(TPC_EP_for_new, Day);
    Hist_TPC_EP_bac_flat->Fill(TPC_EP_bac_new, Day);
    Hist_TPC_EP_full_flat->Fill(TPC_EP_full_new, Day);
}
//////////////////////////////////
void ShiftPhiPOI(int tr)
{
    int index = 0, index2 = 0;
    if(pVz > 0) index2 = (Charge > 0) ? 1 : 2;
    else index2 = (Charge > 0) ? 3 : 4;
    if(Weight_Read && (TPCmean_FF_1->GetEntries() || TPCmean_RF_1->GetEntries()))
    {
        for(int kk = 0; kk < order; kk++)
        {
            if(pVz > 0) index = (Charge > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
            else      index = (Charge > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if(Eta > 0)
            {
                if(Pt < 0.5)
                {
                    PhiMean_cos[kk] = TPCmean_FF_1->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_FF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if(Pt < 1)
                {
                    PhiMean_cos[kk] = TPCmean_FF_2->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_FF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos[kk] = TPCmean_FF_3->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_FF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if(Pt < 0.5)
                {
                    PhiMean_cos[kk] = TPCmean_RF_1->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_RF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if(Pt < 1)
                {
                    PhiMean_cos[kk] = TPCmean_RF_2->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_RF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos[kk] = TPCmean_RF_3->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_RF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    Phi_new[tr] = Phi;   //store the shifted angles
    for(int jj = 0; jj < order; jj++) Phi_new[tr] += -2 * PhiMean_sin[jj] * cos(jj * Phi + Phi) / (jj + 1) + 2 * PhiMean_cos[jj] * sin(jj * Phi + Phi) / (jj + 1);
    if(Phi_new[tr] > PI) Phi_new[tr] -= 2 * PI;
    if(Phi_new[tr] < -PI) Phi_new[tr] += 2 * PI;
    if(Eta > 0)
    {
        if(Pt < 0.5)      Hist_Phi_FF_new_1->Fill(Phi_new[tr], index2, Eweight);
        else if(Pt < 1)   Hist_Phi_FF_new_2->Fill(Phi_new[tr], index2, Eweight);
        else            Hist_Phi_FF_new_3->Fill(Phi_new[tr], index2, Eweight);
    }
    else
    {
        if(Pt < 0.5)      Hist_Phi_RF_new_1->Fill(Phi_new[tr], index2, Eweight);
        else if(Pt < 1)   Hist_Phi_RF_new_2->Fill(Phi_new[tr], index2, Eweight);
        else            Hist_Phi_RF_new_3->Fill(Phi_new[tr], index2, Eweight);
    }
}
//////////////////////////////////////////
//////////////////////////////////
void ShiftPhiPOI_p(int tr)
{
    int index = 0, index2 = 0;
    if(pVz > 0) index2 = (Charge2 > 0) ? 1 : 2;
    else index2 = (Charge2 > 0) ? 3 : 4;
    if(Weight_Read && (TPCmean_FF_1_p->GetEntries() || TPCmean_RF_1_p->GetEntries()))
    {
        for(int kk = 0; kk < order; kk++)
        {
            if(pVz > 0) index = (Charge2 > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
            else      index = (Charge2 > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if(Eta2 > 0)
            {
                if(Pt2 < 0.5)
                {
                    PhiMean_cos_p[kk] = TPCmean_FF_1_p->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_FF_1_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if(Pt2 < 1)
                {
                    PhiMean_cos_p[kk] = TPCmean_FF_2_p->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_FF_2_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos_p[kk] = TPCmean_FF_3_p->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_FF_3_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if(Pt2 < 0.5)
                {
                    PhiMean_cos_p[kk] = TPCmean_RF_1_p->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_RF_1_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if(Pt2 < 1)
                {
                    PhiMean_cos_p[kk] = TPCmean_RF_2_p->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_RF_2_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos_p[kk] = TPCmean_RF_3_p->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_RF_3_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    Phi_pnew[tr] = Phi2;   //store the shifted angles
    for(int jj = 0; jj < order; jj++) Phi_pnew[tr] += -2 * PhiMean_sin_p[jj] * cos(jj * Phi2 + Phi2) / (jj + 1) + 2 * PhiMean_cos_p[jj] * sin(jj * Phi2 + Phi2) / (jj + 1);
    if(Phi_pnew[tr] > PI) Phi_pnew[tr] -= 2 * PI;
    if(Phi_pnew[tr] < -PI) Phi_pnew[tr] += 2 * PI;
    if(Eta2 > 0)
    {
        if(Pt2 < 0.5)      Hist_Phi_FF_new_1_p->Fill(Phi_pnew[tr], index2, Eweight);
        else if(Pt2 < 1)   Hist_Phi_FF_new_2_p->Fill(Phi_pnew[tr], index2, Eweight);
        else            Hist_Phi_FF_new_3_p->Fill(Phi_pnew[tr], index2, Eweight);
    }
    else
    {
        if(Pt2 < 0.5)      Hist_Phi_RF_new_1_p->Fill(Phi_pnew[tr], index2, Eweight);
        else if(Pt2 < 1)   Hist_Phi_RF_new_2_p->Fill(Phi_pnew[tr], index2, Eweight);
        else            Hist_Phi_RF_new_3_p->Fill(Phi_pnew[tr], index2, Eweight);
    }
}
//////////////////////////////////////////
void ShiftPhiAsso(int tr)
{
    int index = 0;
    if(Weight_Read && (TPCmeanAsso_FF_1->GetEntries() || TPCmeanAsso_RF_1->GetEntries()))
    {
        for(int kk = 0; kk < order; kk++)
        {
            if(pVz > 0) index = (ChargeAsso > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
            else      index = (ChargeAsso > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if(EtaAsso > 0)
            {
                if(PtAsso < 0.5)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_FF_1->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_FF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if(PtAsso < 1)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_FF_2->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_FF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_FF_3->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_FF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if(PtAsso < 0.5)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_RF_1->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_RF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if(PtAsso < 1)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_RF_2->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_RF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_RF_3->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_RF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    PhiAsso_new[tr] = PhiAsso;
    for(int jj = 0; jj < order; jj++) PhiAsso_new[tr] += -2 * PhiMeanAsso_sin[jj] * cos(jj * PhiAsso + PhiAsso) / (jj + 1) + 2 * PhiMeanAsso_cos[jj] * sin(jj * PhiAsso + PhiAsso) / (jj + 1);
}
/////////////////////////////////////
void FillPhiAsso()    // shift parameters for Particles of EP
{
    int index = 0;
    for(int kk = 0; kk < order; kk++)
    {
        if(pVz > 0) index = (ChargeAsso > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
        else      index = (ChargeAsso > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if(EtaAsso > 0)
        {
            if(PtAsso < 0.5)
            {
                pTPCmeanPhiAsso_FF_1->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_FF_1->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else if(PtAsso < 1)
            {
                pTPCmeanPhiAsso_FF_2->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_FF_2->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else
            {
                pTPCmeanPhiAsso_FF_3->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_FF_3->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
        }
        if(EtaAsso < 0)
        {
            if(PtAsso < 0.5)
            {
                pTPCmeanPhiAsso_RF_1->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_RF_1->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else if(PtAsso < 1)
            {
                pTPCmeanPhiAsso_RF_2->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_RF_2->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else
            {
                pTPCmeanPhiAsso_RF_3->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_RF_3->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
        }
    }
}
///////////////////////////////////
void FillPhiPOI()   // shift parameters for Particles of Interest
{
    int index = 0, index2 = 0;
    if(pVz > 0) index2 = (Charge > 0) ? 1 : 2;
    else index2 = (Charge > 0) ? 3 : 4;
    for(int kk = 0; kk < order; kk++)
    {
        if(pVz > 0) index = (Charge > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
        else      index = (Charge > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if(Eta > 0)
        {
            if(Pt < 0.5)
            {
                pTPCmeanPhi_FF_1->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_FF_1->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if(kk == 0) Hist_Phi_FF_1->Fill(Phi, index2, Eweight);
            }
            else if(Pt < 1)
            {
                pTPCmeanPhi_FF_2->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_FF_2->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if(kk == 0) Hist_Phi_FF_2->Fill(Phi, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_FF_3->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_FF_3->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if(kk == 0) Hist_Phi_FF_3->Fill(Phi, index2, Eweight);
            }
        }
        if(Eta < 0)
        {
            if(Pt < 0.5)
            {
                pTPCmeanPhi_RF_1->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_RF_1->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if(kk == 0) Hist_Phi_RF_1->Fill(Phi, index2, Eweight);
            }
            else if(Pt < 1)
            {
                pTPCmeanPhi_RF_2->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_RF_2->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if(kk == 0) Hist_Phi_RF_2->Fill(Phi, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_RF_3->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_RF_3->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if(kk == 0) Hist_Phi_RF_3->Fill(Phi, index2, Eweight);
            }
        }
    }
}
/////////////////////////////////
///////////////////////////////////
void FillPhiPOI_p()   // shift parameters for Particles of Interest (Protons)
{
    int index = 0, index2 = 0;
    if(pVz > 0) index2 = (Charge2 > 0) ? 1 : 2;
    else index2 = (Charge2 > 0) ? 3 : 4;
    for(int kk = 0; kk < order; kk++)
    {
        if(pVz > 0) index = (Charge2 > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
        else      index = (Charge2 > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if(Eta2 > 0)
        {
            if(Pt2 < 0.5)
            {
                pTPCmeanPhi_FF_1_p->Fill(index,  Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_FF_1_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if(kk == 0) Hist_Phi_FF_1_p->Fill(Phi2, index2, Eweight);
            }
            else if(Pt2 < 1)
            {
                pTPCmeanPhi_FF_2_p->Fill(index,  Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_FF_2_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if(kk == 0) Hist_Phi_FF_2_p->Fill(Phi2, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_FF_3_p->Fill(index,  Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_FF_3_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if(kk == 0) Hist_Phi_FF_3_p->Fill(Phi2, index2, Eweight);
            }
        }
        if(Eta2 < 0)
        {
            if(Pt2 < 0.5)
            {
                pTPCmeanPhi_RF_1_p->Fill(index,  Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_RF_1_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if(kk == 0) Hist_Phi_RF_1_p->Fill(Phi2, index2, Eweight);
            }
            else if(Pt2 < 1)
            {
                pTPCmeanPhi_RF_2_p->Fill(index,  Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_RF_2_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if(kk == 0) Hist_Phi_RF_2_p->Fill(Phi2, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_RF_3_p->Fill(index,  Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_RF_3_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if(kk == 0) Hist_Phi_RF_3_p->Fill(Phi2, index2, Eweight);
            }
        }
    }
}
/////////////////////////////////
void WriteHistogram(int c, int o, TString JobIDName)
{
    char fname_out[200];
    TString Name2 = "sched";
    Name2.Append(JobIDName);
    if(o != 1) sprintf(fname_out, "cen%d.gamma112_fullEP_eff_pT02_module.root", c);
    if(o == 1) sprintf(fname_out, "cen%d.v2_fullEP_eff_pT02_module.root", c);
    Name2.Append(fname_out);
    TFile *fout = new TFile(Name2, "RECREATE");

    //        gROOT->GetList()->ls();
    TList *list = gROOT->GetList();//GetListOfKeys();
    TIter next(list);
    TKey *key;
    TObject *obj;
    while ((key = (TKey *)next()))
    {
        TString tempStr(key->GetName());
        if (tempStr.Contains("Temp")) continue;
        if(o == 1 && (tempStr.Contains("Parity") || tempStr.Contains("Delta"))) continue;
        obj = gROOT->Get(key->GetName());
        if (!obj) continue;
        if(obj->IsA() == TDirectory::Class())
        {
            delete obj;
            obj = NULL;
            continue;
        }
        obj->Write();
    }

    for(int i = 0; i < 5; i++)
    {
        sprintf(fname_out, "Hist_v2_pt_obs2_caysm_%d", i);
        Hist_v2_pt_obs2_caysm[i]->Write(fname_out, TObject::kOverwrite);
        sprintf(fname_out, "Hist_v2_pt_obs2_p_caysm_%d", i);
        Hist_v2_pt_obs2_p_caysm[i]->Write(fname_out, TObject::kOverwrite);
    }

    num_lam_tree->Write();
    num_proton_tree->Write();
    num_lam_final->Write();
    num_proton_final->Write();
    num_gamma_final->Write();

    if(o != 1){
        for(int l = 0; l < 15; l++){
            pParity_int_obs3_splitpt[l]->Write();
            pParity_int_ss_obs3_splitpt[l]->Write();
            pDelta_int_ss_obs3_splitpt[l]->Write();
        }
    }

	fout->Write();
    fout->Close();
}
//////////////////////////////////////////////////
void WriteWeight(TString OutFileName)
{
    TFile *fWgtNew = new TFile(OutFileName, "UPDATE");
    Hist_netChAsym->Write();
    pTPCmeanPhi_FF_1->Write();
    pTPCmeanPhi_RF_1->Write();
    pTPCmeanPhi_FF_1_p->Write();
    pTPCmeanPhi_RF_1_p->Write();
    pTPCmeanPhiAsso_FF_1->Write();
    pTPCmeanPhiAsso_RF_1->Write();
    pTPCmeanPhi_FF_2->Write();
    pTPCmeanPhi_RF_2->Write();
    pTPCmeanPhi_FF_2_p->Write();
    pTPCmeanPhi_RF_2_p->Write();
    pTPCmeanPhiAsso_FF_2->Write();
    pTPCmeanPhiAsso_RF_2->Write();
    pTPCmeanPhi_FF_3->Write();
    pTPCmeanPhi_RF_3->Write();
    pTPCmeanPhi_FF_3_p->Write();
    pTPCmeanPhi_RF_3_p->Write();
    pTPCmeanPhiAsso_FF_3->Write();
    pTPCmeanPhiAsso_RF_3->Write();
    pTPC_EP_east->Write();
    pTPC_EP_west->Write();
    pTPC_EP_for->Write();
    pTPC_EP_bac->Write();
    pTPC_EP_full->Write();
    //  rc->Write();
    fWgtNew->Close();
}
/////////////////////////////////////////////////////
int ReadWeight(char *InFileName)
{
    TFile *fWgt = new TFile(InFileName, "READ");
    if(!fWgt->IsOpen()) return 0;
    if(fWgt->IsOpen())
    {
        TOF_eff = (TH1D *)fWgt->Get("rc");
        if(TOF_eff && TOF_eff->GetEntries())
        {
            float cont = TOF_eff->GetBinContent(20);
            TOF_eff->Scale(1.25 / cont);
        }
        TPCmean_FF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_1");
        TPCmean_RF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_1");
        TPCmean_FF_1_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_1_p");
        TPCmean_RF_1_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_1_p");
        TPCmeanAsso_FF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_FF_1");
        TPCmeanAsso_RF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_RF_1");
        TPCmean_FF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_2");
        TPCmean_RF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_2");
        TPCmean_FF_2_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_2_p");
        TPCmean_RF_2_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_2_p");
        TPCmeanAsso_FF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_FF_2");
        TPCmeanAsso_RF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_RF_2");
        TPCmean_FF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_3");
        TPCmean_RF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_3");
        TPCmean_FF_3_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_3_p");
        TPCmean_RF_3_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_3_p");
        TPCmeanAsso_FF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_FF_3");
        TPCmeanAsso_RF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_RF_3");
        Read_TPC_EP_full = (TProfile2D *)fWgt->Get("pTPC_EP_full");
        Read_TPC_EP_east = (TProfile2D *)fWgt->Get("pTPC_EP_east");
        Read_TPC_EP_west = (TProfile2D *)fWgt->Get("pTPC_EP_west");
        Read_TPC_EP_for = (TProfile2D *)fWgt->Get("pTPC_EP_for");
        Read_TPC_EP_bac = (TProfile2D *)fWgt->Get("pTPC_EP_bac");
        cout << "Loaded: TPC/BBC/EPD/ZDC EP corrections" << endl;
        TH1D *Read_netChAsym   = (TH1D *)fWgt->Get("Hist_netChAsym");
        MeanNetChargeAsym = Read_netChAsym->GetMean();
        RMSNetChargeAsym = Read_netChAsym->GetRMS();
        StdDevNetChargeAsym = Read_netChAsym->GetStdDev();
        cout << "Loaded: charge asymmetry " << endl;
        cout << "MeanNetChargeAsym = " << MeanNetChargeAsym << endl;
        cout << "StdDevNetChargeAsym = " << StdDevNetChargeAsym << endl;
    }
    return 1;
}
///////////////////////////////////////////////////
void FillCMW()
{
    p_v2_Ach->Fill(net_Nch_Asym, v2_sub);
    if(net_Nch_Asym <= -99) return;
    if(net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(1, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(1, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(1, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(1, Pt);
        }
    }
    else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym)))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(2, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(2, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(2, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(2, Pt);
        }
    }
    else if(net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym)))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(3, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(3, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(3, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(3, Pt);
        }
    }
    else if(net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(4, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(4, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(4, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(4, Pt);
        }
    }
    else
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(5, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(5, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(5, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(5, Pt);
        }
    }
}
//////////////////////////////////////////////////
void FillGamma(int ord)
{
    if(ord == 4)
    {
        pParity_int_ss_oppo_run->Fill(Day2, correlator4, Eweight / eff / eff2);
        pDelta_int_ss_oppo_run->Fill(Day2, correlator3, Eweight / eff / eff2);
    }
    if(ord == 3)
    {
        pParity_int_ss_same_run->Fill(Day2, correlator4, Eweight / eff / eff2);
        pDelta_int_ss_same_run->Fill(Day2, correlator3, Eweight / eff / eff2);
    }
    pParity_int_obs1->Fill(ord, correlator0);
    pParity_int_obs3->Fill(ord, correlator0, Eweight / eff / eff2);
    pParity_int_obs3_splitpt[(int)floor((Pt-0.5)/0.1)]->Fill(ord, correlator0, Eweight / eff / eff2);
    pParity_int_ss_obs1->Fill(ord, correlator4);
    pParity_int_ss_obs3->Fill(ord, correlator4, Eweight / eff / eff2);
    pParity_int_ss_obs3_splitpt[(int)floor((Pt-0.5)/0.1)]->Fill(ord, correlator4, Eweight / eff / eff2);
    pDelta_int_ss_obs1->Fill(ord, correlator3);
    pDelta_int_ss_obs3->Fill(ord, correlator3, Eweight / eff / eff2);
    pDelta_int_ss_obs3_splitpt[(int)floor((Pt-0.5)/0.1)]->Fill(ord, correlator3, Eweight / eff / eff2);
    pParity_eta_ss_obs1->Fill(ord, 0.5 * (Eta + Eta2), correlator4);
    pParity_eta_ss_obs3->Fill(ord, 0.5 * (Eta + Eta2), correlator4, Eweight / eff / eff2);
    pParity_Deta_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4);
    pParity_Deta_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4, Eweight / eff / eff2);
    
    //new graphs
    pParity_pt_ss_obs1->Fill(ord, Pt, correlator0, Eweight / eff / eff2); //Gamma132
    pParity_pt_ss_obs3->Fill(ord, Pt, correlator4, Eweight / eff / eff2); //Gamma112
    pDelta_pt_ss_obs3->Fill(ord, Pt, correlator3, Eweight / eff / eff2);
    
    pParity_Dpt_ss_obs1->Fill(ord, fabs(Pt - Pt2), correlator4);
    pParity_Dpt_ss_obs3->Fill(ord, fabs(Pt - Pt2), correlator4, Eweight / eff / eff2);
    pDelta_eta_ss_obs1->Fill(ord, 0.5 * (Eta + Eta2), correlator3);
    pDelta_eta_ss_obs3->Fill(ord, 0.5 * (Eta + Eta2), correlator3, Eweight / eff / eff2);
    pDelta_Deta_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator3);
    pDelta_Deta_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator3, Eweight / eff / eff2);
    pDelta_pt_ss_obs1->Fill(ord, 0.5 * (Pt + Pt2), correlator3);
    pDelta_Dpt_ss_obs1->Fill(ord, fabs(Pt - Pt2), correlator3);
    pDelta_Dpt_ss_obs3->Fill(ord, fabs(Pt - Pt2), correlator3, Eweight / eff / eff2);
    if(fabs(Eta) < 0.5 && fabs(Eta2) < 0.5)
    {
        pTemp_parity_e->Fill(ord, correlator4e, 1. / eff / eff2);
        pTemp_parity_w->Fill(ord, correlator4w, 1. / eff / eff2);
        pTemp_delta->Fill(ord, correlator3, 1. / eff / eff2);
        if(fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
        {
            pTemp_parity_e_noHBT->Fill(ord, correlator4e, 1. / eff / eff2);
            pTemp_parity_w_noHBT->Fill(ord, correlator4w, 1. / eff / eff2);
            pTemp_delta_noHBT->Fill(ord, correlator3, 1. / eff / eff2);
        }
    }
    if(fabs(Pt - Pt2) > 0.15)
    {
        pParity_Deta_highDpt_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4);
        pParity_Deta_highDpt_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4, Eweight / eff / eff2);
        pDelta_Deta_highDpt_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4);
        pDelta_Deta_highDpt_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4, Eweight / eff / eff2);
    }
    if(fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
    {
        pParity_noHBT_ss_obs1->Fill(ord, correlator4);
        pParity_noHBT_ss_obs3->Fill(ord, correlator4, Eweight / eff / eff2);
        pDelta_noHBT_ss_obs1->Fill(ord, correlator3);
        pDelta_noHBT_ss_obs3->Fill(ord, correlator3, Eweight / eff / eff2);
    }

    if(ord == 3)
    {
        pParity_int_pt_ss_obs3->Fill(Pt, correlator4, Eweight / eff / eff2);
        pParity_int_pt_ssb_obs3->Fill(Pt, correlator0, Eweight / eff / eff2);
    }
    if(ord == 4)
    {
        pParity_int_pt_os_obs3->Fill(Pt, correlator4, Eweight / eff / eff2);
        pParity_int_pt_osb_obs3->Fill(Pt, correlator0, Eweight / eff / eff2);
    }
}

void determine_particlelists()
{

    if(particle1->size() != 0) particle1->clear();
    if(particle2->size() != 0) particle2->clear();

    if(particle_option1 == 0)  // first track is LAMBDA
    {

        for(int part1 = 0 ; part1 < p_lambda->size() ; part1++)
        {

            particle1->push_back(p_lambda->at(part1));
        }
    }
    else if(particle_option1 == 1)
    {

        for(int part1 = 0 ; part1 < p_proton->size() ; part1++)
        {
            p_proton->at(part1).mass = 0.93827;
            particle1->push_back(p_proton->at(part1));
        }
    }


    if(particle_option2 == 0)  // second track is PROTON
    {

        for(int part2 = 0 ; part2 < p_proton->size() ; part2++)
        {
            p_proton->at(part2).mass = 0.93827;
            particle2->push_back(p_proton->at(part2));
        }
    }
}
