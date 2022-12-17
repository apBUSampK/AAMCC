//
// Created by apbus_amp_k on 12.11.22.
//

#include "G4LorentzVector.hh"

#include "MCIniReader.hh"

#include "TClonesArray.h"

MCIniReader::MCIniReader(MCIniReader && rha) noexcept  :
    ftree(rha.ftree), curr(nullptr), curr_st(nullptr), asum(rha.asum), iter(rha.iter), treen(rha.treen) {
    rha.ftree = nullptr;
    ftree->SetBranchAddress("event", &curr);
    ftree->SetBranchAddress("iniState", &curr_st);
}

MCIniReader& MCIniReader::operator=(MCIniReader && rha) noexcept {
    if (this == &rha)
        return *this;
    ftree = rha.ftree;
    rha.ftree = nullptr;
    curr = nullptr;
    rha.curr = nullptr;
    curr_st = nullptr;
    rha.curr_st = nullptr;
    asum = rha.asum;
    iter = rha.iter;
    treen = rha.treen;
    ftree->SetBranchAddress("event", &curr);
    ftree->SetBranchAddress("iniState", &curr_st);
}

MCIniReader::MCIniReader(const std::unique_ptr<TFile>& tfile) : curr(nullptr), curr_st(nullptr), iter(0) {
    tfile->GetObject("events", ftree);          // Legacy weirdness. ptr has to be lvalue thus preventing smart pointer usage
    treen = ftree->GetEntries();
    ftree->SetBranchAddress("event", &curr);
    ftree->SetBranchAddress("iniState", &curr_st);
    URun* run_data;
    tfile->GetObject("run", run_data);          // legacy
    asum = run_data->GetATarg() + run_data->GetAProj();
}

AAMCCinput MCIniReader::operator()() {
//    if (iter >= treen)
//        throw std::out_of_range("No more events in file!");
    ftree->GetEntry(iter);
    iter++;
    AAMCCinput cache;
    auto arr = curr->GetParticleList();
    int size = curr->GetNpa();
    for (int i = 0; i < size; i++) {
        auto particle = (UParticle*) arr->At(i);
        if(particle->GetIndex() >= asum)
            continue;
        aamcc::Nucleon nucl;            //new nucleon to write into cache->nucleons
        switch (particle->GetPdg()) {
            case 2212:                  //proton pdg code
                nucl.isospin = true;
                break;
            case 2112:                  //neutron pdg code
                nucl.isospin = false;
                break;
            default:
                continue;               //break is unnecessary
        }
        nucl.isParticipant = false;
        nucl.x = particle->X();
        nucl.y = particle->Y();         //coords
        nucl.z = particle->Z();
        if(particle->Pz() > 0) {        //proj
            nucl.Nucl = "A";
            cache.FermiMomA_x += particle->Px() * GeV;  //MCIni stores energy/impulses in GeV
            cache.FermiMomA_y += particle->Py() * GeV;
            cache.FermiMomA_z += particle->Pz() * GeV;
        }
        else {                          //targ
            nucl.Nucl = "B";
            cache.FermiMomB_x += particle->Px() * GeV;
            cache.FermiMomB_y += particle->Py() * GeV;
            cache.FermiMomB_z += particle->Pz() * GeV;
        }
        cache.nucleons.push_back(std::move(nucl));     //write the nucleon data
    }
    return cache;
}
