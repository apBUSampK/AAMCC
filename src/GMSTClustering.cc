#include "../include/GMSTClustering.hh"
#include <queue>

bool cd_comp(const GNode& left, const GNode& right) {
    return left.height > right.height;
}

GMSTCluster::GMSTCluster(G4int Z_in, G4int A_in): 
A(A_in), Z(Z_in)
{
};

GMSTCluster::~GMSTCluster(){
};

GMSTClustering::GMSTClustering(){
CritDist = 100;
};

GMSTClustering::GMSTClustering(G4double CD_in, G4double variation, G4double single_silh) : CritDist(CD_in), variation(variation), single_silh(single_silh)
{
};


GMSTClustering::~GMSTClustering() = default;

void GMSTClustering::SetUp(NucleonVector* nucleons_in, G4double ExA, G4double ExB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB) {
    SpecAa = Z = SpecAb = Zb = 0;
    this->ExA = ExA;
    this->ExB = ExB;
    this->boostA = boostA;
    this->boostB = boostB;
    SpecAa = nucleons_in->GetA("A");
    SpecAb = nucleons_in->GetA("B");
    nucleons = NucleonVector();
    nucleons_B = NucleonVector();
    for(int iArray = 0; iArray < nucleons_in->size(); iArray++) {
        Nucleon *nucleon = &(nucleons_in->at(iArray));
        if (nucleon->isParticipant == 0) {
            switch ("A" == nucleon->Nucl ? 0 : 1) {
                case 0:
                    nucleons.push_back(*nucleon);
                    if (nucleon->isospin == 1)
                        Z++;
                    break;
                case 1:
                    nucleons_B.push_back(*nucleon);
                    if (nucleon->isospin == 1)
                        Zb++;
            }
        }
    }
    if(ExA > 0 && ExB > 0) this->SetCDExEn(ExA, SpecAa);
    CritDistA = CritDist;
    if(ExA > 0 && ExB > 0) this->SetCDExEn(ExB, SpecAb);
    g = Graph(SpecAa, SpecAa * (SpecAa - 1) / 2);
    g_B = Graph(SpecAb, SpecAb * (SpecAb - 1) / 2);
    //  making full graphs of nucleons
    for(G4int iArray = 0; iArray < nucleons.size(); iArray++){
        const Nucleon &nucleon = nucleons[iArray];
        for(auto iArray_pairs = iArray + 1; iArray_pairs < nucleons.size(); iArray_pairs++) {
            const Nucleon &nucleon_pair = nucleons[iArray_pairs];
            g.addEdge(iArray, iArray_pairs, std::sqrt(
                    pow(nucleon.GetX() - nucleon_pair.GetX(), 2) + pow(nucleon.GetY() - nucleon_pair.GetY(), 2) +
                    pow(nucleon.GetZ() - nucleon_pair.GetZ(), 2)));
        }
    }
    for(G4int iArray = 0; iArray < nucleons_B.size(); iArray++){
        const Nucleon &nucleon = nucleons_B[iArray];
        for(auto iArray_pairs = iArray + 1; iArray_pairs < nucleons_B.size(); iArray_pairs++) {
            const Nucleon &nucleon_pair = nucleons_B[iArray_pairs];
            g_B.addEdge(iArray, iArray_pairs, std::sqrt(
                    pow(nucleon.GetX() - nucleon_pair.GetX(), 2) + pow(nucleon.GetY() - nucleon_pair.GetY(), 2) +
                    pow(nucleon.GetZ() - nucleon_pair.GetZ(), 2)));
        }
    }
}

void GMSTClustering::SetCDExEn(G4double Ex, G4int NoN) {
   if(Ex/G4double(NoN) < 2.17 * MeV){ CritDist = d0;}
   else{//kappa = std::pow(1/(1+aColRel*std::pow(G4double(SpecAa),-0.333333333333333333)+aSurfRel*std::pow(G4double(SpecAa),-0.666666666666666)),0.3333333333333);
       kappa = 1;
       CritDist = d0*kappa*std::pow(Ex/(eps0*G4double(SpecAa)), 1. / 3. * alphaPow);}
}

std::vector<G4FragmentVector> GMSTClustering::GetClusters(cut) {
    std::vector<std::vector<G4int>> clusters_A, clusters_B;

    //get dendrograms for both nucleons
    GTree treeA = g.AdvancedKruskalMSTDendro();
    GTree treeB = g_B.AdvancedKruskalMSTDendro();
    //get the clustering with corresponding CD:
    for (const auto & iter : treeA.get_cluster(CritDistA))
        clusters_A.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));
    for (const auto & iter : treeB.get_cluster(CritDist))
        clusters_B.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //compile output vectors and calculate momentum
    return CalculateMomentum(CompileVector(clusters_A, clusters_B), ExA, ExB, boostA, boostB);
};

std::vector<G4FragmentVector> GMSTClustering::GetClusters(silhouette) {
    std::vector<std::vector<G4int>> clusters_A, clusters_B;

    //get dendrograms for both nucleons
    GTree treeA = g.AdvancedKruskalMSTDendro();
    GTree treeB = g_B.AdvancedKruskalMSTDendro();

    //get the clustering with the biggest applicable CD:
    std::vector<GNode> current = treeA.get_cluster(CritDistA * (1.0 + variation));
    std::sort(current.begin(), current.end(), cd_comp);
    G4double max_silh = -1, min_dist = kInfinity;
    std::vector<GNode> best = current;
    //calculate mean cluster silhouettes until min applicable CD
    //separate nucleon clusters are given silhouette of single_silh
    if (!current.empty())
        while (current.front().height > CritDistA * (1.0 - (variation > 1 ? 1 : variation))) {
            G4double silh = 0, dist = abs(current.front().height - CritDistA);
            for (auto iter = current.cbegin(); iter != current.cend(); iter++) {
                if (iter->size > 1)
                    for (G4int i = 0; i < iter->size; i++) {
                        G4double inner = 0, outer = -1;
                        for (G4int j = 0; j < iter->size; j++)
                            inner += g.adj[iter->V[i] - 1][iter->V[j] - 1];
                        inner /= (iter->size - 1);
                        for (auto jter = current.cbegin(); jter != current.cend(); jter++)
                            if (iter != jter) {
                                G4double buff = 0;
                                for (G4int j = 0; j < jter->size; j++)
                                    buff += g.adj[iter->V[i] - 1][jter->V[j] - 1];
                                buff /= jter->size;
                                if (outer < 0 || buff < outer)
                                    outer = buff;
                            }
                        silh += ((outer - inner) / max(outer, inner));
                    }
                else
                    silh += single_silh;
            }
            silh /= SpecAa;
            if (silh > max_silh || (silh == max_silh && dist < min_dist)) {
                max_silh = silh;
                best = current;
                min_dist = dist;
            }
            //divide the biggest cluster into two
            current.push_back(*current.front().children.first);
            current.push_back(*current.front().children.second);
            current.erase(current.begin());
            sort(current.begin(), current.end(), cd_comp);
        }
    //conversion to vector<vector<int>>
    for (const auto &iter: best)
        clusters_A.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //repeat for B nucleus

    //get the clustering with the biggest applicable CD:
    current = treeB.get_cluster(CritDist * (1.0 + variation));
    std::sort(current.begin(), current.end(), cd_comp);
    max_silh = -1, min_dist = kInfinity;
    best = current;
    //calculate mean cluster silhouettes until min applicable CD
    //separate nucleon clusters are given silhouette of single_silh
    if (!current.empty())
        while (current.front().height > CritDist * (1.0 - (variation > 1 ? 1 : variation))) {
            G4double silh = 0, dist = abs(current.front().height - CritDist);
            for (auto iter = current.cbegin(); iter != current.cend(); iter++) {
                if (iter->size > 1)
                    for (G4int i = 0; i < iter->size; i++) {
                        G4double inner = 0, outer = -1;
                        for (G4int j = 0; j < iter->size; j++)
                            inner += g_B.adj[iter->V[i] - 1][iter->V[j] - 1];
                        inner /= (iter->size - 1);
                        for (auto jter = current.cbegin(); jter != current.cend(); jter++)
                            if (iter != jter) {
                                G4double buff = 0;
                                for (G4int j = 0; j < jter->size; j++)
                                    buff += g_B.adj[iter->V[i] - 1][jter->V[j] - 1];
                                buff /= jter->size;
                                if (outer < 0 || buff < outer)
                                    outer = buff;
                            }
                        silh += ((outer - inner) / max(outer, inner));
                    }
                else
                    silh += single_silh;
            }
            silh /= SpecAb;
            if (silh > max_silh || (silh == max_silh && dist < min_dist)) {
                max_silh = silh;
                best = current;
                min_dist = dist;
            }
            //divide the biggest cluster into two
            current.push_back(*current.front().children.first);
            current.push_back(*current.front().children.second);
            current.erase(current.begin());
            sort(current.begin(), current.end(), cd_comp);
        }
    //conversion to vector<vector<int>>
    for (const auto &iter: best)
        clusters_B.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //compile output vectors and calculate momentum
    return CalculateMomentum(CompileVector(clusters_A, clusters_B), ExA, ExB, boostA, boostB);
}

std::vector<G4FragmentVector> GMSTClustering::GetClusters(max_alpha) {
    std::vector<std::vector<G4int>> clusters_A, clusters_B;

    //get dendrograms for both nucleons
    GTree treeA = g.AdvancedKruskalMSTDendro();
    GTree treeB = g_B.AdvancedKruskalMSTDendro();

    //get the clustering with the biggest applicable CD:
    std::vector<GNode> current = treeA.get_cluster(CritDistA * (1.0 + variation));
    std::sort(current.begin(), current.end(), cd_comp);
    G4int max_alpha = 0;
    G4double min_dist = kInfinity;
    std::vector<GNode> best = current;
    //find the cluster with the largest alpha particles count
    if (!current.empty())
        while (current.front().height > CritDistA * (1.0 - (variation > 1 ? 1 : variation))) {
            G4int alpha = 0;
            G4double dist = abs(current.front().height - CritDistA);
            for (auto & iter : current) {
                G4int z_count = 0;
                for (G4int i = 0; i < iter.size; i++)
                    if ((nucleons.at(iter.V[i] - 1)).isospin == 1)
                        z_count++;
                if (z_count == 2 && iter.size == 4)
                    alpha++;
            }
            if (alpha > max_alpha || (alpha == max_alpha && dist < min_dist)) {
                max_alpha = alpha;
                best = current;
                min_dist = dist;
            }
            //divide the biggest cluster into two
            current.push_back(*current.front().children.first);
            current.push_back(*current.front().children.second);
            current.erase(current.begin());
            sort(current.begin(), current.end(), cd_comp);
        }
    //conversion to vector<vector<int>>
    for (const auto &iter: best)
        clusters_A.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //repeat for B nucleus

    //get the clustering with the biggest applicable CD:
    current = treeB.get_cluster(CritDist * (1.0 + variation));
    std::sort(current.begin(), current.end(), cd_comp);
    max_alpha = 0;
    min_dist = kInfinity;
    best = current;
    //find the cluster with the largest alpha particles count
    if (!current.empty())
        while (current.front().height > CritDist * (1.0 - (variation > 1 ? 1 : variation))) {
            G4int alpha = 0;
            G4double dist = abs(current.front().height - CritDist);
            for (auto & iter : current) {
                G4int z_count = 0;
                for (G4int i = 0; i < iter.size; i++)
                    if ((nucleons_B.at(iter.V[i] - 1)).isospin == 1)
                        z_count++;
                if (z_count == 2 && iter.size == 4)
                    alpha++;
            }
            if (alpha > max_alpha || (alpha == max_alpha && dist < min_dist)) {
                max_alpha = alpha;
                best = current;
                min_dist = dist;
            }
            //divide the biggest cluster into two
            current.push_back(*current.front().children.first);
            current.push_back(*current.front().children.second);
            current.erase(current.begin());
            sort(current.begin(), current.end(), cd_comp);
        }
    //conversion to vector<vector<int>>
    for (const auto &iter: best)
        clusters_B.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //compile output vectors and calculate momentum
    return CalculateMomentum(CompileVector(clusters_A, clusters_B), ExA, ExB, boostA, boostB);
}

std::vector<G4FragmentVector> GMSTClustering::GetClusters(alpha_destroy) {
    std::vector<std::vector<G4int>> clusters_A, clusters_B;

    //get dendrograms for both nucleons
    GTree treeA = g.AdvancedKruskalMSTDendro();
    GTree treeB = g_B.AdvancedKruskalMSTDendro();

    //get the clustering with the biggest applicable CD:
    std::vector<GNode> current = treeA.get_cluster(CritDistA * (1.0 + variation));
    std::sort(current.begin(), current.end(), cd_comp);
    G4int max_alpha = 0;
    vector<G4double> min_dist(5, kInfinity);
    std::vector<std::vector<GNode>> all(5, std::vector<GNode>());

    //find the cluster with the largest alpha particles count
    if (!current.empty())
        while (current.front().height > CritDistA * (1.0 - (variation > 1 ? 1 : variation))) {
            G4int alpha = 0;
            G4double dist = abs(current.front().height - CritDistA);
            for (auto & iter : current) {
                G4int z_count = 0;
                for (G4int i = 0; i < iter.size; i++)
                    if ((nucleons.at(iter.V[i] - 1)).isospin == 1)
                        z_count++;
                if (z_count == 2 && iter.size == 4)
                    alpha++;
            }
            if (alpha > max_alpha)
                max_alpha = alpha;
            if (all[alpha].empty() || dist < min_dist[alpha]) {
                all[alpha] = current;
                min_dist[alpha] = dist;
            }
            //divide the biggest cluster into two
            current.push_back(*current.front().children.first);
            current.push_back(*current.front().children.second);
            current.erase(current.begin());
            sort(current.begin(), current.end(), cd_comp);
        }
    //destroy clusters
    G4double p_dest = 0.002;
    G4int n_destr = gRandom->Binomial(max_alpha, p_dest);
    G4int n_remain = max_alpha - n_destr;
    while(all[n_remain].empty()) {
        n_remain = (n_remain) ? (n_remain - 1) : 4; //loop over [0; 4]
        if (n_remain == max_alpha - n_destr)
            break;
    }
    //conversion to vector<vector<int>>
    for (const auto &iter: all[n_remain])
        clusters_A.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //repeat for B nucleus

    //get the clustering with the biggest applicable CD:
    current = treeB.get_cluster(CritDist * (1.0 + variation));
    std::sort(current.begin(), current.end(), cd_comp);
    max_alpha = 0;
    min_dist.assign(5, kInfinity);
    all.assign(5, std::vector<GNode>());

    //find the cluster with the largest alpha particles count
    if (!current.empty())
        while (current.front().height > CritDist * (1.0 - (variation > 1 ? 1 : variation))) {
            G4int alpha = 0;
            G4double dist = abs(current.front().height - CritDist);
            for (auto & iter : current) {
                G4int z_count = 0;
                for (G4int i = 0; i < iter.size; i++)
                    if ((nucleons_B.at(iter.V[i] - 1)).isospin == 1)
                        z_count++;
                if (z_count == 2 && iter.size == 4)
                    alpha++;
            }
            if (alpha > max_alpha)
                max_alpha = alpha;
            if (all[alpha].empty() || dist < min_dist[alpha]) {
                all[alpha] = current;
                min_dist[alpha] = dist;
            }
            //divide the biggest cluster into two
            current.push_back(*current.front().children.first);
            current.push_back(*current.front().children.second);
            current.erase(current.begin());
            sort(current.begin(), current.end(), cd_comp);
        }
    //destroy clusters
    n_destr = gRandom->Binomial(max_alpha, p_dest);
    n_remain = max_alpha - n_destr;
    while(all[n_remain].empty()) {
        n_remain = (n_remain) ? (n_remain - 1) : 4; //loop over [0; 4]
        if (n_remain == max_alpha - n_destr)
            break;
    }
    //conversion to vector<vector<int>>
    for (const auto &iter: all[n_remain])
        clusters_B.emplace_back(std::vector<G4int>(iter.V, iter.V + iter.size));

    //compile output vectors and calculate momentum
    return CalculateMomentum(CompileVector(clusters_A, clusters_B), ExA, ExB, boostA, boostB);
}

std::vector<G4FragmentVector> GMSTClustering::CompileVector(const std::vector<std::vector<G4int>>& clusters_A, const std::vector<std::vector<G4int>>& clusters_B) {
    std::vector<G4FragmentVector> outClusters;
    G4FragmentVector output_vector_A;
    G4FragmentVector output_vector_B;

    // Filling a Cluster Vector (side A)
    for(G4int i = 0; i < clusters_A.size(); ++i) {
        G4int Z_clust = 0;
        G4int A_clust = 0;
        for(G4int j = 0; j < clusters_A[i].size(); ++j) {
            Nucleon *nucleon=&(nucleons.at((clusters_A[i])[j] - 1));
            if(nucleon->isospin == 1)
            {
                Z_clust += 1;
            }
            A_clust += 1;
        }

        G4Fragment* frag = new G4Fragment();
        frag->SetA(A_clust);
        frag->SetZ(Z_clust);
        output_vector_A.push_back(frag);
    }

    // Filling a Cluster Vector (side B)
    for(G4int i = 0; i < clusters_B.size(); ++i) {
        G4int Z_clust = 0;
        G4int A_clust = 0;
        for(G4int j = 0; j < clusters_B[i].size(); ++j) {
            Nucleon *nucleon=&(nucleons_B.at((clusters_B[i])[j] - 1));
            if(nucleon->isospin == 1)
            {
                Z_clust += 1;
            }
            A_clust += 1;
        }
        G4Fragment* frag = new G4Fragment();
        frag->SetA(A_clust);
        frag->SetZ(Z_clust);
        output_vector_B.push_back(frag);
    }

    outClusters.push_back(output_vector_A);
    outClusters.push_back(output_vector_B);
    return outClusters;
}

std::vector<G4FragmentVector> GMSTClustering::CalculateMomentum(std::vector<G4FragmentVector> noMomClusters, G4double ExEnA, G4double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB) {
    std::vector<G4FragmentVector> momClusters = noMomClusters;
    std::vector<G4double> MstMassVector_A;
    MstMassVector_A.reserve(noMomClusters.at(0).size());

    std::vector<G4double> MstMassVector_B;
    MstMassVector_B.reserve(noMomClusters.at(1).size());

    G4double SumMassMst = 0;
    G4double SumMassMstEx = 0;
    if(!noMomClusters.at(0).empty()) {

        for (G4int i = 0; i < noMomClusters.at(0).size(); ++i) {
            G4int clfrag_A = (noMomClusters.at(0)[i])->GetA();
            G4int clfrag_Z = (noMomClusters.at(0)[i])->GetZ();
            //Excitation energy computing
            G4double energy = 0;
            if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
                energy = ExEnA * G4double(clfrag_A) / G4double(SpecAa);
            }

            G4double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

            SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

            SumMassMstEx += NuclearMass;

            MstMassVector_A.push_back(NuclearMass);
        }

        G4double PrefragmentMass_A = SumMassMst + ExEnA;

        std::vector<G4LorentzVector *> *momentumVectorA;
        //if is commented as a part of a long history
        if (PrefragmentMass_A < (SumMassMstEx + 1e-5*MeV)) { PrefragmentMass_A += 1e-5*MeV; }
            momentumVectorA = phaseSpaceDecay.Decay(PrefragmentMass_A, MstMassVector_A);
        //}
        for (int I = 0; I < momClusters.at(0).size(); ++I) {
            momClusters.at(0).at(I)->SetMomentum((*momentumVectorA->at(I)).boost(boostA));
        }
        momentumVectorA->clear();
    }
    //side B
    SumMassMst = 0;
    SumMassMstEx = 0;
    if(!noMomClusters.at(1).empty()) {
        for (G4int i = 0; i < noMomClusters.at(1).size(); ++i) {
            G4int clfrag_A = (noMomClusters.at(1)[i])->GetA();
            G4int clfrag_Z = (noMomClusters.at(1)[i])->GetZ();
            //Excitation energy computing
            G4double energy = 0;
            if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
                energy = ExEnB * G4double(clfrag_A) / G4double(SpecAb);
            }

            G4double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

            SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

            SumMassMstEx += NuclearMass;

            MstMassVector_B.push_back(NuclearMass);
        }

        G4double PrefragmentMass_B = SumMassMst + ExEnB;
        std::vector<G4LorentzVector *> *momentumVectorB;
        //if is commented as a part of a history
        if (PrefragmentMass_B < (SumMassMstEx + 1e-5*MeV)) {PrefragmentMass_B += 1e-5*MeV;}
            momentumVectorB = phaseSpaceDecay.Decay(PrefragmentMass_B, MstMassVector_B);
        //}

        for (int I = 0; I < momClusters.at(1).size(); ++I) {
            momClusters.at(1).at(I)->SetMomentum((*momentumVectorB->at(I)).boost(boostB));
        }
        momentumVectorB->clear();
    }

    MstMassVector_A.clear();
    MstMassVector_B.clear();
    noMomClusters.clear();

    return momClusters;
};

Graph::Graph(G4int V, G4int E)
{
    this->V = V;
    this->E = E;
    adj = std::vector<std::vector<G4double>>(V, std::vector<G4double>(V, 0));
}

Graph::Graph()
{
    this->V = 0;
    this->E = 0;
    adj = std::vector<std::vector<G4double>>(V, std::vector<G4double>(V, 0));
}

Graph::~Graph() = default;

void Graph::addEdge(G4int u, G4int v, G4double w)
{
    edges.push_back({ w, {u, v} });
    adj[u][v] = adj[v][u] = w;
}

GTree Graph::AdvancedKruskalMSTDendro()
{
    // Sort edges in increasing order on basis of cost 
    sort(edges.begin(), edges.end());

    // Create a tree handler
    GTree tr(V);

    // Iterate through all sorted edges 
    std::vector< std::pair<G4double, iPair> >::iterator it;
    for (it = edges.begin(); it != edges.end(); it++)
    {
        G4int u = it->second.first;
        G4int v = it->second.second;

        // Check if the selected edge is creating 
        // a cycle or not (Cycle is created if u 
        // and v belong to same set) 
        if (tr.get_node(u) != tr.get_node(v))
        {
            // create new node
            tr.merge(u, v, it->first);
        }
    }

    return tr;
}

GNode::GNode(G4int n) {
    size = 1;
    V = new G4int[1];
    V[0] = n;
    height = 0;
    children = std::make_pair(nullptr, nullptr);
}

GNode::GNode(std::shared_ptr<GNode> first, std::shared_ptr<GNode> second, G4double height) {
    size = first->size + second->size;
    V = new G4int[size];
    for (G4int i = 0; i < size; i++) {
        G4int buff = (i < first->size) ? first->V[i] : second->V[i - first->size];
        V[i] = buff;
    }
    this->height = height;
    children = std::make_pair(first, second);
}

GNode::GNode(const GNode& right) {
    size = right.size;
    V = new G4int[size];
    for (G4int i = 0; i < size; i++)
        V[i] = right.V[i];
    height = right.height;
    if (right.children.first == nullptr)
        children = std::make_pair(std::shared_ptr<GNode>(nullptr),std::shared_ptr<GNode>(nullptr));
    else
        children = std::make_pair(make_shared<GNode>(*right.children.first), make_shared<GNode>(*right.children.second));
}

GNode& GNode::operator=(const GNode& right) {
    size = right.size;
    delete[] V;
    V = new G4int[size];
    for (G4int i = 0; i < size; i++)
        V[i] = right.V[i];
    height = right.height;
    if (right.children.first == nullptr)
        children = std::make_pair(std::shared_ptr<GNode>(nullptr),std::shared_ptr<GNode>(nullptr));
    else
        children = std::make_pair(make_shared<GNode>(*right.children.first), make_shared<GNode>(*right.children.second));
    return *this;
}

GNode::~GNode() {
    delete[] V;
}

GTree::GTree(int size) {
    this->size = size;
    nodes = std::vector<std::shared_ptr<GNode>>(size);
    for (G4int i = 0; i < size; i++)
        nodes[i] = make_shared<GNode>(i + 1);
}

GTree::GTree(const GTree& right) {
    size = right.size;
    nodes = std::vector<std::shared_ptr<GNode>>(size);
    G4bool* cpd = new G4bool[size];
    for (G4int i = 0; i < size; i++)
        cpd[i] = false;
    for (G4int i = 0; i < size; i++) {
        if (!cpd[i]) {
            nodes[i] = make_shared<GNode>(*right.nodes[i]);
            for (G4int j = 0; j < nodes[i]->size; j++) {
                cpd[nodes[i]->V[j]] = true;
                nodes[nodes[i]->V[j]] = nodes[i];
            }
        }
    }
    delete[] cpd;
}

GTree& GTree::operator=(const GTree& right) {
    size = right.size;
    nodes = std::vector<std::shared_ptr<GNode>>(size);
    G4bool *cpd = new G4bool[size];
    for (G4int i = 0; i < size; i++)
        cpd[i] = false;
    for (G4int i = 0; i < size; i++) {
        if (!cpd[i]) {
            nodes[i] = make_shared<GNode>(*right.nodes[i]);
            for (G4int j = 0; j < nodes[i]->size; j++) {
                cpd[nodes[i]->V[j]] = true;
                nodes[nodes[i]->V[j]] = nodes[i];
            }
        }
    }
    delete[] cpd;
    return *this;
}

void GTree::merge(G4int a, G4int b, G4double height) {
    G4int t = nodes[a]->size + nodes[b]->size;
    G4int* chngPtr = new G4int[t];
    for (G4int i = 0; i < nodes[a]->size; i++)
        chngPtr[i] = nodes[a]->V[i] - 1;
    for (G4int i = 0; i < nodes[b]->size; i++)
        chngPtr[i + nodes[a]->size] = nodes[b]->V[i] - 1;
    nodes[a] = make_shared<GNode>(nodes[a], nodes[b], height);
    for (G4int i = 0; i < t; i++)
        nodes[chngPtr[i]] = nodes[a];
    delete[] chngPtr;
}

std::vector<GNode> GTree::get_cluster(G4double CD) {
    if (size < 1)
        return {};
    std::vector<GNode> cluster;
    queue<GNode*> width;
    width.push(nodes[0].get());
    while(!width.empty()) {
        GNode* node = width.front();
        width.pop();
        if (node->height <= CD)
            cluster.push_back(*node);
        else if (node->children.first != nullptr) {
            width.push(node->children.first.get());
            width.push(node->children.second.get());
        }
    }
    return cluster;
}