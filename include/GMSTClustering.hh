#ifndef GMSTClustering_h
#define GMSTClustering_h 1

#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <vector>
#include <memory>
#include "../TGlauber/TGlauberMC.hh"
#include "../TGlauber/TGlauNucleon.hh"
#include "TObject.h"
#include "TVector3.h"
#include "Nucleon.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4FermiPhaseSpaceDecay.hh"
using namespace std;

// Creating shortcut for an integer pair 
typedef  std::pair<G4int, G4int> iPair;

struct GNode {
    std::pair<std::shared_ptr<GNode>, std::shared_ptr<GNode>> children;
    G4double height;
    G4int* V;
    G4int size;

    //constructor from children nodes
    GNode(std::shared_ptr<GNode>, std::shared_ptr<GNode>, G4double);

    //initial nodes constructor
    GNode(G4int);

    GNode(const GNode&);
    GNode& operator=(const GNode&);

    //destructor
    ~GNode();
};


//the tree handler. Isn't intended to be constructed outside the Graph methods.
class GTree {
public:
    //constructor
    GTree(G4int);

    GTree(const GTree&);
    GTree& operator= (const GTree&);

    //merge 2 nodes
    void merge(G4int a, G4int b, G4double height);

    //get clusterization
    std::vector<GNode> get_cluster(G4double CD);

    //get node
    inline const GNode* get_node(G4int a) {return nodes[a].get();}

    //width print (for testing purposes)
    friend ostream& operator<< (ostream&, const GTree&);

private:
    G4int size;
    std::vector<std::shared_ptr<GNode>> nodes;
};

// Structure to represent a graph 
struct Graph
{
	// Vert and edges
    G4int V, E;
    std::vector< std::pair<G4double, iPair> > edges;
    //adjacency matrix
    std::vector<std::vector<G4double>> adj;

    // Constructors 
    Graph(G4int V, G4int E);
    Graph();

    // Destructor
	~Graph();

    // Utility function to add an edge 
    void addEdge(G4int u, G4int v, G4double w);
    // Function to find MST using Kruskal's MST algorithm with hierarchy
    GTree AdvancedKruskalMSTDendro();
};

class GMSTCluster{

	public: 
	GMSTCluster(G4int Z_in, G4int A_in);
	~GMSTCluster();

	public: 
	inline G4int GetZ() {return Z;};
	inline G4int GetA() {return A;};
	inline void SetZ(G4int Z_in) {Z = Z_in;}
	inline void SetA(G4int A_in) {A = A_in;}
	inline void PushBackCoordinateVector(TVector3 vec_in) {coord.push_back(vec_in);}
        inline std::vector<TVector3> GetCoordinates() {return coord;}

	private: 
	G4int Z;
	G4int A;
	std::vector<TVector3> coord;
};

typedef std::vector<GMSTCluster> GMSTClusterVector;

class GMSTClustering{

	public:
    GMSTClustering();
    GMSTClustering(G4double CD_in, G4double SpecAa, G4double SpecAb, G4double single_silh = 0.0, G4double variation = 20.0);
    ~GMSTClustering();

    void SetUp(NucleonVector* nucleons_in, G4double ExA, G4double ExB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB);

    GTree GetTree();
    struct cut{};
    struct silhouette{};
    struct max_alpha{};
    struct alpha_destroy{};
    std::vector<G4FragmentVector> GetClusters(cut);
    std::vector<G4FragmentVector> GetClusters(silhouette);
    std::vector<G4FragmentVector> GetClusters(max_alpha);
    std::vector<G4FragmentVector> GetClusters(alpha_destroy);

	inline void SetCD(G4double CD_in) {CritDist = CD_in; d0 = CD_in;}
	inline G4double GetCD(std::string side) {return (side == "A" ?  CritDistA : CritDist);}
    void SetCDExEn(G4double Ex, G4int NoN);


	private:

    std::vector<G4FragmentVector> CompileVector(const std::vector<std::vector<G4int>>& clusters_A, const std::vector<std::vector<G4int>>& clusters_B);
    std::vector<G4FragmentVector> CalculateMomentum(std::vector<G4FragmentVector> noMomClusters, G4double ExEnA, G4double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB);

    G4double CritDistA;
    G4double CritDist;
    G4double variation;
    G4double single_silh;
	G4double kappa;
	G4double eps0 = 2.17*MeV;
	G4double alphaPow = -1.02;
	G4double d0 = -999;
	G4double aColRel = 5.315; //divided by aVol
	G4double aSurfRel  = 0.017; //divided by aVol
	G4double aVol     = 0.054;
    G4double SpecAa, SpecAb;
    G4int Z, Zb;
    G4double ExA, ExB;

    NucleonVector nucleons, nucleons_B;
    Graph g, g_B;
    CLHEP::Hep3Vector boostA, boostB;


    G4FermiPhaseSpaceDecay phaseSpaceDecay;


};


#endif
