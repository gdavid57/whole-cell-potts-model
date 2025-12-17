#ifndef CELL_STATE_HPP
#define CELL_STATE_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <filesystem>
#include <array>

using namespace std;
using std::ifstream;
using std::ios;
using std::ofstream;

namespace fs = std::filesystem;

// Structure for tRNAs
struct tList
{
    vector<int> index;     // codon index
    int cnum;              // number of charged tRNA molecules
    int tnum;              // number of total tRNA molecules
    int len;               // number of different codons that tRNA can decode
    int aa;                // amino acid coded by the tRNA
    vector<int> dem;       // demande for the tRNA
    vector<double> k;      // k-value
    vector<double> wobble; // wobble of the corresponding codon
};

struct CellStepResult
{
    int produced_proteins;
    int free_ribosomes;
    double physical_time;
};

// CellState class to encapsulate the simulation
class CellState
{
private:
    // Private methods
    double compute_propensity(int j, int i) const;
    void update_kinetic_parameters(double volume);
    void initialize_ribosomes();

    // File reading methods
    bool read_mrna_count(const string &mRNA_data_path);
    bool read_trna_abundances(const string &tRNA_data_path);
    bool read_mrna_sizes(const string &mRNA_data_path);
    bool read_initiation_rates(const string &mRNA_data_path);
    bool read_vmax_km_parameters(const string &tRNA_data_path);
    bool read_codon_trna_index_table(const string &tRNA_data_path);
    bool read_wobble_factors(const string &tRNA_data_path);
    bool read_mrna_sequences(const string &mRNA_data_path);

    // Init method
    void init(const string &tRNA_data_path, const string &mRNA_data_path, double init_cell_volume, int num_ribosomes);

public:
    // Configuration parameters
    double beta;          // beta
    double gamma;         // scaling and finetuning factor for the intitation rates
    double init_charging; // initial charging level of tRNAs
    double radius;        // radius of a ribosome

    // Simulation constants
    int NT;             // number of different tRNAs
    double r;           // constant of proportionality k and ctRNAs
    double transc_rate; // transcription rate
    int Ncod;           // total number of codons
    int Naa;            // total number of amino acids

    // Simulation state variables
    int M;                 // number of mRNAs
    int Nrib;              // number of free ribosomes
    double t;              // physical time
    int t0;                // time step
    double current_volume; // current cell volume

    // Main data structures
    vector<tList> TRNA;
    vector<vector<int>> a;
    vector<vector<int>> rib;
    vector<vector<int>> mRNA;
    vector<vector<double>> rho;
    vector<vector<double>> rho_sites;

    // Variables for kinetic parameters
    vector<double> Vmax;
    vector<double> Km;
    vector<double> chi;
    vector<double> Vt;
    vector<double> init_rates;
    vector<double> alpha;
    vector<int> proteins;
    vector<int> N;    // mRNA lengths
    vector<int> tnum; // number of tRNAs

    // Tables and indices
    vector<std::array<int, 61>> ind_cod; // [2][Ncod]
    vector<std::array<char, 3>> cod;     // [Ncod][3]
    vector<int> den;                     // [Naa]

    // Fine-tuning factors
    vector<double> deltaKm;
    vector<double> deltaVmax;

    // Constructor parameters for pickling
    string trna_data_path;
    string mrna_data_path;
    double init_cell_volume;
    int num_ribosomes;

    // Constructor and destructor
    CellState(const string &tRNA_data_path, const string &mRNA_data_path, double init_cell_volume, int num_ribosomes);

    // Full state constructor for pickle/deserialization
    CellState(double beta, double gamma, double init_charging, double radius,
              int NT, double r, double transc_rate, int Ncod, int Naa,
              int M, int Nrib, double t, int t0, double current_volume,
              const vector<tList> &TRNA,
              const vector<vector<int>> &a,
              const vector<vector<int>> &rib,
              const vector<vector<int>> &mRNA,
              const vector<vector<double>> &rho,
              const vector<vector<double>> &rho_sites,
              const vector<double> &Vmax,
              const vector<double> &Km,
              const vector<double> &chi,
              const vector<double> &Vt,
              const vector<double> &init_rates,
              const vector<double> &alpha,
              const vector<int> &proteins,
              const vector<int> &N,
              const vector<int> &tnum,
              const vector<array<int, 61>> &ind_cod,
              const vector<array<char, 3>> &cod,
              const vector<int> &den,
              const vector<double> &deltaKm,
              const vector<double> &deltaVmax,
              const string &trna_data_path,
              const string &mrna_data_path,
              double init_cell_volume,
              int num_ribosomes);

    ~CellState();

    // Main method
    CellStepResult run_step(double volume, int num_steps);
};

#endif // CELL_STATE_HPP