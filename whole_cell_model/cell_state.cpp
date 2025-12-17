#include "cell_state.hpp"

// For generating random integers and doubles:
std::random_device rd;                          // Non-deterministic random seed
std::mt19937 gen(rd());                         // Mersenne Twister random number generator
std::uniform_real_distribution<> dis(0.0, 1.0); // For r1 and r2 of gillespie simulation

// Constructor
CellState::CellState(const string &tRNA_data_path, const string &mRNA_data_path, double init_cell_volume, int num_ribosomes)
    : trna_data_path(tRNA_data_path), mrna_data_path(mRNA_data_path),
      init_cell_volume(init_cell_volume), num_ribosomes(num_ribosomes)
{
    // Configuration parameters
    beta = 1e07;
    gamma = pow(10, -5);
    init_charging = 0.5;
    radius = 22e-9;

    // Simulation constants
    NT = 40;
    r = 1.;
    transc_rate = 1.;
    Ncod = 61;
    Naa = 20;

    t = 0.0;
    t0 = 0;
    current_volume = 0.0;

    init(tRNA_data_path, mRNA_data_path, init_cell_volume, num_ribosomes);
}

// Full state constructor for pickle/deserialization
CellState::CellState(double beta, double gamma, double init_charging, double radius,
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
                     int num_ribosomes)
    : beta(beta), gamma(gamma), init_charging(init_charging), radius(radius),
      NT(NT), r(r), transc_rate(transc_rate), Ncod(Ncod), Naa(Naa),
      M(M), Nrib(Nrib), t(t), t0(t0), current_volume(current_volume),
      TRNA(TRNA), a(a), rib(rib), mRNA(mRNA), rho(rho), rho_sites(rho_sites),
      Vmax(Vmax), Km(Km), chi(chi), Vt(Vt), init_rates(init_rates), alpha(alpha),
      proteins(proteins), N(N), tnum(tnum), ind_cod(ind_cod), cod(cod), den(den),
      deltaKm(deltaKm), deltaVmax(deltaVmax), trna_data_path(trna_data_path),
      mrna_data_path(mrna_data_path), init_cell_volume(init_cell_volume),
      num_ribosomes(num_ribosomes)
{
    cout << "CellState restored from full state with M = " << M << " mRNAs and " << Nrib << " free ribosomes" << endl;
}

// Destructor
CellState::~CellState()
{
    // vectors automatically handle deallocation
    cout << "CellState() destruction" << endl;
}

// compute_propensity method
double CellState::compute_propensity(int j, int i) const
{
    if (i == -100)
        return alpha[j]; // Return alpha for initiation
    else if (i == -200)
        return beta; // Return beta for termination
    else if (i == -300)
        return 0.;
    else if (i >= 0)
        return TRNA[ind_cod[0][i]].k[ind_cod[1][i]]; // Return k for elongation
    else if ((i < 0) && (i >= -NT))
        return chi[-i - 1]; // Return chi for recharge
    else
        cerr << "i" << '\t' << i << "Problem with propensities" << endl;
    return 0.;
}

// init method
void CellState::init(const string &tRNA_data_path, const string &mRNA_data_path, double init_cell_volume, int num_ribosomes)
{
    // Initial number of ribosomes and storing initial volume
    Nrib = num_ribosomes;
    current_volume = init_cell_volume;

    // Dynamic array allocation
    ind_cod.resize(2);
    cod.resize(61);
    den.resize(Naa);
    deltaKm.resize(NT);
    deltaVmax.resize(NT);

    // Initialize fine-tuning factors
    double default_deltaKm[40] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double default_deltaVmax[40] = {4.5, 4, 0.007, 1, 1, 0.0005, 0.1, 25, 0.005, 0.005, 0.01, 1, 0.006, 0.001, 0.02, 0.008, 0.065, 0.009, 1, 1, 1, 0.01, 0.03, 1, 1, 0.01, 15, 1.7, 1.5, 0.0001, 0.85, 0.85, 7, 6.5, 0.65, 0.08, 0.05, 26, 1, 0.01};

    for (int i = 0; i < NT; i++)
    {
        deltaKm[i] = default_deltaKm[i];
        deltaVmax[i] = default_deltaVmax[i];
    }

    // Reading number of mRNAs
    if (!read_mrna_count(mRNA_data_path))
        return;

    // Reading tRNA abundances
    if (!read_trna_abundances(tRNA_data_path))
        return;

    // Memory allocation for all structures
    N.resize(M);
    init_rates.resize(M);
    alpha.resize(M);
    proteins.resize(M);
    Vmax.resize(NT);
    Km.resize(NT);
    chi.resize(NT);
    Vt.resize(NT);

    TRNA.resize(NT);
    a.resize(M + 1);
    rib.resize(M);
    mRNA.resize(M);
    rho.resize(M);
    rho_sites.resize(M);

    // Initialize protein counters
    for (int i = 0; i < M; i++)
    {
        proteins[i] = 0;
    }

    // Number of ribosomes already initialized in parameters

    // Reading mRNA sizes
    if (!read_mrna_sizes(mRNA_data_path))
        return;

    // Reading initiation rates
    if (!read_initiation_rates(mRNA_data_path))
        return;

    // Reading Vmax/Km parameters
    if (!read_vmax_km_parameters(tRNA_data_path))
        return;

    // Reading codon-tRNA index table
    if (!read_codon_trna_index_table(tRNA_data_path))
        return;

    // Memory reservation for tRNA vectors
    for (int j = 0; j < NT; j++)
    {
        TRNA[j].index.reserve(TRNA[j].len);
        TRNA[j].wobble.reserve(TRNA[j].len);
        TRNA[j].k.reserve(TRNA[j].len);
        TRNA[j].dem.reserve(TRNA[j].len);

        TRNA[j].index.resize(TRNA[j].len);
        TRNA[j].wobble.resize(TRNA[j].len);
        TRNA[j].k.resize(TRNA[j].len);
        TRNA[j].dem.resize(TRNA[j].len);
    }

    // Reading indices and wobble factors
    if (!read_wobble_factors(tRNA_data_path))
        return;

    // Calculate initiation rates adjusted with provided volume
    for (int j = 0; j < M; j++)
    {
        init_rates[j] = (gamma * init_rates[j]) * M_PI * radius * radius * radius * (4 / 3) * (1 / init_cell_volume);
        alpha[j] = init_rates[j] * Nrib;
    }

    // Initialize tRNA state
    for (int j = 0; j < NT; j++)
    {
        TRNA[j].tnum = transc_rate * tnum[j];
        TRNA[j].cnum = ceil(init_charging * TRNA[j].tnum);

        // Calculate kinetic constants k
        const double cnum_factor = r * TRNA[j].cnum;
        for (int i = 0; i < TRNA[j].len; i++)
        {
            TRNA[j].k[i] = (1.0 - TRNA[j].wobble[i]) * cnum_factor;
            TRNA[j].dem[i] = 0;
        }
    }

    // Calculate denominators for Michaelis-Menten
    for (int j = 0; j < Naa; j++)
    {
        den[j] = 0;
    }
    for (int j = 0; j < NT; j++)
    {
        den[TRNA[j].aa] = den[TRNA[j].aa] + TRNA[j].tnum - TRNA[j].cnum;
    }

    // Calculate tRNA recharge rates (chi)
    for (int j = 0; j < NT; j++)
    {
        chi[j] = Vmax[j] * (TRNA[j].tnum - TRNA[j].cnum) / (Km[j] + den[TRNA[j].aa]);
    }

    // Initialize vectors for mRNAs
    for (int i = 0; i < M; i++)
    {
        rib[i].clear();
        rho[i].clear();
        rho_sites[i].clear();
        a[i].clear();
        mRNA[i].clear();

        // Reserve capacity
        rib[i].reserve(N[i] / 10);
        a[i].reserve(N[i] / 10 + 2);
        mRNA[i].reserve(N[i]);
        rho[i].reserve(N[i]);
        rho_sites[i].reserve(N[i]);
    }
    a[M].clear();
    a[M].reserve(NT);

    // Reading mRNA sequences
    if (!read_mrna_sequences(mRNA_data_path))
        return;

    // Initialize ribosomes on mRNAs
    initialize_ribosomes();

    // Initialize action lists
    for (int i = 0; i < M; i++)
    {
        a[i].push_back(-100); // initiation action
    }
    for (int j = 0; j < NT; j++)
    {
        a[M].push_back(-j - 1); // recharge actions
    }

    // Initialize densities
    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N[j]; i++)
        {
            rho[j].push_back(0.);
            rho_sites[j].push_back(0.);
        }
    }

    cout << "CellState initialized with M = " << M << " mRNAs and " << Nrib << " free ribosomes" << endl;
}

// Method to initialize ribosomes
void CellState::initialize_ribosomes()
{
    for (int i = 0; i < M; i++)
    {
        for (int j = N[i] - 6; j > 10; j = j - 38)
        {
            if (Nrib > 2)
            {
                rib[i].push_back(j);
                a[i].push_back(mRNA[i][j]);
                Nrib -= 1;
            }
            if (Nrib <= 1)
            {
                cout << "Warning: Low number of free ribosomes during initialization" << endl;
                break;
            }
        }
    }
}

// File reading methods
bool CellState::read_mrna_count(const string &mRNA_data_path)
{
    ifstream file_M(mRNA_data_path + "/C_numberOfMRNAs.dat");
    if (file_M.is_open())
    {
        file_M >> M;
        file_M.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << mRNA_data_path + "/C_numberOfMRNAs.dat" << endl;
        return false;
    }
}

bool CellState::read_trna_abundances(const string &tRNA_data_path)
{
    tnum.resize(NT);
    ifstream file_tnum(tRNA_data_path + "/abundances_LB.dat");
    if (file_tnum.is_open())
    {
        for (int j = 0; j < NT; j++)
        {
            file_tnum >> tnum[j];
        }
        file_tnum.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << tRNA_data_path + "/abundances_LB.dat" << endl;
        return false;
    }
}

bool CellState::read_mrna_sizes(const string &mRNA_data_path)
{
    ifstream file_sizes(mRNA_data_path + "/sizes_mRNAs.dat");
    if (file_sizes.is_open())
    {
        for (int j = 0; j < M; j++)
        {
            file_sizes >> N[j];
        }
        file_sizes.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << mRNA_data_path + "/sizes_mRNAs.dat" << endl;
        return false;
    }
}

bool CellState::read_initiation_rates(const string &mRNA_data_path)
{
    ifstream file_init(mRNA_data_path + "/initiation_rates.dat");
    if (file_init.is_open())
    {
        for (int j = 0; j < M; j++)
        {
            file_init >> init_rates[j];
        }
        file_init.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << mRNA_data_path + "/initiation_rates.dat" << endl;
        return false;
    }
}

bool CellState::read_vmax_km_parameters(const string &tRNA_data_path)
{
    ifstream file_vmax(tRNA_data_path + "/C_Vmax_Km_Ecoli_format_LB.dat");
    if (file_vmax.is_open())
    {
        for (int j = 0; j < NT; j++)
        {
            file_vmax >> Vt[j] >> Km[j] >> TRNA[j].aa >> TRNA[j].len;
            Km[j] = Km[j] * deltaKm[j];
            Vmax[j] = Vt[j] * deltaVmax[j];
        }
        file_vmax.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << tRNA_data_path + "/C_Vmax_Km_Ecoli_format_LB.dat" << endl;
        return false;
    }
}

bool CellState::read_codon_trna_index_table(const string &tRNA_data_path)
{
    ifstream file_codons(tRNA_data_path + "/C_Codons_indices_table_Ecoli.dat");
    if (file_codons.is_open())
    {
        char dummy[50];
        for (int j = 0; j < Ncod; j++)
        {
            file_codons >> dummy >> dummy >> cod[j][0] >> cod[j][1] >> cod[j][2] >> ind_cod[0][j] >> ind_cod[1][j];
        }
        file_codons.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << tRNA_data_path + "/C_Codons_indices_table_Ecoli.dat" << endl;
        return false;
    }
}

bool CellState::read_wobble_factors(const string &tRNA_data_path)
{
    ifstream file_wobble(tRNA_data_path + "/ind_wobble.dat");
    if (file_wobble.is_open())
    {
        for (int j = 0; j < NT; j++)
        {
            for (int i = 0; i < TRNA[j].len; i++)
            {
                file_wobble >> TRNA[j].index[i] >> TRNA[j].wobble[i];
            }
        }
        file_wobble.close();
        return true;
    }
    else
    {
        cerr << "Can't open file " << tRNA_data_path + "/ind_wobble.dat" << endl;
        return false;
    }
}

bool CellState::read_mrna_sequences(const string &mRNA_data_path)
{
    ifstream file_seq(mRNA_data_path + "/mRNAs_sequences.dat");
    string line;
    int i = -1;
    if (file_seq)
    {
        while (getline(file_seq, line) && (i < M - 1))
        {
            i++;

            for (int tt = 0; tt < N[i] - 1; tt++)
            {
                string str0 = string(1, line.at(3 * tt));
                string str1 = string(1, line.at(3 * tt + 1));
                string str2 = string(1, line.at(3 * tt + 2));
                string str = str0 + str1 + str2;

                for (int j = 0; j < Ncod; j++)
                {
                    string cod0 = string(1, cod[j][0]);
                    string cod1 = string(1, cod[j][1]);
                    string cod2 = string(1, cod[j][2]);
                    string cods = cod0 + cod1 + cod2;

                    if (cods.compare(str) == 0)
                    {
                        mRNA[i].push_back(j);
                        break;
                    }
                }
            }
            mRNA[i].push_back(-200); // termination codon
        }
        file_seq.close();
        return true;
    }
    else
    {
        cerr << "Can't open mRNA sequences file: " << mRNA_data_path + "/mRNAs_sequences.dat" << endl;
        return false;
    }
}

// Method to update kinetic parameters
void CellState::update_kinetic_parameters(double volume)
{
    // Update initiation rates according to volume
    // init_rates already contains the factor (1/initial_volume), so we apply the volume ratio
    for (int i = 0; i < M; i++)
    {
        alpha[i] = init_rates[i] * (current_volume / volume) * Nrib;
    }

    // Update current volume
    current_volume = volume;

    // Recalculate denominators for Michaelis-Menten
    for (int j = 0; j < Naa; j++)
        den[j] = 0;
    for (int j = 0; j < NT; j++)
        den[TRNA[j].aa] = den[TRNA[j].aa] + TRNA[j].tnum - TRNA[j].cnum;

    // Update recharge rates and kinetic constants
    for (int j = 0; j < NT; j++)
    {
        chi[j] = Vmax[j] * (TRNA[j].tnum - TRNA[j].cnum) / (Km[j] + den[TRNA[j].aa]);

        const double cnum_factor = r * TRNA[j].cnum;
        for (int i = 0; i < TRNA[j].len; i++)
        {
            TRNA[j].k[i] = (1.0 - TRNA[j].wobble[i]) * cnum_factor;
        }
    }
}

// run_step method
CellStepResult CellState::run_step(double volume, int num_steps)
{
    int initial_proteins = 0;
    for (int i = 0; i < M; i++)
        initial_proteins += proteins[i];

    // Main Gillespie loop
    for (int step = 0; step < num_steps; step++)
    {
        t0++;

        // Update kinetic parameters
        update_kinetic_parameters(volume);

        // Calculate sum of propensities (a0)
        double a0 = 0.0;
        for (int i = 0; i < M + 1; i++)
        {
            const auto &a_ref = a[i];
            const int a_size = static_cast<int>(a_ref.size());
            for (int j = 0; j < a_size; j++)
            {
                a0 += compute_propensity(i, a_ref[j]);
            }
        }

        // Check tRNA consistency
        for (int i = 0; i < NT; i++)
        {
            if (TRNA[i].cnum < 0)
            {
                cerr << "ERROR: tRNA discharged: " << i << " in step=" << t0 << endl;
            }
        }

        // Generate random numbers
        double r1 = dis(gen);
        double r2 = dis(gen);

        // Calculate time interval
        double tau = 1 / a0 * (log(1 / r1));
        t = t + tau;

        // Select next reaction
        double test = r2 * a0;
        int i = 0;
        int j = 0;

        double suma = compute_propensity(0, a[0][0]);

        while (test > suma)
        {
            const int a_j_size = static_cast<int>(a[j].size());
            if (i < a_j_size - 1)
            {
                i++;
                suma = suma + compute_propensity(j, a[j][i]);
            }
            else
            {
                j++;
                i = 0;
                if (j > M)
                {
                    cerr << "Warning: " << "\t" << "t0 == " << t0 << endl;
                }
                suma = suma + compute_propensity(j, a[j][i]);
            }
        }

        int mu2 = j;
        int mu1 = i;

        // Execute selected reaction
        if (a0 == 0)
        {
            cout << "No possible actions!" << endl;
        }
        else
        {
            if (mu2 < M)
            { // Reaction on mRNA
                // Complete implementation of ribosome reactions
                // [Code identical to original for all reactions]

                if (rib[mu2].size() <= 1)
                {
                    if (mu1 == 0)
                    {
                        if (rib[mu2].size() == 0)
                        { // INITIATION
                            rib[mu2].push_back(0);
                            int position = rib[mu2][mu1];
                            int codon = mRNA[mu2][position];
                            a[mu2][mu1] = codon;
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];

                            TRNA[index].dem[sub_index]++;
                            a[mu2].push_back(-300);
                            Nrib -= 1;
                        }
                        else if (rib[mu2][0] == 9)
                        { // Special ELONGATION 9->10
                            int position = rib[mu2][0];
                            int codon = mRNA[mu2][position];
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];

                            TRNA[index].cnum--;
                            TRNA[index].dem[sub_index]--;

                            rib[mu2][mu1] = 10;
                            position = rib[mu2][mu1];
                            codon = mRNA[mu2][position];

                            index = ind_cod[0][codon];
                            sub_index = ind_cod[1][codon];

                            TRNA[index].dem[sub_index]++;
                            a[mu2][0] = codon;
                            a[mu2][1] = -100;
                        }
                        else if (rib[mu2][0] == N[mu2] - 1)
                        { // TERMINATION
                            rib[mu2].erase(rib[mu2].begin());
                            a[mu2].erase(a[mu2].begin());
                            proteins[mu2]++;
                            Nrib += 1;
                        }
                        else if (rib[mu2][0] == N[mu2] - 2)
                        {
                            int position = rib[mu2][0];
                            int codon = mRNA[mu2][position];
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];

                            TRNA[index].cnum--;
                            TRNA[index].dem[sub_index]--;

                            rib[mu2][0] = N[mu2] - 1;
                            a[mu2][0] = -200;
                        }
                        else
                        { // Standard ELONGATION
                            int position = rib[mu2][0];
                            int codon = mRNA[mu2][position];
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];

                            TRNA[index].cnum--;
                            TRNA[index].dem[sub_index]--;

                            rib[mu2][0]++;
                            position = rib[mu2][0];
                            codon = mRNA[mu2][position];
                            index = ind_cod[0][codon];
                            sub_index = ind_cod[1][codon];
                            TRNA[index].dem[sub_index]++;
                            a[mu2][0] = codon;
                        }
                    }
                    else
                    { // INITIATION with ribosome already present
                        rib[mu2].push_back(0);
                        int position = rib[mu2][mu1];
                        int codon = mRNA[mu2][position];
                        int index = ind_cod[0][codon];
                        int sub_index = ind_cod[1][codon];

                        Nrib -= 1;

                        if (rib[mu2][0] > 10)
                        {
                            TRNA[index].dem[sub_index]++;
                            a[mu2][1] = codon;
                        }
                        else
                        {
                            a[mu2][1] = -300;
                        }
                        a[mu2].push_back(-300);
                    }
                }
                else
                { // Case with more than one ribosome
                    if (mu1 == 0)
                    { // Rightmost RIBOSOME
                        if (rib[mu2][0] == N[mu2] - 1)
                        { // TERMINATION
                            Nrib += 1;

                            if (rib[mu2][1] == N[mu2] - 11)
                            {
                                int position = rib[mu2][1];
                                int codon = mRNA[mu2][position];
                                int index = ind_cod[0][codon];
                                int sub_index = ind_cod[1][codon];
                                a[mu2][1] = codon;
                                TRNA[index].dem[sub_index]++;
                            }
                            rib[mu2].erase(rib[mu2].begin());
                            a[mu2].erase(a[mu2].begin());
                            proteins[mu2]++;
                        }
                        else if (rib[mu2][0] == N[mu2] - 2)
                        {
                            int position = rib[mu2][0];
                            int codon = mRNA[mu2][position];
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];

                            TRNA[index].cnum--;
                            TRNA[index].dem[sub_index]--;

                            rib[mu2][0]++;
                            a[mu2][0] = -200;

                            if (rib[mu2][1] == N[mu2] - 12)
                            {
                                position = rib[mu2][1];
                                codon = mRNA[mu2][position];
                                index = ind_cod[0][codon];
                                sub_index = ind_cod[1][codon];
                                a[mu2][1] = codon;
                                TRNA[index].dem[sub_index]++;
                            }
                        }
                        else
                        { // Standard ELONGATION
                            int position = rib[mu2][0];
                            int codon = mRNA[mu2][position];
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];

                            TRNA[index].cnum--;
                            TRNA[index].dem[sub_index]--;

                            rib[mu2][0]++;
                            position = rib[mu2][0];
                            codon = mRNA[mu2][position];
                            index = ind_cod[0][codon];
                            sub_index = ind_cod[1][codon];
                            a[mu2][0] = codon;
                            TRNA[index].dem[sub_index]++;

                            if (rib[mu2][1] == (rib[mu2][0] - 11))
                            {
                                position = rib[mu2][1];
                                codon = mRNA[mu2][position];
                                index = ind_cod[0][codon];
                                sub_index = ind_cod[1][codon];
                                a[mu2][1] = codon;
                                TRNA[index].dem[sub_index]++;
                            }
                        }
                    }
                    else if (mu1 == static_cast<int>(a[mu2].size()) - 1)
                    { // Leftmost RIBOSOME - INITIATION
                        Nrib -= 1;
                        rib[mu2].push_back(0);
                        if (rib[mu2][rib[mu2].size() - 2] > 10)
                        {
                            int position = rib[mu2][rib[mu2].size() - 1];
                            int codon = mRNA[mu2][position];
                            int index = ind_cod[0][codon];
                            int sub_index = ind_cod[1][codon];
                            a[mu2][mu1] = codon;
                            TRNA[index].dem[sub_index]++;
                        }
                        else
                        {
                            a[mu2][mu1] = -300;
                        }
                        a[mu2].push_back(-300);
                    }
                    else if (mu1 == static_cast<int>(a[mu2].size()) - 2 && rib[mu2][static_cast<int>(rib[mu2].size()) - 1] == 9)
                    {
                        // Ribosome hop from site 9 to site 10
                        int position = rib[mu2][rib[mu2].size() - 1];
                        int codon = mRNA[mu2][position];
                        int index = ind_cod[0][codon];
                        int sub_index = ind_cod[1][codon];

                        TRNA[index].cnum--;
                        TRNA[index].dem[sub_index]--;

                        rib[mu2][rib[mu2].size() - 1]++;
                        a[mu2][a[mu2].size() - 1] = -100;

                        if (rib[mu2][rib[mu2].size() - 2] > 20)
                        {
                            position = rib[mu2][rib[mu2].size() - 1];
                            codon = mRNA[mu2][position];
                            index = ind_cod[0][codon];
                            sub_index = ind_cod[1][codon];
                            TRNA[index].dem[sub_index]++;
                            a[mu2][mu1] = codon;
                        }
                        else
                        {
                            a[mu2][mu1] = -300;
                        }
                    }
                    else
                    { // RIBOSOME in the middle - ELONGATION
                        int position = rib[mu2][mu1];
                        int codon = mRNA[mu2][position];
                        int index = ind_cod[0][codon];
                        int sub_index = ind_cod[1][codon];

                        TRNA[index].cnum--;
                        TRNA[index].dem[sub_index]--;

                        rib[mu2][mu1]++;

                        if (rib[mu2][mu1 - 1] != rib[mu2][mu1] + 10)
                        {
                            position = rib[mu2][mu1];
                            codon = mRNA[mu2][position];
                            index = ind_cod[0][codon];
                            sub_index = ind_cod[1][codon];
                            a[mu2][mu1] = codon;
                            TRNA[index].dem[sub_index]++;
                        }
                        else
                        {
                            a[mu2][mu1] = -300;
                        }

                        if (mu1 < static_cast<int>(rib[mu2].size()) - 1)
                        {
                            if (rib[mu2][mu1 + 1] == rib[mu2][mu1] - 11)
                            {
                                position = rib[mu2][mu1 + 1];
                                codon = mRNA[mu2][position];
                                index = ind_cod[0][codon];
                                sub_index = ind_cod[1][codon];
                                a[mu2][mu1 + 1] = codon;
                                TRNA[index].dem[sub_index]++;
                            }
                        }
                    }
                }
            }
            else
            { // Recharge reaction
                int index = -a[mu2][mu1] - 1;
                TRNA[index].cnum++;
            }
        }
    }

    // Calculate total number of proteins produced during this step
    int final_proteins = 0;
    for (int i = 0; i < M; i++)
        final_proteins += proteins[i];
    int produced_proteins = final_proteins - initial_proteins;

    CellStepResult result;
    result.produced_proteins = produced_proteins;
    result.free_ribosomes = Nrib;
    result.physical_time = t;

    return result;
}