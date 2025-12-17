#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "cell_state.hpp"

namespace py = pybind11;

PYBIND11_MODULE(py_cell_state, m)
{
    m.doc() = "Python bindings for CellState whole-cell model library";

    // Bind the CellStepResult struct
    py::class_<CellStepResult>(m, "CellStepResult")
        .def_readwrite("produced_proteins", &CellStepResult::produced_proteins)
        .def_readwrite("free_ribosomes", &CellStepResult::free_ribosomes)
        .def_readwrite("physical_time", &CellStepResult::physical_time)
        .def("__repr__", [](const CellStepResult &r)
             { return "CellStepResult(produced_proteins=" + std::to_string(r.produced_proteins) +
                      ", free_ribosomes=" + std::to_string(r.free_ribosomes) +
                      ", physical_time=" + std::to_string(r.physical_time) + ")"; })
        .def(py::pickle(
            [](const CellStepResult &r)
            {
                return py::make_tuple(r.produced_proteins, r.free_ribosomes, r.physical_time);
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state!");
                CellStepResult r;
                r.produced_proteins = t[0].cast<int>();
                r.free_ribosomes = t[1].cast<int>();
                r.physical_time = t[2].cast<double>();
                return r;
            }));

    // Bind the tList struct first
    py::class_<tList>(m, "tList")
        .def_readwrite("index", &tList::index)
        .def_readwrite("cnum", &tList::cnum)
        .def_readwrite("tnum", &tList::tnum)
        .def_readwrite("len", &tList::len)
        .def_readwrite("aa", &tList::aa)
        .def_readwrite("dem", &tList::dem)
        .def_readwrite("k", &tList::k)
        .def_readwrite("wobble", &tList::wobble)
        .def(py::pickle(
            [](const tList &t)
            {
                return py::make_tuple(t.index, t.cnum, t.tnum, t.len, t.aa, t.dem, t.k, t.wobble);
            },
            [](py::tuple tuple)
            {
                if (tuple.size() != 8)
                    throw std::runtime_error("Invalid state!");
                tList t;
                t.index = tuple[0].cast<std::vector<int>>();
                t.cnum = tuple[1].cast<int>();
                t.tnum = tuple[2].cast<int>();
                t.len = tuple[3].cast<int>();
                t.aa = tuple[4].cast<int>();
                t.dem = tuple[5].cast<std::vector<int>>();
                t.k = tuple[6].cast<std::vector<double>>();
                t.wobble = tuple[7].cast<std::vector<double>>();
                return t;
            }));

    // Bind the main CellState class
    py::class_<CellState>(m, "CellState")
        .def(py::init<const std::string &, const std::string &, double, int>(),
             py::arg("trna_data_path"), py::arg("mrna_data_path"),
             py::arg("init_cell_volume"), py::arg("num_ribosomes"),
             "Initialize CellState with data paths, initial cell volume, and number of ribosomes")
        .def(py::init<double, double, double, double,
                      int, double, double, int, int,
                      int, int, double, int, double,
                      const std::vector<tList> &,
                      const std::vector<std::vector<int>> &,
                      const std::vector<std::vector<int>> &,
                      const std::vector<std::vector<int>> &,
                      const std::vector<std::vector<double>> &,
                      const std::vector<std::vector<double>> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::vector<int> &,
                      const std::vector<int> &,
                      const std::vector<int> &,
                      const std::vector<std::array<int, 61>> &,
                      const std::vector<std::array<char, 3>> &,
                      const std::vector<int> &,
                      const std::vector<double> &,
                      const std::vector<double> &,
                      const std::string &,
                      const std::string &,
                      double,
                      int>(),
             py::arg("beta"), py::arg("gamma"), py::arg("init_charging"), py::arg("radius"),
             py::arg("NT"), py::arg("r"), py::arg("transc_rate"), py::arg("Ncod"), py::arg("Naa"),
             py::arg("M"), py::arg("Nrib"), py::arg("t"), py::arg("t0"), py::arg("current_volume"),
             py::arg("TRNA"), py::arg("a"), py::arg("rib"), py::arg("mRNA"),
             py::arg("rho"), py::arg("rho_sites"), py::arg("Vmax"), py::arg("Km"),
             py::arg("chi"), py::arg("Vt"), py::arg("init_rates"), py::arg("alpha"),
             py::arg("proteins"), py::arg("N"), py::arg("tnum"), py::arg("ind_cod"),
             py::arg("cod"), py::arg("den"), py::arg("deltaKm"), py::arg("deltaVmax"),
             py::arg("trna_data_path"), py::arg("mrna_data_path"),
             py::arg("init_cell_volume"), py::arg("num_ribosomes"),
             "Initialize CellState from full state (for deserialization)")
        .def("run_step", &CellState::run_step,
             py::arg("volume"), py::arg("num_steps"),
             "Run simulation step with given volume and number of steps")
        // Expose all public attributes for serialization
        .def_readwrite("beta", &CellState::beta)
        .def_readwrite("gamma", &CellState::gamma)
        .def_readwrite("init_charging", &CellState::init_charging)
        .def_readwrite("radius", &CellState::radius)
        .def_readwrite("NT", &CellState::NT)
        .def_readwrite("r", &CellState::r)
        .def_readwrite("transc_rate", &CellState::transc_rate)
        .def_readwrite("Ncod", &CellState::Ncod)
        .def_readwrite("Naa", &CellState::Naa)
        .def_readwrite("M", &CellState::M)
        .def_readwrite("Nrib", &CellState::Nrib)
        .def_readwrite("t", &CellState::t)
        .def_readwrite("t0", &CellState::t0)
        .def_readwrite("current_volume", &CellState::current_volume)
        .def_readwrite("TRNA", &CellState::TRNA)
        .def_readwrite("a", &CellState::a)
        .def_readwrite("rib", &CellState::rib)
        .def_readwrite("mRNA", &CellState::mRNA)
        .def_readwrite("rho", &CellState::rho)
        .def_readwrite("rho_sites", &CellState::rho_sites)
        .def_readwrite("Vmax", &CellState::Vmax)
        .def_readwrite("Km", &CellState::Km)
        .def_readwrite("chi", &CellState::chi)
        .def_readwrite("Vt", &CellState::Vt)
        .def_readwrite("init_rates", &CellState::init_rates)
        .def_readwrite("alpha", &CellState::alpha)
        .def_readwrite("proteins", &CellState::proteins)
        .def_readwrite("N", &CellState::N)
        .def_readwrite("tnum", &CellState::tnum)
        .def_readwrite("ind_cod", &CellState::ind_cod)
        .def_readwrite("cod", &CellState::cod)
        .def_readwrite("den", &CellState::den)
        .def_readwrite("deltaKm", &CellState::deltaKm)
        .def_readwrite("deltaVmax", &CellState::deltaVmax)
        .def_readwrite("trna_data_path", &CellState::trna_data_path)
        .def_readwrite("mrna_data_path", &CellState::mrna_data_path)
        .def_readwrite("init_cell_volume", &CellState::init_cell_volume)
        .def_readwrite("num_ribosomes", &CellState::num_ribosomes)
        .def(py::pickle(
            [](const CellState &cs)
            {
                return py::make_tuple(
                    cs.beta, cs.gamma, cs.init_charging, cs.radius,
                    cs.NT, cs.r, cs.transc_rate, cs.Ncod, cs.Naa,
                    cs.M, cs.Nrib, cs.t, cs.t0, cs.current_volume,
                    cs.TRNA, cs.a, cs.rib, cs.mRNA, cs.rho, cs.rho_sites,
                    cs.Vmax, cs.Km, cs.chi, cs.Vt, cs.init_rates, cs.alpha,
                    cs.proteins, cs.N, cs.tnum, cs.ind_cod, cs.cod, cs.den,
                    cs.deltaKm, cs.deltaVmax, cs.trna_data_path, cs.mrna_data_path,
                    cs.init_cell_volume, cs.num_ribosomes);
            },
            [](py::tuple t)
            {
                if (t.size() != 38)
                    throw std::runtime_error("Invalid state!");
                return CellState(
                    t[0].cast<double>(),                            // beta
                    t[1].cast<double>(),                            // gamma
                    t[2].cast<double>(),                            // init_charging
                    t[3].cast<double>(),                            // radius
                    t[4].cast<int>(),                               // NT
                    t[5].cast<double>(),                            // r
                    t[6].cast<double>(),                            // transc_rate
                    t[7].cast<int>(),                               // Ncod
                    t[8].cast<int>(),                               // Naa
                    t[9].cast<int>(),                               // M
                    t[10].cast<int>(),                              // Nrib
                    t[11].cast<double>(),                           // t
                    t[12].cast<int>(),                              // t0
                    t[13].cast<double>(),                           // current_volume
                    t[14].cast<std::vector<tList>>(),               // TRNA
                    t[15].cast<std::vector<std::vector<int>>>(),    // a
                    t[16].cast<std::vector<std::vector<int>>>(),    // rib
                    t[17].cast<std::vector<std::vector<int>>>(),    // mRNA
                    t[18].cast<std::vector<std::vector<double>>>(), // rho
                    t[19].cast<std::vector<std::vector<double>>>(), // rho_sites
                    t[20].cast<std::vector<double>>(),              // Vmax
                    t[21].cast<std::vector<double>>(),              // Km
                    t[22].cast<std::vector<double>>(),              // chi
                    t[23].cast<std::vector<double>>(),              // Vt
                    t[24].cast<std::vector<double>>(),              // init_rates
                    t[25].cast<std::vector<double>>(),              // alpha
                    t[26].cast<std::vector<int>>(),                 // proteins
                    t[27].cast<std::vector<int>>(),                 // N
                    t[28].cast<std::vector<int>>(),                 // tnum
                    t[29].cast<std::vector<std::array<int, 61>>>(), // ind_cod
                    t[30].cast<std::vector<std::array<char, 3>>>(), // cod
                    t[31].cast<std::vector<int>>(),                 // den
                    t[32].cast<std::vector<double>>(),              // deltaKm
                    t[33].cast<std::vector<double>>(),              // deltaVmax
                    t[34].cast<std::string>(),                      // trna_data_path
                    t[35].cast<std::string>(),                      // mrna_data_path
                    t[36].cast<double>(),                           // init_cell_volume
                    t[37].cast<int>()                               // num_ribosomes
                );
            }));
}