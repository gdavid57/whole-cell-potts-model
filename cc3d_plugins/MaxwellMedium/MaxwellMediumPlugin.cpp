#include <CompuCell3D/CC3D.h>

using namespace CompuCell3D;

#include "MaxwellMediumPlugin.h"
#include <cmath>
#include <algorithm>

MaxwellMediumPlugin::MaxwellMediumPlugin() :
    pUtils(0),
    lockPtr(0),
    xmlData(0),
    cellFieldG(0),
    boundaryStrategy(0),
    N(0), sumX(0.0), sumY(0.0), sumZ(0.0),
    sumX2(0.0), sumY2(0.0), sumZ2(0.0),
    sumXY(0.0), sumXZ(0.0), sumYZ(0.0),
    lastRelaxationMCS(0),
    theta_snapshot(0.0),
    lambdaElastic(0.0), refArea(1.0), refVolume(1.0),
    tau_r(0.0), theta_relax(0.0), is2D(true),
    cx0(0.0), cy0(0.0), cz0(0.0), refCenterInitialized(false) {}

MaxwellMediumPlugin::~MaxwellMediumPlugin() {
    pUtils->destroyLock(lockPtr);
    delete lockPtr;
    lockPtr = 0;
}

void MaxwellMediumPlugin::init(Simulator *simulator, CC3DXMLElement *_xmlData) {
    xmlData = _xmlData;
    sim = simulator;
    potts = simulator->getPotts();
    cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();
    pUtils = sim->getParallelUtils();
    lockPtr = new ParallelUtilsOpenMP::OpenMPLock_t;
    pUtils->initLock(lockPtr);

    update(xmlData, true);
    potts->registerEnergyFunctionWithName(this, "MaxwellMedium");
    simulator->registerSteerableObject(this);
    potts->registerCellGChangeWatcher(this);

    // Initialise accumulators
    N = 0;
    sumX = sumY = sumZ = 0.0;
    sumX2 = sumY2 = sumZ2 = 0.0;
    sumXY = sumXZ = sumYZ = 0.0;
    
    Dim3D fieldDim = cellFieldG->getDim();
    is2D = (fieldDim.x == 1 || fieldDim.y == 1 || fieldDim.z == 1);
    if(is2D){
        refArea = (refArea <= 0.0) ? static_cast<double>(fieldDim.x * fieldDim.y * fieldDim.z) / std::max({(long)fieldDim.x, (long)fieldDim.y, (long)fieldDim.z}) : refArea;
    } else {
        refVolume = (refVolume <= 0.0) ? static_cast<double>(fieldDim.x * fieldDim.y * fieldDim.z) : refVolume;
    }

    for (int x = 0; x < fieldDim.x; ++x) {
        for (int y = 0; y < fieldDim.y; ++y) {
            for (int z = 0; z < fieldDim.z; ++z) {
                Point3D pt(x, y, z);
                CellG *cell = cellFieldG->get(pt);
                if (!cell) continue;
                Coordinates3D<double> coord = BoundaryStrategy::getInstance()->calculatePointCoordinates(pt);
                N += 1;
                sumX += coord.x;
                sumY += coord.y;
                sumZ += coord.z;
                sumX2 += coord.x * coord.x;
                sumY2 += coord.y * coord.y;
                sumZ2 += coord.z * coord.z;
                sumXY += coord.x * coord.y;
                sumXZ += coord.x * coord.z;
                sumYZ += coord.y * coord.z;
            }
        }
    }

    // Initialize relaxation state
    lastRelaxationMCS = sim->getStep();
    theta_snapshot = getThetaGlobal();
    // Fix reference center to initial center of mass
    ensureRefCenterInitialized();
}

// Calculate global strain from accumulated sums
double MaxwellMediumPlugin::computeThetaGlobal(long long n, double sx, double sy, double sz,
                                          double sxx, double syy, double szz,
                                          double sxy, double sxz, double syz) {
    if (n <= 1) return 0.0;

    // Moments with respect to the fixed initial center (cx0, cy0, cz0)
    ensureRefCenterInitialized();
    double invn = 1.0 / static_cast<double>(n);

    if (is2D) {
        double Gxx = (sxx * invn) - 2.0 * cx0 * (sx * invn) + cx0 * cx0;
        double Gyy = (syy * invn) - 2.0 * cy0 * (sy * invn) + cy0 * cy0;
        double Gxy = (sxy * invn) - cx0 * (sy * invn) - cy0 * (sx * invn) + cx0 * cy0;
        double trace = Gxx + Gyy;
        double diff = Gxx - Gyy;
        double term = std::sqrt(std::max(0.0, diff * diff + 4.0 * Gxy * Gxy));
        double lambda1 = 0.5 * (trace + term);
        double lambda2 = 0.5 * (trace - term);
        lambda1 = std::max(0.0, lambda1);
        lambda2 = std::max(0.0, lambda2);
        double R1 = std::sqrt(lambda1);
        double R2 = std::sqrt(lambda2);
        double areaGyr = 4.0 * M_PI * R1 * R2;
        return areaGyr / refArea;
    } else {
        double Gxx = (sxx * invn) - 2.0 * cx0 * (sx * invn) + cx0 * cx0;
        double Gyy = (syy * invn) - 2.0 * cy0 * (sy * invn) + cy0 * cy0;
        double Gzz = (szz * invn) - 2.0 * cz0 * (sz * invn) + cz0 * cz0;
        double Gxy = (sxy * invn) - cx0 * (sy * invn) - cy0 * (sx * invn) + cx0 * cy0;
        double Gxz = (sxz * invn) - cx0 * (sz * invn) - cz0 * (sx * invn) + cx0 * cz0;
        double Gyz = (syz * invn) - cy0 * (sz * invn) - cz0 * (sy * invn) + cy0 * cz0;
        double detG = Gxx * (Gyy * Gzz - Gyz * Gyz) - Gxy * (Gxy * Gzz - Gxz * Gyz) + Gxz * (Gxy * Gyz - Gyy * Gxz);
        if (detG < 0.0) detG = 0.0;
        double Rprod = std::sqrt(detG);
        double volumeGyr = (32.0 / 3.0) * M_PI * Rprod;
        return volumeGyr / refVolume;
    }
}

inline void MaxwellMediumPlugin::ensureRefCenterInitialized(){
    if (refCenterInitialized) return;
    if (N <= 0) return; // wait until cells exist to set the reference center
    cx0 = sumX / static_cast<double>(N);
    cy0 = sumY / static_cast<double>(N);
    cz0 = sumZ / static_cast<double>(N);
    refCenterInitialized = true;
}

// Calculate energy from strain
double MaxwellMediumPlugin::energyFromTheta(double theta_global) {
    double theta_e = theta_global - theta_relax;
    double ref = is2D ? refArea : refVolume;
    return 0.5 * lambdaElastic * ref * theta_e * theta_e;
}

inline double MaxwellMediumPlugin::_get_non_nan_energy(double energy) {
    if (energy != energy) return 0.0; // NaN check
    return energy;
}

double MaxwellMediumPlugin::changeEnergy(const Point3D &pt, const CellG *newCell, const CellG *oldCell) {
    if (newCell == oldCell) return 0.0;

    // Apply exact exponential relaxation once per MCS
    advanceRelaxationIfNewMCS();

    double theta_global_init = computeThetaGlobal(N, sumX, sumY, sumZ, sumX2, sumY2, sumZ2, sumXY, sumXZ, sumYZ);
    double Hinit = energyFromTheta(theta_global_init);

    long long nTrial = N;
    double sxTrial = sumX, syTrial = sumY, szTrial = sumZ;
    double sxxTrial = sumX2, syyTrial = sumY2, szzTrial = sumZ2;
    double sxyTrial = sumXY, sxzTrial = sumXZ, syzTrial = sumYZ;

    Coordinates3D<double> coord = boundaryStrategy->calculatePointCoordinates(pt);

    if (oldCell) {
        nTrial -= 1;
        sxTrial -= coord.x; syTrial -= coord.y; szTrial -= coord.z;
        sxxTrial -= coord.x * coord.x; syyTrial -= coord.y * coord.y; szzTrial -= coord.z * coord.z;
        sxyTrial -= coord.x * coord.y; sxzTrial -= coord.x * coord.z; syzTrial -= coord.y * coord.z;
    }
    if (newCell) {
        nTrial += 1;
        sxTrial += coord.x; syTrial += coord.y; szTrial += coord.z;
        sxxTrial += coord.x * coord.x; syyTrial += coord.y * coord.y; szzTrial += coord.z * coord.z;
        sxyTrial += coord.x * coord.y; sxzTrial += coord.x * coord.z; syzTrial += coord.y * coord.z;
    }

    double theta_global_final = computeThetaGlobal(nTrial, sxTrial, syTrial, szTrial, sxxTrial, syyTrial, szzTrial, sxyTrial, sxzTrial, syzTrial);
    double Hfinal = energyFromTheta(theta_global_final);

    return _get_non_nan_energy(Hfinal - Hinit);
}

void MaxwellMediumPlugin::update(CC3DXMLElement *_xmlData, bool _fullInitFlag) {
    automaton = potts->getAutomaton();
    ASSERT_OR_THROW("CELL TYPE PLUGIN WAS NOT PROPERLY INITIALIZED", automaton);

    if (_xmlData) {
        if (_xmlData->getFirstElement("LambdaElastic"))
            lambdaElastic = _xmlData->getFirstElement("LambdaElastic")->getDouble();
        if (_xmlData->getFirstElement("ReferenceArea"))
            refArea = _xmlData->getFirstElement("ReferenceArea")->getDouble();
        if (_xmlData->getFirstElement("ReferenceVolume"))
            refVolume = _xmlData->getFirstElement("ReferenceVolume")->getDouble();
        if (_xmlData->getFirstElement("TauR"))
            tau_r = _xmlData->getFirstElement("TauR")->getDouble();
    }

    Dim3D fieldDimLocal = potts->getCellFieldG()->getDim();
    is2D = (fieldDimLocal.x == 1 || fieldDimLocal.y == 1 || fieldDimLocal.z == 1);
    if (is2D) {
        if (refArea <= 0.0)
            refArea = static_cast<double>(fieldDimLocal.x * fieldDimLocal.y * fieldDimLocal.z) /
                      static_cast<double>(std::max({(long)fieldDimLocal.x, (long)fieldDimLocal.y, (long)fieldDimLocal.z}));
    } else {
        if (refVolume <= 0.0)
            refVolume = static_cast<double>(fieldDimLocal.x * fieldDimLocal.y * fieldDimLocal.z);
    }

    boundaryStrategy = BoundaryStrategy::getInstance();
}

std::string MaxwellMediumPlugin::toString() {
    return "MaxwellMedium";
}

std::string MaxwellMediumPlugin::steerableName() {
    return toString();
}

void MaxwellMediumPlugin::field3DChange(const Point3D &pt, CellG *newCell, CellG *oldCell) {
    Coordinates3D<double> coord = boundaryStrategy->calculatePointCoordinates(pt);

    if (oldCell) {
        --N;
        sumX -= coord.x; sumY -= coord.y; sumZ -= coord.z;
        sumX2 -= coord.x * coord.x; sumY2 -= coord.y * coord.y; sumZ2 -= coord.z * coord.z;
        sumXY -= coord.x * coord.y; sumXZ -= coord.x * coord.z; sumYZ -= coord.y * coord.z;
    }

    if (newCell) {
        ++N;
        sumX += coord.x; sumY += coord.y; sumZ += coord.z;
        sumX2 += coord.x * coord.x; sumY2 += coord.y * coord.y; sumZ2 += coord.z * coord.z;
        sumXY += coord.x * coord.y; sumXZ += coord.x * coord.z; sumYZ += coord.y * coord.z;
    }
    
    
}

double MaxwellMediumPlugin::getThetaGlobal() {
    return computeThetaGlobal(N, sumX, sumY, sumZ, sumX2, sumY2, sumZ2, sumXY, sumXZ, sumYZ);
}

void MaxwellMediumPlugin::advanceRelaxationIfNewMCS() {
    if (tau_r <= 0.0) return;
    long long currentMCS = sim->getStep();
    if (currentMCS <= lastRelaxationMCS) return;

    long long deltaMCS = currentMCS - lastRelaxationMCS;
    double dt = static_cast<double>(deltaMCS); // assume 1 time unit per MCS

    // Closed-form solution with piecewise-constant theta on [t_k, t_{k+1})
    // theta_relax(t+dt) = e^{-dt/tau} theta_relax(t) + (1 - e^{-dt/tau}) theta(t)
    double theta_now = theta_snapshot; // value from previous MCS snapshot
    double factor = std::exp(-dt / tau_r);
    theta_relax = factor * theta_relax + (1.0 - factor) * theta_now;

    // Prepare for next interval: snapshot current theta
    theta_snapshot = getThetaGlobal();
    lastRelaxationMCS = currentMCS;
}
