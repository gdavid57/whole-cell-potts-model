#ifndef MAXWELLMEDIUMPLUGIN_H
#define MAXWELLMEDIUMPLUGIN_H

#include <CompuCell3D/CC3D.h>

#include "MaxwellMediumDLLSpecifier.h"

class CC3DXMLElement;

namespace CompuCell3D {

    class Simulator;
    class Potts3D;
    class Automaton;
    class BoundaryStrategy;
    class ParallelUtilsOpenMP;

    template <class T> class Field3D;
    template <class T> class WatchableField3D;

    class MAXWELLMEDIUM_EXPORT MaxwellMediumPlugin : public Plugin, public EnergyFunction, public CellGChangeWatcher {

    private:
        // Accumulated sums for gyration tensor
        long long N;
        double sumX, sumY, sumZ;
        double sumX2, sumY2, sumZ2;
        double sumXY, sumXZ, sumYZ;

        // Relaxation state (per-MCS exact update)
        long long lastRelaxationMCS;
        double theta_snapshot;

        // Physical and reference parameters
        double lambdaElastic;
        double refArea;
        double refVolume;
        double tau_r;         // Maxwell relaxation time
        double theta_relax;   // Relaxed strain
        bool is2D;

        // Fixed reference center (initial center of mass) for non-translation-invariant strain
        double cx0, cy0, cz0;
        bool refCenterInitialized;

        // Internal utility functions
        inline double _get_non_nan_energy(double energy);
        double computeThetaGlobal(long long n, double sx, double sy, double sz, double sxx, double syy, double szz, double sxy, double sxz, double syz);
        double energyFromTheta(double theta_global);
        void advanceRelaxationIfNewMCS();
        inline void ensureRefCenterInitialized();

        CC3DXMLElement *xmlData;
        Potts3D *potts;
        Simulator *sim;
        ParallelUtilsOpenMP *pUtils;
        ParallelUtilsOpenMP::OpenMPLock_t *lockPtr;
        Automaton *automaton;
        BoundaryStrategy *boundaryStrategy;
        WatchableField3D<CellG *> *cellFieldG;

    public:
        MaxwellMediumPlugin();
        virtual ~MaxwellMediumPlugin();

        // Energy function interface
        virtual double changeEnergy(const Point3D &pt, const CellG *newCell, const CellG *oldCell);

        // Lattice change watcher
        virtual void field3DChange(const Point3D &pt, CellG *newCell, CellG *oldCell) override;

        virtual void init(Simulator *simulator, CC3DXMLElement *_xmlData = 0);
        
        virtual void update(CC3DXMLElement *_xmlData, bool _fullInitFlag = false);
        virtual std::string steerableName();
        virtual std::string toString();

        // Python accessors
        virtual void setLambdaElastic(double lambdaElasticParam) { this->lambdaElastic = lambdaElasticParam; }
        virtual double getLambdaElastic() { return this->lambdaElastic; }
        virtual double getRefArea() { return this->refArea; }
        virtual double getRefVolume() { return this->refVolume; }
        virtual double getTauR() { return this->tau_r; }
        virtual double getThetaRelax() { return this->theta_relax; }
        virtual double getThetaGlobal();

    };

};

#endif
