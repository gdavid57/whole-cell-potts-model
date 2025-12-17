

#ifndef BACTERIASHAPEBISPLUGIN_H
#define BACTERIASHAPEBISPLUGIN_H

#include <CompuCell3D/CC3D.h>

#include "BacteriaShapeBisData.h"
#include "BacteriaShapeBisDLLSpecifier.h"
#include <vector>
#include <tuple>

class CC3DXMLElement;

namespace CompuCell3D {

    class Simulator;
    class Potts3D;
    class Automaton;
    class BoundaryStrategy;
    class ParallelUtilsOpenMP;

    template <class T> class Field3D;
    template <class T> class WatchableField3D;

    class BACTERIASHAPEBIS_EXPORT  BacteriaShapeBisPlugin : public Plugin ,public EnergyFunction ,public CellGChangeWatcher  {

    private:    

        ExtraMembersGroupAccessor<BacteriaShapeBisData> bacteriaShapeBisDataAccessor;                
        CC3DXMLElement *xmlData;        
        Potts3D *potts;
        Simulator *sim;
        ParallelUtilsOpenMP *pUtils;            
        ParallelUtilsOpenMP::OpenMPLock_t *lockPtr;        
        Automaton *automaton;
        BoundaryStrategy *boundaryStrategy;
        WatchableField3D<CellG *> *cellFieldG;
        bool is2D;

    public:
        BacteriaShapeBisPlugin();
        virtual ~BacteriaShapeBisPlugin();
        ExtraMembersGroupAccessor<BacteriaShapeBisData> * getBacteriaShapeBisDataAccessorPtr(){return & bacteriaShapeBisDataAccessor;}                
        
        double energyPenalty;
        double lambdaShape;

        // Amplitudes and criteria
        double rotationAmplitude;
        double translationAmplitude; // constant translation amplitude along orientation
        double rotationCriterion;    // cadence for orientation updates
        bool allowBackward;          // keep for compatibility (sign on translation)
        
        //Energy function interface
        virtual double changeEnergy(const Point3D &pt, const CellG *newCell, const CellG *oldCell);        
        
        virtual void setMajorAxisLength(CellG *cell, double _majorAxisLength);
        virtual void setMinorAxisLength(CellG *cell, double _minorAxisLength);
        virtual void setTheta(CellG *cell, double _theta);
        virtual double getTheta(CellG *cell);
        virtual double getLambdaShape();

        // 3D orientation getters (angles kept for reporting)
        virtual double getPhi(CellG *cell);

        // 3D axis vector control (primary orientation used by energy)
        virtual void setAxis(CellG *cell, double ax, double ay, double az);
        
        // Individual axis component getters to avoid SWIG std::vector issues
        virtual double getAxisX(CellG *cell);
        virtual double getAxisY(CellG *cell);
        virtual double getAxisZ(CellG *cell);

        // Callback triggered each time a lattice pixel changes ownership
        virtual void field3DChange(const Point3D &pt, CellG *newCell, CellG *oldCell);

        // Behavior updater (only rotation on rotationCriterion cadence)
        void updateBehavior(CellG *cell);
        // Orientation helper also refreshes lambda vectors immediately
        void updateOrientation(CellG *cell);

        // Shape testing functions
        virtual int isPointInBacillus(double xCOM, double yCOM, double px, double py, double L, double l, double angle);
        virtual int isPointInBacillus3D_Axis(double xCOM, double yCOM, double zCOM, double px, double py, double pz,
                                            double majorAxis, double minorAxis, double ax, double ay, double az);
        
        // Helper functions for 3D
        virtual void computePrincipalAxisFromInertia(double Ixx, double Iyy, double Izz,
                                                     double Ixy, double Ixz, double Iyz,
                                                     double &ax, double &ay, double &az);

        virtual void init(Simulator *simulator, CC3DXMLElement *_xmlData=0);
        virtual void extraInit(Simulator *simulator);
        //Steerrable interface

        virtual void update(CC3DXMLElement *_xmlData, bool _fullInitFlag=false);
        virtual std::string steerableName();

        virtual std::string toString();

    };

};

#endif
