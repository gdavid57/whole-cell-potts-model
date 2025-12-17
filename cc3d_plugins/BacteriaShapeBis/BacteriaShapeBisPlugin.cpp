#include <CompuCell3D/CC3D.h>        

using namespace CompuCell3D;

#include "BacteriaShapeBisPlugin.h"
#include <cmath>
#include <random>

BacteriaShapeBisPlugin::BacteriaShapeBisPlugin():
pUtils(0),
lockPtr(0),
xmlData(0) ,
cellFieldG(0),
boundaryStrategy(0),
is2D(true)
{}

BacteriaShapeBisPlugin::~BacteriaShapeBisPlugin() {

    pUtils->destroyLock(lockPtr);

    delete lockPtr;

    lockPtr=0;

}

void BacteriaShapeBisPlugin::init(Simulator *simulator, CC3DXMLElement *_xmlData) {

    xmlData=_xmlData;
    sim=simulator;
    potts=simulator->getPotts();
    cellFieldG = (WatchableField3D<CellG *> *)potts->getCellFieldG();

    pUtils=sim->getParallelUtils();

    lockPtr=new ParallelUtilsOpenMP::OpenMPLock_t;

    pUtils->initLock(lockPtr); 

    update(xmlData,true);

    potts->getCellFactoryGroupPtr()->registerClass(&bacteriaShapeBisDataAccessor);

    potts->registerEnergyFunctionWithName(this,"BacteriaShapeBis");

    // register callback for lattice pixel ownership changes
    potts->registerCellGChangeWatcher(this);

    simulator->registerSteerableObject(this);

    // Detect dimensionality following MaxwellMedium pattern
    Dim3D fieldDim = cellFieldG->getDim();
    is2D = (fieldDim.x == 1 || fieldDim.y == 1 || fieldDim.z == 1);

    // Initialize constant translation lambda vectors from current orientations
    {
        CellInventory *cellInventoryPtr = &potts->getCellInventory();
        for (CellInventory::cellInventoryIterator cInvItr = cellInventoryPtr->cellInventoryBegin();
             cInvItr != cellInventoryPtr->cellInventoryEnd(); ++cInvItr){
            CellG *cell = cellInventoryPtr->getCell(cInvItr);
            if (!cell) continue;
            auto data = bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr);
            double ax, ay, az;
            if (is2D){
                ax = std::cos(data->theta); ay = std::sin(data->theta); az = 0.0;
            } else {
                ax = data->axis_x; ay = data->axis_y; az = data->axis_z;
                double n = std::sqrt(ax*ax + ay*ay + az*az) + 1e-12; ax/=n; ay/=n; az/=n;
            }
            // Directly use translationAmplitude (signed) to determine force vector
            cell->lambdaVecX = translationAmplitude * ax;
            cell->lambdaVecY = translationAmplitude * ay;
            cell->lambdaVecZ = translationAmplitude * az;
        }
    }

}

void BacteriaShapeBisPlugin::extraInit(Simulator *simulator){
}

void BacteriaShapeBisPlugin::update(CC3DXMLElement *_xmlData, bool _fullInitFlag){

    //PARSE XML IN THIS FUNCTION
    //For more information on XML parser function please see CC3D code or lookup XML utils API
    automaton = potts->getAutomaton();
    ASSERT_OR_THROW("CELL TYPE PLUGIN WAS NOT PROPERLY INITIALIZED YET. MAKE SURE THIS IS THE FIRST PLUGIN THAT YOU SET", automaton)
   set<unsigned char> cellTypesSet;

    // Current parameters only (no legacy)
    CC3DXMLElement * rotAmpElem = xmlData->getFirstElement("RotationAmplitude");
    rotationAmplitude = rotAmpElem ? rotAmpElem->getDouble() : 0.0;

    CC3DXMLElement * transAmpElem = xmlData->getFirstElement("TranslationAmplitude");
    translationAmplitude = transAmpElem ? transAmpElem->getDouble() : 0.0;

    CC3DXMLElement * rotCriElem = xmlData->getFirstElement("RotationCriterion");
    rotationCriterion = rotCriElem ? rotCriElem->getDouble() : 1000000.0;

    CC3DXMLElement * lamShaElem = xmlData->getFirstElement("LambdaShape");
    lambdaShape = lamShaElem ? lamShaElem->getDouble() : 10.0;

    // Keep allowBackward for optional sign inversion of constant translation
    CC3DXMLElement * allBacElem = xmlData->getFirstElement("AllowBackward");
    allowBackward = allBacElem ? allBacElem->getBool() : false;

    boundaryStrategy=BoundaryStrategy::getInstance();

}

// Simplified changeEnergy: orientation is updated elsewhere (field3DChange)
double BacteriaShapeBisPlugin::changeEnergy(const Point3D &pt,const CellG *newCell,const CellG *oldCell) {    

    double energy = 0;
    
    if (oldCell){
        auto data = bacteriaShapeBisDataAccessor.get(oldCell->extraAttribPtr);
        int V_pt_target;
        
        if (is2D) {
            // 2D case: use existing 2D function
            V_pt_target = isPointInBacillus(oldCell->xCOM, oldCell->yCOM, pt.x, pt.y, 
                                          data->majorAxisLength, data->minorAxisLength, data->theta);
        } else {
            // 3D case: use axis-based function
            V_pt_target = isPointInBacillus3D_Axis(oldCell->xCOM, oldCell->yCOM, oldCell->zCOM,
                                                   pt.x, pt.y, pt.z,
                                                   data->majorAxisLength, data->minorAxisLength,
                                                   data->axis_x, data->axis_y, data->axis_z);
        }
        energy += lambdaShape * (2.0 * V_pt_target - 1.0);
    }

    if(newCell){
        auto data = bacteriaShapeBisDataAccessor.get(newCell->extraAttribPtr);
        int V_pt_target;
        
        if (is2D) {
            // 2D case: use existing 2D function
            V_pt_target = isPointInBacillus(newCell->xCOM, newCell->yCOM, pt.x, pt.y, 
                                          data->majorAxisLength, data->minorAxisLength, data->theta);
        } else {
            // 3D case: use axis-based function
            V_pt_target = isPointInBacillus3D_Axis(newCell->xCOM, newCell->yCOM, newCell->zCOM,
                                                   pt.x, pt.y, pt.z,
                                                   data->majorAxisLength, data->minorAxisLength,
                                                   data->axis_x, data->axis_y, data->axis_z);
        }
        energy += lambdaShape * (1.0 - 2.0 * V_pt_target);
    }

    return energy; 
}

// ----------------------------------------------------------------------------------
// Behavior update: called on each pixel change
void BacteriaShapeBisPlugin::field3DChange(const Point3D &pt, CellG *newCell, CellG *oldCell){
    updateBehavior(newCell);
    updateBehavior(oldCell);
}

void BacteriaShapeBisPlugin::updateBehavior(CellG *cell){
    if(!cell) return;
    static thread_local std::mt19937 rng{std::random_device{}()};
    auto data = bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr);
    data->rotationStep += 1;
    if(data->rotationStep < rotationCriterion) return;
    data->rotationStep = 0;
    updateOrientation(cell);
}

void BacteriaShapeBisPlugin::updateOrientation(CellG *cell){
    if(!cell) return;
    static thread_local std::mt19937 rng{std::random_device{}()};
    auto data = bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr);
    if (is2D){
        double Ixx = cell->iXX;
        double Iyy = cell->iYY;
        double Ixy = cell->iXY;
        double thetaGeometric = 0.5 * atan2(-2.0 * Ixy, Iyy - Ixx);
        std::uniform_real_distribution<double> noiseDist(-1.0, 1.0);
        double randomNoise = noiseDist(rng) * rotationAmplitude;
        data->theta = thetaGeometric + randomNoise;
    } else {
        double ax_new, ay_new, az_new;
        computePrincipalAxisFromInertia(cell->iXX, cell->iYY, cell->iZZ,
                                        cell->iXY, cell->iXZ, cell->iYZ,
                                        ax_new, ay_new, az_new);

        double ax_old = data->axis_x, ay_old = data->axis_y, az_old = data->axis_z;
        if (ax_old*ax_new + ay_old*ay_new + az_old*az_new < 0.0) {
            ax_new = -ax_new; ay_new = -ay_new; az_new = -az_new;
        }

        const double alpha = 0.0;
        double ax = alpha*ax_old + (1.0-alpha)*ax_new;
        double ay = alpha*ay_old + (1.0-alpha)*ay_new;
        double az = alpha*az_old + (1.0-alpha)*az_new;
        double n = std::sqrt(ax*ax + ay*ay + az*az) + 1e-12;
        ax /= n; ay /= n; az /= n;

        if (rotationAmplitude > 0.0) {
            std::uniform_real_distribution<double> uni(-1.0, 1.0);
            double eps = uni(rng) * rotationAmplitude;
            double rx = uni(rng), ry = uni(rng), rz = uni(rng);
            double dot = rx*ax + ry*ay + rz*az;
            rx -= dot*ax; ry -= dot*ay; rz -= dot*az;
            double rn = std::sqrt(rx*rx + ry*ry + rz*rz) + 1e-12;
            rx /= rn; ry /= rn; rz /= rn;
            double c = std::cos(eps), s = std::sin(eps);
            double ax2 = c*ax + s*rx;
            double ay2 = c*ay + s*ry;
            double az2 = c*az + s*rz;
            double n2 = std::sqrt(ax2*ax2 + ay2*ay2 + az2*az2) + 1e-12;
            ax = ax2 / n2; ay = ay2 / n2; az = az2 / n2;
        }

        data->axis_x = ax; data->axis_y = ay; data->axis_z = az;
        data->theta = std::atan2(ay, ax);
        data->phi = std::atan2(az, std::sqrt(ax*ax + ay*ay));
    }
    // After updating orientation, immediately refresh lambda vectors for constant translation
    double ax, ay, az;
    if (is2D){
        ax = std::cos(data->theta); ay = std::sin(data->theta); az = 0.0;
    } else {
        ax = data->axis_x; ay = data->axis_y; az = data->axis_z;
        double n = std::sqrt(ax*ax + ay*ay + az*az) + 1e-12; ax/=n; ay/=n; az/=n;
    }
    // Directly use translationAmplitude (signed)
    cell->lambdaVecX = translationAmplitude * ax;
    cell->lambdaVecY = translationAmplitude * ay;
    cell->lambdaVecZ = translationAmplitude * az;
}

// removed per-MCS mobility updater

int BacteriaShapeBisPlugin::isPointInBacillus(double xCOM, double yCOM, double px, double py, 
                                                  double majorAxis, double minorAxis, double theta) {
    // Rotation inverse du point dans le repère du bacille
    double x_rel = (px - xCOM) * cos(-theta) - (py - yCOM) * sin(-theta);
    double y_rel = (px - xCOM) * sin(-theta) + (py - yCOM) * cos(-theta);

    // Longueur des demi-axes
    double halfL = majorAxis / 2.0;
    double halfl = minorAxis / 2.0;

    // Vérification dans le rectangle central
    if (fabs(x_rel) <= halfL - halfl && fabs(y_rel) <= halfl) {
        return 1;
    }

    // Vérification dans les demi-sphères
    double dx = fabs(x_rel) - (halfL - halfl);
    double dy = y_rel;
    if (dx * dx + dy * dy <= halfl * halfl) {
        return 1; 
    }

    // Cas extérieur
    return 0;
}

int BacteriaShapeBisPlugin::isPointInBacillus3D_Axis(double xCOM, double yCOM, double zCOM,
                                                    double px, double py, double pz,
                                                    double majorAxis, double minorAxis,
                                                    double ax, double ay, double az) {
    // Relative position
    double rx = px - xCOM, ry = py - yCOM, rz = pz - zCOM;

    // Projection along the axis
    double s = rx*ax + ry*ay + rz*az;

    // Radial squared distance to the axis
    double r2 = rx*rx + ry*ry + rz*rz - s*s;

    double halfL = majorAxis / 2.0;
    double radius = minorAxis / 2.0;

    if (std::fabs(s) <= halfL - radius) {
        return (r2 <= radius*radius) ? 1 : 0;
    }

    double ds = std::fabs(s) - (halfL - radius);
    return (ds*ds + r2 <= radius*radius) ? 1 : 0;
}

void BacteriaShapeBisPlugin::computePrincipalAxisFromInertia(double Ixx, double Iyy, double Izz,
                                                             double Ixy, double Ixz, double Iyz,
                                                             double &ax, double &ay, double &az) {
    // Symmetric 3x3 Jacobi eigen decomposition (few sweeps)
    double A[3][3] = { {Ixx, Ixy, Ixz}, {Ixy, Iyy, Iyz}, {Ixz, Iyz, Izz} };
    double V[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };

    auto rotate = [&](int p, int q){
        if (A[p][q] == 0.0) return;
        double app = A[p][p], aqq = A[q][q], apq = A[p][q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0 ? 1.0 : -1.0) / (std::fabs(tau) + std::sqrt(1.0 + tau*tau));
        double c = 1.0 / std::sqrt(1.0 + t*t);
        double s = t * c;

        double App = app - t * apq;
        double Aqq = aqq + t * apq;
        A[p][p] = App; A[q][q] = Aqq; A[p][q] = A[q][p] = 0.0;

        for (int r=0; r<3; ++r) if (r != p && r != q) {
            double arp = A[r][p], arq = A[r][q];
            A[r][p] = A[p][r] = c*arp - s*arq;
            A[r][q] = A[q][r] = c*arq + s*arp;
        }
        for (int r=0; r<3; ++r) {
            double vrp = V[r][p], vrq = V[r][q];
            V[r][p] = c*vrp - s*vrq;
            V[r][q] = c*vrq + s*vrp;
        }
    };

    for (int it=0; it<8; ++it) {
        int p=0, q=1; double m = std::fabs(A[0][1]);
        if (std::fabs(A[0][2]) > m) { p=0; q=2; m = std::fabs(A[0][2]); }
        if (std::fabs(A[1][2]) > m) { p=1; q=2; }
        if (m < 1e-12) break;
        rotate(p,q);
    }

    int k = 0;
    if (A[1][1] < A[k][k]) k = 1;
    if (A[2][2] < A[k][k]) k = 2;

    ax = V[0][k]; ay = V[1][k]; az = V[2][k];
    double n = std::sqrt(ax*ax + ay*ay + az*az) + 1e-12;
    ax /= n; ay /= n; az /= n;
}

std::string BacteriaShapeBisPlugin::toString(){
    return "BacteriaShapeBis";
}

std::string BacteriaShapeBisPlugin::steerableName(){
    return toString();
}

void BacteriaShapeBisPlugin::setMajorAxisLength(CellG *cell, double _majorAxisLength){
    bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->majorAxisLength = _majorAxisLength;
}

void BacteriaShapeBisPlugin::setMinorAxisLength(CellG *cell, double _minorAxisLength){
    bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->minorAxisLength = _minorAxisLength;
}

void BacteriaShapeBisPlugin::setTheta(CellG *cell, double _theta){
    // Only meaningful in 2D; in 3D orientation is controlled by axis vector
    if (is2D) {
        bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->theta = _theta;
    }
}

double BacteriaShapeBisPlugin::getTheta(CellG *cell){
    return bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->theta;
}

double BacteriaShapeBisPlugin::getLambdaShape(){
    return this->lambdaShape;
}

double BacteriaShapeBisPlugin::getPhi(CellG *cell){
    return bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->phi;
}


void BacteriaShapeBisPlugin::setAxis(CellG *cell, double ax, double ay, double az){
    auto data = bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr);
    // normalize to keep unit vector
    double n = std::sqrt(ax*ax + ay*ay + az*az) + 1e-12;
    ax /= n; ay /= n; az /= n;
    
    data->axis_x = ax; 
    data->axis_y = ay; 
    data->axis_z = az;
    
    // keep angles in sync for reporting convenience
    data->theta = std::atan2(ay, ax);
    data->phi   = std::atan2(az, std::sqrt(ax*ax + ay*ay));

    // Update lambdaVec based on new axis and translation amplitude
    // Directly use translationAmplitude (signed)
    cell->lambdaVecX = translationAmplitude * ax;
    cell->lambdaVecY = translationAmplitude * ay;
    cell->lambdaVecZ = translationAmplitude * az;
}

double BacteriaShapeBisPlugin::getAxisX(CellG *cell){
    return bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->axis_x;
}

double BacteriaShapeBisPlugin::getAxisY(CellG *cell){
    return bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->axis_y;
}

double BacteriaShapeBisPlugin::getAxisZ(CellG *cell){
    return bacteriaShapeBisDataAccessor.get(cell->extraAttribPtr)->axis_z;
}
