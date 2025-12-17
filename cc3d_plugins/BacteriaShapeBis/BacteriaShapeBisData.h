

#ifndef BACTERIASHAPEBISPATA_H
#define BACTERIASHAPEBISPATA_H

#include <vector>
#include <random>
#include <cmath>
#include "BacteriaShapeBisDLLSpecifier.h"


namespace CompuCell3D {

   class BACTERIASHAPEBIS_EXPORT BacteriaShapeBisData{

      public:

         BacteriaShapeBisData(){
            // Initialize orientation angles with random values between -π/2 and π/2
            static thread_local std::mt19937 rng{std::random_device{}()};
            std::uniform_real_distribution<double> dist(-M_PI/2.0, M_PI/2.0);
            theta = dist(rng);
            phi = dist(rng);
            // Initialize principal axis from angles (stable 3D control)
            axis_x = std::cos(phi) * std::cos(theta);
            axis_y = std::cos(phi) * std::sin(theta);
            axis_z = std::sin(phi);
            double n = std::sqrt(axis_x*axis_x + axis_y*axis_y + axis_z*axis_z) + 1e-12;
            axis_x /= n; axis_y /= n; axis_z /= n;
         };
         ~BacteriaShapeBisData(){};
         
         // Orientation parameters (spherical angles for axis)
         double theta = 0.0;  // azimuth around z-axis (2D and 3D)
         double phi = 0.0;    // elevation from xy-plane (3D only)
         
         // Principal axis (unit vector) for 3D capsule orientation
         double axis_x = 1.0;
         double axis_y = 0.0;
         double axis_z = 0.0;
         
        // Shape parameters
        double rotationStep = 0;
         double majorAxisLength = 20.0;  // Length along principal axis
         double minorAxisLength = 10.0;  // Diameter perpendicular to principal axis


   };
};

#endif

