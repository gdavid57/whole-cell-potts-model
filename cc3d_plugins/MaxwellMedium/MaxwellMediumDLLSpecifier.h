

#ifndef MAXWELLMEDIUM_EXPORT_H

#define MAXWELLMEDIUM_EXPORT_H

    #if defined(_WIN32)

      #ifdef MaxwellMediumShared_EXPORTS

          #define MAXWELLMEDIUM_EXPORT __declspec(dllexport)

          #define MAXWELLMEDIUM_EXPIMP_TEMPLATE

      #else

          #define MAXWELLMEDIUM_EXPORT __declspec(dllimport)

          #define MAXWELLMEDIUM_EXPIMP_TEMPLATE extern

      #endif

    #else

         #define MAXWELLMEDIUM_EXPORT

         #define MAXWELLMEDIUM_EXPIMP_TEMPLATE

    #endif

#endif

