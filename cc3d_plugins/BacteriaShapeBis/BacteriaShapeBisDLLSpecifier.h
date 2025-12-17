

#ifndef BACTERIASHAPEBIS_EXPORT_H

#define BACTERIASHAPEBIS_EXPORT_H

    #if defined(_WIN32)

      #ifdef BacteriaShapeBisShared_EXPORTS

          #define BACTERIASHAPEBIS_EXPORT __declspec(dllexport)

          #define BACTERIASHAPEBIS_EXPIMP_TEMPLATE

      #else

          #define BACTERIASHAPEBIS_EXPORT __declspec(dllimport)

          #define BACTERIASHAPEBIS_EXPIMP_TEMPLATE extern

      #endif

    #else

         #define BACTERIASHAPEBIS_EXPORT

         #define BACTERIASHAPEBIS_EXPIMP_TEMPLATE

    #endif

#endif

