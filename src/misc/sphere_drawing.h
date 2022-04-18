#ifndef CGL_UTIL_SPHEREDRAWING_H
#define CGL_UTIL_SPHEREDRAWING_H

#include <vector>

#include "CGL/CGL.h"

namespace CGL {
namespace Misc {
  
class SphereMesh {
public:
  // Supply the desired number of vertices
  SphereMesh(int num_lat = 40, int num_lon = 40);
};


}; // namespace Misc
} // namespace CGL

#endif // CGL_UTIL_SPHEREDRAWING_H
