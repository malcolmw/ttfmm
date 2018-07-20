#include <deque>

#define PI 3.141592653589793
#define EARTH_RADIUS 6371.0

typedef struct GeoGrid {
  double ***value;
  int ***frozen;
  int nlat, nlon, nz; // Number of grid nodes
  double lat0, lon0, z0; // Origin node of grid
  double dlat, dlon, dz; // Distance between nodes
} GeoGrid ;

typedef struct GridNode {
  int ilat, ilon, iz; // Number of grid nodes
  double value;
} GridNode ;

void hpsort(unsigned long n, std::deque<GridNode> &ra);

