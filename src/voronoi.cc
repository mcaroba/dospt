#include "voro++/install/include/voro++/voro++.hh"
using namespace voro;

int i;

extern "C" {

void voronoi_volumes_(int *n, double L[], double x[], double y[], double z[], double volumes[]) {

  container con(0.,L[0],0.,L[1],0.,L[2],1,1,1,true,true,true,8);

  i = 0;
  while (i < *n) {
    con.put(i, x[i], y[i], z[i]);
    i += 1;
  }

  voronoicell v;

  i = 0;
  while (i < *n){
    con.compute_cell(v,0,i);
    volumes[i] = v.volume();
    i += 1;
  }

}
}
