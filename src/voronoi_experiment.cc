#include "voro++/install/include/voro++/voro++.hh"
using namespace voro;

int i;



//***********************************************************************
// Code for the dodecahedron envelope for vacuum calculations
// This is copied from Chris Rycroft's "irregular packing example code":
//
// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));
// Create a wall class that, whenever called, will replace the Voronoi cell
// with a prescribed shape, in this case a dodecahedron
class wall_initial_shape : public wall {
        public:
                wall_initial_shape(double R) {

                        // Create a dodecahedron
                        v.init(-R,R,-R,R,-R,R);
                        v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
                        v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
                        v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
                        v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
                };
                bool point_inside(double x,double y,double z) {return true;}
                bool cut_cell(voronoicell &c,double x,double y,double z) {

                        // Set the cell to be equal to the dodecahedron
                        c=v;
                        return true;
                }
                bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {

                        // Set the cell to be equal to the dodecahedron
                        c=v;
                        return true;
                }
        private:
                voronoicell v;
};
//***********************************************************************





extern "C" {

void voronoi_volumes_(int *n, double L[], double x[], double y[], double z[], \
                      bool vacuum, double R[], double volumes[]) {

  container con(0.,L[0],0.,L[1],0.,L[2],1,1,1,true,true,true,8);

  // Create the "initial shape" wall class and add it to the container
  i = 0;
  while (i < *n) {
    wis = wall_initial_shape(R[i]);
    con.add_wall(wis);
    i += 1
  }

  // Add the particles to the container
  i = 0;
  while (i < *n) {
    con.put(i, x[i], y[i], z[i]);
    i += 1;
  }

//  voronoicell v;

  i = 0;
  while (i < *n){
    con.compute_cell(v,0,i);
    volumes[i] = v.volume();
    i += 1;
  }

}
}
