#include "voro++/install/include/voro++/voro++.hh"
using namespace voro;
#include <iostream>
using namespace std;

int i;
int j;
int i_vacuum;

//***********************************************************************
// Code for the dodecahedron envelope for vacuum calculations
// This is copied from Chris Rycroft's "irregular packing example code":
//
// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));
double Rh;
// Create a wall class that, whenever called, will replace the Voronoi cell
// with a prescribed shape, in this case a dodecahedron
class wall_initial_shape : public wall {
        public:
                wall_initial_shape(double R) {
                        // Create a dodecahedron of the same volume as a sphere of radius R
                        Rh = 2. / sqrt(Phi*Phi+1.) * 0.91045318140924225 * R;
                        v.init(-2*Rh,2*Rh,-2*Rh,2*Rh,-2*Rh,2*Rh);
                        v.plane(0,Phi*Rh,1*Rh);v.plane(0,-Phi*Rh,1*Rh);v.plane(0,Phi*Rh,-1*Rh);
                        v.plane(0,-Phi*Rh,-1*Rh);v.plane(1*Rh,0,Phi*Rh);v.plane(-1*Rh,0,Phi*Rh);
                        v.plane(1*Rh,0,-Phi*Rh);v.plane(-1*Rh,0,-Phi*Rh);v.plane(Phi*Rh,1*Rh,0);
                        v.plane(-Phi*Rh,1*Rh,0);v.plane(Phi*Rh,-1*Rh,0);v.plane(-Phi*Rh,-1*Rh,0);
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
                      int *vacuum, double R[], double volumes[]) {

  container con(0.,L[0],0.,L[1],0.,L[2],1,1,1,true,true,true,8);

  i_vacuum = *vacuum;


// This works if there is a vacuum in the simulation box (specified by the user)
// This is slower since we reset the container walls for each atom to account for
// atom-specific radii
  if ( i_vacuum == 1 ) {
    // Add the particles to the container
    i = 0;
    while (i < *n) {
      wall_initial_shape wis(R[i]);
      con.add_wall(wis);
      j = 0;
      while (j < *n) {
        con.put(j, x[j], y[j], z[j]);
        j += 1;
      }
      voronoicell v;
      con.compute_cell(v,0,i);
      volumes[i] = v.volume();
      i += 1;
    }

  } else {
    // Add the particles to the container
    i = 0;
    while (i < *n) {
      con.put(i, x[i], y[i], z[i]);
      i += 1;
    }
    voronoicell v;
    i = 0;
    while (i < *n) {
      con.compute_cell(v,0,i);
      volumes[i] = v.volume();
      i += 1;
    }
  }



}
}
