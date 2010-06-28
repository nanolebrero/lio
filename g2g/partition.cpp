/* includes */
#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cmath>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;
using namespace G2G;

void g2g_compute_group_weights(PointGroup& group);

/********************
 * PointGroup
 ********************/
void PointGroup::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

void PointGroup::compute_weights(void) {
  g2g_compute_group_weights(*this);
}

/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : PointGroup(), atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : PointGroup(), atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/
Cube::Cube(void) : PointGroup() { }

/************************************************************
 * Construct partition
 ************************************************************/

/* global variables */
list<PointGroup> final_partition;

/* methods */
void regenerate_partition(void)
{
	cout << "<============ G2G Partition (" << fortran_vars.grid_type << ")============>" << endl;

	/* determina el exponente minimo para cada tipo de atomo */
	cout << "determining minimum exponents" << endl;	
	vector<double> min_exps(120, numeric_limits<double>::max());	// uno por elemento de la tabla periodica

	for (uint i = 0; i < fortran_vars.m; i++) {
		uint contractions = fortran_vars.contractions(i);
		uint nuc = fortran_vars.nucleii(i) - 1;
		uint nuc_type = fortran_vars.atom_types(nuc);
		for (uint j = 0; j < contractions; j++) {
			min_exps[nuc_type] = min(min_exps[nuc_type], fortran_vars.a_values(i, j));
		}
	}

	vector<double> min_exps_func(fortran_vars.m, numeric_limits<double>::max());	// uno por funcion
	for (uint i = 0; i < fortran_vars.m; i++) {
		uint contractions = fortran_vars.contractions(i);
		for (uint j = 0; j < contractions; j++) {
			min_exps_func[i] = min(min_exps_func[i], fortran_vars.a_values(i, j));
		}
	}
	
	/* permite encontrar el prisma conteniendo el sistema */
	cout << "determining x0 and x1" << endl;
	double3 x0 = make_double3(0,0,0);
	double3 x1 = make_double3(0,0,0);
	for (uint atom = 0; atom < fortran_vars.atoms; atom++) {
		double3 atom_position(fortran_vars.atom_positions(atom));

		uint atom_type = fortran_vars.atom_types(atom);
		double max_radius = sqrt(max_function_exponent / min_exps[atom_type]);
		_DBG(cout << "tipo: " << atom_type << " " << min_exps[atom_type] << " radio: " << max_radius << endl);
		if (atom == 0) { x0 = x1 = atom_position; }
		else {
			if ((atom_position.x - max_radius) < x0.x) x0.x = (atom_position.x - max_radius);
			if ((atom_position.y - max_radius) < x0.y) x0.y = (atom_position.y - max_radius);
			if ((atom_position.z - max_radius) < x0.z) x0.z = (atom_position.z - max_radius);

			if ((atom_position.x + max_radius) > x1.x) x1.x = (atom_position.x + max_radius);
			if ((atom_position.y + max_radius) > x1.y) x1.y = (atom_position.y + max_radius);
			if ((atom_position.z + max_radius) > x1.z) x1.z = (atom_position.z + max_radius);
		}		
	}

	/* el prisma tiene vertices (x,y) */
	cout << "x0 " << x0 << endl;  // vertice inferior, izquierdo y mas lejano
	cout << "x1 " << x1 << endl;  // vertice superior, derecho y mas cercano

  /* la particion en cubos */
	uint3 prism_size = ceil_uint3((x1 - x0) / little_cube_size);
	cout << "prism size: " << prism_size.x << " " << prism_size.y << " " << prism_size.z << endl;

	vector< vector < vector< Cube > > > prism(prism_size.x,
					vector< vector < Cube > >(prism_size.y,
									vector < Cube >(prism_size.z)));

  /* initialize spheres */
  vector<Sphere> sphere_array(fortran_vars.atoms);
  for (uint atom = 0; atom < fortran_vars.atoms; atom++) {
    uint atom_shells = fortran_vars.shells(atom);
    uint included_shells = (uint)ceil(sphere_radius * atom_shells);
		
		double radius;
		if (included_shells == 0) radius = 0;
		else {
			double x = cos((M_PI / (atom_shells + 1)) * (atom_shells - included_shells + 1));
	    double rm = fortran_vars.rm(atom);
	  	radius = rm * (1.0 + x) / (1.0 - x);
		}
		
		_DBG(cout << "esfera incluye " << included_shells << " capas de " << atom_shells << " (radio: " << radius << ")" << endl);
    sphere_array[atom] = Sphere(atom, radius);
  }

	cout << "precomputing distances..." << endl;
	for (uint i = 0; i < fortran_vars.atoms; i++) {
		const double3& atom_i_position(fortran_vars.atom_positions(i));
		double nearest_neighbor_dist = numeric_limits<double>::max();
		
		double sphere_i_radius = sphere_array[i].radius;
		
		for (uint j = 0; j < fortran_vars.atoms; j++) {
			const double3& atom_j_position(fortran_vars.atom_positions(j));
			double dist = length(atom_i_position - atom_j_position);
			fortran_vars.atom_atom_dists(i, j) = dist;
      //_DBG(cout << "distancia atomo " << i << " -> " << j << " : " << dist << endl);
			if (i != j) {
				nearest_neighbor_dist = min(nearest_neighbor_dist, dist);
        if (dist <= sphere_i_radius || dist <= sphere_array[j].radius) { throw runtime_error("other atom contained in sphere"); }
				if ((sphere_i_radius + sphere_array[j].radius) >= dist) { throw runtime_error("Overlapping sphere radius!"); }
      }
		}
		fortran_vars.nearest_neighbor_dists(i) = nearest_neighbor_dist;
	}
	
	cout << "computing points and assigning to cubes/spheres..." << endl;

	uint puntos_totales = 0;

  uint puntos_finales = 0;
  uint funciones_finales = 0;

  uint costo = 0;
	
	final_partition.clear();
	
	Timer t_total;
	t_total.start();
	/* computa las posiciones de los puntos (y los guarda) */
	for (uint atom = 0; atom < fortran_vars.atoms; atom++) {		
		uint atom_shells = fortran_vars.shells(atom);
		
		const double3& atom_position(fortran_vars.atom_positions(atom));

		double t0 = M_PI / (atom_shells + 1);
		double rm = fortran_vars.rm(atom);

		for (uint shell = 0; shell < atom_shells; shell++) {			
			double t1 = t0 * (shell + 1);
			double x = cos(t1);
			double w = t0 * abs(sin(t1));
			double r1 = rm * (1.0 + x) / (1.0 - x);
      double wrad = w * (r1 * r1) * rm * 2.0 / ((1.0 - x) * (1.0 - x));
			
			for (uint point = 0; point < (uint)fortran_vars.grid_size; point++) {
				puntos_totales++;
				
				double3 rel_point_position = make_double3(fortran_vars.e(point,0), fortran_vars.e(point,1), fortran_vars.e(point,2));
				double3 point_position = atom_position + rel_point_position * r1;

				bool inside_prism = ((x0.x <= point_position.x && point_position.x <= x1.x) &&
														(x0.y <= point_position.y && point_position.y <= x1.y) &&
														(x0.z <= point_position.z && point_position.z <= x1.z));

				if (inside_prism) {
					double point_weight = wrad * fortran_vars.wang(point); // integration weight
          Point point_object(atom, shell, point, point_position, point_weight);

          /* assign to sphere? */
          uint included_shells = (uint)ceil(sphere_radius * atom_shells);
          if (shell >= (atom_shells - included_shells)) {
            Sphere& sphere = sphere_array[atom];
            sphere.add_point(point_object);
          }
          /* or insert into corresponding cube */
          else {
            uint3 cube_coord = floor_uint3((point_position - x0) / little_cube_size);
            prism[cube_coord.x][cube_coord.y][cube_coord.z].add_point(point_object);
          }
				}
        else {
          //cout << "point outside prism: " << point_position.x << " " << point_position.y << " " << point_position.z << endl;
        }
			}
		}		
	}
	t_total.stop();
  
	cout << "Time: " << t_total << endl;
	cout << "Grilla original: " << puntos_totales << " puntos, " << puntos_totales * fortran_vars.m << " funciones" << endl;

  uint nco_m = 0;
  uint m_m = 0;

  puntos_finales = 0;
	cout << "filling cube parameters and adding to partition..." << endl;
	t_total.start();
	for (uint i = 0; i < prism_size.x; i++) {
		for (uint j = 0; j < prism_size.y; j++) {
			for (uint k = 0; k < prism_size.z; k++) {
				Cube& cube = prism[i][j][k];

				double3 cube_coord_abs = x0 + make_uint3(i,j,k) * little_cube_size;

        cube.assign_significative_functions(cube_coord_abs, min_exps_func);
				if (cube.functions.empty()) { /*cout << "cube with " << cube.number_of_points << " points has no functions" << endl;*/ continue;  }
        if (cube.number_of_points < min_points_per_cube) { /*cout << "not enough points" << endl;*/ continue; }

        final_partition.push_back(cube);

        PointGroup& group = final_partition.back();

        group.compute_weights();
        if (group.number_of_points < min_points_per_cube) { /*final_partition.pop_back(); cout << "not enough points" << endl;*/ continue; }

        // para hacer histogramas
#ifdef HISTOGRAM
        cout << "[" << fortran_vars.grid_type << "] cubo: (" << i << "," << j << "," << k << "): " << group.number_of_points << " puntos; " <<
          group.total_functions() << " funciones, vecinos: " << group.nucleii.size() << endl;
#endif

        puntos_finales += group.number_of_points;
				funciones_finales += group.number_of_points * group.total_functions();
        costo += group.number_of_points * (group.total_functions() + group.total_functions());
        //cout << "cubo: funcion x punto: " << cube.total_functions() / (double)cube.number_of_points << endl;

        nco_m += group.total_functions() * fortran_vars.nco;
        m_m += group.total_functions() * group.total_functions();
			}
		}
	}

	if (sphere_radius > 0) {
	  cout << "filling sphere parameters and adding to partition..." << endl;
	  for (uint i = 0; i < fortran_vars.atoms; i++) {
			Sphere& sphere = sphere_array[i];
	
			assert(sphere.number_of_points > 0);
			
			sphere.assign_significative_functions(min_exps_func);	
			assert(sphere.total_functions() > 0);

      if (sphere.number_of_points < min_points_per_cube) { cout << "not enough points" << endl; continue; }
      final_partition.push_front(sphere);

      PointGroup& group = final_partition.front();
      assert(group.number_of_points != 0);
	    group.compute_weights();
      t_total.sync();

      if (group.number_of_points < min_points_per_cube) { final_partition.pop_front(); cout << "not enough points" << endl; continue; }

#ifdef HISTOGRAM
			cout << "sphere: " << group.number_of_points << " puntos, " << group.total_functions() <<
        " funciones | funcion x punto: " << group.total_functions() / (double)group.number_of_points <<
        " vecinos: " << group.nucleii.size() << endl;
#endif

      assert(group.number_of_points > 0);

	    puntos_finales += group.number_of_points;
	    funciones_finales += group.number_of_points * group.total_functions();
      costo += group.number_of_points * (group.total_functions() + group.total_functions());

      nco_m += group.total_functions() * fortran_vars.nco;
      m_m += group.total_functions() * group.total_functions();
	  }
	}
  t_total.sync();
	t_total.stop();

	cout << "total: " << t_total << endl;
	cout << "Grilla final: " << puntos_finales << " puntos (recordar que los de peso 0 se tiran), " << funciones_finales << " funciones" << endl ;
  cout << "Costo: " << costo << endl;
  cout << "NCOxM: " << nco_m << " MxM: " << m_m << endl;
  cout << "Particion final: " << final_partition.size() << " grupos" << endl;
}
