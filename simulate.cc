#include "cross_sections.cc"
#include "grid_2d.cc"
#include "grid_3d.cc"
#include "material.cc"
#include "proton_beam.cc"
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Call " << argv[0] << " config-path" << std::endl;
    return 1;
  }
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));

  std::vector<Atom> atoms;
  int a, z;
  std::string name;
  std::ifstream file;
  file.open("./Materials/atoms.txt");
  std::string line, token;
  while (getline(file, line)) {
    std::stringstream iss;
    iss << line;
    getline(iss, token, ' ');
    name = token.c_str();
    getline(iss, token, ' ');
    z = atoi(token.c_str());
    getline(iss, token, ' ');
    a = atoi(token.c_str());
    if (name == "hydrogen") {
      Atom tmp(a, z, "./Splines/" + name + "_el_ruth_rate.txt",
               "./Splines/" + name + "_el_ruth_angle.txt");
      atoms.push_back(tmp);
    } else {
      Atom tmp(a, z, "./Splines/" + name + "_el_ruth_rate.txt",
               "./Splines/" + name + "_ne_rate.txt",
               "./Splines/" + name + "_el_ruth_angle.txt",
               "./Splines/" + name + "_ne_energyangle_cdf.txt");
      atoms.push_back(tmp);
    }
  }
  file.close();

  std::vector<std::string> material_names;
  file.open("./Materials/materials.txt");
  while (getline(file, line)) {
    material_names.push_back(line);
  }
  file.close();

  std::vector<Material> materials(material_names.size());
  for (unsigned int i = 0; i < material_names.size(); i++) {
    materials[i].read_material("./Materials/" + material_names[i], atoms);
  }

  double nozzle_radius, e0, initial_x_sd;
  std::vector<double> x(3, 0), w(2, 0);
  w[0] = M_PI / 2;

  libconfig::Config cfg;
  cfg.readFile(argv[1]);
  cfg.lookupValue("nozzle_radius", nozzle_radius);
  cfg.lookupValue("initial_x_sd", initial_x_sd);

  double initial_e_mean, initial_e_sd;
  cfg.lookupValue("initial_e_mean", initial_e_mean);
  cfg.lookupValue("initial_e_sd", initial_e_sd);

  double dt, absorption_e, air_gap;
  int nrep;
  cfg.lookupValue("step_size", dt);
  cfg.lookupValue("absorption_energy", absorption_e);
  cfg.lookupValue("replicates", nrep);
  cfg.lookupValue("air_gap", air_gap);

  double grid_dx;
  cfg.lookupValue("grid_dx", grid_dx);

  int output_1d, output_2d, output_3d, output_2d_slice;
  std::string path_1d, path_2d, path_3d, path_slice;
  cfg.lookupValue("output_1d", output_1d);
  cfg.lookupValue("output_2d", output_2d);
  cfg.lookupValue("output_3d", output_3d);
  cfg.lookupValue("output_2d_slice", output_2d_slice);
  if (output_1d == 1) {
    cfg.lookupValue("out_path_1d", path_1d);
  }
  if (output_2d == 1) {
    cfg.lookupValue("out_path_2d", path_2d);
  }
  if (output_3d == 1) {
    cfg.lookupValue("out_path_3d", path_3d);
  }
  if (output_2d_slice == 1) {
    cfg.lookupValue("out_path_2d_slice", path_slice);
  }

  proton_path p(initial_e_mean + 3 * initial_e_sd, dt, air_gap, absorption_e,
                materials);
  int n = p.energy.size();
  // We need all grids to be initialised for code to compile
  // but make unneeded ones tiny
  int tmp = n * grid_dx / dt;
  if (output_1d == 0 && output_2d == 0) {
    tmp = 1;
  }
  Grid_2d g2d(tmp, grid_dx);
  tmp = n * grid_dx / dt;
  if (output_3d == 0) {
    tmp = 1;
  }
  Grid_3d g3d(tmp, grid_dx);
  tmp = n * grid_dx / dt;
  if (output_2d_slice == 0) {
    tmp = 1;
  }
  Grid_2d g2d_slice(tmp, grid_dx);

  int len;
  for (int i = 0; i < nrep; i++) {
    e0 = initial_e_mean + gsl_ran_gaussian_ziggurat(gen, initial_e_sd);
    do {
      x[1] = gsl_ran_gaussian_ziggurat(gen, initial_x_sd);
      x[2] = gsl_ran_gaussian_ziggurat(gen, initial_x_sd);
    } while (x[1] * x[1] + x[2] * x[2] > nozzle_radius * nozzle_radius);
    p.reset(e0, x, w);
    len = p.simulate(dt, absorption_e, air_gap, gen, materials);
    if (output_1d == 1 || output_2d == 1) {
      g2d.add(p.x, p.s, len);
    }
    if (output_3d == 1) {
      g3d.add(p.x, p.s, len);
    }
    if (output_2d_slice == 1) {
      g2d_slice.add_to_slice(p.x, p.s, len);
    }
  }
  if (output_1d == 1) {
    g2d.print_1d(path_1d);
  }
  if (output_2d == 1) {
    g2d.print(path_2d);
  }
  if (output_3d == 1) {
    g3d.print(path_3d);
  }
  if (output_2d_slice == 1) {
    g2d_slice.print(path_slice);
  }
  gsl_rng_free(gen);
  return 1;
}
