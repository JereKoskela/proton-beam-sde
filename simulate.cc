#include "cross_sections.cc"
#include "grid.cc"
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
  double Rutherford_cutoff, Backscatter_cutoff;
  libconfig::Config cfg;
  cfg.readFile(argv[1]);
  cfg.lookupValue("Rutherford_cutoff", Rutherford_cutoff);
  cfg.lookupValue("Backscatter_cutoff", Backscatter_cutoff);
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
      Atom tmp(a, z, "./Splines/" + name + "_el_ruth_cross_sec.txt",
               Rutherford_cutoff, Backscatter_cutoff);
      atoms.push_back(tmp);
    } else {
      Atom tmp(a, z, "./Splines/" + name + "_ne_rate.txt",
               "./Splines/" + name + "_el_ruth_cross_sec.txt",
               "./Splines/" + name + "_ne_energyangle_cdf.txt",
               Rutherford_cutoff);
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
  std::string out_path;
  cfg.lookupValue("out_path", out_path);

  proton_path p(initial_e_mean + 3 * initial_e_sd, dt, air_gap, absorption_e,
                materials);
  int n = p.energy.size();
  Grid grid(n * grid_dx / dt, grid_dx);

  int len;
  for (int i = 0; i < nrep; i++) {
    e0 = initial_e_mean + gsl_ran_gaussian_ziggurat(gen, initial_e_sd);
    do {
      x[1] = gsl_ran_gaussian_ziggurat(gen, initial_x_sd);
      x[2] = gsl_ran_gaussian_ziggurat(gen, initial_x_sd);
    } while (x[1] * x[1] + x[2] * x[2] > nozzle_radius * nozzle_radius);
    p.reset(e0, x, w);
    len = p.simulate(dt, absorption_e, air_gap, gen, materials);
    grid.add(p.x, p.s, len);
  }
  grid.print(out_path);
  gsl_rng_free(gen);
  return 1;
}
