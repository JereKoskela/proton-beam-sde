#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef CS
#define CS

struct CS_1d {

  CS_1d(const std::string filename) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    std::stringstream iss;

    getline(file, line);
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }

    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      rate.push_back(atof(token.c_str()));
    }
  }

  CS_1d(const CS_1d &other) : energy(other.energy), rate(other.rate) {}

  CS_1d() : energy(), rate() {}

  double evaluate(const double e) const {
    double ret = 0;
    int r;
    if (energy.size() > 0) {
      if (e <= energy[0]) {
        ret = rate[0];
      } else if (e >= energy.back()) {
        ret = rate.back();
      } else {
        r = std::distance(energy.begin(),
                          std::lower_bound(energy.begin(), energy.end(), e));
        ret = ((energy[r] - e) * rate[r - 1] + (e - energy[r - 1]) * rate[r]) /
              (energy[r] - energy[r - 1]);
      }
    }
    return ret;
  }

  std::vector<double> energy, rate;
};

struct CS_3d {

  CS_3d(const std::string filename) : energy(), exit_energy(), cdf(), rvalue() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_energy.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss5;
      iss5 << line;
      while (getline(iss5, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      rvalue.push_back(tmp_vec);
    }
    file.close();
  }

  CS_3d(const CS_3d &other)
      : energy(other.energy), exit_energy(other.exit_energy), cdf(other.cdf),
        rvalue(other.rvalue) {}

  CS_3d() : energy(), exit_energy(), cdf(), rvalue() {}

  void sample_from_energy_index(const double energy_index, const double u,
                                double &out_energy_cm,
                                double &out_rvalue) const {
    double diff = 0;
    int density_index =
        std::distance(cdf[energy_index].begin(),
                      std::lower_bound(cdf[energy_index].begin(),
                                       cdf[energy_index].end(), u));
    if (density_index == 0 || density_index == int(cdf[energy_index].size())) {
      density_index =
          std::min(density_index, int(cdf[energy_index].size()) - 1);
      out_energy_cm = exit_energy[energy_index][density_index];
      out_rvalue = rvalue[energy_index][density_index];
    } else {
      diff = (u - cdf[energy_index][density_index - 1]) /
             (cdf[energy_index][density_index] -
              cdf[energy_index][density_index - 1]);
      out_energy_cm = exit_energy[energy_index][density_index] * diff +
                      exit_energy[energy_index][density_index - 1] * (1 - diff);
      out_rvalue = rvalue[energy_index][density_index] * diff +
                   rvalue[energy_index][density_index - 1] * (1 - diff);
    }
    return;
  }

  void sample(const double e, double &r, double &out_e_cm, gsl_rng *gen) const {
    int energy_index = std::distance(
        energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    double u = gsl_rng_uniform(gen);
    double out_energy_cm;
    double out_energy_cm_2;
    double out_rvalue;
    double out_rvalue_2;
    double diff;
    if (energy_index == 0 || energy_index == int(energy.size())) {
      energy_index = std::min(energy_index, int(energy.size()) - 1);
      sample_from_energy_index(energy_index, u, out_energy_cm, out_rvalue);
    } else {
      sample_from_energy_index(energy_index, u, out_energy_cm, out_rvalue);
      sample_from_energy_index(energy_index - 1, u, out_energy_cm_2, out_rvalue_2);
      diff = (e - energy[energy_index - 1]) /
             (energy[energy_index] - energy[energy_index - 1]);
      out_energy_cm = out_energy_cm * diff + out_energy_cm_2 * (1 - diff);
      out_rvalue = out_rvalue * diff + out_rvalue_2 * (1 - diff);
    }
    r = out_rvalue;
    out_e_cm = out_energy_cm;
    return;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_energy, cdf, rvalue;
};

struct CS_2d {

  CS_2d(const std::string filename) : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_angle.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
  }

  CS_2d(const CS_2d &other)
      : energy(other.energy), exit_angle(other.exit_angle), cdf(other.cdf) {}

  CS_2d() : energy(), exit_angle(), cdf() {}

  double sample_from_energy_index(const double energy_index,
                                  const double u) const {
    double ret = 0;
    double diff = 0;
    int density_index =
        std::distance(cdf[energy_index].begin(),
                      std::lower_bound(cdf[energy_index].begin(),
                                       cdf[energy_index].end(), u));
    if (density_index == 0 || density_index == int(cdf[energy_index].size())) {
      density_index =
          std::min(density_index, int(cdf[energy_index].size()) - 1);
      ret = exit_angle[energy_index][density_index];
    } else {
      diff = (u - cdf[energy_index][density_index - 1]) /
             (cdf[energy_index][density_index] -
              cdf[energy_index][density_index - 1]);
      ret = exit_angle[energy_index][density_index] * diff +
            exit_angle[energy_index][density_index - 1] * (1 - diff);
    }
    return ret;
  }

  double sample(const double e, gsl_rng *gen) const {
    int energy_index = std::distance(
        energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    double u = gsl_rng_uniform(gen);
    double out_angle_cm;
    double out_angle_cm_2;
    double diff;
    if (energy_index == 0 || energy_index == int(energy.size())) {
      energy_index = std::min(energy_index, int(energy.size()) - 1);
      out_angle_cm = sample_from_energy_index(energy_index, u);
    } else {
      out_angle_cm = sample_from_energy_index(energy_index, u);
      out_angle_cm_2 = sample_from_energy_index(energy_index - 1, u);
      diff = (e - energy[energy_index - 1]) /
             (energy[energy_index] - energy[energy_index - 1]);
      out_angle_cm = out_angle_cm * diff + out_angle_cm_2 * (1 - diff);
    }
    return out_angle_cm;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_angle, cdf;
};

#endif
