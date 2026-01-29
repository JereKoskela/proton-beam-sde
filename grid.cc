#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef GRID
#define GRID

struct Grid {

  Grid(const int n, const double dx_)
      : dx(dx_), x_pp(n), x_pm(n), x_mp(n), x_mm(n) {
    std::vector<double> tmp(1, 0);
    std::vector<std::vector<double>> tmp2(1, tmp);
    for (int i = 0; i < n; i++) {
      x_pp[i] = tmp2;
      x_pm[i] = tmp2;
      x_mp[i] = tmp2;
      x_mm[i] = tmp2;
    }
  }

  void print(const std::string filename) {
    std::ofstream outfile(filename);
    for (unsigned int ix = 0; ix < x_pp.size(); ix++) {
      for (unsigned int j = 0; j < x_pp[ix].size(); j++) {
        for (unsigned int k = 0; k < x_pp[ix][j].size(); k++) {
          if (x_pp[ix][j][k] > 0) {
            outfile << dx * ix << " " << dx * j << " " << dx * k << " "
                    << x_pp[ix][j][k] << std::endl;
          }
        }
      }
      for (unsigned int j = 0; j < x_pm[ix].size(); j++) {
        for (unsigned int k = 0; k < x_pm[ix][j].size(); k++) {
          if (x_pm[ix][j][k] > 0) {
            outfile << dx * ix << " " << dx * j << " " << -dx * (k + 1) << " "
                    << x_pm[ix][j][k] << std::endl;
          }
        }
      }
      for (unsigned int j = 0; j < x_mp[ix].size(); j++) {
        for (unsigned int k = 0; k < x_mp[ix][j].size(); k++) {
          if (x_mp[ix][j][k] > 0) {
            outfile << dx * ix << " " << -dx * (j + 1) << " " << dx * k << " "
                    << x_mp[ix][j][k] << std::endl;
          }
        }
      }
      for (unsigned int j = 0; j < x_mm[ix].size(); j++) {
        for (unsigned int k = 0; k < x_mm[ix][j].size(); k++) {
          if (x_mm[ix][j][k] > 0) {
            outfile << dx * ix << " " << -dx * (j + 1) << " " << -dx * (k + 1)
                    << " " << x_mm[ix][j][k] << std::endl;
          }
        }
      }
    }
    outfile.close();
    return;
  }

  void add(const std::vector<std::vector<double>> &y,
           const std::vector<double> &s, const int len) {
    int ix;
    unsigned int iy, iz;
    std::vector<double> tmp(1, 0);
    std::vector<std::vector<double>> tmp2(1, tmp);
    for (int i = 1; i < len; i++) {
      const double x = 0.5 * (y[i][0] + y[i - 1][0]);
      const double yy = 0.5 * (y[i][1] + y[i - 1][1]);
      const double zz = 0.5 * (y[i][2] + y[i - 1][2]);
      ix = floor(x / dx);
      iy = floor(fabs(yy) / dx);
      iz = floor(fabs(zz) / dx);
      if (ix >= 0) {
        if (ix >= int(x_pp.size())) {
          x_pp.resize(ix + 1, tmp2);
          x_pm.resize(ix + 1, tmp2);
          x_mp.resize(ix + 1, tmp2);
          x_mm.resize(ix + 1, tmp2);
        }
        if (yy >= 0 && zz >= 0) {
          if (iy >= x_pp[ix].size()) {
            x_pp[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_pp[ix][iy].size()) {
            x_pp[ix][iy].resize(iz + 1, 0);
          }
          x_pp[ix][iy][iz] += s[i];
        } else if (yy >= 0 && zz < 0) {
          if (iy >= x_pm[ix].size()) {
            x_pm[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_pm[ix][iy].size()) {
            x_pm[ix][iy].resize(iz + 1, 0);
          }
          x_pm[ix][iy][iz] += s[i];
        } else if (yy < 0 && zz >= 0) {
          if (iy >= x_mp[ix].size()) {
            x_mp[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_mp[ix][iy].size()) {
            x_mp[ix][iy].resize(iz + 1, 0);
          }
          x_mp[ix][iy][iz] += s[i];
        } else {
          if (iy >= x_mm[ix].size()) {
            x_mm[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_mm[ix][iy].size()) {
            x_mm[ix][iy].resize(iz + 1, 0);
          }
          x_mm[ix][iy][iz] += s[i];
        }
      }
    }
    return;
  }

  double dx;
  std::vector<std::vector<std::vector<double>>> x_pp, x_pm, x_mp, x_mm;
};

#endif
