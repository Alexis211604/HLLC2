#ifndef FLUX_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include "operation.h"

void compute_star(std::vector<double> U_l, std::vector<double> U_r, double & p_star, double & u_star, double gamma, int x_or_y);
void vect_S(std::vector<double> U_l, std::vector<double> U_r, double & S_l, double & S_m, double & S_r, double gamma, int x_or_y);
std::vector<double> flux_rusanov(std::vector<double> U_l, std::vector<double> U_r, double gamma, int x_or_y);
std::vector<double> flux_HLL(std::vector<double> U_l, std::vector<double> U_r, double gamma, int x_or_y);
std::vector<double> flux_HLLC(std::vector<double> U_l, std::vector<double> U_r, double gamma, int x_or_y);

#define FLUX_H
#endif