#ifndef TIMESCHEME_H

#include <iostream>
#include <vector>
#include "flux.h"
#include "operation.h"

void Euler(std::vector<double> &uij, std::vector<double> u_mimj, std::vector<double> u_imj, std::vector<double> u_pimj, std::vector<double> u_mij, std::vector<double> u_ij, std::vector<double> u_pij, std::vector<double> u_mipj, std::vector<double> u_ipj, std::vector<double> u_pipj, std::vector<double> u_ppij, std::vector<double> u_mmij, std::vector<double> u_ippj, std::vector<double> u_immj, double dt, double dx, double dy, double gamma, int choix_schema, int choix_ordre);
void RK2_SSP(std::vector<double> &uij, std::vector<double> u_mimj, std::vector<double> u_imj, std::vector<double> u_pimj, std::vector<double> u_mij, std::vector<double> u_ij, std::vector<double> u_pij, std::vector<double> u_mipj, std::vector<double> u_ipj, std::vector<double> u_pipj, std::vector<double> u_ppij, std::vector<double> u_mmij, std::vector<double> u_ippj, std::vector<double> u_immj, double dt, double dx, double dy, double gamma, int choix_schema, int choix_ordre);
void RK4(std::vector<double> &uij, std::vector<double> u_mimj, std::vector<double> u_imj, std::vector<double> u_pimj, std::vector<double> u_mij, std::vector<double> u_ij, std::vector<double> u_pij, std::vector<double> u_mipj, std::vector<double> u_ipj, std::vector<double> u_pipj, std::vector<double> u_ppij, std::vector<double> u_mmij, std::vector<double> u_ippj, std::vector<double> u_immj, double dt, double dx, double dy, double gamma, int choix_schema, int choix_ordre);


#define TIMESCHEME_H
#endif