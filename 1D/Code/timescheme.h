#ifndef TIMESCHEME_H

#include <iostream>
#include <vector>
#include "flux.h"
#include "operation.h"

void Euler(std::vector<double> &ui, std::vector<double> ummi_old, std::vector<double> umi_old, std::vector<double> ui_old, std::vector<double> upi_old, std::vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre);
void RK2_SSP(std::vector<double> &ui, std::vector<double> ummi_old, std::vector<double> umi_old, std::vector<double> ui_old, std::vector<double> upi_old, std::vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre);
void RK2_SSP_demi(std::vector<double> &ui, std::vector<double> ummi_old, std::vector<double> umi_old, std::vector<double> ui_old, std::vector<double> upi_old, std::vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre);
void RK4(std::vector<double> &ui, std::vector<double> ummi_old, std::vector<double> umi_old, std::vector<double> ui_old, std::vector<double> upi_old, std::vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre);


#define TIMESCHEME_H
#endif