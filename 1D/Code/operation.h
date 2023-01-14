#ifndef OPERATION_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

void lecture(int & n, double & tf, std::string & test, double & disc, double & rho_l, double & u_l, double & p_l, double & rho_r, double & u_r, double & p_r, int & choix_schema, int & choix_ordre);
std::vector<double> build_U(double rho, double u, double E);
void U0(std::vector<std::vector<double>> & U, double rho_l, double u0_l, double E_l, double rho_r, double u0_r, double E_r, double dx, double x0, double disc, int n);
std::vector<double> F(std::vector<double> Ui, double p, double E);
double pression(std::vector<double> Ui, double gamma);
double minmod(double sigma_centre, double sigma_droite, double sigma_gauche);
double phi(double theta);
std::vector<double> ft(double dx, std::vector<double> ummi_old, std::vector<double> umi_old, std::vector<double> ui_old, std::vector<double> upi_old, std::vector<double> uppi_old, double gamma, int choix_schema, int choix_ordre);
void error(std::vector<std::vector<double>> U, std::string result, int n, double dx, double gamma, double b, int fichier, std::string test);

#define OPERATION_H
#endif