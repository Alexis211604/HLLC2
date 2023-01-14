#ifndef OPERATION_H

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>

void lecture(int & n, double & tf, std::string & test, double & disc_x, double & disc_y, std::vector<double> & rho, std::vector<double> & u, std::vector<double> & v, std::vector<double> & p, int & choix_schema, int & choix_ordre);
void U0(std::vector<std::vector<double>> & U_rho, std::vector<std::vector<double>> & U_rhou, std::vector<std::vector<double>> & U_rhov, std::vector<std::vector<double>> & U_E, std::vector<double> rho, std::vector<double> u, std::vector<double> v, std::vector<double> E, double dx, double x0, double dy, double y0, double disc_x, double disc_y, int n);
std::vector<double> Fx(std::vector<double> Ui, double p, double E);
std::vector<double> Fy(std::vector<double> Ui, double p, double E);
double pression(std::vector<double> Ui, double gamma);
double phi(double theta);
void Grand_gamma(std::vector<double> & gamma10, std::vector<double> & gamma01, double dx, std::vector<double> u_1, std::vector<double> u_2, std::vector<double> u_3, std::vector<double> u_4, std::vector<double> u_5, std::vector<double> u_6, std::vector<double> u_7, std::vector<double> u_8, std::vector<double> u_9);
double minmod(double sigma_centre, double sigma_droite, double sigma_gauche);
std::vector<double> ft(double dx, double dy, std::vector<double> u_mimj, std::vector<double> u_imj, std::vector<double> u_pimj, std::vector<double> u_mij, std::vector<double> u_ij, std::vector<double> u_pij, std::vector<double> u_mipj, std::vector<double> u_ipj, std::vector<double> u_pipj, std::vector<double> u_ppij, std::vector<double> u_mmij, std::vector<double> u_ippj, std::vector<double> u_immj, double gamma, int choix_schema, int choix_ordre);
void error(std::vector<std::vector<double>> U_rho, std::vector<std::vector<double>> U_rhou, std::vector<std::vector<double>> U_rhov, std::vector<std::vector<double>> U_E, std::string result, int n, double dx, double dy, double gamma, double b, int fichier, std::string test);



#define OPERATION_H
#endif