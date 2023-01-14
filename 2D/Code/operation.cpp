#ifndef OPERATION_CPP

#include "operation.h"
#include "flux.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string> 

using namespace std;

// Lecture des variables 
void lecture(int & n, double & tf, string & test, double & disc_x, double & disc_y, vector<double> & rho, vector<double> & u, vector<double> & v, vector<double> & p, int & choix_schema, int & choix_ordre){

    string ignore;  

    std::ifstream monFlux("Initialisation/initialisation.data", std::ios::in);

    if (monFlux){

        monFlux >> ignore; 
        monFlux >> choix_schema; 
        monFlux >> ignore; 
        monFlux >> choix_ordre;
        monFlux >> ignore; 
        monFlux >> n; 
        monFlux >> ignore; 
        monFlux >> test; 

        monFlux.close();
    }

    std::ifstream monFlux1("Initialisation/Conf_" + test + ".data", std::ios::in);

    if (monFlux1){
        monFlux1 >> ignore; 
        monFlux1 >> tf; 
        monFlux1 >> ignore; 
        monFlux1 >> disc_x;
        monFlux1 >> ignore; 
        monFlux1 >> disc_y; 
        monFlux1 >> ignore; 
        monFlux1 >> rho[0]; 
        monFlux1 >> ignore; 
        monFlux1 >> u[0]; 
        monFlux1 >> ignore; 
        monFlux1 >> v[0]; 
        monFlux1 >> ignore; 
        monFlux1 >> p[0];
        monFlux1 >> ignore; 
        monFlux1 >> rho[1]; 
        monFlux1 >> ignore; 
        monFlux1 >> u[1]; 
        monFlux1 >> ignore; 
        monFlux1 >> v[1]; 
        monFlux1 >> ignore; 
        monFlux1 >> p[1];
        monFlux1 >> ignore; 
        monFlux1 >> rho[2]; 
        monFlux1 >> ignore; 
        monFlux1 >> u[2]; 
        monFlux1 >> ignore; 
        monFlux1 >> v[2]; 
        monFlux1 >> ignore; 
        monFlux1 >> p[2];
        monFlux1 >> ignore; 
        monFlux1 >> rho[3]; 
        monFlux1 >> ignore; 
        monFlux1 >> u[3]; 
        monFlux1 >> ignore; 
        monFlux1 >> v[3]; 
        monFlux1 >> ignore; 
        monFlux1 >> p[3];
        
        monFlux1.close();
    } 
}

// Initialisation de U (matrice; ligne = 3 composantes; colonne = n uj)
void U0(vector<vector<double>> & U_rho, vector<vector<double>> & U_rhou, vector<vector<double>> & U_rhov, vector<vector<double>> & U_E, vector<double> rho, vector<double> u, vector<double> v, vector<double> E, double dx, double x0, double dy, double y0, double disc_x, double disc_y, int n){
    
    double xi, yi; // Maillage cartésien uniforme : yj = yi

    for (int i = 0; i < n; i++)
    {
        xi = x0 + i*dx;

        for (int j = 0; j < n; j++)
        {
            yi = y0 + j*dy;
            if ((xi > disc_x) && (yi > disc_y))
            {
                U_rho[i][j] = rho[0];
                U_rhou[i][j] = rho[0]*u[0];
                U_rhov[i][j] = rho[0]*v[0];
                U_E[i][j] = E[0];
            }
            else if ((xi <= disc_x) && (yi > disc_y))
            {
                U_rho[i][j] = rho[1];
                U_rhou[i][j] = rho[1]*u[1];
                U_rhov[i][j] = rho[1]*v[1];
                U_E[i][j] = E[1];
                
            }
            else if ((xi <= disc_x) && (yi <= disc_y))
            {
                U_rho[i][j] = rho[2];
                U_rhou[i][j] = rho[2]*u[2];
                U_rhov[i][j] = rho[2]*v[2];
                U_E[i][j] = E[2];
            }
        
            else
            {
                U_rho[i][j] = rho[3];
                U_rhou[i][j] = rho[3]*u[3];
                U_rhov[i][j] = rho[3]*v[3];
                U_E[i][j] = E[3];
            }
        }
        
    }
}


// Creation du vecteur Fx
vector<double> Fx(vector<double> Ui, double p, double E){

    vector<double> Fxi(4,0); 

    Fxi[0] = Ui[1];
    Fxi[1] = pow(Ui[1],2)/Ui[0] + p;
    Fxi[2] = Ui[1]*Ui[2]/Ui[0];
    Fxi[3] = (Ui[1]/Ui[0])*(E+p);

    return Fxi;
}

// Creation du vecteur Fy
vector<double> Fy(vector<double> Ui, double p, double E){

    vector<double> Fyi(4,0); 

    Fyi[0] = Ui[2];
    Fyi[1] = Ui[1]*Ui[2]/Ui[0];
    Fyi[2] = pow(Ui[2],2)/Ui[0] + p;
    Fyi[3] = (Ui[2]/Ui[0])*(E+p);

    return Fyi;
}


// Calcul de P
double pression(vector<double> Ui, double gamma){
    double p;
    p = (Ui[3] - (Ui[0]*(pow(Ui[1]/Ui[0],2)+pow(Ui[2]/Ui[0],2)))/2 ) * (gamma - 1.);
    return p;
}

// limiter
double phi(double theta){

    double theta_lim1, theta_lim;

    // Superbee
    theta_lim1 = max(0., min(1., 2*theta));
    theta_lim = max(theta_lim1, min(2., theta));
    // Van Leer
    // theta_lim = (theta + abs(theta))/(1+theta);
    // Minmod
    // if (theta < 0)
    // {
    //     theta_lim = 0;
    // }
    // else if (1 < abs(theta))
    // {
    //     theta_lim = 1.;
    // }
    // else{
    //     theta_lim = theta;
    // }
    

    return theta_lim;
}

void Grand_gamma(vector<double> & gamma10, vector<double> & gamma01, double dx, vector<double> u_1, vector<double> u_2, vector<double> u_3, vector<double> u_4, vector<double> u_5, vector<double> u_6, vector<double> u_7, vector<double> u_8, vector<double> u_9){

    for (int k = 0; k < u_1.size(); k++)
    {
        gamma10[k] = (1./(6*dx)) * (u_3[k] - u_1[k] + u_6[k] - u_4[k] + u_9[k] - u_7[k]);
        gamma01[k] = (1./(6*dx)) * (u_7[k] - u_1[k] + u_8[k] - u_2[k] + u_9[k] - u_3[k]);
    }
}

double minmod(double sigma_centre, double sigma_droite, double sigma_gauche){

    double theta_lim1, theta_lim;

    // Minmod
    if ((sigma_centre < 0) && (sigma_droite < 0) && (sigma_gauche < 0))
    {
        theta_lim1 = max(sigma_centre, sigma_droite);
        theta_lim = max(theta_lim1, sigma_gauche);
    }
    else if ((sigma_centre > 0) && (sigma_droite > 0) && (sigma_gauche > 0))
    {
        theta_lim1 = min(sigma_centre, sigma_droite);
        theta_lim = min(theta_lim1, sigma_gauche);
    }
    else{
        theta_lim = 0.;
    }
    
    return theta_lim;
}

vector<double> ft(double dx, double dy, vector<double> u_mimj, vector<double> u_imj, vector<double> u_pimj, vector<double> u_mij, vector<double> u_ij, vector<double> u_pij, vector<double> u_mipj, vector<double> u_ipj, vector<double> u_pipj, vector<double> u_ppij, vector<double> u_mmij, vector<double> u_ippj, vector<double> u_immj, double gamma, int choix_schema, int choix_ordre){

    vector<double> res(u_ij.size(),0);
    vector<double> flux_mdemi_x(u_ij.size(),0), flux_pdemi_x(u_ij.size(),0);
    vector<double> flux_mdemi_y(u_ij.size(),0), flux_pdemi_y(u_ij.size(),0);

    vector<double> pij_xpi(u_ij.size(),0), ppij_xpi(u_ij.size(),0), pij_xmi(u_ij.size(),0), pmij_xmi(u_ij.size(),0), pij_ypj(u_ij.size(),0), pipj_ypj(u_ij.size(),0), pimj_ymj(u_ij.size(),0), pij_ymj(u_ij.size(),0);
    vector<double> gamma10(u_ij.size(),0), gamma01(u_ij.size(),0);

    vector<double> sigma_centre_xi(4,0), sigma_droite_xi(4,0), sigma_gauche_xi(4,0), sigma_centre_xmi(4,0), sigma_droite_xmi(4,0), sigma_gauche_xmi(4,0), sigma_centre_xpi(4,0), sigma_droite_xpi(4,0), sigma_gauche_xpi(4,0);
    vector<double> sigma_centre_yj(4,0), sigma_haut_yj(4,0), sigma_bas_yj(4,0), sigma_centre_ymj(4,0), sigma_haut_ymj(4,0), sigma_bas_ymj(4,0), sigma_centre_ypj(4,0), sigma_haut_ypj(4,0), sigma_bas_ypj(4,0);
    vector<double> pix(4,0), pmix(4,0), ppix(4,0), pjy(4,0), pmjy(4,0), ppjy(4,0);
    vector<double> theta_pi_x(4,0), theta_i_x(4,0), theta_mi_x(4,0), theta_pj_y(4,0), theta_j_y(4,0), theta_mj_y(4,0);


    // Grand_gamma(gamma10, gamma01, dx, u_mimj, u_imj, u_pimj, u_mij, u_ij, u_pij, u_mipj, u_ipj, u_pipj);
    for (int k = 0; k < u_ij.size(); k++)
    {
        // pi x
        sigma_centre_xi[k] = (u_pij[k] - u_mij[k])/(2*dx);
        sigma_droite_xi[k] = (u_pij[k] - u_ij[k])/(dx);
        sigma_gauche_xi[k] = (u_ij[k] - u_mij[k])/(dx);

        // pi+1 x
        sigma_centre_xpi[k] = (u_ppij[k] - u_ij[k])/(2*dx);
        sigma_droite_xpi[k] = (u_ppij[k] - u_pij[k])/(dx);
        sigma_gauche_xpi[k] = (u_pij[k] - u_ij[k])/(dx);

        // pi-1 y
        sigma_centre_xmi[k] = (u_ij[k] - u_mmij[k])/(2*dx);
        sigma_droite_xmi[k] = 2*(u_ij[k] - u_mij[k])/(dx);
        sigma_gauche_xmi[k] = 2*(u_mij[k] - u_mmij[k])/(dx);

        // pj
        sigma_centre_yj[k] = (u_ipj[k] - u_imj[k])/(2*dx);
        sigma_haut_yj[k] = (u_ipj[k] - u_ij[k])/(dx);
        sigma_bas_yj[k] = (u_ij[k] - u_imj[k])/(dx);

        // pj+1 y
        sigma_centre_ypj[k] = (u_ippj[k] - u_ij[k])/(2*dx);
        sigma_haut_ypj[k] = (u_ippj[k] - u_ipj[k])/(dx);
        sigma_bas_ypj[k] = (u_ipj[k] - u_ij[k])/(dx);

        // pj-1 y
        sigma_centre_ymj[k] = (u_ij[k] - u_immj[k])/(2*dx);
        sigma_haut_ymj[k] = 2*(u_ij[k] - u_imj[k])/(dx);
        sigma_bas_ymj[k] = 2*(u_imj[k] - u_immj[k])/(dx);
    }

    // // Calcul des theta
    for (int k = 0; k < u_ij.size(); k++)
    {
    // MUSCL
        theta_mi_x[k] = min(1.,(u_mij[k] - u_mmij[k])/(u_ij[k] - u_mij[k]));
        theta_i_x[k] = min(1.,(u_ij[k] - u_mij[k])/(u_pij[k] - u_ij[k]));
        theta_pi_x[k] = min(1.,(u_pij[k] - u_ij[k])/(u_ppij[k] - u_pij[k]));
        theta_mj_y[k] = min(1.,(u_imj[k] - u_immj[k])/(u_ij[k] - u_imj[k]));
        theta_j_y[k] = min(1.,(u_ij[k] - u_imj[k])/(u_ipj[k] - u_ij[k]));
        theta_pj_y[k] = min(1.,(u_ipj[k] - u_ij[k])/(u_ippj[k] - u_ipj[k]));
    }

    // Calcul des coeff (gamma)
    for (int k = 0; k < u_ij.size(); k++)
    {
        // pix[k] = minmod(sigma_centre_xi[k], phi(theta_i_x[k]) * sigma_droite_xi[k], phi(theta_i_x[k]) * sigma_gauche_xi[k]);
        // ppix[k] = minmod(sigma_centre_xpi[k], phi(theta_pi_x[k]) * sigma_droite_xpi[k], phi(theta_pi_x[k]) * sigma_gauche_xpi[k]);
        // pmix[k] = minmod(sigma_centre_xmi[k], phi(theta_mi_x[k]) * sigma_droite_xmi[k], phi(theta_mi_x[k]) * sigma_gauche_xmi[k]);
        // pjy[k] = minmod(sigma_centre_yj[k], phi(theta_j_y[k]) * sigma_haut_yj[k], phi(theta_j_y[k]) * sigma_bas_yj[k]);
        // ppjy[k] = minmod(sigma_centre_ypj[k], phi(theta_pj_y[k]) * sigma_haut_ypj[k], phi(theta_pj_y[k]) * sigma_bas_ypj[k]);
        // pmjy[k] = minmod(sigma_centre_ymj[k], phi(theta_mj_y[k]) * sigma_haut_ymj[k], phi(theta_mj_y[k]) * sigma_bas_ymj[k]);
        pix[k] = minmod(sigma_centre_xi[k], sigma_droite_xi[k], sigma_gauche_xi[k]);
        ppix[k] = minmod(sigma_centre_xpi[k], sigma_droite_xpi[k], sigma_gauche_xpi[k]);
        pmix[k] = minmod(sigma_centre_xmi[k], sigma_droite_xmi[k], sigma_gauche_xmi[k]);
        pjy[k] = minmod(sigma_centre_yj[k], sigma_haut_yj[k], sigma_bas_yj[k]);
        ppjy[k] = minmod(sigma_centre_ypj[k], sigma_haut_ypj[k], sigma_bas_ypj[k]);
        pmjy[k] = minmod(sigma_centre_ymj[k], sigma_haut_ymj[k], sigma_bas_ymj[k]);
    }

    for (int k = 0; k < u_ij.size(); k++)
    {
        // pij_xpi[k] = u_ij[k] + gamma10[k] * (dx/2.);
        // ppij_xpi[k] = u_pij[k] - gamma10[k] * (dx/2.);
        // pij_xmi[k] = u_ij[k] - gamma10[k] * (dx/2.);
        // pmij_xmi[k] = u_mij[k] + gamma10[k] * (dx/2.);
        // pij_ypj[k] = u_ij[k] + gamma01[k] * (dx/2.);
        // pipj_ypj[k] = u_ipj[k] - gamma01[k] * (dx/2.);
        // pimj_ymj[k] = u_imj[k] + gamma01[k] * (dx/2.);
        // pij_ymj[k] = u_ij[k] - gamma01[k] * (dx/2.);
        pij_xpi[k] = u_ij[k] + pix[k] * (dx/2.);
        ppij_xpi[k] = u_pij[k] - ppix[k] * (dx/2.);
        pij_xmi[k] = u_ij[k] - pix[k] * (dx/2.);
        pmij_xmi[k] = u_mij[k] + pmix[k] * (dx/2.);
        pij_ypj[k] = u_ij[k] + pjy[k] * (dx/2.);
        pipj_ypj[k] = u_ipj[k] - ppjy[k] * (dx/2.);
        pimj_ymj[k] = u_imj[k] + pmjy[k] * (dx/2.);
        pij_ymj[k] = u_ij[k] - pjy[k] * (dx/2.);
    }
    

    if ((choix_schema == 1)  && (choix_ordre == 1))
    {
        flux_pdemi_x = flux_rusanov(u_ij, u_pij, gamma, 1);
        flux_mdemi_x = flux_rusanov(u_mij, u_ij, gamma, 1);
        flux_pdemi_y = flux_rusanov(u_ij, u_ipj, gamma, 2);
        flux_mdemi_y = flux_rusanov(u_imj, u_ij, gamma, 2);
    }
    else if ((choix_schema == 2)  && (choix_ordre == 1))
    {
        flux_pdemi_x = flux_HLL(u_ij, u_pij, gamma, 1);
        flux_mdemi_x = flux_HLL(u_mij, u_ij, gamma, 1);
        flux_pdemi_y = flux_HLL(u_ij, u_ipj, gamma, 2);
        flux_mdemi_y = flux_HLL(u_imj, u_ij, gamma, 2);
    }
    else if ((choix_schema == 3)  && (choix_ordre == 1))
    {
        flux_pdemi_x = flux_HLLC(u_ij, u_pij, gamma, 1);
        flux_mdemi_x = flux_HLLC(u_mij, u_ij, gamma, 1);
        flux_pdemi_y = flux_HLLC(u_ij, u_ipj, gamma, 2);
        flux_mdemi_y = flux_HLLC(u_imj, u_ij, gamma, 2);
    }
    else if ((choix_schema == 1)  && (choix_ordre == 2))
    {
        flux_pdemi_x = flux_rusanov(pij_xpi, ppij_xpi, gamma, 1);
        flux_mdemi_x = flux_rusanov(pmij_xmi, pij_xmi, gamma, 1);
        flux_pdemi_y = flux_rusanov(pij_ypj, pipj_ypj, gamma, 2);
        flux_mdemi_y = flux_rusanov(pimj_ymj, pij_ymj, gamma, 2);
    }
    else if ((choix_schema == 2)  && (choix_ordre == 2))
    {
        flux_pdemi_x = flux_HLL(pij_xpi, ppij_xpi, gamma, 1);
        flux_mdemi_x = flux_HLL(pmij_xmi, pij_xmi, gamma, 1);
        flux_pdemi_y = flux_HLL(pij_ypj, pipj_ypj, gamma, 2);
        flux_mdemi_y = flux_HLL(pimj_ymj, pij_ymj, gamma, 2);
    }
    else if ((choix_schema == 3)  && (choix_ordre == 2))
    {
        flux_pdemi_x = flux_HLLC(pij_xpi, ppij_xpi, gamma, 1);
        flux_mdemi_x = flux_HLLC(pmij_xmi, pij_xmi, gamma, 1);
        flux_pdemi_y = flux_HLLC(pij_ypj, pipj_ypj, gamma, 2);
        flux_mdemi_y = flux_HLLC(pimj_ymj, pij_ymj, gamma, 2);
    }
    else
    {
        cout << "Pb dans ft, on prend rusanov ordre 1" << endl;
        flux_pdemi_x = flux_rusanov(u_ij, u_pij, gamma, 1);
        flux_mdemi_x = flux_rusanov(u_mij, u_ij, gamma, 1);
        flux_pdemi_y = flux_rusanov(u_ij, u_ipj, gamma, 2);
        flux_mdemi_y = flux_rusanov(u_imj, u_ij, gamma, 2);
    }

    for (int k = 0; k < u_ij.size(); k++)
    {
        res[k] = - (1./dx) * (flux_pdemi_x[k] - flux_mdemi_x[k]) - (1./dy) * (flux_pdemi_y[k] - flux_mdemi_y[k]);
    }

    return res;
}

// Calcul erreur
void error(vector<vector<double>> U_rho, vector<vector<double>> U_rhou, vector<vector<double>> U_rhov, vector<vector<double>> U_E, string result, int n, double dx, double dy, double gamma, double b, int fichier, string test){

    string ignore;
    double delta_x(0), delta_y(0), density(0), vx(0), vy(0), pressure(0), energy(0), e_int(0);
    double error_rho(0), error_u(0), error_v(0), error_p(0), error_e(0);
    int kx, ky;
    double xi, yi, x0(0), y0(0);
    vector<double> U_int(4,0);

    std::ifstream monFlux("Resultats/Conf_" + test + "/tf_HLLC_fin_1000.txt", std::ios::in);

    for (int i = 1; i < n/fichier-1; i++)
    // for (int i = 0; i < n; i++)
    {
        kx=i*fichier;
        xi = x0 + kx*delta_x;
        for (int j = 1; j < n/fichier-1; j++)
        // for (int j = 0; j < n; j++)
        {
            ky=j*fichier;
            yi = y0 + ky*delta_y;

            if (monFlux){
            
            monFlux >> delta_x;
            monFlux >> delta_y;
            monFlux >> density;
            monFlux >> vx;
            monFlux >> vy;
            monFlux >> pressure;
            monFlux >> energy;
            }
        // cout << delta_x << " " << delta_y << " " << density << " " << vx << " " << vy << " " << pressure << " " << energy << endl;
        U_int[0] = U_rho[kx][ky];
        U_int[1] = U_rhou[kx][ky];
        U_int[2] = U_rhov[kx][ky];
        U_int[3] = U_E[kx][ky];

        error_rho += pow(U_rho[kx][ky] - density, 2);
        error_u += pow(U_rhou[kx][ky]/U_rho[kx][ky] - vx, 2);
        error_v += pow(U_rhov[kx][ky]/U_rho[kx][ky] - vy, 2);
        error_p += pow(pression(U_int, gamma) - pressure, 2);
        e_int = (pression(U_int, gamma)*(1-b*U_int[0]))/((gamma-1)*U_int[0]);
        error_e += pow(e_int - energy, 2);

        }
    }

    ofstream mon_flux1;
    mon_flux1.open("Resultats/Conf_" + test + "/erreur_" + result, ios::out);
    mon_flux1 << dx*dy << " " << sqrt(error_rho*dx*dy) << " " << sqrt(error_u*dx*dy) << " " << sqrt(error_v*dx*dy) << " " << sqrt(error_p*dx*dy) << " " << sqrt(error_e*dx*dy) << endl;
    mon_flux1.close();
    
    monFlux.close();

}

#define OPERATION_CPP
#endif