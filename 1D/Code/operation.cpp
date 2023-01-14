#ifndef OPERATION_CPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string> 
#include "flux.h"
#include "operation.h"

using namespace std;

// Lecture des variables
void lecture(int & n, double & tf, string & test, double & disc, double & rho_l, double & u_l, double & p_l, double & rho_r, double & u_r, double & p_r, int & choix_schema, int & choix_ordre){

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

    std::ifstream monFlux1("Initialisation/test_" + test + ".data", std::ios::in);

    if (monFlux1){

        monFlux1 >> ignore;
        monFlux1 >> tf;
        monFlux1 >> ignore;
        monFlux1 >> disc;
        monFlux1 >> ignore;
        monFlux1 >> rho_l;
        monFlux1 >> ignore;
        monFlux1 >> u_l;
        monFlux1 >> ignore;
        monFlux1 >> p_l;
        monFlux1 >> ignore;
        monFlux1 >> rho_r;
        monFlux1 >> ignore;
        monFlux1 >> u_r;
        monFlux1 >> ignore;
        monFlux1 >> p_r;

        monFlux1.close();
    }
}



// Creation du vecteur U
vector<double> build_U(double rho, double u, double E){

    vector<double> Ui(3,0);

    Ui[0] = rho;
    Ui[1] = rho*u;
    Ui[2] = E;

    return Ui;
}

// Initialisation de U (matrice; ligne = 3 composantes; colonne = n uj)
void U0(vector<vector<double>> & U, double rho_l, double u0_l, double E_l, double rho_r, double u0_r, double E_r, double dx, double x0, double disc, int n){

    U.resize(n,vector<double> (3));
    vector<double> U_l(3,0), U_r(3,0);

    double xi;

    U_l = build_U(rho_l, u0_l, E_l);
    U_r = build_U(rho_r, u0_r, E_r);

    for (int i = 0; i < n; i++)
    {
        xi = x0 + i*dx;

        if (xi < disc)
        {
            for (int k = 0; k < 3; k++)
            {
                U[i][k] = U_l[k];
            }

        }
        else
        {
            for (int k = 0; k < 3; k++)
            {
                U[i][k] = U_r[k];
            }
        }
    }
}


// Creation du vecteur F
vector<double> F(vector<double> Ui, double p, double E){

    vector<double> Fi(3,0);

    Fi[0] = Ui[1];
    Fi[1] = pow(Ui[1],2)/Ui[0] + p;
    Fi[2] = (Ui[1]/Ui[0])*(E+p);

    return Fi;
}


// Calcul de P
double pression(vector<double> Ui, double gamma){
    double p;
    p = (Ui[2] - (Ui[0]*pow(Ui[1]/Ui[0],2))/2 ) * (gamma - 1.);
    return p;
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


vector<double> ft(double dx, vector<double> ummi_old, vector<double> umi_old, vector<double> ui_old, vector<double> upi_old, vector<double> uppi_old, double gamma, int choix_schema, int choix_ordre){

    vector<double> res(3,0);
    vector<double> flux_mdemi(3,0), flux_pdemi(3,0);
    vector<double> pi(3,0), pmi(3,0), ppi(3,0);
    vector<double> sigma_centre_i(3,0), sigma_droite_i(3,0), sigma_gauche_i(3,0), sigma_centre_mi(3,0), sigma_droite_mi(3,0), sigma_gauche_mi(3,0), sigma_centre_pi(3,0), sigma_droite_pi(3,0), sigma_gauche_pi(3,0);
    vector<double> ui_ppi(3,0), ui_mpi(3,0), upi_mppi(3,0), umi_ppmi(3,0);
    vector<double> theta_pi(3,0), theta_i(3,0), theta_mi(3,0);

    for (int k = 0; k < ui_old.size(); k++)
    {
        // pi
        sigma_centre_i[k] = (upi_old[k] - umi_old[k])/(2*dx);
        sigma_droite_i[k] = (upi_old[k] - ui_old[k])/(dx);
        sigma_gauche_i[k] = (ui_old[k] - umi_old[k])/(dx);

        // pi+1
        sigma_centre_pi[k] = (uppi_old[k] - ui_old[k])/(2*dx);
        sigma_droite_pi[k] = (uppi_old[k] - upi_old[k])/(dx);
        sigma_gauche_pi[k] = (upi_old[k] - ui_old[k])/(dx);

        // pi-1
        sigma_centre_mi[k] = (ui_old[k] - ummi_old[k])/(2*dx);
        sigma_droite_mi[k] = (ui_old[k] - umi_old[k])/(dx);
        sigma_gauche_mi[k] = (umi_old[k] - ummi_old[k])/(dx);
    }

    // // Calcul des theta
    for (int k = 0; k < ui_old.size(); k++)
    {
    // MUSCL
        theta_mi[k] = min(1.,(umi_old[k] - ummi_old[k])/(ui_old[k] - umi_old[k]));
        theta_i[k] = min(1.,(ui_old[k] - umi_old[k])/(upi_old[k] - ui_old[k]));
        theta_pi[k] = min(1.,(upi_old[k] - ui_old[k])/(uppi_old[k] - upi_old[k]));
        // if ((ui_old[k] - umi_old[k]) < 1e-16) //theta_i
        // {
        //     theta_mi[k] = 0.; //0.35;
        // }
        // else{
        //     theta_mi[k] = min(1.,(umi_old[k] - ummi_old[k])/(ui_old[k] - umi_old[k]));
        // }
                    
        // if ((upi_old[k] - ui_old[k]) < 1e-16) //theta_i+1
        // {
        //     theta_i[k] = 0.;
        // }
        // else{
        //     theta_i[k] = min(1.,(ui_old[k] - umi_old[k])/(upi_old[k] - ui_old[k]));
        // }
        // if ((uppi_old[k] - upi_old[k]) < 1e-16) //theta_i-1
        // {
        //     theta_pi[k] = 0.;
        // }
        // else{
        //     theta_pi[k] = min(1.,(upi_old[k] - ui_old[k])/(uppi_old[k] - upi_old[k]));
        // }
    }

     // Minimiser les oscillations (les phi sont pris en compte ici, on pourrait les mettre plus bas)
    for (int k = 0; k < ui_old.size(); k++)
    {
        pi[k] = minmod(sigma_centre_i[k], phi(theta_i[k]) * sigma_droite_i[k], phi(theta_i[k]) * sigma_gauche_i[k]);
        ppi[k] = minmod(sigma_centre_pi[k], phi(theta_pi[k]) * sigma_droite_pi[k], phi(theta_pi[k]) * sigma_gauche_pi[k]);
        pmi[k] = minmod(sigma_centre_mi[k], sigma_droite_mi[k] * phi(theta_mi[k]), sigma_gauche_mi[k] * phi(theta_mi[k]));
        // pi[k] = minmod(sigma_centre_i[k], sigma_droite_i[k], sigma_gauche_i[k]);
        // ppi[k] = minmod(sigma_centre_pi[k], sigma_droite_pi[k], sigma_gauche_pi[k]);
        // pmi[k] = minmod(sigma_centre_mi[k], sigma_droite_mi[k], sigma_gauche_mi[k] * phi(theta_mi[k]));
    }

    // Calcul ui+pi+1
    for (int k = 0; k < ui_old.size(); k++)
    {
        // ui_ppi[k] = ui_old[k] + phi(theta_i[k])*(dx/2)*pi[k];
        // ui_mpi[k] = ui_old[k] - phi(theta_i[k])*(dx/2)*pi[k];
        // upi_mppi[k] = upi_old[k] - phi(theta_pi[k])*(dx/2)*ppi[k];
        // umi_ppmi[k] = umi_old[k] + phi(theta_mi[k])*(dx/2)*pmi[k];
        ui_ppi[k] = ui_old[k] + (dx/2)*pi[k];
        ui_mpi[k] = ui_old[k] - (dx/2)*pi[k];
        upi_mppi[k] = upi_old[k] - (dx/2)*ppi[k];
        umi_ppmi[k] = umi_old[k] + (dx/2)*pmi[k];
    }
    

    if ((choix_schema == 1) && (choix_ordre == 2))
    {
        // for (int k= 0; k < 3; k++)
        // {
        //     flux_pdemi[k] = flux_rusanov(ui_old, upi_old, gamma)[k] + phi(theta_i[k]) * (flux_rusanov(ui_ppi, upi_mppi, gamma)[k] - flux_rusanov(ui_old, upi_old, gamma)[k]);
        //     flux_mdemi[k] = flux_rusanov(umi_old, ui_old, gamma)[k] + phi(theta_mi[k]) * (flux_rusanov(umi_ppmi, ui_mpi, gamma)[k] - flux_rusanov(umi_old, ui_old, gamma)[k]);
        // }
        flux_pdemi = flux_rusanov(ui_ppi, upi_mppi, gamma);
        flux_mdemi = flux_rusanov(umi_ppmi, ui_mpi, gamma);
    }
    else if ((choix_schema == 2) && (choix_ordre == 2))
    {
        flux_pdemi = flux_HLL(ui_ppi, upi_mppi, gamma);
        flux_mdemi = flux_HLL(umi_ppmi, ui_mpi, gamma);
    }
    else if ((choix_schema == 3) && (choix_ordre == 2))
    {
        flux_pdemi = flux_HLLC(ui_ppi, upi_mppi, gamma);
        flux_mdemi = flux_HLLC(umi_ppmi, ui_mpi, gamma);
    }
    else if ((choix_schema == 1) && (choix_ordre == 1))
    {
        flux_pdemi = flux_rusanov(ui_old, upi_old, gamma);
        flux_mdemi = flux_rusanov(umi_old, ui_old, gamma);
    }
    else if ((choix_schema == 2) && (choix_ordre == 1))
    {
        flux_pdemi = flux_HLL(ui_old, upi_old, gamma);
        flux_mdemi = flux_HLL(umi_old, ui_old, gamma);
    }
    else if ((choix_schema == 3) && (choix_ordre == 1))
    {
        flux_pdemi = flux_HLLC(ui_old, upi_old, gamma);
        flux_mdemi = flux_HLLC(umi_old, ui_old, gamma);
    }
    else
    {
        cout << "Pb dans les choix dans ft, on prend Rusanov ordre 1" << endl;
        flux_pdemi = flux_rusanov(ui_old, upi_old, gamma);
        flux_mdemi = flux_rusanov(umi_old, ui_old, gamma);
    }

    for (int k = 0; k < ui_old.size(); k++)
    {
        res[k] = - (1./dx) * (flux_pdemi[k] - flux_mdemi[k]);
    }

    return res;
}

// Calcul erreur
void error(vector<vector<double>> U, string result, int n, double dx, double gamma, double b, int fichier, string test){

    string ignore;
    double delta_x(0), density(0), velocity(0), pressure(0), energy(0), e_int(0);
    double error_rho(0), error_u(0), error_p(0), error_e(0);
    int k;

    std::ifstream monFlux("Resultats/Test_" + test + "/tf_Rusanov_test_2_10000.txt", std::ios::in);

    for (int i = 0; i < 100; i++)
    {
        k=i*fichier;

        if (monFlux){
            
            monFlux >> delta_x;
            monFlux >> density;
            monFlux >> velocity;
            monFlux >> pressure;
            monFlux >> energy;
        }

        // cout << delta_x << " " << density << " " << velocity << " " << pressure << " " << energy << endl;

        error_rho += pow(U[k][0] - density, 2);
        error_u += pow(U[k][1]/U[k][0] - velocity, 2);
        error_p += pow(pression(U[k], gamma) - pressure, 2);
        e_int = (pression(U[k], gamma)*(1-b*U[k][0]))/((gamma-1)*U[k][0]);
        error_e += pow(e_int - energy, 2);

    }

    monFlux.close();
    
    ofstream mon_flux1;
    mon_flux1.open("Resultats/Test_" + test + "/erreur_" + result, ios::app);
    mon_flux1 << dx << " " << sqrt(error_rho*dx) << " " << sqrt(error_u*dx) << " " << sqrt(error_p*dx) << " " << sqrt(error_e*dx) << endl;
    mon_flux1.close();
    
}

#define OPERATION_CPP
#endif