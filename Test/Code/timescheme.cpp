#ifndef TIMESCHEME_CPP

#include "timescheme.h"
#include "flux.h"
#include "operation.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string> 

using namespace std;

// Euler 
void Euler(vector<double> &uij, vector<double> u_mimj, vector<double> u_imj, vector<double> u_pimj, vector<double> u_mij, vector<double> u_ij, vector<double> u_pij, vector<double> u_mipj, vector<double> u_ipj, vector<double> u_pipj, vector<double> u_ppij, vector<double> u_mmij, vector<double> u_ippj, vector<double> u_immj, double dt, double dx, double dy, double gamma, int choix_schema, int choix_ordre){
    
    vector<double> F_time(u_ij.size(),0);
    
    F_time = ft(dx, dy, u_mimj, u_imj, u_pimj, u_mij, u_ij, u_pij, u_mipj, u_ipj, u_pipj, u_ppij, u_mmij, u_ippj, u_immj, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < u_ij.size(); k++)
    {
        uij[k] = u_ij[k] + dt*F_time[k];
    }
}

void RK2_SSP(vector<double> &uij, vector<double> u_mimj, vector<double> u_imj, vector<double> u_pimj, vector<double> u_mij, vector<double> u_ij, vector<double> u_pij, vector<double> u_mipj, vector<double> u_ipj, vector<double> u_pipj, vector<double> u_ppij, vector<double> u_mmij, vector<double> u_ippj, vector<double> u_immj, double dt, double dx, double dy, double gamma, int choix_schema, int choix_ordre){

    vector<double> k1(u_ij.size(), 0), k2(u_ij.size(), 0);
    vector<double> u_star_mimj(u_ij.size(), 0), u_star_imj(u_ij.size(), 0), u_star_pimj(u_ij.size(), 0), u_star_mij(u_ij.size(), 0), u_star_ij(u_ij.size(), 0), u_star_pij(u_ij.size(), 0), u_star_mipj(u_ij.size(), 0), u_star_ipj(u_ij.size(), 0), u_star_pipj(u_ij.size(), 0), u_star_ppij(u_ij.size(), 0), u_star_mmij(u_ij.size(), 0), u_star_ippj(u_ij.size(), 0), u_star_immj(u_ij.size(), 0);
    double u_int(0.);

    k1 = ft(dx, dy, u_mimj, u_imj, u_pimj, u_mij, u_ij, u_pij, u_mipj, u_ipj, u_pipj, u_ppij, u_mmij, u_ippj, u_immj, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < u_ij.size(); k++)
    {
        u_star_mimj[k] = u_mimj[k] + dt*k1[k];
        u_star_imj[k] = u_imj[k] + dt*k1[k];
        u_star_pimj[k] = u_pimj[k] + dt*k1[k];
        u_star_mij[k] = u_mij[k] + dt*k1[k];
        u_star_ij[k] = u_ij[k] + dt*k1[k];
        u_star_pij[k] = u_pij[k] + dt*k1[k];
        u_star_mipj[k] = u_mipj[k] + dt*k1[k];
        u_star_ipj[k] = u_ipj[k] + dt*k1[k];
        u_star_pipj[k] = u_pipj[k] + dt*k1[k];
        u_star_ppij[k] = u_star_ppij[k] + dt*k1[k];
        u_star_mmij[k] = u_star_mmij[k] + dt*k1[k];
        u_star_ippj[k] = u_star_ippj[k] + dt*k1[k];
        u_star_immj[k] = u_star_immj[k] + dt*k1[k];

    }

    k2  = ft(dx, dy, u_star_mimj, u_star_imj, u_star_pimj, u_star_mij, u_star_ij, u_star_pij, u_star_mipj, u_star_ipj, u_star_pipj, u_star_ppij, u_star_mmij, u_star_ippj, u_star_immj, gamma, choix_schema, choix_schema);
    
    for (int k = 0; k < u_ij.size(); k++)
    {
        u_int = 0.5*u_ij[k] + 0.5*(u_star_ij[k] + dt*k2[k]);
        if ((u_int > 1e-16) || (u_int < 1e16))
        {
            uij[k] = 0.5*u_ij[k] + 0.5*(u_star_ij[k] + dt*k2[k]);
        }
        else{
            u_ij[k] = 0.;
            
        }
        
        // uij[k] = 0.5*dt*k1[k] + u_ij[k] + 0.5*dt*k2[k];
    }
}

void RK4(vector<double> &uij, vector<double> u_mimj, vector<double> u_imj, vector<double> u_pimj, vector<double> u_mij, vector<double> u_ij, vector<double> u_pij, vector<double> u_mipj, vector<double> u_ipj, vector<double> u_pipj, vector<double> u_ppij, vector<double> u_mmij, vector<double> u_ippj, vector<double> u_immj, double dt, double dx, double dy, double gamma, int choix_schema, int choix_ordre){
    // ui = ui_old - (dt/dx) * (flux_pdemi - flux_mdemi);

    vector<double> k1(u_ij.size(), 0),k2(u_ij.size(), 0),k3(u_ij.size(), 0),k4(u_ij.size(), 0),ktot(u_ij.size(), 0);

    k1 = ft(dx, dy, u_mimj, u_imj, u_pimj, u_mij, u_ij, u_pij, u_mipj, u_ipj, u_pipj, u_ppij, u_mmij, u_ippj, u_immj, gamma, choix_schema, choix_ordre);

    //////k2
    vector<double> u_mimj_k2(u_ij.size(), 0), u_imj_k2(u_ij.size(), 0), u_pimj_k2(u_ij.size(), 0), u_mij_k2(u_ij.size(), 0), u_ij_k2(u_ij.size(), 0), u_pij_k2(u_ij.size(), 0), u_mipj_k2(u_ij.size(), 0), u_ipj_k2(u_ij.size(), 0), u_pipj_k2(u_ij.size(), 0), u_ppij_k2(u_ij.size(), 0), u_mmij_k2(u_ij.size(), 0), u_ippj_k2(u_ij.size(), 0), u_immj_k2(u_ij.size(), 0);
    for (int k = 0; k < u_ij.size(); k++)
    {
        u_mimj_k2[k] = u_mimj[k] + dt*k1[k]/2;
        u_imj_k2[k] = u_imj[k] + dt*k1[k]/2;
        u_pimj_k2[k] = u_pimj[k] + dt*k1[k]/2;
        u_mij_k2[k] = u_mij[k] + dt*k1[k]/2;
        u_ij_k2[k] = u_ij[k] + dt*k1[k]/2;
        u_pij_k2[k] = u_pij[k] + dt*k1[k]/2;
        u_mipj_k2[k] = u_mipj[k] + dt*k1[k]/2;
        u_ipj_k2[k] = u_ipj[k] + dt*k1[k]/2;
        u_pipj_k2[k] = u_pipj[k] + dt*k1[k]/2;
        u_ppij_k2[k] = u_ppij[k] + dt*k1[k]/2;
        u_mmij_k2[k] = u_mmij[k] + dt*k1[k]/2;
        u_ippj_k2[k] = u_ippj[k] + dt*k1[k]/2;
        u_immj_k2[k] = u_immj[k] + dt*k1[k]/2;
    }
    k2  = ft(dx, dy, u_mimj_k2, u_imj_k2, u_pimj_k2, u_mij_k2, u_ij_k2, u_pij_k2, u_mipj_k2, u_ipj_k2, u_pipj_k2, u_ppij_k2, u_mmij_k2, u_ippj_k2, u_immj_k2, gamma, choix_schema, choix_ordre);

    //////k3
    vector<double> u_mimj_k3(u_ij.size(), 0), u_imj_k3(u_ij.size(), 0), u_pimj_k3(u_ij.size(), 0), u_mij_k3(u_ij.size(), 0), u_ij_k3(u_ij.size(), 0), u_pij_k3(u_ij.size(), 0), u_mipj_k3(u_ij.size(), 0), u_ipj_k3(u_ij.size(), 0), u_pipj_k3(u_ij.size(), 0), u_ppij_k3(u_ij.size(), 0), u_mmij_k3(u_ij.size(), 0), u_ippj_k3(u_ij.size(), 0), u_immj_k3(u_ij.size(), 0);
    for (int k = 0; k < u_ij.size(); k++)
    {
        u_mimj_k3[k] = u_mimj[k] + dt*k2[k]/2;
        u_imj_k3[k] = u_imj[k] + dt*k2[k]/2;
        u_pimj_k3[k] = u_pimj[k] + dt*k2[k]/2;
        u_mij_k3[k] = u_mij[k] + dt*k2[k]/2;
        u_ij_k3[k] = u_ij[k] + dt*k2[k]/2;
        u_pij_k3[k] = u_pij[k] + dt*k2[k]/2;
        u_mipj_k3[k] = u_mipj[k] + dt*k2[k]/2;
        u_ipj_k3[k] = u_ipj[k] + dt*k2[k]/2;
        u_pipj_k3[k] = u_pipj[k] + dt*k2[k]/2;
        u_ppij_k3[k] = u_ppij[k] + dt*k2[k]/2;
        u_mmij_k3[k] = u_mmij[k] + dt*k2[k]/2;
        u_ippj_k3[k] = u_ippj[k] + dt*k2[k]/2;
        u_immj_k3[k] = u_immj[k] + dt*k2[k]/2;
    }
    k3  = ft(dx, dy, u_mimj_k3, u_imj_k3, u_pimj_k3, u_mij_k3, u_ij_k3, u_pij_k3, u_mipj_k3, u_ipj_k3, u_pipj_k3, u_ppij_k3, u_mmij_k3, u_ippj_k3, u_immj_k3, gamma, choix_schema, choix_ordre);

    //////k4
    vector<double> u_mimj_k4(u_ij.size(), 0), u_imj_k4(u_ij.size(), 0), u_pimj_k4(u_ij.size(), 0), u_mij_k4(u_ij.size(), 0), u_ij_k4(u_ij.size(), 0), u_pij_k4(u_ij.size(), 0), u_mipj_k4(u_ij.size(), 0), u_ipj_k4(u_ij.size(), 0), u_pipj_k4(u_ij.size(), 0), u_ppij_k4(u_ij.size(), 0), u_mmij_k4(u_ij.size(), 0), u_ippj_k4(u_ij.size(), 0), u_immj_k4(u_ij.size(), 0);
    for (int k = 0; k < u_ij.size(); k++)
    {
        u_mimj_k4[k] = u_mimj[k] + dt*k3[k];
        u_imj_k4[k] = u_imj[k] + dt*k3[k];
        u_pimj_k4[k] = u_pimj[k] + dt*k3[k];
        u_mij_k4[k] = u_mij[k] + dt*k3[k];
        u_ij_k4[k] = u_ij[k] + dt*k3[k];
        u_pij_k4[k] = u_pij[k] + dt*k3[k];
        u_mipj_k4[k] = u_mipj[k] + dt*k3[k];
        u_ipj_k4[k] = u_ipj[k] + dt*k3[k];
        u_pipj_k4[k] = u_pipj[k] + dt*k3[k];
        u_ppij_k4[k] = u_ppij[k] + dt*k3[k];
        u_mmij_k4[k] = u_mmij[k] + dt*k3[k];
        u_ippj_k4[k] = u_ippj[k] + dt*k3[k];
        u_immj_k4[k] = u_immj[k] + dt*k3[k];
    }
    k4  = ft(dx, dy, u_mimj_k4, u_imj_k4, u_pimj_k4, u_mij_k4, u_ij_k4, u_pij_k4, u_mipj_k4, u_ipj_k4, u_pipj_k4, u_ppij_k4, u_mmij_k4, u_ippj_k4, u_immj_k4, gamma, choix_schema, choix_ordre);

    // Ktot
    for (int k = 0; k < u_ij.size(); k++)
    {
        ktot[k]= k1[k] + 2*k2[k] + 2*k3[k] + k4[k];
        uij[k] = u_ij[k] + (dt/6)*ktot[k];
        // cout << "ktot " << ktot[k] << endl;
    }
}


#define TIMESCHEME_CPP
#endif