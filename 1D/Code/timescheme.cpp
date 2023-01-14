#ifndef TIMESCHEME_CPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string> 
#include "timescheme.h"
#include "flux.h"
#include "operation.h"

using namespace std;

void Euler(vector<double> &ui, vector<double> ummi_old, vector<double> umi_old, vector<double> ui_old, vector<double> upi_old, vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre){
    vector<double> F_time(3,0);
    
    F_time = ft(dx, ummi_old, umi_old, ui_old, upi_old, uppi_old, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < ui_old.size(); k++)
    {
        ui[k] = ui_old[k] + dt*F_time[k];
        // Mood
        // if (((ui[k] > umi_old[k]) && (ui[k] > upi_old[k])) || ((ui[k] < umi_old[k]) && (ui[k] < upi_old[k])))
        // {
        //     ui[k] = ui_old[k] + dt*ft(dx, ummi_old, umi_old, ui_old, upi_old, uppi_old, gamma, choix_schema, 1)[k];
        // }
    }

    // if ((ui[0] < 0) || (pression(ui, gamma) < 0) || ((pression(ui, gamma))/((gamma-1)*ui[0]) < 0) )
    // {
    //     cout << "Mood " << endl;
    //     for (int k = 0; k < ui_old.size(); k++)
    //     {
    //         ui[k] = ui_old[k] + dt*ft(dx, ummi_old, umi_old, ui_old, upi_old, uppi_old, gamma, choix_schema, 1)[k];
    //     }
    // }
    

}

void RK2_SSP(vector<double> &ui, vector<double> ummi_old, vector<double> umi_old, vector<double> ui_old, vector<double> upi_old, vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre){

    vector<double> k1(ui_old.size(), 0), k2(ui_old.size(), 0), u_star_i(ui_old.size(), 0), u_star_mi(ui_old.size(), 0), u_star_pi(ui_old.size(), 0), u_star_mmi(ui_old.size(), 0), u_star_ppi(ui_old.size(), 0);

    k1 = ft(dx, ummi_old, umi_old, ui_old, upi_old, uppi_old, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < ui_old.size(); k++)
    {
        u_star_i[k] = ui_old[k] + dt*k1[k];
        u_star_mi[k] = umi_old[k] + dt*k1[k];
        u_star_pi[k] = upi_old[k] + dt*k1[k];
        u_star_mmi[k] = ummi_old[k] + dt*k1[k];
        u_star_ppi[k] = uppi_old[k] + dt*k1[k];
        // cout << "ktot " << ktot[k] << endl;
    }

    k2  = ft(dx, u_star_mmi, u_star_mi, u_star_i, u_star_pi, u_star_ppi, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < ui_old.size(); k++)
    {
        ui[k] = 0.5*ui_old[k] + 0.5*(u_star_i[k] + dt*k2[k]);
    }
}

void RK2_SSP_demi(vector<double> &ui, vector<double> ummi_old, vector<double> umi_old, vector<double> ui_old, vector<double> upi_old, vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre){
    
    vector<double> k1(ui_old.size(), 0), k2(ui_old.size(), 0), u_star_i(ui_old.size(), 0), u_star_mi(ui_old.size(), 0), u_star_pi(ui_old.size(), 0), u_star_mmi(ui_old.size(), 0), u_star_ppi(ui_old.size(), 0);

    k1 = ft(dx, ummi_old, umi_old, ui_old, upi_old, uppi_old, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < ui_old.size(); k++)
    {
        u_star_i[k] = ui_old[k] + (2./3.)*dt*k1[k];
        u_star_mi[k] = umi_old[k] + (2./3.)*dt*k1[k];
        u_star_pi[k] = upi_old[k] + (2./3.)*dt*k1[k];
        u_star_mmi[k] = ummi_old[k] + (2./3.)*dt*k1[k];
        u_star_ppi[k] = uppi_old[k] + (2./3.)*dt*k1[k];
        // cout << "ktot " << ktot[k] << endl;
    }

    k2  = ft(dx, u_star_mmi, u_star_mi, u_star_i, u_star_pi, u_star_ppi, gamma, choix_schema, choix_ordre);

    for (int k = 0; k < ui_old.size(); k++)
    {
        ui[k] = (5./8.)*ui_old[k] + (3./8.)*u_star_i[k] + 0.75*(dt*k2[k]);
    }
}

void RK4(vector<double> &ui, vector<double> ummi_old, vector<double> umi_old, vector<double> ui_old, vector<double> upi_old, vector<double> uppi_old, double dt, double dx, double gamma, int choix_schema, int choix_ordre){
    // ui = ui_old - (dt/dx) * (flux_pdemi - flux_mdemi);

    vector<double> k1(ui_old.size(), 0),k2(ui_old.size(), 0),k3(ui_old.size(), 0),k4(ui_old.size(), 0),ktot(ui_old.size(), 0);

    k1 = ft(dx, ummi_old, umi_old, ui_old, upi_old, uppi_old, gamma, choix_schema, choix_ordre);

    //////k2
    vector<double> u_i_k2(ui_old.size(), 0), u_pi_k2(ui_old.size(), 0), u_mi_k2(ui_old.size(), 0), u_ppi_k2(ui_old.size(), 0), u_mmi_k2(ui_old.size(), 0);
    for (int k = 0; k < ui_old.size(); k++)
    {
        u_i_k2[k] = ui_old[k] +dt*k1[k]/2.;
        u_pi_k2[k] = upi_old[k] +dt*k1[k]/2.;
        u_mi_k2[k] = umi_old[k] +dt*k1[k]/2.;
        u_ppi_k2[k] = uppi_old[k] +dt*k1[k]/2.;
        u_mmi_k2[k] = ummi_old[k] +dt*k1[k]/2.;
    }
    k2 = ft(dx, u_mmi_k2, u_mi_k2, u_i_k2, u_pi_k2, u_ppi_k2, gamma, choix_schema, choix_ordre);

    //////k3
    vector<double> u_i_k3(ui_old.size(), 0), u_pi_k3(ui_old.size(), 0), u_mi_k3(ui_old.size(), 0), u_ppi_k3(ui_old.size(), 0), u_mmi_k3(ui_old.size(), 0);
    for (int k = 0; k < ui_old.size(); k++)
    {
        u_i_k3[k] = ui_old[k] +dt*k2[k]/2.;
        u_pi_k3[k] = upi_old[k] +dt*k2[k]/2.;
        u_mi_k3[k] = umi_old[k] +dt*k2[k]/2.;
        u_mmi_k3[k] = ummi_old[k] +dt*k2[k]/2.;
        u_ppi_k3[k] = uppi_old[k] +dt*k2[k]/2.;
    }
    k3 = ft(dx, u_mmi_k3, u_mi_k3, u_i_k3, u_pi_k3, u_ppi_k3, gamma, choix_schema, choix_ordre);

    //////k4
    vector<double> u_i_k4(ui_old.size(), 0), u_pi_k4(ui_old.size(), 0), u_mi_k4(ui_old.size(), 0), u_ppi_k4(ui_old.size(), 0), u_mmi_k4(ui_old.size(), 0);
    for (int k = 0; k < ui_old.size(); k++)
    {
        u_i_k4[k] = ui_old[k] +dt*k3[k];
        u_pi_k4[k] = upi_old[k] +dt*k3[k];
        u_mi_k4[k] = umi_old[k] +dt*k3[k];
        u_mmi_k4[k] = ummi_old[k] +dt*k3[k];
        u_ppi_k4[k] = uppi_old[k] +dt*k3[k];
    }
    k4 = ft(dx, u_mmi_k4, u_mi_k4, u_i_k4, u_pi_k4, u_ppi_k4, gamma, choix_schema, choix_ordre);

    // Ktot
    for (int k = 0; k < ui_old.size(); k++)
    {
        ktot[k]= k1[k] + 2*k2[k] + 2*k3[k] + k4[k];
        ui[k] = ui_old[k] + (dt/6)*ktot[k];
        // cout << "ktot " << ktot[k] << endl;
    }
}


#define TIMESCHEME_CPP
#endif