#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "timescheme.h"
#include "flux.h"
#include "operation.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////
int main(){

    // Initialisation variables physiques
    double rho_l(0), rho_r(0), u_l(0), u_r(0), p_l(0), p_r(0);
    // Initialisation domaine
    double Nx, x0, xi, delta_x, CFL, tf(0), disc(0);
    int n(0);
    string test, result;
    int choix_schema, choix_ordre;

    lecture(n, tf, test, disc, rho_l, u_l, p_l, rho_r, u_r, p_r, choix_schema, choix_ordre);

    // Domaine et pas d'espace
    CFL = 0.8;
    Nx = 1.;
    delta_x = Nx/n;
    x0 =0.;

    // Initialisation variables pour E
    double b, gamma, E_l, E_r, eint_l, eint_r;
    // double E_m, eint_m;

    b = 0.;
    gamma = 1.4;

    eint_l = (p_l*(1-b*rho_l))/((gamma-1)*rho_l);
    eint_r = (p_r*(1-b*rho_r))/((gamma-1)*rho_r);
    // eint_m = (0.01*(1-b*rho_r))/((gamma-1)*rho_r);

    E_l = rho_l*pow(u_l,2)/2. + rho_l*eint_l;
    E_r = rho_r*pow(u_r,2)/2. + rho_r*eint_r;
    // E_m = rho_r*pow(u_r,2)/2. + rho_r*eint_m;


    // Initialisation des variables utiles
    int kmax = 100;
    double delta_t = 0.00001, t = 0.;
    vector<double> flux_pdemi(3,0), flux_mdemi(3,0);
    vector<vector<double>> U_old(n,vector<double> (3)), U;

    // Initialisation de U
    U0(U, rho_l, u_l, E_l, rho_r, u_r, E_r, delta_x, x0, disc, n);

    // choix_schema de la méthode de résolution
    // int choix_schema, choix_ordre;
    // cout << "Choisir la méthode de résolution : 1 Rusanov, 2 HLL et 3 HLLC" << endl;
    // cin >> choix_schema;

    // cout << "Choisir l'ordre de la méthode de résolution : 1 ou 2" << endl;
    // cin >> choix_ordre;

    if ((choix_schema == 1) && (choix_ordre == 1))
    {
        result = "Rusanov_ordre_1.dat";
    }
    else if ((choix_schema == 2) && (choix_ordre == 1))
    {
        result = "HLL_ordre_1.dat";
    }
    else if ((choix_schema == 3) && (choix_ordre == 1))
    {
        result = "HLLC_ordre_1.dat";
    }
    else if ((choix_schema == 1) && (choix_ordre == 2))
    {
        result = "Rusanov_ordre_2.dat";
    }
    else if ((choix_schema == 2) && (choix_ordre == 2))
    {
        result = "HLL_ordre_2.dat";
    }
    else if ((choix_schema == 3) && (choix_ordre == 2))
    {
        result = "HLLC_ordre_2.dat";
    }
    else
    {
        cout << "Ce n'est pas le bon choix_schema, par défaut on prend Rusanov" << endl;
        result = "Rusanov_" + test + ".dat";
    }


    ///////// TEST S/////////////////
    // double S_l, S_m, S_r;
    // vect_S(U[90], U[90], S_l, S_m, S_r, gamma);
    // cout << S_l << " " << S_m << " " << S_r << endl;

    // Sauvegarde de la solution initiale
    // vector<double> p_i(n,0), e_i(n,0);
    // ofstream mon_flux;
    // mon_flux.open("Resultats/t0_" + result, ios::out);
    // for (int i = 0; i < n; i++)
    // {
    //     xi = x0 + i*delta_x;
    //     p_i[i] = pression(U[i], gamma);
    //     e_i[i] = (p_i[i]*(1-b*U[i][0]))/((gamma-1)*U[i][0]);
    //     mon_flux << xi << " " << U[i][0] << " " << U[i][1]/U[i][0] << " " << p_i[i] << " " << e_i[i] << endl;
    // }

    // // Sauvegarde de la solution finale
    int k(0), fichier(1);
    // fichier= n/100;
    


    // mon_flux.close();
    int l(1);
    // Itération pour avoir les Ui
    while (t <= tf)
    {
        vector<double> p_f(n/fichier,0), e_f(n/fichier,0);
        ofstream mon_flux2; // mon_flux2.open("Resultats/tf_" + result, ios::out);
        mon_flux2.open("Resultats/Test_" + test + "/tf_" + to_string(l) + "_" + result, ios::out);
        double lambda_max = -1;
        // Boucle pour garder les valeurs de Un
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                U_old[i][j] = U[i][j];
            }
            if (U[i][1]/U[i][0]<0)
            {
                lambda_max = max(lambda_max, abs(U[i][1]/U[i][0]) + sqrt(gamma*pression(U[i], gamma)/U[i][0]));
            }
            else{
                lambda_max = max(lambda_max, U[i][1]/U[i][0] + sqrt(gamma*pression(U[i], gamma)/U[i][0]));
            }
        }
        // cout  << lambda_max << endl;

        for (int i = 2; i < n-2; i++)
        {
            Euler(U[i], U_old[i-2], U_old[i-1], U_old[i], U_old[i+1], U_old[i+2], delta_t, delta_x, gamma, choix_schema, choix_ordre);
            // RK4(U[i], U_old[i-2], U_old[i-1], U_old[i], U_old[i+1], U_old[i+2], delta_t, delta_x, gamma, choix_schema, choix_ordre);
            // RK2_SSP(U[i], U_old[i-2], U_old[i-1], U_old[i], U_old[i+1], U_old[i+2], delta_t, delta_x, gamma, choix_schema, choix_ordre);
            k=i*fichier;
            xi = x0 + k*delta_x;
            p_f[i] = pression(U[k], gamma);
            e_f[i] = (p_f[i]*(1-b*U[k][0]))/((gamma-1)*U[k][0]);
            mon_flux2 << xi << " " << U[k][0] << " " << U[k][1]/U[k][0] << " " << p_f[i] << " " << e_f[i] << endl;

        }

        
        delta_t = delta_x*CFL/(2*lambda_max);
        t = t + delta_t;
        l++;
        cout <<"Iteration : " << l-1 << ", Time = " << t << endl;


        mon_flux2.close();
    }

    // // // Sauvegarde de la solution finale
    // int k, fichier(1);
    fichier= n/100;
    vector<double> p_f(n/fichier,0), e_f(n/fichier,0);
    ofstream mon_flux2;
    
    // mon_flux2.open("Resultats/tf_" + result, ios::out);
    mon_flux2.open("Resultats/Test_" + test + "/tf_" + result, ios::out);
    for (int i = 0; i < n/fichier; i++)
    {
        k=i*fichier;
        xi = x0 + k*delta_x;
        p_f[i] = pression(U[k], gamma);
        e_f[i] = (p_f[i]*(1-b*U[k][0]))/((gamma-1)*U[k][0]);
        mon_flux2 << xi << " " << U[k][0] << " " << U[k][1]/U[k][0] << " " << p_f[i] << " " << e_f[i] << endl;
    }

    mon_flux2.close();

    error(U, result, n, delta_x, gamma, b, fichier, test);


    return 0;
}
