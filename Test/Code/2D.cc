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
    vector<double> rho_1(4,0), u0_1(4,0), v0_1(4,0), p_1(4,0);
    // Initialisation domaine
    double Nx, Ny, x0_1, y0_1, xi, yi, delta_x, delta_y, CFL, tf(0), disc_x(0), disc_y(0);
    double tf_simu(0);
    int n(0);
    string test, result; 
    int choix_schema, choix_ordre;

    lecture(n, tf, test, disc_x, disc_y, rho_1, u0_1, v0_1, p_1, choix_schema, choix_ordre);

    // Domaine et pas d'espace
    CFL = 0.475;
    Nx = 10.;
    Ny = 10.;
    delta_x = Nx/n; 
    delta_y = Ny/n; 
    x0_1 =0.;
    y0_1 =0.;
    int nn = n*n;

    // Initialisation variables pour E
    double b, gamma;
    vector<double> E(4,0), eint(4,0);
    
    b = 0.;
    gamma = 1.4;

    double pi;
    pi= 4.0*atan(1.0);
    
    // for (int k = 0; k < E.size(); k++)
    // {
    //     eint[k] = (p[k]*(1-b*rho[k]))/((gamma-1)*rho[k]);
    //     E[k] = rho[k]*(pow(u0[k],2)+pow(v0[k],2))/2. + rho[k]*eint[k];
    // }

       // Initialisation des variables utiles
    int kmax = 100;
    double delta_t = 0.00001, t = 0.;
    vector<double> flux_pdemi_x(4,0), flux_mdemi_x(4,0), flux_pdemi_y(4,0), flux_mdemi_y(4,0);
    vector<vector<double>> U_rho_old(n,vector<double> (n)), U_rhou_old(n,vector<double> (n)), U_rhov_old(n,vector<double> (n)), U_E_old(n,vector<double> (n));
    vector<vector<double>> U_rho(n,vector<double> (n)), U_rhou(n,vector<double> (n)), U_rhov(n,vector<double> (n)), U_E(n,vector<double> (n));


  double u_inf = 0.5;
  double v_inf = 0.0;
  double beta = u_inf;
  double x0, y0;
  /* Initial solution */
  x0 = 5.0, y0 = 5.0;
    for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
        xi = i*delta_x;
        yi = j*delta_y;
        double rx, ry;
        rx = (xi - x0);
        ry = (yi - y0);
        if (rx < -5)      { rx += 10; }
        else if (rx > 5)  { rx -= 10; }
        double rsq = rx*rx + ry*ry;
        double rho, u, v, P;
        double du, dv;
        rho = pow(1.0 - ((gamma-1.0)*beta*beta)/(8.0*gamma*pi*pi) * exp(1.0-rsq), 1.0/(gamma-1.0));
        P   = pow(rho,gamma);
        du  = - beta/(2.0*pi) * exp(0.5*(1.0-rsq)) * ry;
        dv  =   beta/(2.0*pi) * exp(0.5*(1.0-rsq)) * rx;
        u   = u_inf + du;
        v   = v_inf + dv;
        U_rho[i][j] = rho;
        U_rhou[i][j] = rho*u;
        U_rhov[i][j] = rho*v;
        U_E[i][j] = P/(gamma-1.0) + 0.5*rho*(u*u+v*v);
      }
    }

    
 
    // Initialisation de U
    //U0(U_rho, U_rhou, U_rhov, U_E, rho, u0, v0, E, delta_x, x0, delta_y, y0, disc_x, disc_y, n);    

    if ((choix_schema == 1) && (choix_ordre == 1))
    {
        result = "Rusanov_ordre_1_" + test + ".txt";
    }
    else if ((choix_schema == 2) && (choix_ordre == 1))
    {
        result = "HLL_ordre_1_" + test + ".txt";
    }
    else if ((choix_schema == 3) && (choix_ordre == 1))
    {
        result = "HLLC_ordre_1_" + test + ".txt";
    }
    else if ((choix_schema == 1) && (choix_ordre == 2))
    {
        result = "Rusanov_ordre_2_" + test + ".txt";
    }
    else if ((choix_schema == 2) && (choix_ordre == 2))
    {
        result = "HLL_ordre_2_" + test + ".txt";
    }
    else if ((choix_schema == 3) && (choix_ordre == 2))
    {
        result = "HLLC_ordre_2_" + test + ".txt";
    }
    else
    {
        cout << "Ce n'est pas le bon choix, par défaut on prend Rusanov" << endl;
        result = "Rusanov_" + test + ".txt";
    }
    

    ///////// TEST S/////////////////
    // double S_l, S_m, S_r;
    // vect_S(U[90], U[90], S_l, S_m, S_r, gamma);
    // cout << S_l << " " << S_m << " " << S_r << endl;

    // Sauvegarde de la solution initiale
    vector<vector<double>> p_i(n,vector<double> (n)), e_i(n,vector<double> (n));
    vector<double> U_int(4,0), U_mimj(4,0), U_imj(4,0), U_pimj(4,0), U_mij(4,0), U_ij(4,0), U_pij(4,0), U_mipj(4,0), U_ipj(4,0),U_pipj(4,0), U_ppij(4,0), U_mmij(4,0), U_ippj(4,0), U_immj(4,0);

    // ofstream mon_flux;
    // mon_flux.open("Resultats/t0_" + result, ios::out);
    // for (int i = 0; i < n; i++)
    // {
    //     xi = x0_1 + i*delta_x;
    //     for (int j = 0; j < n; j++)
    //     {
    //         U_int[0] = U_rho[i][j];
    //         U_int[1] = U_rhou[i][j];
    //         U_int[2] = U_rhov[i][j];
    //         U_int[3] = U_E[i][j];

    //         yi = y0_1 + j*delta_y;
    //         p_i[i][j] = pression(U_int, gamma);
    //         e_i[i][j] = (p_i[i][j]*(1-b*U_int[0]))/((gamma-1)*U_int[0]);
    //         mon_flux << xi << " " << yi << " " << U_int[0] << " " << U_int[1]/U_int[0] << " " << p_i[i][j] << " " << e_i[i][j] << endl;
    //     }
        
    // }
    // mon_flux.close();

    // // Itération pour avoir les Ui
    while (t <= tf)
    {
        double lambda_max_x = -1;
        double lambda_max_y = -1;
        // Boucle pour garder les valeurs de Un
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                U_rho_old[i][j] = U_rho[i][j];
                U_rhou_old[i][j] = U_rhou[i][j];
                U_rhov_old[i][j] = U_rhov[i][j];
                U_E_old[i][j] = U_E[i][j];

                U_int[0] = U_rho[i][j];
                U_int[1] = U_rhou[i][j];
                U_int[2] = U_rhov[i][j];
                U_int[3] = U_E[i][j];

                lambda_max_x = max(lambda_max_x, abs(U_rhou[i][j]/U_rho[i][j]) + sqrt(gamma*pression(U_int, gamma)/U_rho[i][j]));
                lambda_max_y = max(lambda_max_y, abs(U_rhov[i][j]/U_rho[i][j]) + sqrt(gamma*pression(U_int, gamma)/U_rho[i][j]));
                
                // if (U_rhou[i][j]/U_rho[i][j]<0)
                // {
                //     lambda_max_x = max(lambda_max_x, abs(U_rhou[i][j]/U_rho[i][j]) + sqrt(gamma*pression(U_int, gamma)/U_rho[i][j]));
                //     lambda_max_y = max(lambda_max_y, abs(U_rho[i][j]/U_rho[i][j]) + sqrt(gamma*pression(U_int, gamma)/U_rho[i][j]));
                // }
                // else{
                //     lambda_max_x = max(lambda_max_x, U_rhou[i][j]/U_rho[i][j] + sqrt(gamma*pression(U_int, gamma)/U_rho[i][j]));
                //     lambda_max_y = max(lambda_max_y, U_rhov[i][j]/U_rho[i][j] + sqrt(gamma*pression(U_int, gamma)/U_rho[i][j]));
                // }
            
            }
            
        }

        // Calcul de U_ij
        for (int i = 2; i < n-2; i++)
        {
            for (int j = 2; j < n-2; j++)
            {
                // U_mimj
                U_mimj[0] = U_rho_old[i-1][j-1];
                U_mimj[1] = U_rhou_old[i-1][j-1];
                U_mimj[2] = U_rhov_old[i-1][j-1];
                U_mimj[3] = U_E_old[i-1][j-1];

                // U_imj
                U_imj[0] = U_rho_old[i][j-1];
                U_imj[1] = U_rhou_old[i][j-1];
                U_imj[2] = U_rhov_old[i][j-1];
                U_imj[3] = U_E_old[i][j-1];

                // U_pimj
                U_pimj[0] = U_rho_old[i+1][j-1];
                U_pimj[1] = U_rhou_old[i+1][j-1];
                U_pimj[2] = U_rhov_old[i+1][j-1];
                U_pimj[3] = U_E_old[i+1][j-1];

                // U_mij
                U_mij[0] = U_rho_old[i-1][j];
                U_mij[1] = U_rhou_old[i-1][j];
                U_mij[2] = U_rhov_old[i-1][j];
                U_mij[3] = U_E_old[i-1][j];

                // U_ij
                U_ij[0] = U_rho_old[i][j];
                U_ij[1] = U_rhou_old[i][j];
                U_ij[2] = U_rhov_old[i][j];
                U_ij[3] = U_E_old[i][j];

                // U_pij
                U_pij[0] = U_rho_old[i+1][j];
                U_pij[1] = U_rhou_old[i+1][j];
                U_pij[2] = U_rhov_old[i+1][j];
                U_pij[3] = U_E_old[i+1][j];

                // U_mipj
                U_mipj[0] = U_rho_old[i-1][j+1];
                U_mipj[1] = U_rhou_old[i-1][j+1];
                U_mipj[2] = U_rhov_old[i-1][j+1];
                U_mipj[3] = U_E_old[i-1][j+1];

                // U_ipj
                U_ipj[0] = U_rho_old[i][j+1];
                U_ipj[1] = U_rhou_old[i][j+1];
                U_ipj[2] = U_rhov_old[i][j+1];
                U_ipj[3] = U_E_old[i][j+1];

                // U_pipj
                U_pipj[0] = U_rho_old[i+1][j+1];
                U_pipj[1] = U_rhou_old[i+1][j+1];
                U_pipj[2] = U_rhov_old[i+1][j+1];
                U_pipj[3] = U_E_old[i+1][j+1];

                // U_ppij
                U_ppij[0] = U_rho_old[i+2][j];
                U_ppij[1] = U_rhou_old[i+2][j];
                U_ppij[2] = U_rhov_old[i+2][j];
                U_ppij[3] = U_E_old[i+2][j];
                // U_mmij
                U_mmij[0] = U_rho_old[i-2][j];
                U_mmij[1] = U_rhou_old[i-2][j];
                U_mmij[2] = U_rhov_old[i-2][j];
                U_mmij[3] = U_E_old[i-2][j];
                // U_ippj
                U_ippj[0] = U_rho_old[i][j+2];
                U_ippj[1] = U_rhou_old[i][j+2];
                U_ippj[2] = U_rhov_old[i][j+2];
                U_ippj[3] = U_E_old[i][j+2];
                // U_ipj
                U_immj[0] = U_rho_old[i][j-2];
                U_immj[1] = U_rhou_old[i][j-2];
                U_immj[2] = U_rhov_old[i][j-2];
                U_immj[3] = U_E_old[i][j-2];

                // Euler(U_ij, U_mimj, U_imj, U_pimj, U_mij, U_ij, U_pij, U_mipj, U_ipj, U_pipj, U_ppij, U_mmij, U_ippj, U_immj, delta_t, delta_x, delta_y, gamma, choix_schema, choix_ordre);
                // RK2_SSP(U_ij, U_mimj, U_imj, U_pimj, U_mij, U_ij, U_pij, U_mipj, U_ipj, U_pipj, U_ppij, U_mmij, U_ippj, U_immj, delta_t, delta_x, delta_y, gamma, choix_schema, choix_ordre);
                RK4(U_ij, U_mimj, U_imj, U_pimj, U_mij, U_ij, U_pij, U_mipj, U_ipj, U_pipj, U_ppij, U_mmij, U_ippj, U_immj, delta_t, delta_x, delta_y, gamma, choix_schema, choix_ordre);

                    
                U_rho[i][j] = U_ij[0];
                U_rhou[i][j] = U_ij[1];
                U_rhov[i][j] = U_ij[2];
                U_E[i][j] = U_ij[3];
            }    
        }
  
        double lambda_max = max(lambda_max_x, lambda_max_y);
        delta_t = CFL*delta_x/(2*lambda_max);
        t = t + delta_t;
        cout << "t = " << t << endl;   
    }
    tf_simu = t;

    // // Sauvegarde de la solution finale
    int kx, ky, fichier(1);
    // fichier= n/100;
    vector<vector<double>> p_f(n,vector<double> (n)), e_f(n,vector<double> (n));
    ofstream mon_flux2;
    mon_flux2.open("Resultats/tf_" + result, ios::out);
    for (int i = 2; i < n/fichier-2; i++)
    // for (int i = 0; i < n; i++)
    {
        kx=i*fichier;
        xi = x0_1 + kx*delta_x;
        for (int j = 2; j < n/fichier-2; j++)
        // for (int j = 0; j < n; j++)
        {
            mon_flux2 << " " << endl;
            ky=j*fichier;
            yi = y0_1 + ky*delta_y;

            U_int[0] = U_rho[kx][ky];
            U_int[1] = U_rhou[kx][ky];
            U_int[2] = U_rhov[kx][ky];
            U_int[3] = U_E[kx][ky];

            p_f[kx][ky] = pression(U_int, gamma);
            e_f[kx][ky] = (p_f[kx][ky]*(1-b*U_int[0]))/((gamma-1)*U_int[0]);
            mon_flux2 << xi << " " << yi << " " << U_int[0] << " " << U_int[1]/U_int[0] << " " << U_int[2]/U_int[0] << " " << p_f[kx][ky] << " " << e_f[kx][ky] << endl;
        }
    }
    mon_flux2.close();


    /* Exact solution */
    vector<vector<double>> U_rho_exact(n,vector<double> (n)), U_rhou_exact(n,vector<double> (n)), U_rhov_exact(n,vector<double> (n)), U_E_exact(n,vector<double> (n));
    x0 = 5.0+tf_simu*u_inf, y0 = 5.0+tf_simu*v_inf;
    if (x0 > 10) x0 -= 10; //periodic domain
   // printf("Final time: %lf, Vortex center: %lf, %lf\n",tff,x0,y0);
    for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
        xi = i*delta_x;
        yi = j*delta_y;
        double rx, ry;
        rx = (xi - x0);
        ry = (yi - y0);
        if (rx < -5)      { rx += 10; }
        else if (rx > 5)  { rx -= 10; }
      double rsq = rx*rx + ry*ry;
      double rho, u, v, P;
      double du, dv;
      rho = pow(1.0 - ((gamma-1.0)*beta*beta)/(8.0*gamma*pi*pi) * exp(1.0-rsq), 1.0/(gamma-1.0));
      P   = pow(rho,gamma);
      du  = - beta/(2.0*pi) * exp(0.5*(1.0-rsq)) * ry;
      dv  =   beta/(2.0*pi) * exp(0.5*(1.0-rsq)) * rx;
      u   = u_inf + du;
      v   = v_inf + dv;
      U_rho_exact[i][j]= rho;
      U_rhou_exact[i][j] = rho*u;
      U_rhov_exact[i][j] = rho*v;
      U_E_exact[i][j] = P/(gamma-1.0) + 0.5*rho*(u*u+v*v);
      }
    }
    vector<vector<double>> p_f_exact(n,vector<double> (n)), e_f_exact(n,vector<double> (n));
    ofstream mon_flux3;
    mon_flux3.open("Resultats/Exact", ios::out);
    for (int i = 2; i < n/fichier-2; i++)
    // for (int i = 0; i < n; i++)
    {
        kx=i*fichier;
        xi = x0_1 + kx*delta_x;
        for (int j = 2; j < n/fichier-2; j++)
        // for (int j = 0; j < n; j++)
        {
            ky=j*fichier;
            yi = y0_1 + ky*delta_y;

            U_int[0] = U_rho_exact[kx][ky];
            U_int[1] = U_rhou_exact[kx][ky];
            U_int[2] = U_rhov_exact[kx][ky];
            U_int[3] = U_E_exact[kx][ky];

            p_f_exact[kx][ky] = pression(U_int, gamma);
            e_f_exact[kx][ky] = (p_f[kx][ky]*(1-b*U_int[0]))/((gamma-1)*U_int[0]);
            mon_flux3 << xi << " " << yi << " " << U_int[0] << " " << U_int[1]/U_int[0] << " " << U_int[2]/U_int[0] << " " << p_f[kx][ky] << " " << e_f[kx][ky] << endl;
        }
    }
    mon_flux3.close();

    error(U_rho, U_rhou, U_rhov, U_E, U_rho_exact, U_rhou_exact, U_rhov_exact, U_E_exact, result, n, delta_x, delta_y, gamma, b, fichier);
    

    return 0;
}