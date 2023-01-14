#ifndef FLUX_CPP

#include "timescheme.h"
#include "flux.h"
#include "operation.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string> 
#include "flux.h"
#include "operation.h"

using namespace std;

void compute_star(vector<double> U_l, vector<double> U_r, double & p_star, double & u_star, double gamma){

    p_star = 1;
    u_star = 1;

    double rho_bar, a_bar, a_l, a_r, p_l, p_r;

    p_l = pression(U_l, gamma);
    p_r = pression(U_r, gamma);

    a_l = sqrt(gamma*p_l/U_l[0]);
    a_r = sqrt(gamma*p_r/U_r[0]);


    rho_bar = (U_l[0] + U_r[0])/2.;
    a_bar = (a_l + a_r)/2.;

    p_star = (p_l + p_r)/2. - (U_r[1]/U_r[0] - U_l[1]/U_l[0])*rho_bar*a_bar/2.;
    u_star = (U_l[1]/U_l[0] + U_r[1]/U_r[0])/2. - (p_r - p_l)/(2.*rho_bar*a_bar);

}

void vect_S(vector<double> U_l, vector<double> U_r, double & S_l, double & S_m, double & S_r, double gamma){

    double a_l, a_r, p_star, u_star;
    double q_l, q_r, H_l, H_r;

    compute_star(U_l, U_r, p_star, u_star, gamma);

    H_l = p_star/pression(U_l, gamma);
    H_r = p_star/pression(U_r, gamma);

    a_l = sqrt(gamma*pression(U_l, gamma)/U_l[0]);
    a_r = sqrt(gamma*pression(U_r, gamma)/U_r[0]);

    if (H_l <= 1.)
    {
        q_l = 1.;
    }
    else{
        q_l = sqrt(1+((gamma+1)/(2*gamma))*(H_l-1.));
    }

    if (H_r <= 1.)
    {
        q_r = 1.;
    }
    else{
        q_r = sqrt(1+((gamma+1)/(2*gamma))*(H_r-1.));
    }


    S_l = U_l[1]/U_l[0] - a_l*q_l;
    S_m = u_star;
    S_r = U_r[1]/U_r[0] + a_r*q_r;
}

// Creation du flux Rusanov
vector<double> flux_rusanov(vector<double> U_l, vector<double> U_r, double gamma){

    vector<double> flux_pdemi(3,0), F_l(3,0), F_r(3,0);
    double p_l = pression(U_l, gamma);
    double p_r = pression(U_r, gamma);
    double vp;

    // Choisir en vp (valeurs propres max entre droite et gauche) et val_p (max de toutes les val p)
    if ((U_l[1]/U_l[0]<0) && (U_r[1]/U_r[0]<0)){
        vp = max(abs(U_l[1]/U_l[0]) + sqrt(gamma*pression(U_l, gamma)/U_l[0]), abs(U_r[1]/U_r[0]) + sqrt(gamma*pression(U_r, gamma)/U_r[0]));
    }
    else{
        vp = max(U_l[1]/U_l[0] + sqrt(gamma*pression(U_l, gamma)/U_l[0]), U_r[1]/U_r[0] + sqrt(gamma*pression(U_r, gamma)/U_r[0]));
    }
    
    
    F_l = F(U_l, p_l, U_l[2]);
    F_r = F(U_r, p_r, U_r[2]);

    for (int i = 0; i < U_l.size(); i++)
    {
        flux_pdemi[i] = (F_l[i]+F_r[i])/2. - (vp/2.)*(U_r[i]-U_l[i]);
    }

    return flux_pdemi;
}

// Creation du flux
vector<double> flux_HLL(vector<double> U_l, vector<double> U_r, double gamma){

    vector<double> flux_pdemi(3,0), F_l(3,0), F_m(3,0), F_r(3,0);
    double S_l, S_m, S_r, p_l, p_r;

    vect_S(U_l, U_r, S_l, S_m, S_r, gamma);

    p_l = pression(U_l, gamma);
    p_r = pression(U_r, gamma);

    F_l = F(U_l, p_l, U_l[2]);
    F_r = F(U_r, p_r, U_r[2]);

    for (int k = 0; k < U_l.size(); k++)
    {
        F_m[k] = (S_r*F_l[k] - S_l*F_r[k] + S_l*S_r*(U_r[k] - U_l[k]))/(S_r - S_l);
    }

    // cout << "Sl = " << S_l << " Sm = " << S_m <<  " et Sr = " << S_r << endl;

    if (S_l >= 0)
    {
        for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_l[k];
        }
    }
    else if ((S_l <= 0) && (S_r>=0))
    {
       for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_m[k];
        }
    }
    else if (S_r <= 0)
    {
        for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_r[k];
        }
    }
    else{
        cout << "Problème sur les S_k" << endl;
    }

    return flux_pdemi;
}



// Creation du flux
vector<double> flux_HLLC(vector<double> U_l, vector<double> U_r, double gamma){

    vector<double> flux_pdemi(3,0), F_l(3,0), F_r(3,0);
    vector<double> U_lstar(3,0), U_rstar(3,0), F_lstar(3,0), F_rstar(3,0);

    double S_l, S_star, S_r, p_l, p_r;

    // On recupere les pressions
    p_l = pression(U_l, gamma);
    p_r = pression(U_r, gamma);

    // On recupere S_k
    double S_m;
    vect_S(U_l, U_r, S_l, S_m, S_r, gamma);
    // cout << "-----------------------------------------------------------------------" << endl;
    // cout << "Sm = " << S_m << endl;
    // S_l = min(U_l[1]/U_l[0] - sqrt(gamma*p_l/U_l[0]), U_r[1]/U_r[0] - sqrt(gamma*p_l/U_r[0]));
    // S_r = max(U_l[1]/U_l[0] + sqrt(gamma*p_l/U_r[0]), U_r[1]/U_r[0] + sqrt(gamma*p_r/U_r[0]));

    // Calcul de S_star
    double S_star_num, S_star_den;

    S_star_num = p_r - p_l + U_l[1]*(S_l - U_l[1]/U_l[0]) - U_r[1]*(S_r - U_r[1]/U_r[0]);
    S_star_den = U_l[0]*(S_l - U_l[1]/U_l[0]) - U_r[0]*(S_r - U_r[1]/U_r[0]);

    S_star = S_star_num/S_star_den;
    // cout << "S_star = " << S_star << endl;

    // On recupere F_r et F_l
    F_l = F(U_l, p_l, U_l[2]);
    F_r = F(U_r, p_r, U_r[2]);

    // Creation des U_star
    double U_l_3 = U_l[2]/U_l[0] + (S_star - U_l[1]/U_l[0]) * (S_star + p_l/(U_l[0] * (S_l - U_l[1]/U_l[0])));

    U_lstar[0] = U_l[0]*((S_l - U_l[1]/U_l[0])/(S_l - S_star));
    U_lstar[1] = U_l[0]*((S_l - U_l[1]/U_l[0])/(S_l - S_star))*S_star;
    U_lstar[2] = U_l[0]*((S_l - U_l[1]/U_l[0])/(S_l - S_star)) * U_l_3;

    double U_r_3 = U_r[2]/U_r[0] + (S_star - U_r[1]/U_r[0]) * (S_star + (p_r/(U_r[0] * (S_r - U_r[1]/U_r[0]))));

    U_rstar[0] = U_r[0]*((S_r - U_r[1]/U_r[0])/(S_r - S_star));
    U_rstar[1] = U_r[0]*((S_r - U_r[1]/U_r[0])/(S_r - S_star))*S_star;
    U_rstar[2] = U_r[0]*((S_r - U_r[1]/U_r[0])/(S_r - S_star)) * U_r_3;

    // cout << "Sl = " << S_l << " S* = " << S_star << " et Sr = " << S_r << endl;

    // Calcul des F_kstar
    for (int k = 0; k < U_l.size(); k++)
    {
        F_lstar[k] = F_l[k] + S_l*(U_lstar[k] - U_l[k]);
        F_rstar[k] = F_r[k] + S_r*(U_rstar[k] - U_r[k]);
    }

    ///////////////////////////////// Variation de HLLC
    // vector<double> D_star(3,0);
    // D_star[0] = 0.;
    // D_star[1] = 1.;
    // D_star[2] = S_star;

    // Variation 1
    // for (int k = 0; k < U_l.size(); k++)
    // {
    //     F_lstar[k] = (S_star*(S_l*U_l[k] - F_l[k]) + S_l*(p_l + U_l[0]*(S_l - U_l[1]/U_l[0])*(S_star - U_l[1]/U_l[0]))*D_star[k])/(S_l - S_star);
    //     F_rstar[k] = (S_star*(S_r*U_r[k] - F_r[k]) + S_r*(p_r + U_r[0]*(S_r - U_r[1]/U_r[0])*(S_star - U_r[1]/U_r[0]))*D_star[k])/(S_r - S_star);
    // }

    // Variation 2
    // double p_LR;
    // p_LR = (p_l + p_r + U_l[0]*(S_l - U_l[1]/U_l[0])*(S_star - U_l[1]/U_l[0]) + U_r[0]*(S_r - U_r[1]/U_r[0])*(S_star - U_r[1]/U_r[0]))/2;

    // for (int k = 0; k < U_l.size(); k++)
    // {
    //     F_lstar[k] = (S_star*(S_l*U_l[k] - F_l[k]) + S_l*p_LR*D_star[k])/(S_l - S_star);
    //     F_rstar[k] = (S_star*(S_r*U_r[k] - F_r[k]) + S_r*p_LR*D_star[k])/(S_r - S_star);
    // }

    // Calcul de flux suivant les differents cas
    if (S_l >= 0)
    {
        for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_l[k];
        }
    }
    else if ((S_l < 0) && (S_star>0))
    {
       for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_lstar[k];
        }
    }
    else if ((S_star <= 0) && (S_r>0))
    {
       for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_rstar[k];
        }
    }
    else if (S_r <= 0)
    {
        for (int k = 0; k < U_l.size(); k++)
        {
            flux_pdemi[k] = F_r[k];
        }
    }
    else{
        cout << "Problème sur les S_k" << endl;
    }

    return flux_pdemi;
}

#define FLUX_CPP
#endif