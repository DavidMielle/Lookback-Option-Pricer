#include <iostream>
#include <cmath>
#include "LookbackOptionDavid.h"

int main()
{
    // Paramètres d'entrée
    double S0 = 50;           // Prix initial du sous-jacent
    double r = 0.1;           // Taux d'intérêt sans risque
    double sigma = 0.4;        // Volatilité
    double T = 0.25;            // Maturité
    int numSimulations = 10000; // Nombre de simulations Monte Carlo

    // Création de l'objet option de type Call
    LookbackOption callOption(S0, r, sigma, T, LookbackOption::OptionType::Call, numSimulations);

    // Simulation du prix de l'option Call
    double callPrice = callOption.calculatePrice();
    std::cout << "Prix de l'option Lookback Call David: " << callPrice << std::endl;

    // Création de l'objet option de type Put
    LookbackOption putOption(S0, r, sigma, T, LookbackOption::OptionType::Put, numSimulations);

    // Simulation du prix de l'option Put
    double putPrice = putOption.calculatePrice();
    std::cout << "Prix de l'option Lookback Put David: " << putPrice << std::endl;

    //grecques call
    /*double deltaCall = callOption.calculateDelta();
    std::cout << "delta call : " << deltaCall << std::endl;

    double gammaCall = callOption.calculateGamma();
    std::cout << "gamma call : " << gammaCall << std::endl;

    double thetaCall = callOption.calculateTheta();
    std::cout << "theta call : " << thetaCall << std::endl;

    double VegaCall = callOption.calculateVega();
    std::cout << "Vega call : " << VegaCall << std::endl;

    double rhoCall = callOption.calculateRho();
    std::cout << "rho call : " << rhoCall << std::endl;*/




    return 0;
}
