#include <iostream>
#include <cmath>
#include "LookbackOption.h"

int main() {
    // Paramètres d'entrée
    double S0 = 100;       // Prix initial du sous-jacent
    double r = 0.05;       // Taux d'intérêt sans risque
    double sigma = 0.2;    // Volatilité
    double T = 1.0;        // Maturité
    int numSimulations = 10000; // Nombre de simulations Monte Carlo

    // Création de l'objet option de type Call
    LookbackOption callOption(S0, r, sigma, T, LookbackOption::OptionType::Call, numSimulations);

    // Simulation du prix de l'option Call
    double callPrice = callOption.calculatePrice();
    std::cout << "Prix de l'option Lookback Call : " << callPrice << std::endl;

    // Création de l'objet option de type Put
    LookbackOption putOption(S0, r, sigma, T, LookbackOption::OptionType::Put, numSimulations);

    // Simulation du prix de l'option Put
    double putPrice = putOption.calculatePrice();
    std::cout << "Prix de l'option Lookback Put : " << putPrice << std::endl;

    return 0;
}
