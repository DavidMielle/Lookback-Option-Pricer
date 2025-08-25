//
// Created by Antoine on 30/10/2023.
//

#ifndef PROJETINFOM2QF_LOOKBACKOPTION_H
#define PROJETINFOM2QF_LOOKBACKOPTION_H

#include <vector>
#include <random>

class LookbackOption {
public:
    enum class OptionType { Call, Put };

    // Constructeur
    LookbackOption(double S0, double r, double sigma, double T, OptionType type, int numSimulations);

    // Méthodes pour calculer le prix et les grecques
    double calculatePrice();
    double calculateDelta();
    double calculateGamma();
    double calculateTheta();
    double calculateRho();
    double calculateVega();

private:
    double S0;            // Prix initial du sous-jacent
    double r;             // Taux d'intérêt sans risque
    double sigma;         // Volatilité
    double T;             // Maturité
    OptionType type;      // Type d'option (Call ou Put)
    int numSimulations;   // Nombre de simulations pour Monte Carlo

    std::vector<double> simulatedPrices;  // Stocke les prix simulés du sous-jacent

    // Méthodes privées pour aider aux calculs
    double simulatePath();
    double payoff(const std::vector<double>& prices);
    double blackScholesPrice(double S, double K);
    static double normalRandom();  // Génère un nombre aléatoire suivant une distribution normale
};


#endif //PROJETINFOM2QF_LOOKBACKOPTION_H
