#include "LookbackOption.h"
#include <cmath>
#include <algorithm>
#include <iostream>

LookbackOption::LookbackOption(double S0, double r, double sigma, double T, OptionType type, int numSimulations)
        : S0(S0), r(r), sigma(sigma), T(T), type(type), numSimulations(numSimulations) {
    simulatedPrices.resize(numSimulations + 1);
    simulatedPrices[0] = S0;
}

double LookbackOption::simulatePath() {
    double dt = T / numSimulations;
    for (int i = 1; i <= numSimulations; ++i) {
        simulatedPrices[i] = simulatedPrices[i - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * normalRandom());
    }
    return simulatedPrices.back();
}

double LookbackOption::payoff(const std::vector<double>& prices) {
    double Smin = *std::min_element(prices.begin(), prices.end());
    double Smax = *std::max_element(prices.begin(), prices.end());
    if (type == OptionType::Call) {
        return std::max(prices.back() - Smin, 0.0);
    } else {
        return std::max(Smax - prices.back(), 0.0);
    }
}

double LookbackOption::blackScholesPrice(double S, double K) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    double N_d1 = 0.5 * (1.0 + erf(d1 / sqrt(2.0)));
    double N_d2 = 0.5 * (1.0 + erf(d2 / sqrt(2.0)));

    if (type == OptionType::Call) {
        return S * N_d1 - K * exp(-r * T) * N_d2;
    } else {
        return K * exp(-r * T) * (1 - N_d2) - S * (1 - N_d1);
    }
}


double LookbackOption::normalRandom() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(gen);
}

double LookbackOption::calculatePrice() {
    double sum = 0.0;
    for (int i = 0; i < numSimulations; ++i) {
        simulatePath();
        sum += payoff(simulatedPrices);
    }
    return (exp(-r * T) * sum) / numSimulations;
}

// J'utilise la Méthode des différences finies pour calculer les parametres.
// Est ce que c'est optimal ? Le code met pas mal de temps à run,
// peut-être que ça vient de là ? à optimiser ?

double LookbackOption::calculateDelta() {
    double dS = 0.01 * S0;  // Petite variation du prix du sous-jacent S0 vers le haut
    double priceUp = calculatePrice();
    S0 -= 2 * dS; // Petite variation du prix du sous-jacent S0 vers le bas
    double priceDown = calculatePrice();
    S0 += dS;  // Remettre S0 à sa valeur originale
    return (priceUp - priceDown) / (2 * dS);
}

double LookbackOption::calculateGamma() {
    double dS = 0.01 * S0;
    double priceUp = calculatePrice();
    double originalPrice = calculatePrice();
    S0 -= 2 * dS;
    double priceDown = calculatePrice();
    S0 += dS;
    return (priceUp - 2 * originalPrice + priceDown) / (dS * dS);
}

double LookbackOption::calculateTheta() {
    double dt = 0.01 * T;
    T -= dt;
    double priceLater = calculatePrice();
    T += dt;
    return (priceLater - calculatePrice()) / dt;
}

double LookbackOption::calculateRho() {
    double dr = 0.01 * r;
    r += dr;
    double priceUp = calculatePrice();
    r -= 2 * dr;
    double priceDown = calculatePrice();
    r += dr;
    return (priceUp - priceDown) / (2 * dr);
}

double LookbackOption::calculateVega() {
    double dSigma = 0.01 * sigma;
    sigma += dSigma;
    double priceUp = calculatePrice();
    sigma -= 2 * dSigma;
    double priceDown = calculatePrice();
    sigma += dSigma;
    return (priceUp - priceDown) / (2 * dSigma);
}


