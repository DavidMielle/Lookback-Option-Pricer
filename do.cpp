#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

// Fonction pour calculer le minimum d'un vecteur entre deux indices
double calculateMin(const std::vector<double>& vec, std::size_t start, std::size_t end) {
    if (start < 0 || end >= vec.size() || start > end) {
        // Gestion d'erreurs ou traitement selon tes besoins
        return 0.0;  // Valeur par défaut à ajuster
    }

    // Utilisation de std::min_element pour trouver le minimum dans le vecteur
    auto minElement = std::min_element(vec.begin() + start, vec.begin() + end + 1);

    return *minElement;
}

// Fonction pour calculer le maximum d'un vecteur entre deux indices
double calculateMax(const std::vector<double>& vec, std::size_t start, std::size_t end) {
    if (start < 0 || end >= vec.size() || start > end) {
        // Gestion d'erreurs ou traitement selon tes besoins
        return 0.0;  // Valeur par défaut à ajuster
    }

    // Utilisation de std::max_element pour trouver le maximum dans le vecteur
    auto maxElement = std::max_element(vec.begin() + start, vec.begin() + end + 1);

    return *maxElement;
}

// Fonction pour calculer le mouvement brownien
double calculateBrownianMotion() {
    // Initialiser le générateur de nombres aléatoires
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<double> distribution(0.0, 1.0);

    // Retourner un échantillon de la distribution normale standard
    return distribution(gen);
}

// Fonction pour calculer le prix d'un call européen avec le modèle de Black-Scholes
double blackScholesCall(double S, double K, double r, double sigma, double T) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    double N_d1 = 0.5 * (1.0 + erf(d1 / sqrt(2.0)));
    double N_d2 = 0.5 * (1.0 + erf(d2 / sqrt(2.0)));

    double callPrice = S * N_d1 - K * exp(-r * T) * N_d2;

    return callPrice;
}

double blackScholesPut(double S, double K, double r, double sigma, double T){
    return blackScholesCall(S, K, r, sigma, T) - S + K * exp(-r * T);
}
// Fonction principale pour calculer les options lookback avec le modèle de Black-Scholes
void calculateLookbackOptionBS(double S0, double r, double sigma, double T, int numSimulations) {
    // Paramètres de simulation
    double dt = T / numSimulations;
    double Smin;
    double Smax;
    double dS = 1;
    double dr = 0.1;      // Petite variation pour le taux d'intérêt
    double dSigma = 0.1;  // Petite variation pour la volatilité

    // Vecteur pour stocker les prix du sous-jacent pendant la simulation
    std::vector<double> prices(numSimulations + 1);
    prices[0] = S0;

    //init
    Smin = S0;
    Smax = S0;

    // Boucle de simulation Monte Carlo
    for (int i = 1; i <= numSimulations; ++i) {
        // Simulation du mouvement brownien géométrique
        prices[i] = prices[i - 1] * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * calculateBrownianMotion());

        // Mettre à jour Smin et Smax
        Smin = calculateMin(prices, 0, i);
        Smax = calculateMax(prices, 0, i);

        // std::cout << "Smin : " << Smin << std::endl;
        // std::cout << "Smax : " << Smax << std::endl;
        // std::cout << "prices[" << i << "] : " << prices[i] << std::endl;
        // À faire : Accumuler les résultats pour les statistiques (si nécessaire)
    }

    // Calcul du prix de l'option lookback avec le modèle de Black-Scholes
    double callLookbackPayoff = std::max(prices.back() - Smin, 0.0);
    double putLookbackPayoff = std::max(Smax - prices.back(), 0.0);

    // Calcul du prix de l'option lookback européenne avec le modèle de Black-Scholes
    double callLookbackPrice = blackScholesCall(S0, callLookbackPayoff, r, sigma, T);
    double putLookbackPrice = blackScholesPut(S0, putLookbackPayoff, r, sigma, T);

    // Affichage des résu ltats
    // std::cout << "Call Lookback Payoff: " << callLookbackPayoff << std::endl;
    // std::cout << "Put Lookback Payoff: " << putLookbackPayoff << std::endl;
    std::cout << "Call Lookback Price: " << callLookbackPrice << std::endl;
    std::cout << "Put Lookback Price: " << putLookbackPrice << std::endl;

    // Calcul des grecques par différences finies
    double delta_call = (blackScholesCall(prices.back() + dS, callLookbackPayoff, r, sigma, T) -
                    blackScholesCall(prices.back() - dS, callLookbackPayoff, r, sigma, T)) / (2.0 * dS);

    double gamma_call = (blackScholesCall(prices.back() + dS, callLookbackPayoff, r, sigma, T) -
                    2.0 * blackScholesCall(prices.back(), callLookbackPayoff, r, sigma, T) +
                    blackScholesCall(prices.back() - dS, callLookbackPayoff, r, sigma, T)) / (dS * dS);

    double theta_call = (blackScholesCall(prices.back(), callLookbackPayoff, r, sigma, T - dt) -
                    blackScholesCall(prices.back(), callLookbackPayoff, r, sigma, T)) / dt;

    double rho_call = (blackScholesCall(prices.back(), callLookbackPayoff, r + dr, sigma, T) -
                    blackScholesCall(prices.back(), callLookbackPayoff, r - dr, sigma, T)) / (2.0 * dr);

    double vega_call = (blackScholesCall(prices.back(), callLookbackPayoff, r, sigma + dSigma, T) -
                    blackScholesCall(prices.back(), callLookbackPayoff, r, sigma - dSigma, T)) / (2.0 * dSigma);



    double delta_put = (blackScholesPut(prices.back() + dS, putLookbackPayoff, r, sigma, T) -
                    blackScholesPut(prices.back() - dS, putLookbackPayoff, r, sigma, T)) / (2.0 * dS);

    double gamma_put = (blackScholesPut(prices.back() + dS, putLookbackPayoff, r, sigma, T) -
                    2.0 * blackScholesPut(prices.back(), putLookbackPayoff, r, sigma, T) +
                    blackScholesPut(prices.back() - dS, putLookbackPayoff, r, sigma, T)) / (dS * dS);

    double theta_put = (blackScholesPut(prices.back(), putLookbackPayoff, r, sigma, T - dt) -
                    blackScholesPut(prices.back(), putLookbackPayoff, r, sigma, T)) / dt;

    double rho_put = (blackScholesPut(prices.back(), putLookbackPayoff, r + dr, sigma, T) -
                    blackScholesPut(prices.back(), putLookbackPayoff, r - dr, sigma, T)) / (2.0 * dr);

    double vega_put = (blackScholesPut(prices.back(), putLookbackPayoff, r, sigma + dSigma, T) -
                    blackScholesPut(prices.back(), putLookbackPayoff, r, sigma - dSigma, T)) / (2.0 * dSigma);

    // Affichage des résultats ou stockage dans un fichier Excel
    std::cout << "Delta_call : " << delta_call << std::endl;
    std::cout << "Gamma_call : " << gamma_call << std::endl;
    std::cout << "Theta_call : " << theta_call << std::endl;
    std::cout << "Rho_call : " << rho_call << std::endl;
    std::cout << "Vega_call : " << vega_call << std::endl;
    std::cout << "Delta_put : " << delta_put << std::endl;
    std::cout << "Gamma_put : " << gamma_put << std::endl;
    std::cout << "Theta_put : " << theta_put << std::endl;
    std::cout << "Rho_put : " << rho_put << std::endl;
    std::cout << "Vega_put : " << vega_put << std::endl;
}


int main() {
    // Entrées
    double S0, r, sigma, T;
    // double q = 0;
    S0 = 100;
    r = 0.1;
    sigma = 0.4;
    T = 0.25;
    int numSimulations;

    // Calcul de la CDF pour la distribution normale standard
    

    // À faire : Demander d'autres entrées nécessaires
    // compute_terms_put(S0, S0, r, q, sigma, T); //valeurs égales à celles du Hull

    // compute_initial_price_put(S0, S0, r, q, sigma, T); //pareil donc c'est good

    // Calculer les options lookback
    std::cout << "Entrez le nombre de simulations : ";
    std::cin >> numSimulations;
    calculateLookbackOptionBS(S0, r, sigma, T, numSimulations);

    return 0;
}







void compute_terms_call(double S0, double Smin, double r, double q, double sigma, double T){
    double a1, a2, a3, Y1;
    a1 = (log(S0/Smin) + (r - q + sigma*sigma/2)*T)/(sigma*sqrt(T));
    a2 = a1 - sigma*sqrt(T);
    a3 = (log(S0/Smin) + (q - r + sigma*sigma/2)*T)/(sigma*sqrt(T));
    Y1 = (-2*(r - q -sigma*sigma/2)*log(S0/Smin))/(sigma*sigma);

    std::cout << "a1 : " << a1 << ", a2 : " << a2 << ", a3 : " << a3 << ", Y1 : " << Y1 << std::endl; 
}

void compute_terms_put(double S0, double Smax, double r, double q, double sigma, double T){
    double b1, b2, b3, Y2;
    b1 = (log(Smax/S0) + (q - r + sigma*sigma/2)*T)/(sigma*sqrt(T));
    b2 = b1 - sigma*sqrt(T);
    b3 = (log(Smax/S0) + (r - q - sigma*sigma/2)*T)/(sigma*sqrt(T));
    Y2 = (2*(r - q -sigma*sigma/2)*log(Smax/S0))/(sigma*sigma);

    std::cout << "b1 : " << b1 << ", b2 : " << b2 << ", b3 : " << b3 << ", Y2 : " << Y2 << std::endl; 
}

double cdf(double variable){
    double cdf = 0.5 * (1.0 + std::erf(variable / std::sqrt(2.0)));
    return cdf;
}

double compute_initial_price_put(double S0, double Smax, double r, double q, double sigma, double T){
    double b1, b2, b3, Y2;
    double initPrice;
    b1 = (log(Smax/S0) + (q - r + sigma*sigma/2)*T)/(sigma*sqrt(T));
    b2 = b1 - sigma*sqrt(T);
    b3 = (log(Smax/S0) + (r - q - sigma*sigma/2)*T)/(sigma*sqrt(T));
    Y2 = (2*(r - q -sigma*sigma/2)*log(Smax/S0))/(sigma*sigma);


    std::cout << "test avant" << std::endl;
    initPrice = Smax*exp(-r*T)*(cdf(b1) - exp(Y2)*cdf(-b3)*(sigma*sigma)/(2*(r-q))) + S0*exp(-q*T)*(sigma*sigma)/(2*(r-q))*cdf(-b2) - S0*exp(-q*T)*cdf(b2); 
    std::cout << "test apres" << std::endl;
    std::cout << initPrice << std::endl;
    return initPrice;
}

