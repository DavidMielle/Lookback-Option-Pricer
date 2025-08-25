--------------Tuto DLL sur VS Community 2022--------------------
1) Créer un nouveau projet : https://cdn.discordapp.com/attachments/950478077442543656/1167816304946851941/image.png?ex=654f80b0&is=653d0bb0&hm=b5e34f44f5c795779f011b7a835bcfa6b377488da5a6eeefbd41c5850b70e95d&
2) Dans un fichier quelconque, rendre les fonctions C++ disponibles en implémentant :

```cpp

extern "C" {
    __declspec(dllexport) double lookback_option_price(double T, bool isPut, double S0, double r, double sigma, int num_simulations, double dt) {
        return super appel à ta fonction;
    }
}

```
3) Côté VBA, implémenter :

```vba

Public Declare PtrSafe Function lookback_option_price Lib _
"chemin_dll" (ByVal T As Double, ByVal bool As Boolean, ByVal S0 As Double, ByVal r As Double, ByVal sigma As Double, ByVal num_simulations As Long, ByVal dt As Double) As Double

Function LookbackOptionPriceWrapper(T As Double, bool As Boolean, S0 As Double, r As Double, sigma As Double, num_simulations As Long, dt As Double) As Double
    LookbackOptionPriceWrapper = lookback_option_price(T, bool, S0, r, sigma, num_simulations, dt)
End Function

```
