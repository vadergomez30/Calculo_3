#pragma once
#include <bits/stdc++.h>
using namespace std;
 
// ─────────────────────────────────────────────
//  Funciones disponibles
// ─────────────────────────────────────────────
 
//  f1(x) = x^4 * sqrt(3 + 2x^2) / 3
double integ_f1(double x) {
    return (x * x * x * x * sqrt(3.0 + 2.0 * x * x)) / 3.0;
}
 
//  f2(x) = x^5 / (x^2 + 4)^(1/5)
double integ_f2(double x) {
    return (x * x * x * x * x) / pow(x * x + 4.0, 1.0 / 5.0);
}
 
// ─────────────────────────────────────────────
//  Regla del trapecio compuesta
//  Divide [a,b] en 2^k subintervalos
// ─────────────────────────────────────────────
double trapecio(double (*f)(double), double a, double b, int k) {
    int n = 1 << k;          // n = 2^k
    double h = (b - a) / n;
    double suma = f(a) + f(b);
    for (int i = 1; i < n; i++)
        suma += 2.0 * f(a + i * h);
    return suma * h / 2.0;
}
 
// ─────────────────────────────────────────────
//  Extrapolacion de Romberg
//
//  Construye la tabla R donde:
//    R[k][0] = trapecio con 2^k subintervalos
//    R[k][j] = extrapolacion de Richardson de orden j
//
//  Criterio de parada: |R[k][k] - R[k-1][k-1]| < tol
//  tol se deriva de los digitos de precision pedidos: tol = 0.5 * 10^(-digitos)
// ─────────────────────────────────────────────
void romberg(double (*f)(double), double a, double b, int digitos) {
 
    double tol = 0.5 * pow(10.0, -digitos);
    int max_nivel = 20;   // limite de seguridad
 
    // Tabla de Romberg (se expande dinamicamente)
    vector<vector<double>> R(max_nivel + 1, vector<double>(max_nivel + 1, 0.0));
 
    cout << fixed << setprecision(digitos + 2);
    cout << "\nTabla de Romberg:\n";
    cout << setw(6) << "k";
    for (int j = 0; j <= 6; j++) cout << setw(18) << "O(h^" << 2*(j+1) << ")";
    cout << "\n";
 
    R[0][0] = trapecio(f, a, b, 0);
    cout << setw(6) << 0 << setw(18) << R[0][0] << "\n";
 
    for (int k = 1; k <= max_nivel; k++) {
        // Columna 0: trapecio con 2^k subintervalos
        R[k][0] = trapecio(f, a, b, k);
 
        // Extrapolaciones de Richardson
        for (int j = 1; j <= k; j++) {
            double factor = pow(4.0, j);
            R[k][j] = (factor * R[k][j-1] - R[k-1][j-1]) / (factor - 1.0);
        }
 
        // Imprimir fila
        cout << setw(6) << k;
        for (int j = 0; j <= k && j <= 6; j++)
            cout << setw(18) << R[k][j];
        cout << "\n";
 
        // Criterio de convergencia sobre la diagonal
        if (k >= 1 && fabs(R[k][k] - R[k-1][k-1]) < tol) {
            cout << "\nConvergio en k = " << k << "\n";
            cout << "Resultado: " << R[k][k] << "\n";
            cout << "Error estimado: " << fabs(R[k][k] - R[k-1][k-1]) << "\n";
            return;
        }
    }
 
    // Si no convergio, reportar la mejor estimacion
    cout << "\nNo convergio en " << max_nivel << " niveles.\n";
    cout << "Mejor estimacion: " << R[max_nivel][max_nivel] << "\n";
}
 
// ─────────────────────────────────────────────
//  Menu principal de integracion
// ─────────────────────────────────────────────
void integracion() {
    while (true) {
        cout << "     INTEGRACION NUMERICA - ROMBERG\n";
        cout << "Seleccione la funcion a integrar:\n\n";
        cout << "  1. f(x) = x^4 * sqrt(3 + 2x^2) / 3\n";
        cout << "  2. f(x) = x^5 / (x^2 + 4)^(1/5)\n";
        cout << "  3. Salir\n\n";
        cout << "Opcion: ";
 
        int opc;
        cin >> opc;
        if (cin.fail()) {
            cin.clear(); cin.ignore(10000, '\n');
            cout << "Opcion invalida.\n";
            continue;
        }
        if (opc == 3) {
            cout << "Saliendo.\n";
            return;
        }
        if (opc != 1 && opc != 2) {
            cout << "Opcion invalida.\n";
            continue;
        }
 
        double a, b;
        cout << "Ingrese el limite inferior de integracion (a): ";
        cin >> a;
        if (cin.fail()) {
            cin.clear(); cin.ignore(10000, '\n');
            cout << "Valor invalido.\n"; continue;
        }
 
        cout << "Ingrese el limite superior de integracion (b): ";
        cin >> b;
        if (cin.fail()) {
            cin.clear(); cin.ignore(10000, '\n');
            cout << "Valor invalido.\n"; continue;
        }
 
        if (a >= b) {
            cout << "Error: el limite inferior debe ser menor que el superior.\n";
            continue;
        }
 
        int digitos;
        cout << "Ingrese los digitos de precision deseados (ej. 6): ";
        cin >> digitos;
        if (cin.fail() || digitos < 1 || digitos > 15) {
            cin.clear(); cin.ignore(10000, '\n');
            cout << "Valor invalido. Ingrese un entero entre 1 y 15.\n"; continue;
        }
 
        cout << "\nIntegrando";
        if (opc == 1) cout << " f(x) = x^4 * sqrt(3 + 2x^2) / 3";
        else          cout << " f(x) = x^5 / (x^2 + 4)^(1/5)";
        cout << "  en [" << a << ", " << b << "]\n";
        cout << "Tolerancia: 0.5e-" << digitos << "\n";
 
        if (opc == 1) romberg(integ_f1, a, b, digitos);
        else          romberg(integ_f2, a, b, digitos);
 
        cout << "\nDesea integrar otra funcion? (s/n): ";
        char cont; cin >> cont;
        if (cont == 'n' || cont == 'N') {
            cout << "Saliendo.\n";
            return;
        }
    }
}
 