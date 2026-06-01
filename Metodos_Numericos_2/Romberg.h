#pragma once
#include <bits/stdc++.h>
using namespace std;
//  f1(x) = x^4 * sqrt(3 + 2x^2) / 3
double integ_f1(double x) {
    return (x * x * x * x * sqrt(3.0 + 2.0 * x * x)) / 3.0;
}

//  f2(x) = x^5 / (x^2 + 4)^(1/5)
double integ_f2(double x) {
    return (x * x * x * x * x) / pow(x * x + 4.0, 1.0 / 5.0);
}

double trapecio(double (*f)(double), double a, double b, int k) {
    int n = 1 << k;          // n = 2^k
    double h = (b - a) / n;
    double suma = f(a) + f(b);
    for (int i = 1; i < n; i++)
        suma += 2.0 * f(a + i * h);
    return suma * h / 2.0;
}

void romberg(double (*f)(double), double a, double b, int decimales, int filas) {

    // filas = numero de niveles k pedidos por el usuario (minimo 3 para llegar a O(h^6))
    if (filas < 3) filas = 3;

    // La tabla solo llega hasta la columna 2 (O(h^6))
    const int MAX_COL = 3;
    vector<vector<double>> R(filas, vector<double>(MAX_COL, 0.0));

    // Llenar columna 0: trapecio puro
    for (int k = 0; k < filas; k++)
        R[k][0] = trapecio(f, a, b, k);

    // Extrapolacion columna 1: O(h^4)
    // R[k][1] = (4^1 * R[k][0] - R[k-1][0]) / (4^1 - 1)
    for (int k = 1; k < filas; k++)
        R[k][1] = (4.0 * R[k][0] - R[k-1][0]) / 3.0;

    // Extrapolacion columna 2: O(h^6)
    // R[k][2] = (4^2 * R[k][1] - R[k-1][1]) / (4^2 - 1)
    for (int k = 2; k < filas; k++)
        R[k][2] = (16.0 * R[k][1] - R[k-1][1]) / 15.0;

    // ── Imprimir tabla ──
    cout << fixed << setprecision(decimales);
    cout << "\nTabla de Romberg:\n\n";
    cout << setw(6) << "k"
         << setw(decimales + 10) << "O(h^2)  [Trapecio]"
         << setw(decimales + 10) << "O(h^4)"
         << setw(decimales + 10) << "O(h^6)"
         << "\n";
    cout << string(6 + 3*(decimales + 10), '-') << "\n";

    for (int k = 0; k < filas; k++) {
        cout << setw(6) << k;
        cout << setw(decimales + 10) << R[k][0];
        if (k >= 1) cout << setw(decimales + 10) << R[k][1];
        else        cout << setw(decimales + 10) << " ";
        if (k >= 2) cout << setw(decimales + 10) << R[k][2];
        else        cout << setw(decimales + 10) << " ";
        cout << "\n";
    }

    cout << "\nResultado final (O(h^6)): " << R[filas-1][2] << "\n";
}
void integracion() {
    while (true) {
        cout << "\n======================================\n";
        cout << "  INTEGRACION NUMERICA - ROMBERG\n";
        cout << "======================================\n";
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
        if (opc == 3) { cout << "Saliendo.\n"; return; }
        if (opc != 1 && opc != 2) { cout << "Opcion invalida.\n"; continue; }

        double a, b;
        cout << "Ingrese el limite inferior (a): ";
        cin >> a;
        if (cin.fail()) { cin.clear(); cin.ignore(10000,'\n'); cout << "Valor invalido.\n"; continue; }

        cout << "Ingrese el limite superior (b): ";
        cin >> b;
        if (cin.fail()) { cin.clear(); cin.ignore(10000,'\n'); cout << "Valor invalido.\n"; continue; }

        if (a >= b) { cout << "Error: a debe ser menor que b.\n"; continue; }

        int decimales;
        cout << "Ingrese los digitos de precision (decimales a mostrar, ej. 6): ";
        cin >> decimales;
        if (cin.fail() || decimales < 1 || decimales > 15) {
            cin.clear(); cin.ignore(10000,'\n');
            cout << "Valor invalido. Ingrese un entero entre 1 y 15.\n"; continue;
        }

        int filas;
        cout << "Ingrese el numero de filas de la tabla (minimo 3 para llegar a O(h^6)): ";
        cin >> filas;
        if (cin.fail() || filas < 3) {
            cin.clear(); cin.ignore(10000,'\n');
            cout << "Se usaran 3 filas (minimo).\n";
            filas = 3;
        }

        cout << "\nIntegrando";
        if (opc == 1) cout << " f(x) = x^4 * sqrt(3+2x^2)/3";
        else          cout << " f(x) = x^5 / (x^2+4)^(1/5)";
        cout << "  en [" << a << ", " << b << "]\n";

        if (opc == 1) romberg(integ_f1, a, b, decimales, filas);
        else          romberg(integ_f2, a, b, decimales, filas);

        cout << "\nDesea integrar otra funcion? (s/n): ";
        char cont; cin >> cont;
        if (cont == 'n' || cont == 'N') { cout << "Saliendo.\n"; return; }
    }
}
