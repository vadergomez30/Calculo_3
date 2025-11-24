#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

int main() {
    int n;
    cout << "Metodo de Cholesky para Ax = b\n";
    cout << "Introduce el numero de variables (dimension n): ";
    cin >> n;

    // MATRIZ A y vector b
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    cout << "\nIntroduce la matriz A (" << n << "x" << n << "):\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            cout << "A[" << i+1 << "][" << j+1 << "]: ";
            cin >> A[i][j];
        }

    cout << "\nIntroduce el vector b:\n";
    for (int i = 0; i < n; i++) {
        cout << "b[" << i+1 << "]: ";
        cin >> b[i];
    }

    // MATRIZ L (triangular inferior)
    vector<vector<double>> L(n, vector<double>(n, 0.0));

    // FACTORIZACIÓN DE CHOLESKY
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {

            double suma = 0.0;
            for (int k = 0; k < j; k++)
                suma += L[i][k] * L[j][k];

            if (i == j) {
                double val = A[i][i] - suma;
                if (val <= 0.0) {
                    cerr << "ERROR: La matriz A no es definida positiva.\n";
                    return 1;
                }
                L[i][j] = sqrt(val);
            } else {
                L[i][j] = (A[i][j] - suma) / L[j][j];
            }
        }
    }

    cout << fixed << setprecision(7);

    // IMPRIMIR L
    cout << "\nMatriz L (triangular inferior):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(12) << L[i][j];
        cout << "\n";
    }

    // MATRIZ U = L^T
    vector<vector<double>> U(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U[i][j] = L[j][i];  // transposición

    cout << "\nMatriz U = L^T (triangular superior):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(12) << U[i][j];
        cout << "\n";
    }

    // SUSTITUCIÓN HACIA ADELANTE: L * c = b
    vector<double> c(n);
    for (int i = 0; i < n; i++) {
        double suma = 0.0;
        for (int k = 0; k < i; k++)
            suma += L[i][k] * c[k];
        c[i] = (b[i] - suma) / L[i][i];
    }

    cout << "\nVector c (solucion de Lc = b):\n";
    for (int i = 0; i < n; i++)
        cout << "c[" << i+1 << "] = " << c[i] << "\n";

    // SUSTITUCIÓN HACIA ATRÁS: U * x = c
    vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        double suma = 0.0;
        for (int k = i+1; k < n; k++)
            suma += U[i][k] * x[k];
        x[i] = (c[i] - suma) / U[i][i];
    }

    cout << "\nSolucion x del sistema Ax = b:\n";
    for (int i = 0; i < n; i++)
        cout << "x[" << i+1 << "] = " << x[i] << "\n";

    return 0;
}
