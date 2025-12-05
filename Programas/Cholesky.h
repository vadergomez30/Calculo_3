#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib> 
using namespace std;
void imprimirMatriz(const vector<vector<double>>& M, const string& nombre) {
    int n = M.size();
    cout << "\n" << nombre << ":\n";
    cout << fixed << setprecision(7);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(12) << M[i][j];
        cout << "\n";
    }
}

void imprimirVector(const vector<double>& v, const string& nombre) {
    cout << "\n" << nombre << ":\n";
    for (int i = 0; i < (int)v.size(); i++) {
        cout << nombre << "[" << i+1 << "] = " << v[i] << "\n";
    }
}
bool cholesky(const vector<vector<double>>& A,
              const vector<double>& b,
              vector<vector<double>>& L,
                            vector<vector<double>>& U,
              vector<double>& xChol) {
    int n = A.size();
    L.assign(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {

            double suma = 0.0;
            for (int k = 0; k < j; k++)
                suma += L[i][k] * L[j][k];

            if (i == j) {
                double val = A[i][i] - suma;
                if (val <= 0.0) {
                    cerr << "ERROR: La matriz A no es definida positiva (Cholesky).\n";
                    system("pause");
                    return false;
                }
                L[i][j] = sqrt(val);
            } else {
                L[i][j] = (A[i][j] - suma) / L[j][j];
            }
            imprimirMatriz(L, "Matriz L (Cholesky - parcial)");
        }
    }
    U.assign(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U[i][j] = L[j][i];

    imprimirMatriz(U, "Matriz U = L^T (Cholesky)");
    vector<double> c(n);
    for (int i = 0; i < n; i++) {
        double suma = 0.0;
        for (int k = 0; k < i; k++)
            suma += L[i][k] * c[k];
        c[i] = (b[i] - suma) / L[i][i];
    }

     
    xChol.assign(n, 0.0);
    for (int i = n-1; i >= 0; i--) {
        double suma = 0.0;
        for (int k = i+1; k < n; k++)
            suma += U[i][k] * xChol[k];
        xChol[i] = (c[i] - suma) / U[i][i];
    }

    return true;
}
double sumaLU(const vector<vector<double>>& L,
              const vector<vector<double>>& U,
              int i, int j, int k) {
    if (k < 0) return 0.0;
    return L[i][k] * U[k][j] + sumaLU(L, U, i, j, k - 1);
}

double sumaLc(const vector<vector<double>>& L,
              const vector<double>& c,
              int i, int k) {
    if (k < 0) return 0.0;
    return L[i][k] * c[k] + sumaLc(L, c, i, k - 1);
}

double sumaUx(const vector<vector<double>>& U,
              const vector<double>& x,
              int i, int k) {
    if (k <= i) return 0.0;
    return U[i][k] * x[k] + sumaUx(U, x, i, k - 1);
}

bool doolittleRec(int p,
                  const vector<vector<double>>& A,
                  vector<vector<double>>& L,
                  vector<vector<double>>& U,
                  int n) {
    if (p == n) return true;
    for (int j = p; j < n; ++j) {
        double suma = sumaLU(L, U, p, j, p - 1);
        U[p][j] = A[p][j] - suma;
        imprimirMatriz(U, "Matriz U (Doolittle - parcial)");
    }

    if (fabs(U[p][p]) < 1e-12) {
        cerr << "ERROR: pivote casi cero en Doolittle U[" << p << "][" << p
             << "]. La matriz puede requerir pivoteo.\n";
        system("pause");
        return false;
    }
    for (int i = p + 1; i < n; ++i) {
        double suma = sumaLU(L, U, i, p, p - 1);
        L[i][p] = (A[i][p] - suma) / U[p][p];
        imprimirMatriz(L, "Matriz L (Doolittle - parcial)");
    }

    return doolittleRec(p + 1, A, L, U, n);
}

void forwardSubRec(const vector<vector<double>>& L,
                   const vector<double>& b,
                   vector<double>& c,
                   int i, int n) {
    if (i == n) return;
    double suma = sumaLc(L, c, i, i - 1);
    c[i] = (b[i] - suma) / L[i][i];
    forwardSubRec(L, b, c, i + 1, n);
}

void backSubRec(const vector<vector<double>>& U,
                const vector<double>& c,
                vector<double>& x,
                int i) {
    if (i < 0) return;
    int n = x.size();
    double suma = sumaUx(U, x, i, n - 1);
    x[i] = (c[i] - suma) / U[i][i];
    backSubRec(U, c, x, i - 1);
}

bool doolittle(const vector<vector<double>>& A,
               const vector<double>& b,
               vector<vector<double>>& L,
               vector<vector<double>>& U,
               vector<double>& xDol) {
    int n = A.size();
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
        L[i][i] = 1.0;

    if (!doolittleRec(0, A, L, U, n))
        return false;

    imprimirMatriz(L, "Matriz L (Doolittle - final)");
    imprimirMatriz(U, "Matriz U (Doolittle - final)");

    vector<double> c(n, 0.0);
    forwardSubRec(L, b, c, 0, n);

    xDol.assign(n, 0.0);
    backSubRec(U, c, xDol, n - 1);

    return true;
}
bool croutTridiagonal(const vector<vector<double>>& A,
                      const vector<double>& b,
                      vector<vector<double>>& L,
                      vector<vector<double>>& U,
                      vector<double>& xCrout) {
    int n = A.size();
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i)
        U[i][i] = 1.0;
    if (n == 1) {
        L[0][0] = A[0][0];
        imprimirMatriz(L, "Matriz L (Crout - parcial)");
        imprimirMatriz(U, "Matriz U (Crout - parcial)");
    } else {
        L[0][0] = A[0][0];
        if (fabs(L[0][0]) < 1e-12) {
            cerr << "ERROR: L[1][1] = 0 en Crout.\n";
            system("pause");
            return false;
        }
        imprimirMatriz(L, "Matriz L (Crout - parcial)");

        U[0][1] = A[0][1] / L[0][0];   
        imprimirMatriz(U, "Matriz U (Crout - parcial)");
        for (int i = 1; i <= n-2; ++i) {
            L[i][i-1] = A[i][i-1];
            imprimirMatriz(L, "Matriz L (Crout - parcial)");
            L[i][i] = A[i][i] - L[i][i-1] * U[i-1][i];
            if (fabs(L[i][i]) < 1e-12) {
                cerr << "ERROR: L[" << i+1 << "][" << i+1 << "] = 0 en Crout.\n";
                system("pause");
                return false;
            }
            imprimirMatriz(L, "Matriz L (Crout - parcial)");
            U[i][i+1] = A[i][i+1] / L[i][i];
            imprimirMatriz(U, "Matriz U (Crout - parcial)");
        }

        int i = n - 1;
        L[i][i-1] = A[i][i-1];
        imprimirMatriz(L, "Matriz L (Crout - parcial)");

        L[i][i] = A[i][i] - L[i][i-1] * U[i-1][i];
        if (fabs(L[i][i]) < 1e-12) {
            cerr << "ERROR: L[" << i+1 << "][" << i+1 << "] = 0 en Crout.\n";
            system("pause");
            return false;
        }
        imprimirMatriz(L, "Matriz L (Crout - parcial)");
    }
    imprimirMatriz(L, "Matriz L (Crout - final)");
    imprimirMatriz(U, "Matriz U (Crout - final)");
    vector<double> c(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double suma = 0.0;
        for (int j = 0; j < i; ++j)
            suma += L[i][j] * c[j];
        c[i] = (b[i] - suma) / L[i][i];
    }
    xCrout.assign(n, 0.0);
    for (int i = n-1; i >= 0; --i) {
        double suma = 0.0;
        for (int j = i+1; j < n; ++j)
            suma += U[i][j] * xCrout[j];
        xCrout[i] = (c[i] - suma) / U[i][i]; 
    }

    return true;
}

void menuLU() {
    cout << "Metodos para resolver Ax = b\n";
    cout << "1) Cholesky\n";
    cout << "2) Doolittle\n";
    cout << "3) Crout\n";
    cout << "Elige una opcion: ";
    int opcion;
    cin >> opcion;

    switch (opcion) {
        case 1: {
            cout << "\n=== Metodo de Cholesky ===\n";
            int n;
            cout << "Introduce el numero de variables (dimension n): ";
            cin >> n;

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

            vector<vector<double>> L, U;
            vector<double> xChol;
            if (cholesky(A, b, L, U, xChol)) {
                imprimirVector(xChol, "x (solucion por Cholesky)");
            }
            system("pause");
            break;
        }
        case 2: {
            cout << "\n=== Metodo de Doolittle (recursivo) ===\n";
            int n;
            cout << "Introduce el numero de variables (dimension n): ";
            cin >> n;

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

            vector<vector<double>> L, U;
            vector<double> xDol;
            if (doolittle(A, b, L, U, xDol)) {
                imprimirVector(xDol, "x (solucion por Doolittle)");
            }
            system("pause");
            break;
        }
        case 3: {
            cout << "\n=== Metodo de Crout (tridiagonal) ===\n";
            int n;
            cout << "Introduce el numero de variables (dimension n): ";
            cin >> n;

            vector<vector<double>> A(n, vector<double>(n));
            vector<double> b(n);

            cout << "\nIntroduce la matriz A tridiagonal (" << n << "x" << n << "):\n";
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

            vector<vector<double>> L, U;
            vector<double> xCrout;
            if (croutTridiagonal(A, b, L, U, xCrout)) {
                imprimirVector(xCrout, "x (solucion por Crout)");
            }
            system("pause");
            break;
        }
        default:
            cout << "Opcion no valida.\n";
            system("pause");
    }

    return;
}
