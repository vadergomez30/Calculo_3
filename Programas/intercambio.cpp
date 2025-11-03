#include <bits/stdc++.h>
using namespace std;

/*
  Método de Intercambio (Pivoteo Total) para resolver Ax = b
  - Busca el mayor |A[i][j]| en submatriz k..n-1
  - Intercambia filas y columnas
  - Elimina hacia adelante
  - Sustitución hacia atrás en U*y = b~
  - Deshace permutación de columnas: x[perm[col]] = y[col]
*/

struct SolveResult {
    bool ok;
    vector<double> x;
    string msg;
};

static const double EPS = 1e-12;

SolveResult resolverIntercambio(vector<vector<double>> A, vector<double> b) {
    const int n = (int)A.size();
    if (n == 0 || (int)A[0].size() != n || (int)b.size() != n) {
        return {false, {}, "Dimensiones inválidas de A o b."};
    }

    // perm[col] = índice ORIGINAL de la variable que ahora está en la columna 'col'
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);

    for (int k = 0; k < n; ++k) {
        // 1) Buscar pivote máximo en submatriz A[k..n-1][k..n-1]
        int piv_i = k, piv_j = k;
        double max_val = 0.0;
        for (int i = k; i < n; ++i) {
            for (int j = k; j < n; ++j) {
                double val = fabs(A[i][j]);
                if (val > max_val) {
                    max_val = val;
                    piv_i = i; piv_j = j;
                }
            }
        }

        if (max_val < EPS) {
            return {false, {}, "La matriz es (casi) singular (pivote ~ 0)."};
        }

        // 2) Intercambios fila y columna
        if (piv_i != k) {
            swap(A[piv_i], A[k]);
            swap(b[piv_i], b[k]);
        }
        if (piv_j != k) {
            for (int i = 0; i < n; ++i) swap(A[i][piv_j], A[i][k]);
            swap(perm[piv_j], perm[k]); // MUY IMPORTANTE
        }

        // 3) Eliminación hacia adelante
        const double akk = A[k][k];
        if (fabs(akk) < EPS) return {false, {}, "Pivote cero tras permutar (inestable)."};
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / akk;
            A[i][k] = 0.0;
            for (int j = k + 1; j < n; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    // 4) Sustitución hacia atrás: U * y = b~
    vector<double> y(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[i][j] * y[j];
        if (fabs(A[i][i]) < EPS) return {false, {}, "Pivote cero en sustitución hacia atrás."};
        y[i] = s / A[i][i];
    }

    // 5) Deshacer permutación de columnas
    vector<double> x(n, 0.0);
    for (int col = 0; col < n; ++col) {
        int orig = perm[col];   // índice original de la variable en esta columna
        x[orig] = y[col];
    }

    return {true, x, "OK"};
}

int main() {
    ios::sync_with_stdio(false);

    int n;
    cout << "Ingrese el tamaño n: ";
    if (!(cin >> n) || n <= 0) {
        cerr << "n inválido.\n";
        return 0;
    }

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    cout << "\nIntroduce la matriz A (" << n << "x" << n << ") elemento por elemento:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << "A[" << i+1 << "][" << j+1 << "] = ";
            cin >> A[i][j];
        }
    }

    cout << "\nIntroduce el vector b (" << n << " valores):\n";
    for (int i = 0; i < n; ++i) {
        cout << "b[" << i+1 << "] = ";
        cin >> b[i];
    }

    SolveResult res = resolverIntercambio(A, b);
    if (!res.ok) {
        cerr << "\nError: " << res.msg << "\n";
        return 0;
    }

    cout << fixed << setprecision(10);
    cout << "\nSolución x (orden original de variables):\n";
    for (int i = 0; i < n; ++i)
        cout << "x" << (i + 1) << " = " << res.x[i] << "\n";
    return 0;
    
    cin.tie(nullptr);
}
