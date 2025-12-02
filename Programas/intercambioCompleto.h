#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdlib>  // <-- para system("pause")
using std::cin;
using std::cout;
using std::endl;
using std::fixed;
using std::setprecision;
using std::string;
using std::vector;

using Mat = vector<vector<double>>;

static const double EPS1 = 1e-12;

// ---------- Utilidades de impresión ----------
void printTable(const Mat& A,
                const vector<int>& rowLab,   // fila i representa x_{rowLab[i]}
                const vector<int>& colLab,   // col j representa b_{colLab[j]}
                const string& title,
                int width = 12, int prec = 6)
{
    cout << "\n" << title << "\n";
    cout << fixed << setprecision(prec);

    // Encabezado: columnas b?
    cout << std::setw(width) << "" << " |";
    for (size_t j = 0; j < A.size(); ++j) {
        cout << std::setw(width) << ("b" + std::to_string(colLab[j]));
    }
    cout << "\n";

    // Filas: x?
    for (size_t i = 0; i < A.size(); ++i) {
        cout << std::setw(width - 2) << ("x" + std::to_string(rowLab[i])) << " |";
        for (size_t j = 0; j < A.size(); ++j) {
            double y = (std::fabs(A[i][j]) < 1e-15 ? 0.0 : A[i][j]);
            cout << std::setw(width) << y;
        }
        cout << "\n";
    }
}

// Reacomodo final a A^{-1} (filas x1..xn, columnas b1..bn)
Mat rebuildInverse(const Mat& A,
                   const vector<int>& rowLab,
                   const vector<int>& colLab)
{
    int n = (int)A.size();
    Mat Ainv(n, vector<double>(n, 0.0));

    // Para cada xr (fila)
    for (int r = 1; r <= n; ++r) {
        int fila = -1;
        for (int i = 0; i < n; ++i) if (rowLab[i] == r) { fila = i; break; }

        // Para cada bs (columna)
        for (int s = 1; s <= n; ++s) {
            int col = -1;
            for (int j = 0; j < n; ++j) if (colLab[j] == s) { col = j; break; }
            Ainv[r-1][s-1] = A[fila][col];
        }
    }
    return Ainv;
}

vector<double> matVec(const Mat& M, const vector<double>& v)
{
    int n = (int)M.size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            y[i] += M[i][j] * v[j];
    return y;
}

// ---------- Método de Intercambio ----------
struct SolveResult {
    bool ok;
    Mat Ainv;
    vector<double> x;
    string msg;
};

SolveResult metodoIntercambio(Mat A, const vector<double>& b_in)
{
    int n = (int)A.size();
    if (n == 0 || (int)A[0].size() != n || (int)b_in.size() != n)
        return {false, {}, {}, "Dimensiones inválidas."};

    // Etiquetas iniciales
    vector<int> rowLab(n), colLab(n);
    for (int i = 0; i < n; ++i) { rowLab[i] = i + 1; colLab[i] = i + 1; }

    // n iteraciones
    for (int it = 1; it <= n; ++it) {

        // --- Paso 1: seleccionar pivote (l,k) = argmax |A_ij|
        int l = 0, k = 0;
        double best = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (std::fabs(A[i][j]) > best) {
                    best = std::fabs(A[i][j]);
                    l = i; k = j;
                }

        double alk = A[l][k];
        if (std::fabs(alk) < EPS1)
            return {false, {}, {}, "Pivote casi cero: matriz singular o mal condicionada."};

        cout << "\nIteracion " << it
             << ": entra x" << colLab[k]
             << "  y  sale b" << rowLab[l] << "\n";

        // --- Paso 2: a_ij = a_ij - (a_ik * a_lj)/a_lk    (i != l, j != k)
        for (int i = 0; i < n; ++i) if (i != l) {
            for (int j = 0; j < n; ++j) if (j != k) {
                A[i][j] = A[i][j] - (A[i][k] * A[l][j]) / alk;
            }
        }

        // --- Paso 3: columna pivote (i != l): a_ik = a_ik / a_lk
        for (int i = 0; i < n; ++i) if (i != l) {
            A[i][k] = A[i][k] / alk;
        }

        // --- Paso 4: renglón pivote (j != k): a_lj = - a_lj / a_lk
        for (int j = 0; j < n; ++j) if (j != k) {
            A[l][j] = - A[l][j] / alk;
        }

        // --- Paso 5: pivote: a_lk = 1 / a_lk
        A[l][k] = 1.0 / alk;

        // --- Paso 6: actualizar etiquetas (intercambio x_k <-> b_l)
        int tmp = rowLab[l];          // era b_l
        rowLab[l] = colLab[k];        // ahora la fila representa x_k
        colLab[k] = tmp;              // la columna representa b_l

        // Mostrar tabla de la iteración
        printTable(A, rowLab, colLab, "Estado despues de la iteracion " + std::to_string(it), 12, 6);
    }

    // --- Paso 7: reacomodo a A^{-1} (filas x1..xn, columnas b1..bn)
    Mat Ainv = rebuildInverse(A, rowLab, colLab);

    // --- Pasos 8-9: x = A^{-1} b
    vector<double> x = matVec(Ainv, b_in);

    return {true, Ainv, x, ""};
}

// ---------- Programa principal del método ----------
void Intercambio(){
    //std::ios::sync_with_stdio(false);

    int n;
    cout << "Ingrese n (tamano de A): ";
    if(!(cin >> n) || n <= 0){
        cout << "n invalido.\n";
        system("pause");
        return;
    }

    Mat A(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);

    cout << "Ingresa la matriz A (" << n << "x" << n << ") por filas:\n";
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            cout << "A["<<i+1<<","<<j+1<<"] = ";
            cin >> A[i][j];
        }
    }
    cout << "Ingresa el vector b (" << n << " componentes):\n";
    for(int i=0;i<n;++i){
        cout << "b["<<i+1<<"] = ";
        cin >> b[i];
    }

    auto res = metodoIntercambio(A, b);
    if(!res.ok){
        cout << "\nError: " << res.msg << "\n";
        system("pause");
        return;
    }

    // Mostrar A^{-1} reacomodada y solucion
    cout << "\nA^{-1} (filas x1..xn, columnas b1..bn):\n";
    cout << fixed << setprecision(6);
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            double y = (std::fabs(res.Ainv[i][j]) < 1e-15 ? 0.0 : res.Ainv[i][j]);
            cout << std::setw(12) << y;
        }
        cout << "\n";
    }

    cout << "\nSolucion x = A^{-1} b :\n";
    for(int i=0;i<n;++i){
        double y = (std::fabs(res.x[i]) < 1e-15 ? 0.0 : res.x[i]);
        cout << "x" << (i+1) << " = " << y << "\n";
    }
    cout << endl;

    system("pause");
    return;
}
