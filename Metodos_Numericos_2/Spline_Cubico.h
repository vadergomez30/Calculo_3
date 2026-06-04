#include <bits/stdc++.h>
using namespace std;


typedef vector<vector<double>> Mat;

string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == string::npos) ? "" : s.substr(a, b - a + 1);
}
bool leerDesdeCSV(const string& nombreArchivo, vector<vector<double>>& puntos) {
    ifstream archivo(nombreArchivo);
    if (!archivo.is_open()) {
        cout << "Error: No se pudo abrir el archivo '" << nombreArchivo << "'.\n";
        return false;
    }

    puntos.clear();
    string linea;
    int numLinea = 0;
    while (getline(archivo, linea)) {
        numLinea++;
        linea = trim(linea);
        if (linea.empty()) continue;

    
        char sep = (linea.find(';') != string::npos) ? ';' : ',';

        
        stringstream ss(linea);
        string token1, token2;
        if (!getline(ss, token1, sep) || !getline(ss, token2, sep)) {
            cout << "Advertencia: línea " << numLinea << " mal formada, se omite.\n";
            continue;
        }
        token1 = trim(token1);
        token2 = trim(token2);

    
        try {
            double x = stod(token1);
            double y = stod(token2);
            puntos.push_back({x, y});
        } catch (...) {
            if (numLinea == 1)
                cout << "Info: Se omitió la primera línea (posible encabezado).\n";
            else
                cout << "Advertencia: línea " << numLinea << " no numérica, se omite.\n";
        }
    }

    if (puntos.size() < 2) {
        cout << "Error: Se necesitan al menos 2 puntos válidos en el CSV.\n";
        return false;
    }
    return true;
}

bool leerManual(vector<vector<double>>& puntos) {
    int n;
    cout << "Ingrese la cantidad de puntos: ";
    cin >> n;
    if (cin.fail() || n < 2) {
        cout << "Error: Ingrese un número entero >= 2.\n";
        cin.clear();
        cin.ignore(10000, '\n');
        return false;
    }

    puntos.assign(n, vector<double>(2));
    cout << "Ingrese los puntos (x y) separados por espacio:\n";
    for (int i = 0; i < n; i++) {
        cout << "  Punto " << i << ": ";
        cin >> puntos[i][0] >> puntos[i][1];
        if (cin.fail()) {
            cout << "Error: Ingrese solo números.\n";
            cin.clear();
            cin.ignore(10000, '\n');
            return false;
        }
    }

    cout << "\nLos puntos ingresados son:\n";
    for (int i = 0; i < n; i++)
        cout << "  " << i << ": (" << puntos[i][0] << ", " << puntos[i][1] << ")\n";

    cout << "¿Son correctos? (s/n): ";
    char resp;
    cin >> resp;
    if (resp == 'n' || resp == 'N') {
        cout << "Ingresa el índice a corregir: ";
        int idx;
        cin >> idx;
        if (idx < 0 || idx >= n) {
            cout << "Índice fuera de rango.\n";
            return false;
        }
        cout << "Nuevo valor de x: "; cin >> puntos[idx][0];
        cout << "Nuevo valor de y: "; cin >> puntos[idx][1];
        if (cin.fail()) {
            cout << "Error: valor no numérico.\n";
            cin.clear(); cin.ignore(10000, '\n');
            return false;
        }
    }
    return true;
}

//  Cálculo del spline cúbico natural
//
//  Convención:  S_i(x) = a_i*(x-x_i)^3 + b_i*(x-x_i)^2 + c_i*(x-x_i) + d_i
//
//  Donde S_i'' (x_i) = M_i  (momentos, lo que aquí se resuelve)
//  Condición natural: M_0 = M_{n-1} = 0
//
//  Sistema tridiagonal (nodos interiores i = 1..n-2):
//    h_{i-1}*M_{i-1} + 2*(h_{i-1}+h_i)*M_i + h_i*M_{i+1} = 6*(f_i - f_{i-1})
//
//  Coeficientes:
//    a_i = (M_{i+1} - M_i) / (6*h_i)
//    b_i = M_i / 2
//    c_i = f_i - h_i*(M_{i+1} + 2*M_i) / 6
//    d_i = y_i
//
void calcularSpline(const vector<vector<double>>& mat) {
    int n = mat.size();

    for (int i = 0; i < n - 1; i++) {
        if (mat[i+1][0] <= mat[i][0]) {
            cout << "Error: los valores de x deben ser estrictamente crecientes.\n";
            cout << "  x[" << i << "] = " << mat[i][0]
                 << "  x[" << i+1 << "] = " << mat[i+1][0] << "\n";
            return;
        }
    }

    //Paso 1: h_i y diferencias divididas f_i ──
    vector<double> h(n-1), f(n-1);
    for (int i = 0; i < n-1; i++) {
        h[i] = mat[i+1][0] - mat[i][0];
        f[i] = (mat[i+1][1] - mat[i][1]) / h[i];
    }

    // aso 2: Sistema tridiagonal para momentos interiores ──
    //  m = número de momentos interiores = n-2
    //  fila i del sistema corresponde al nodo interior i+1
    int m = n - 2;

    if (m == 0) {
        // Solo 2 puntos: spline lineal (los momentos son todos 0)
        cout << "\nSolo 2 puntos: el spline es lineal.\n";
        cout << "S1(x) = " << f[0] << "*(x - " << mat[0][0] << ") + " << mat[0][1] << "\n";
        return;
    }

    Mat A(m, vector<double>(m, 0));
    vector<double> B(m);

    for (int i = 0; i < m; i++) {
        // Nodo interior real: k = i+1  (i va de 0 a m-1)
        // h[i]   = h_{k-1}  (intervalo izquierdo del nodo k)
        // h[i+1] = h_{k}    (intervalo derecho  del nodo k)
        A[i][i] = 2.0 * (h[i] + h[i+1]);          // diagonal principal
        if (i > 0)   A[i][i-1] = h[i];             // sub-diagonal: h_{k-1}
        if (i < m-1) A[i][i+1] = h[i+1];           // super-diagonal: h_{k}
        B[i] = 6.0 * (f[i+1] - f[i]);              // lado derecho
    }



    // ── Paso 3: Resolver A*M = B ──
    Mat B_mat(m, vector<double>(1));
    for (int i = 0; i < m; i++) B_mat[i][0] = B[i];

    Mat A_inv = inversa(A);
    if (A_inv.empty()) {
        cout << "Error: La matriz A no es invertible.\n";
        return;
    }

    Mat X = mult(A_inv, B_mat);

    // ── Paso 4: Armar vector completo de momentos M (con frontera natural) ──
    vector<double> M(n, 0.0);   // M[0] = M[n-1] = 0  (frontera natural)
    for (int i = 0; i < m; i++) M[i+1] = X[i][0];

    // ── Paso 5: Coeficientes de cada tramo ──
    vector<double> a(n-1), b(n-1), c(n-1), d(n-1);
    for (int i = 0; i < n-1; i++) {
        a[i] = (M[i+1] - M[i]) / (6.0 * h[i]);
        b[i] = M[i] / 2.0;
        c[i] = f[i] - h[i] * (M[i+1] + 2.0 * M[i]) / 6.0;
        d[i] = mat[i][1];
    }

    // ── Tabla compacta ──
    cout << fixed << setprecision(5);
    int W = 11; // ancho de cada columna

    cout << "\n=== TABLA SPLINE ===\n";
    cout << setw(4)  << "i"
         << setw(W)  << "xi"
         << setw(W)  << "yi"
         << setw(W)  << "hi"
         << setw(W)  << "fi"
         << setw(W)  << "Si"
         << setw(W)  << "ai"
         << setw(W)  << "bi"
         << setw(W)  << "ci"
         << setw(W)  << "di"
         << "\n";
    cout << string(4 + W*9, '-') << "\n";

    for (int i = 0; i < n-1; i++) {
        cout << setw(4)  << i
             << setw(W)  << mat[i][0]
             << setw(W)  << mat[i][1]
             << setw(W)  << h[i]
             << setw(W)  << f[i]
             << setw(W)  << M[i]
             << setw(W)  << a[i]
             << setw(W)  << b[i]
             << setw(W)  << c[i]
             << setw(W)  << d[i]
             << "\n";
    }
    cout << string(4 + W*9, '-') << "\n";

    // ── Polinomios por tramo ──
    cout << "\n=== POLINOMIOS SPLINE ===\n";
    for (int i = 0; i < n-1; i++) {
        cout << "S" << i << "(x) = "
             << a[i] << "*(x - " << mat[i][0] << ")^3 + "
             << b[i] << "*(x - " << mat[i][0] << ")^2 + "
             << c[i] << "*(x - " << mat[i][0] << ") + "
             << d[i]
             << "   [" << mat[i][0] << ", " << mat[i+1][0] << "]\n";
    }

    // ── Paso 6 (opcional): Evaluar el spline en un punto ──
    cout << "\n¿Desea evaluar el spline en un punto? (s/n): ";
    char ev; cin >> ev;
    while (ev == 's' || ev == 'S') {
        cout << "  Ingrese x: ";
        double xq; cin >> xq;
        if (cin.fail()) { cin.clear(); cin.ignore(10000,'\n'); break; }

        // Buscar tramo
        int seg = -1;
        for (int i = 0; i < n-1; i++) {
            if (xq >= mat[i][0] - 1e-12 && xq <= mat[i+1][0] + 1e-12) {
                seg = i; break;
            }
        }
        if (seg == -1) {
            cout << "  x = " << xq << " está fuera del rango ["
                 << mat[0][0] << ", " << mat[n-1][0] << "].\n";
        } else {
            double dx = xq - mat[seg][0];
            double val = ((a[seg]*dx + b[seg])*dx + c[seg])*dx + d[seg];
            cout << "  S" << seg+1 << "(" << xq << ") = " << val << "\n";
        }
        cout << "¿Evaluar en otro punto? (s/n): ";
        cin >> ev;
    }
}

//  Menú principal
void splineCubico() {
    char continuar = 's';
    while (continuar == 's' || continuar == 'S') {
        cout << "       SPLINE CUBICO NATURAL\n";
        cout << "¿Como desea ingresar los datos?\n";
        cout << "  1. Manualmente\n";
        cout << "  2. Desde archivo CSV\n";
        cout << "Opcion: ";
        int opcion; cin >> opcion;
        if (cin.fail()) {
            cin.clear(); cin.ignore(10000, '\n');
            cout << "Opción inválida.\n"; continue;
        }

        vector<vector<double>> puntos;
        bool ok = false;

        if (opcion == 1) {
            ok = leerManual(puntos);
        } else if (opcion == 2) {
            cout << "Ingrese la ruta del archivo CSV: ";
            string ruta; cin >> ruta;
            ok = leerDesdeCSV(ruta, puntos);
            if (ok) {
                cout << "\nPuntos leídos del CSV:\n";
                for (int i = 0; i < (int)puntos.size(); i++)
                    cout << "  " << i << ": ("
                         << puntos[i][0] << ", " << puntos[i][1] << ")\n";
            }
        } else {
            cout << "Opción inválida.\n"; continue;
        }

        if (ok) calcularSpline(puntos);

        cout << "\n¿Desea ingresar otro conjunto de puntos? (s/n): ";
        cin >> continuar;
    }
}