#include <iostream>
#include <vector>
#include <cmath>    
#include <iomanip>  // setw, setprecision
#include <cstdlib>  // system("pause")
using namespace std;

vector<vector<double>> getSubmatrix(const vector<vector<double>>& mat, int row, int col) {
    int n = mat.size();
    vector<vector<double>> submat(n - 1, vector<double>(n - 1));
    int subi = 0;
    for (int i = 0; i < n; i++) {
        if (i == row) continue;
        int subj = 0;
        for (int j = 0; j < n; j++) {
            if (j == col) continue;
            submat[subi][subj] = mat[i][j];
            subj++;
        }
        subi++;
    }
    return submat;
}

double determinant(const vector<vector<double>>& mat) {
    int n = mat.size();
    
    if (n == 1) return mat[0][0];
    if (n == 2) return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]; // caso base 2x2

    double det = 0;
    for (int j = 0; j < n; j++) {
        vector<vector<double>> submat = getSubmatrix(mat, 0, j);
        det += ( (j % 2 == 0 ? 1 : -1) * mat[0][j] * determinant(submat) );
    }
    return det;
}

void Gauss() {
    int n,i,j,k;
    cout << "Ingrese el tamano de la matriz cuadrada (n x n): ";
    cin >> n;

    vector<vector<double>> mat(n, vector<double>(n)); //Matriz
    vector<double> b(n); // vector de resultados
    vector<double> Res(n); //Vector de respuestas

    cout << "Ingrese los elementos de la matriz fila por fila:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "Elemento [" << i+1 << "][" << j+1 << "]: ";
            cin >> mat[i][j];
        }
    }
    //Mostrar Matriz
    cout << "\nLa matriz ingresada es:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "El determinante es: " << determinant(mat) << endl;

    //Método
    if (determinant(mat)==0){
        cout << "No es posible aplicar el metodo...\n";
        system("pause");
        return;
    }

    //Vector de términos independientes
    cout << "\nIngrese el vector de terminos independientes b (" << n << "x1):\n";
    for (i = 0; i < n; i++) {
        cout << "b[" << i+1 << "]: ";
        cin >> b[i];
    }
    cout << "\nVector b:\n";
    for (i = 0; i < n; i++) {
        cout << b[i] << endl;
    }

    // --- Eliminación Gaussiana ---
    for (k = 0; k < n - 1; k++) {
        double piv = mat[k][k];
        if (fabs(piv) < 1e-12) {
            cout << "\nPivote casi cero en la fila " << k+1 << ". No se puede continuar.\n";
            system("pause");
            return;
        }
        for (i = k + 1; i < n; i++) {
            double factor = mat[i][k] / piv;
            for (j = k; j < n; j++) {
                mat[i][j] -= factor * mat[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    //Revisar nuevas matrices
    cout << "\n\t--Nuevas Matrices--\n";

    cout << "\nLa matriz transformada:\n"; 
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            //Para evitar errores de redondeo:
            if (fabs(mat[i][j]) < 1e-10){  // Si es muy cercano a 0
                cout << 0 << "\t";
            }
            else{
                cout << mat[i][j] << "\t";
            }
        }
        cout << endl;
    }

    cout << "\nVector b (transformado):\n";
    for (i = 0; i < n; i++) {
        cout << b[i] << endl;
    }

    // --- Sustitución regresiva ---
    for (i = n - 1; i >= 0; i--) {
        double suma = 0;
        for (j = i + 1; j < n; j++) {
            suma += mat[i][j] * Res[j];
        }
        Res[i] = (b[i] - suma) / mat[i][i];
    }
    cout << "\nSoluciones del sistema:\n";
    for (i = 0; i < n; i++) {
        cout << "x[" << i+1 << "] = " << Res[i] << endl;
    }

    system("pause");
}

void recursiva(vector<vector<double>>& mat, vector<double>&ind,int k, int n){ //Crea la matriz identidad mediante gauss-jordan
    if(k==n){
        cout<<"\nLos valores para x son: ";
        for(int i=0; i<n; i++){
            cout<<"\nx"<<i+1<<" = "<<ind[i]<<" ";
        }
        cout << "\n";
        return;
    }
    cout<<"\niteracion: "<<k+1<<'\n';
    double factor;
    factor = mat[k][k];
    for(int i=k; i<n; i++){
        mat[k][i]/=factor;
    }
    ind[k]/=factor;
    for (int i = k + 1; i < n; ++i) {
        factor = mat[i][k];
        for (int j = k; j < n; ++j) {
            mat[i][j] -= factor * mat[k][j];
        }
        ind[i] -= factor * ind[k];
    }
    for (int p = k ; p >= 0; p--) { 
        for (int i = p - 1; i >= 0; i--) {
            factor = mat[i][p];
            for (int j = 0; j < n; j++) {
                mat[i][j] -= factor * mat[p][j];
            }
            ind[i] -= factor * ind[p];
        }
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<mat[i][j]<<" ";
        }
        cout<<'\n';
    }
    cout<<'\n';
    for(int i=0; i<n; i++){
        cout<<ind[i]<<" ";
    }
    cout<<'\n';
    recursiva(mat,ind,k+1, n);
}


void GaussJordan() {
    cout<<"Metodo de Gauss-Jordan\n";
    cout<<"\nIngresa el tamaño de tu matriz: \n";
    int n; cin>>n;
    vector<vector<double>> mat(n, vector<double>(n));
    vector<double>ind(n);
    cout<<"\nIngrese los elementos de su matriz de coeficientes: \n";
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cin>>mat[i][j];
        }
    }   
    
    if(determinant(mat) == 0){
        cout<<"La matriz no tiene solucion unica\n";
        system("pause");
        return ;
    }

    cout<<"\nIngrese los elementos de su vector de terminos independientes: \n";
    for(int i=0; i<n; i++){
        cin>>ind[i];
    }

    for(int i=0; i<n; i++){
        if(mat[i][i] == 0.0){
            cout<<"La matriz contiene un 0 en la diagonal\n";
            system("pause");
            return ;
        }
    }
    int k=0;
    recursiva(mat,ind,k, n);
    cout<<'\n';

    system("pause");
    return ;
}

vector<vector<double>> identidad2(int n) {
    vector<vector<double>>pos(n,vector<double>(n));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j)pos[i][j]=1;
            else pos[i][j]=0;
        }
    }
    return pos;
}

void resultado(vector<vector<double>>& iden, vector<double>& ind, int n){
    vector<double> solucion(n,0);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            solucion[i]+=ind[j]*iden[i][j];
        }
    }
    cout<<"\nLa solucion del sistema de ecuaciones es: \n";
    for(int i=0; i<n; i++){
        cout<<"x"<<i+1<<" = "<<solucion[i]<<'\n';
    }
    cout<<'\n';
}

void recursiva2(vector<vector<double>>& mat, vector<vector<double>>& iden, int k, vector<double>& ind, int n){
    if(k == n){
        resultado(iden,ind,n);
        return;
    }

    cout << "\nIteracion: " << k+1 << '\n';

    double factor = mat[k][k];

    for(int j = 0; j < n; j++) {
        mat[k][j] /= factor;
        iden[k][j] /= factor;
    }

    for(int i = 0; i < n; i++) {
        if(i == k) continue;
        factor = mat[i][k];
        for(int j = 0; j < n; j++) {
            mat[i][j] -= factor * mat[k][j];
            iden[i][j] -= factor * iden[k][j];
        }
    }
    cout<<"\nMatriz de coeficientes\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << mat[i][j] << " ";
        }
        cout<<'\n';
    }
    cout<<"\nMatriz inversa\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << iden[i][j] << " ";
        }
        cout << '\n';
    }

    recursiva2(mat, iden, k+1, ind, n);
}

void Inversa() {
    cout<<"Metodo de Inversion de Matrices\n";
    cout<<"\nIngrese el tamano de la matriz: ";
    int n; cin>>n;
    vector<vector<double>> mat(n, vector<double>(n));
    vector<vector<double>> iden = identidad2(n);
    vector<double>ind(n);
    cout<<"\nIngrese los elementos de su matriz de coeficientes: \n";
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cin>>mat[i][j];
        }
    }   
    
    if(determinant(mat) == 0){
        cout<<"\nLa matriz no tiene solucion unica\n";
        system("pause");
        return ;
    }

    cout<<"\nIngrese los elementos de su vector de terminos independientes: \n";
    for(int i=0; i<n; i++){
        cin>>ind[i];
    }

    for(int i=0; i<n; i++){
        if(mat[i][i] == 0.0){
            cout<<"\nLa matriz contiene un 0 en la diagonal\n";
            system("pause");
            return ;
        }
    }
    int k=0;
    recursiva2(mat,iden,k,ind,n);
    cout<<'\n';

    system("pause");
    return ;
}

void M_relajacion(){
    int n5;
    cout << "Tamano de la matriz: ";
    cin >> n5;

    vector<vector<double>> A(n5, vector<double>(n5));
    vector<double> b(n5);

    for (int i = 0; i < n5; i++) {
        for (int j = 0; j < n5; j++) {
            cout << "A[" << i+1 << "][" << j+1 << "]: ";
            cin >> A[i][j];
        }
    }

    for (int i = 0; i < n5; i++) {
        cout << "b[" << i+1 << "]: ";
        cin >> b[i];
    }

    bool dd = true;
    cout << "\nDominancia diagonal:\n";
    for (int i = 0; i < n5; i++) {
        double suma = 0.0;
        for (int j = 0; j < n5; j++) {
            if (i != j) suma += fabs(A[i][j]);
        }
        bool fila = fabs(A[i][i]) >= suma;
        if (!fila) dd = false;
        cout << "Fila " << i+1 << ": " << (fila ? "verdadero" : "falso") << "\n";
    }

    double determinante = determinant(A);
    cout << "\nDeterminante de A: " << determinante << "\n";

    if (!dd) {
        cout << "\nLa matriz no cumple dominancia diagonal. No se puede realizar el metodo.\n";
        system("pause");
        return;
    }

    double tol;
    cout << "Tolerancia: ";
    cin >> tol;

    vector<vector<double>> C(n5, vector<double>(n5));
    vector<double> D(n5);
    for (int i = 0; i < n5; i++) {
        double piv = A[i][i];
        for (int j = 0; j < n5; j++) {
            if (i == j) C[i][j] = -1.0;
            else C[i][j] = -A[i][j] / piv;
        }
        D[i] = b[i] / piv;
    }

    cout << fixed << setprecision(8);
    cout << "\nMatriz C:\n";
    for (int i = 0; i < n5; i++) {
        for (int j = 0; j < n5; j++) cout << setw(12) << C[i][j];
        cout << "\n";
    }

    cout << "\nVector D:\n";
    for (int i = 0; i < n5; i++) cout << setw(12) << D[i];
    cout << "\n";

    vector<double> X(n5, 0.0), CX(n5, 0.0), R(n5, 0.0), Xnew(n5,0.0);
    bool detener = false;
    int iter = 0;
    int iter_max = 20000;
    const double eqeps = 1e-12;

    while (!detener && iter < iter_max) {
        cout << "\nIteracion " << iter << ":\n\n";

        cout << "Vector X:\n";
        for (int i = 0; i < n5; i++) cout << setw(12) << X[i];
        cout << "\n\n";

        for (int i = 0; i < n5; i++) {
            double s = 0.0;
            for (int j = 0; j < n5; j++) s += C[i][j] * X[j];
            CX[i] = s;
        }

        cout << "Vector CX:\n";
        for (int i = 0; i < n5; i++) cout << setw(12) << CX[i];
        cout << "\n\n";

        for (int i = 0; i < n5; i++) R[i] = CX[i] + D[i];

        cout << "Vector R:\n";
        for (int i = 0; i < n5; i++) cout << setw(12) << R[i];
        cout << "\n\n";

        double maxv = 0.0;
        for (int i = 0; i < n5; i++) if (fabs(R[i]) > maxv) maxv = fabs(R[i]);
        cout << "MAXIMO: " << maxv << "\n\n";

        if (maxv < tol) {
            detener = true;
            cout << "DETENER: SI\n";
        } else {
            cout << "DETENER: NO\n";
        }

        for (int i = 0; i < n5; i++) {
            if (fabs(fabs(R[i]) - maxv) < eqeps) Xnew[i] = X[i] + R[i];
            else Xnew[i] = X[i];
        }

        X = Xnew;
        iter++;
    }

    if (iter >= iter_max) cout << "\nEl metodo supero el numero maximo de iteraciones\n";

    cout << "\nSolucion final X:\n";
    for (int i = 0; i < n5; i++) cout << setw(12) << X[i];
    cout << "\n";

    system("pause");
    return;
}
