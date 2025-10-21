#include <iostream>
#include <vector>

using namespace std;

// FunciÃ³n para obtener la submatriz eliminando fila `row` y columna `col`
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

// FunciÃ³n recursiva para calcular el determinante
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

    //MÃ©todo
    if (determinant(mat)==0){
        cout << "No es posible aplicar el metodo...\n";
    }

    //Vector de tÃ©rminos independientes
    cout << "\nIngrese el vector de terminos independientes b (" << n << "x1):\n";
    for (i = 0; i < n; i++) {
        cout << "b[" << i+1 << "]: ";
        cin >> b[i];
    }
    cout << "\nVector b:\n";
    for (i = 0; i < n; i++) {
        cout << b[i] << endl;
    }




    //FÃ³rmulas para la TriangularizaciÃ³n del sistema de ecuaciones de la matriz aumentada
    /*
    for ( i = 0; i < n-1; i++){
        for (j = i+1; j < n; j++){
            b[j]=b[j]-(mat[j][i]/mat[i][i])*b[j];
        }
        for (k = n; k < i; k++){
            if (k==i){
                mat[j][k]=0;
                if (j!=i){
                    mat[j][k]=mat[j][k]-(mat[i][k]*mat[j][i])/mat[i][i];
                }
                
            }
        }

    }
    */


    // --- EliminaciÃ³n Gaussiana ---
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            double factor = mat[i][k] / mat[k][k];
            for (j = k; j < n; j++) {
                mat[i][j] -= factor * mat[k][j];
            }
            b[i] -= factor * b[k];
        }
    }


    //Revisar nuevas matrices
    cout << "\n\t--Nuevas MAtrices--\n";

    cout << "\nLa matriz transformada:\n"; 
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            //Para evitar errores de redondeo:
            if (abs(mat[i][j]) < 1e-10){  // Si es muy cercano a 0
            cout << 0 << "\t";
            }
            else{
            cout << mat[i][j] << "\t";
            }
            //cout << mat[i][j] << "\t";
        }
        cout << endl;
    }


    cout << "\nVector b (transformado):\n";
    for (i = 0; i < n; i++) {
        cout << b[i] << endl;
    }

    // --- SustituciÃ³n regresiva ---
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


}
void recursiva(vector<vector<double>>& mat, vector<double>&ind,int k, int n){
    if(k==n){
        cout<<"\nLos valores para x son: ";
        for(int i=0; i<n; i++){
            cout<<"\nx"<<i+1<<" = "<<ind[i]<<" ";
        }
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
        return ;
    }

    cout<<"\nIngrese los elementos de su vector de terminos independientes: \n";
    for(int i=0; i<n; i++){
        cin>>ind[i];
    }

    for(int i=0; i<n; i++){
        if(mat[i][i] == 0.0){
            cout<<"La matriz contiene un 0 en la diagonal\n";
            return ;
        }
    }
    int k=0;
    recursiva(mat,ind,k, n);
    cout<<'\n';


    return ;
}