#include <iostream>
#include <vector>
#include "Ope_Matrices.h"

using namespace std;

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
void recursiva(vector<vector<double>>& mat, vector<double>&ind,int k, int n){ //Crea la matriz identidad mediante gauss-jordan
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
vector<vector<double>> identidad(int n) {
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
    vector<vector<double>> iden = identidad(n);
    vector<double>ind(n);
    cout<<"\nIngrese los elementos de su matriz de coeficientes: \n";
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cin>>mat[i][j];
        }
    }   
    
    if(determinant(mat) == 0){
        cout<<"\nLa matriz no tiene solucion unica\n";
        return ;
    }

    cout<<"\nIngrese los elementos de su vector de terminos independientes: \n";
    for(int i=0; i<n; i++){
        cin>>ind[i];
    }

    for(int i=0; i<n; i++){
        if(mat[i][i] == 0.0){
            cout<<"\nLa matriz contiene un 0 en la diagonal\n";
            return ;
        }
    }
    int k=0;
    recursiva2(mat,iden,k,ind,n);
    cout<<'\n';


    return ;
}