#include <bits/stdc++.h>
using namespace std;
#include "Ope_Matrices.h"
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