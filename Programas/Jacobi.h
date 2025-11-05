#include <bits/stdc++.h>
using namespace std; 
using Mat = vector<vector<double>>;

bool supuestos(const Mat& A){
    int n=A.size();
    for(int i=0;i<n;i++){
        double diag= fabs(A[i][i]);
        double suma=0.0;
        for(int j=0;j<n;j++){
            if(j==i) continue;
            suma+=fabs(A[i][j]);
        }
        if(diag < suma) return false;
    }
    return (determinant(A)==0)? false: true;
}

bool detener(const vector<double> &x_ant, const vector<double> &x_nue, const double &error){
    for(size_t i=0; i<x_ant.size(); i++){
        if(fabs(x_nue[i]-x_ant[i])>error){
            return false;
        }
    }
    return true;
}

void iterativoJacobi(const Mat&A, const vector<double>&b, vector<double>&x_ant, vector<double>&x_nue, const double error, int &count){
    for(size_t i=0; i<A.size(); i++){
        double res=b[i];
        for (size_t j=0; j<A[0].size(); j++){
            if(j==i) continue;
            res-=A[i][j]*x_ant[j];
        }
        x_nue[i]=1/A[i][i] * res;
    }
    cout<<"\nIteracion "<<++count<<endl;
    cout<<"Solucion X: ";
    for(auto i: x_nue){
        cout<<i<<" ";
    } 
    if(detener(x_ant, x_nue, error)){
        cout<<"\nConvergencia alcanzada despues de "<<count<<" iteraciones."<<endl;
        return;
    } else {
        x_ant=x_nue;
        iterativoJacobi(A, b, x_ant, x_nue, error, count);
    }
}

void Jacobi(){
    int n, count=0;
    double error;
    cout << "Ingrese el tamano n de la matriz (n x n): "<<endl;
    cin >> n;
    cout << "Ingrese el error tolerable: "<<endl;
    cin >> error;
    Mat A(n, vector<double>(n,0.0));
    cout << "Ingrese la matriz A ("<<n<<"x"<<n<<"), por filas:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cin>>A[i][j];
        }
    }
    vector<double> b(n,0.0);
    cout << "Ingrese el vector b ("<<n<<" valores):\n";
    for(int i=0;i<n;i++){
        cin>>b[i];
    }
    if(!supuestos(A)){
        cout<<"La matriz no cumple los supuestos necesarios (dominancia diagonal y no singularidad).\n";
        return;
    }
    vector<double> x_ant(n,0.0);
    vector<double> x_nue(n,0.0);
    iterativoJacobi(A, b, x_ant, x_nue, error, count);
}