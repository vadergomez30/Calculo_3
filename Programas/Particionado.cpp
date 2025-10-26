#include <iostream>
#include <vector>
#include "Ope_Matrices.h"

using namespace std;

vector<vector<double>> getD(vector<vector<double>> mat,int P,vector<vector<double>> A22inv){
    vector<vector<double>>A12(P,vector<double>(3));
    for(int i=0; i<P; i++){
        for(int j=0; j<3; j++){
            A12[i][j]=mat[i][j+P];
        }
    }
    vector<vector<double>>A22(3,vector<double>(3));
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            A22[i][j]=mat[i+P][j+P];
        }
    }
    return mMUlt(A12, A22inv);
}
vector<vector<double>> getE(vector<vector<double>> mat,int P,vector<vector<double>> A22inv){
    vector<vector<double>>A21(3,vector<double>(P));
    for(int i=0; i<3; i++){
        for(int j=0; j<P; j++){
            A21[i][j]=mat[i+P][j];
        }
    }
    return mMUlt(A22inv,A21);
}
vector<vector<double>> getC(vector<vector<double>>&mat,vector<vector<double>> D,int P){
    vector<vector<double>>A21(3,vector<double>(P));
    for(int i=0; i<3; i++){
        for(int j=0; j<P; j++){
            A21[i][j]=mat[i+P][j];
        }
    }
    vector<vector<double>>A11(P,vector<double>(P));
    for(int i=0; i<P; i++){
        for(int j=0; j<P; j++){
            A11[i][j]=mat[i][j];
        }
    }

    vector<vector<double>> M = mMUlt(D, A21);

    vector<vector<double>> C(P, vector<double>(P, 0.0));
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < P; ++j) {
            C[i][j] = A11[i][j] - M[i][j];
        }
    }
    return C;
}
vector<vector<double>> getF(vector<vector<double>> mat,vector<vector<double>> D,vector<vector<double>> C,vector<vector<double>> E,vector<vector<double>>A22inv ,int &P){
    vector<vector<double>>idenC=identidad(P);
    vector<vector<double>>ECinv=mMUlt(E,inversa(C,idenC,0,C.size()));
    vector<vector<double>>ECinvD=mMUlt(ECinv,D);

    vector<vector<double>>F(3,(vector<double>(3)));

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = A22inv[i][j] + ECinvD[i][j];
        }
    }
    return F;

}
void resultado(vector<vector<double>>D,vector<vector<double>>E,vector<vector<double>>C,vector<vector<double>>F,int n, int P){
    vector<vector<double>>idenC=identidad(P);
    vector<vector<double>>Cinv=inversa(C,idenC,0,C.size());
    vector<vector<double>>CinvD=mMUlt(Cinv,D);
    vector<vector<double>>ECinv=mMUlt(E,Cinv);

    vector<vector<double>>res(n,vector<double>(n));
    for(int i=0; i<P;i++){
        for (int j=0; j<P; j++){
            res[i][j]=Cinv[i][j];
        }
    }
    for(int i=0; i<P;i++){
        for (int j=0; j<3; j++){
            res[i][j+P]=CinvD[i][j]*-1;
        }
    }
    for(int i=0; i<3;i++){
        for (int j=0; j<P; j++){
            res[i+P][j]=ECinv[i][j]*-1;
        }
    }
    for(int i=0; i<3;i++){
        for (int j=0; j<3; j++){
            res[i+P][j+P]=F[i][j];
        }
    }
    cout<<"\nMatriz inversa\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (abs(res[i][j]) < 1e-10){  // Si es muy cercano a 0
                cout << 0 << "\t";
            }else{
                cout << res[i][j] << " ";
            }
        }
        cout << '\n';
    }
}


int main() {
    cout << "Metodo de Particionado de Matrices\n";
    cout<<"\nIngrese el tamano de la matriz: ";
    int n, P; 
    cin>>n;
    P=n-3;
    if (n<4){
        cout<<"\nLa cantidad de variables no son suficientes para aplicar este metodo.\n";
        return 0;
    }
    vector<vector<double>> mat(n, vector<double>(n));
    vector<double>ind(n);
    cout<<"\nIngrese los elementos de su matriz de coeficientes: \n";
    for(int i=0; i<n; i++){
        cout<< "\nFila "<<i+1<<": \n";
        for(int j=0; j<n; j++){
            cin>>mat[i][j];
        }
    }
    if (determinant(mat) == 0){
        cout<<"\nLa matriz no tiene solucion unica, no se puede aplicar el metodo.\n";
        return 0;
    }
    cout<<"\nIngrese los elementos de su vector independiente: \n";
    for(int i=0; i<n; i++){
        cout<<"Independiente "<<i+1<<": ";
        cin>>ind[i];
    }
    for(int i=0; i<n; i++){
        if(mat[i][i] == 0.0){
            cout<<"\nLa matriz contiene un 0 en la diagonal\n";
            return 0;
        }
    }
    vector<vector<double>>iden=identidad(3);
    vector<vector<double>>A22(3,vector<double>(3));
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            A22[i][j]=mat[i+P][j+P];
        }
    }
    if (determinant(A22) == 0){
        cout<<"\nLa matriz no tiene inversa, no se puede aplicar el metodo.\n";
        return 0;
    }
    vector<vector<double>>A22inv=inversa(A22,iden,0,3);

    vector<vector<double>>D= getD(mat,P,A22inv);
    vector<vector<double>>E= getE(mat,P,A22inv);
    vector<vector<double>>C= getC(mat,D,P);
    vector<vector<double>>F= getF(mat,D,C,E,A22inv,P);
    resultado(D,E,C,F,n,P);
    
}