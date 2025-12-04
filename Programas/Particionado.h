#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>  
using namespace std;

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

vector<vector<double>> inversa(vector<vector<double>>& mat, vector<vector<double>>& iden, int k, int n) {
    if (k >= n) return iden;
    double factor = mat[k][k];
    for (int j = 0; j < n; j++) {
        mat[k][j] /= factor;
        iden[k][j] /= factor;
    }

    for (int i = 0; i < n; i++) {
        if (i == k) continue;
        factor = mat[i][k];
        for (int j = 0; j < n; j++) {
            mat[i][j] -= factor * mat[k][j];
            iden[i][j] -= factor * iden[k][j];
        }
    }

    return inversa(mat, iden, k + 1, n);
}

vector<vector<double>> mMUlt(const vector<vector<double>>& A, const vector<vector<double>>& B){
    int n=A.size();
    int m=B[0].size();
    vector<vector<double>> C(n,vector<double>(m,0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            for (int k = 0; k < (int)B.size(); k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<vector<double>> getD(vector<vector<double>> mat,int P,vector<vector<double>> A22inv){
    int k=mat.size()-P;
    vector<vector<double>>A12(P,vector<double>(k));
    for(int i=0; i<P; i++){
        for(int j=0; j<k; j++){
            A12[i][j]=mat[i][j+P];
        }
    }
    return mMUlt(A12, A22inv);
}

vector<vector<double>> getE(vector<vector<double>> mat,int P,vector<vector<double>> A22inv){
    int k=mat.size()-P;
    vector<vector<double>>A21(k,vector<double>(P));
    for(int i=0; i<k; i++){
        for(int j=0; j<P; j++){
            A21[i][j]=mat[i+P][j];
        }
    }
    return mMUlt(A22inv,A21);
}

vector<vector<double>> getC(vector<vector<double>>&mat,vector<vector<double>> D,int P){
    int k=mat.size()-P;
    vector<vector<double>>A21(k,vector<double>(P));
    for(int i=0; i<k; i++){
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
    int k=mat.size()-P;
    vector<vector<double>>idenC=identidad(P);
    vector<vector<double>>ECinv=mMUlt(E,inversa(C,idenC,0,C.size()));
    vector<vector<double>>ECinvD=mMUlt(ECinv,D);

    vector<vector<double>>F(k,(vector<double>(k,0)));

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            F[i][j] = A22inv[i][j] + ECinvD[i][j];
        }
    }
    return F;

}

void resultado(vector<vector<double>>D,vector<vector<double>>E,vector<vector<double>>C,vector<vector<double>>F,vector<double>ind,int n, int P){
    int k=n-P;
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
        for (int j=0; j<k; j++){
            res[i][j+P]=CinvD[i][j]*-1;
        }
    }
    for(int i=0; i<k;i++){
        for (int j=0; j<P; j++){
            res[i+P][j]=ECinv[i][j]*-1;
        }
    }
    for(int i=0; i<k;i++){
        for (int j=0; j<k; j++){
            res[i+P][j+P]=F[i][j];
        }
    }
    cout<<"\nMatriz inversa\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (abs(res[i][j]) < 1e-10){  
                cout << 0 << "\t";
            }else{
                cout << res[i][j] << " ";
            }
        }
        cout << '\n';
    }

    vector<double> solucion(n,0);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            solucion[i]+=ind[j]*res[i][j];
        }
    }
    cout<<"\nLa solucion del sistema de ecuaciones es: \n";
    for(int i=0; i<n; i++){
        cout<<"x"<<i+1<<" = "<<solucion[i]<<'\n';
    }
    cout<<'\n';
}

void printMatrix(const vector<vector<double>>& M, const string& name = "") {
    if (!name.empty()) cout << "\n" << name << ":\n";
    if (M.empty()) { cout << "(Vacia)\n"; return; }
    cout.setf(ios::fixed);
    cout << setprecision(6);
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M[i].size(); ++j) {
            double v = (std::fabs(M[i][j]) < 1e-10) ? 0.0 : M[i][j];
            cout << v << "\t";
        }
        cout << '\n';
    }
    cout << '\n';
}
void Particionado() {
    cout << "Metodo de Particionado de Matrices\n";
    cout<<"\nIngrese el tamano de la matriz: ";
    int n, P; 
    cin>>n;

    if (!cin || n <= 0) {
        cout << "\nTamano de matriz invalido.\n";
        system("pause");
        return;
    }

    cout<<"\nIngrese el tamaño de la particion (2<= P <= "<<n-2<<"):";
    cin>>P;

    if (!cin || P < 2 || P > n-2) {
        cout << "\nLa particion P no es valida.\n";
        system("pause");
        return;
    }

    int k=n-P;
    if (n<4){
        cout<<"\nLa cantidad de variables no son suficientes para aplicar este metodo.\n";
        system("pause");
        return ;
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
        system("pause");
        return ;
    }

    cout<<"\nIngrese los elementos de su vector independiente: \n";
    for(int i=0; i<n; i++){
        cout<<"Independiente "<<i+1<<": ";
        cin>>ind[i];
    }

    for(int i=0; i<n; i++){
        if(mat[i][i] == 0.0){
            cout<<"\nLa matriz contiene un 0 en la diagonal\n";
            system("pause");
            return ;
        }
    }

    vector<vector<double>>iden=identidad(k);
    vector<vector<double>>A22(k,vector<double>(k));
    for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
            A22[i][j]=mat[i+P][j+P];
        }
    }

    if (determinant(A22) == 0){
        cout<<"\nLa submatriz A22 no tiene inversa, no se puede aplicar el metodo.\n";
        system("pause");
        return ;
    }

    vector<vector<double>>A22inv = inversa(A22, iden, 0, k); 
    if (A22inv.empty()) {
        cout << "\nError: no se pudo invertir A22 (posible singularidad)\n";
        system("pause");
        return;
    }

    printMatrix(A22inv, "A22 inversa");

    vector<vector<double>>D= getD(mat,P,A22inv);
    printMatrix(D, "D");

    vector<vector<double>>E= getE(mat,P,A22inv);
    printMatrix(E, "E");

    vector<vector<double>>C= getC(mat,D,P);
    printMatrix(C, "C");

    vector<vector<double>>F= getF(mat,D,C,E,A22inv,P);
    printMatrix(F, "F");

    resultado(D,E,C,F,ind,n,P);

    system("pause");
    return;
}
