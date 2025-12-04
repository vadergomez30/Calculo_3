#include <bits/stdc++.h>
#include <cstdlib>
using namespace std; 

using Mat = vector<vector<double>>;
static const double EPS = 1e-12;
void printVec(const vector<double>& v, const string& title){
    cout << "\n" << title << ":\n";
    cout.setf(std::ios::fixed); 
    cout << setprecision(8);
    for(double x: v){
        double y = (fabs(x)<1e-12? 0.0: x);
        cout << y << "\n";
    }
}

void printMat(const Mat& A, const string& title){
    cout << "\n" << title << ":\n";
    cout.setf(std::ios::fixed); 
    cout << setprecision(8);
    for(const auto& r: A){
        for(double x: r){
            double y = (fabs(x)<1e-12? 0.0: x);
            cout << y << "\t";
        }
        cout << "\n";
    }
}

Mat mult(const Mat& A, const Mat& B){
    int n = A.size(), p = A[0].size(), m = B[0].size();
    Mat C(n, vector<double>(m,0.0));
    for(int i=0;i<n;i++){
        for(int k=0;k<p;k++){
            double aik = A[i][k];
            if (fabs(aik)<1e-18) continue;
            for(int j=0;j<m;j++){
                C[i][j] += aik * B[k][j];
            }
        }
    }
    return C;
}

vector<double> multMatVec(const Mat& A, const vector<double>& x){
    int n = A.size(), m = A[0].size();
    vector<double> y(n,0.0);
    for(int i=0;i<n;i++){
        double s=0.0;
        for(int j=0;j<m;j++) s += A[i][j]*x[j];
        y[i] = s;
    }
    return y;
}

Mat resta(const Mat& A, const Mat& B){
    int n = A.size(), m = A[0].size();
    Mat C(n, vector<double>(m,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

Mat invGaussJordan(const Mat& A){
    int n = A.size();
    Mat M = A, Inv = identidad2(n);
    for(int col=0; col<n; col++){
        int piv = col;
        double best = fabs(M[piv][col]);
        for(int r=col+1;r<n;r++){
            double v = fabs(M[r][col]);
            if(v > best){ best = v; piv = r; }
        }
        if(best < EPS) 
            throw runtime_error("Matriz singular en inversion (pivote ~ 0).");

        if(piv != col){ 
            swap(M[piv], M[col]); 
            swap(Inv[piv], Inv[col]); 
        }

        double d = M[col][col];
        for(int j=0;j<n;j++){ 
            M[col][j]  /= d; 
            Inv[col][j]/= d; 
        }

        for(int i=0;i<n;i++){
            if(i==col) continue;
            double f = M[i][col];
            if(fabs(f)<EPS) continue;
            for(int j=0;j<n;j++){
                M[i][j]  -= f*M[col][j];
                Inv[i][j]-= f*Inv[col][j];
            }
        }
    }
    return Inv;
}
void partirA(const Mat& A, int p, Mat& A11, Mat& A12, Mat& A21, Mat& A22){
    int n=A.size();
    A11.assign(p,      vector<double>(p));
    A12.assign(p,      vector<double>(n-p));
    A21.assign(n-p,    vector<double>(p));
    A22.assign(n-p,    vector<double>(n-p));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i<p && j<p)       A11[i][j]     = A[i][j];
            else if(i<p && j>=p) A12[i][j-p]   = A[i][j];
            else if(i>=p && j<p) A21[i-p][j]   = A[i][j];
            else                 A22[i-p][j-p] = A[i][j];
        }
    }
}

void partirb(const vector<double>& b, int p, vector<double>& b1, vector<double>& b2){
    int n=b.size();
    b1.assign(p,0.0);
    b2.assign(n-p,0.0);
    for(int i=0;i<p;i++)   b1[i]     = b[i];
    for(int i=p;i<n;i++)   b2[i-p]   = b[i];
}
void GaussParticionado(){
    int n;
    cout << "Hola" << endl;
    cout << "Ingrese n (n>=2): ";
    if(!(cin>>n) || n<2){ 
        cerr<<"n invalido.\n"; 
        system("pause");
        return; 
    }

    Mat A(n, vector<double>(n,0.0));
    cout << "Ingrese la matriz A ("<<n<<"x"<<n<<"), por filas:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(!(cin>>A[i][j])){ 
                cerr<<"Entrada invalida.\n"; 
                system("pause");
                return; 
            }
        }
    }

    vector<double> b(n,0.0);
    cout << "Ingrese el vector b ("<<n<<" valores):\n";
    for(int i=0;i<n;i++){
        if(!(cin>>b[i])){ 
            cerr<<"Entrada invalida.\n"; 
            system("pause");
            return; 
        }
    }

    printMat(A, "Matriz A ingresada");
    printVec(b, "Vector b ingresado");

    double detA = determinant(A);
    cout.setf(std::ios::fixed); 
    cout << setprecision(10);
    cout << "\nDeterminante de A: " << detA << "\n";
    if(fabs(detA)<EPS){
        cout << "El sistema no tiene solucion unica (det(A)=0). Fin.\n";
        system("pause");
        return;
    }

    int p;
    cout << "\nIngrese la particion p (tamano de A11, 1.."<<(n-1)<<"): ";
    if(!(cin>>p) || p<=0 || p>=n){ 
        cerr<<"Particion p invalida.\n"; 
        system("pause");
        return; 
    }
    Mat A11,A12,A21,A22;
    vector<double> b1,b2;
    partirA(A,p,A11,A12,A21,A22);
    partirb(b,p,b1,b2);

    printMat(A11, "A11");
    printMat(A12, "A12");
    printMat(A21, "A21");
    printMat(A22, "A22");
    printVec(b1,  "b1");
    printVec(b2,  "b2");

    try{
        Mat A11_inv = invGaussJordan(A11);
        Mat A12_p   = mult(A11_inv, A12);             
        vector<double> b1_p = multMatVec(A11_inv, b1); 

        printMat(A11_inv, "A11^{-1}");
        printMat(A12_p,   "A12' = A11^{-1} * A12");
        printVec(b1_p,    "b1'  = A11^{-1} * b1");

        Mat A21A12p = mult(A21, A12_p);
        Mat A22_p   = resta(A22, A21A12p);                  
        vector<double> A21b1p = multMatVec(A21, b1_p);
        vector<double> b2_p(b2.size(),0.0);
        for(size_t i=0;i<b2.size();i++) 
            b2_p[i] = b2[i] - A21b1p[i];

        printMat(A22_p, "A22' = A22 - A21 * A12'");
        printVec(b2_p,  "b2'  = b2 - A21 * b1'");
        double detA11  = determinant(A11);
        double detA22p = determinant(A22_p);
        cout << "\nDet(A11): " << detA11 
             << "   Det(A22'): " << detA22p
             << "   Chequeo det(A) ~ det(A11)*det(A22') = " << detA11*detA22p << "\n";

       
        Mat A22_p_inv = invGaussJordan(A22_p);
        vector<double> b2_pp = multMatVec(A22_p_inv, b2_p);

        printMat(A22_p_inv, "(A22')^{-1}");
        printVec(b2_pp,     "b2'' = (A22')^{-1} * b2'");

       
        vector<double> A12pb2pp = multMatVec(A12_p, b2_pp);
        vector<double> b1_pp(b1_p.size(),0.0);
        for(size_t i=0;i<b1_p.size();i++) 
            b1_pp[i] = b1_p[i] - A12pb2pp[i];

        printVec(b1_pp, "b1'' = b1' - A12' * b2''");

       
        vector<double> x(n,0.0);
        for(int i=0;i<p;i++)   x[i]   = b1_pp[i];
        for(int i=p;i<n;i++)   x[i]   = b2_pp[i-p];

        printVec(x, "Solucion x (Gauss-Jordan particionado)");

       
        vector<double> Ax = multMatVec(A, x);
        printVec(Ax, "A*x (debe aproximar b)");

        system("pause");
    }
    catch(const exception& e){
        cerr << "\nError: " << e.what() << "\n";
        cerr << "Revisa tu particion p (A11 debe ser invertible) o el condicionamiento numerico.\n";
        system("pause");
        return;
    }

    return;
}
