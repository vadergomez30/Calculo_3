#include <bits/stdc++.h>
using namespace std;

#define db double
using Mat = vector<vector<db>>;
static const double EPS = 1e-12;

db broyden1_funcion1(db x, db y) { return x*x + x*y + 2*y*y - 5; }
db broyden1_funcion2(db x, db y) { return 5*y - 2*x*y*y + 3; }

db broyden2_funcion1(db x, db y) { return x*x-3*y*y-10; }
db broyden2_funcion2(db x, db y) { return 2*y*y-3*x*y+1; }

db broyden3_funcion1(db x, db y, db z) { return 3*x*x+y*y-8*y+2*z*z-5; }
db broyden3_funcion2(db x, db y, db z) { return -2*x*x-12*x+y*y-3*z*z+10; }
db broyden3_funcion3(db x, db y, db z) { return -x+2*y+5*z; }

db broyden4_funcion1(db x, db y, db z) { return x*x+y*y-2*z+3; }
db broyden4_funcion2(db x, db y, db z) { return x+y+z-5; }
db broyden4_funcion3(db x, db y, db z) { return x*x-y*y+z*z-9; }    

Mat inversa(const Mat& A) {
    int n = A.size();
    Mat M = A; 
    Mat I(n, vector<db>(n, 0));

    for(int i = 0; i < n; i++)
        I[i][i] = 1;

    for(int i = 0; i < n; i++) {

        if (fabs(M[i][i]) < EPS)
            return Mat(); 

        db pivote = M[i][i];

        for(int j = 0; j < n; j++) {
            M[i][j] /= pivote;
            I[i][j] /= pivote;
        }

        for(int k = 0; k < n; k++) {
            if(k != i) {
                db factor = M[k][i];
                for(int j = 0; j < n; j++) {
                    M[k][j] -= factor * M[i][j];
                    I[k][j] -= factor * I[i][j];
                }
            }
        }
    }
    return I;
}

Mat mult(const Mat& A, const Mat& B){
    int n = A.size(), p = A[0].size(), m = B[0].size();
    Mat C(n, vector<db>(m,0.0));
    for(int i=0;i<n;i++)
        for(int k=0;k<p;k++)
            for(int j=0;j<m;j++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

Mat suma(const Mat& A, const Mat& B){
    int n = A.size(), m = A[0].size();
    Mat C(n, vector<db>(m,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

Mat resta(const Mat& A, const Mat& B){
    int n = A.size(), m = A[0].size();
    Mat C(n, vector<db>(m,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

Mat transpuesta(const Mat& A){
    int n = A.size(), m = A[0].size();
    Mat At(m, vector<db>(n,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            At[j][i] = A[i][j];
    return At;
}

void broyden(Mat x1, Mat F1, Mat J_inv, double error, Mat x0, Mat F0, int opc, int max_iter){

    int iter = 1;
    int n = x1.size();

    while(true){
        if(iter > max_iter){
            cout << "No se encontro solucion en " << max_iter << " iteraciones.\n";
            return;
        }
        if(n == 2){
        double normF = max(fabs(F1[0][0]), fabs(F1[1][0]));
        if(normF < error) break;
        }
        else if(n == 3){
        double normF = max({fabs(F1[0][0]), fabs(F1[1][0]), fabs(F1[2][0])});
        if(normF < error) break;
        }   
        Mat dx = resta(x1, x0);
        Mat dF = resta(F1, F0);

        Mat temp = mult(transpuesta(dx), mult(J_inv, dF));
        double denom = temp[0][0];

        Mat numerador = mult(
                            resta(dx, mult(J_inv, dF)),
                            mult(transpuesta(dx), J_inv)
                        );

        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                numerador[i][j] /= denom;

        Mat J_new = suma(J_inv, numerador);


        J_inv = J_new;
        x0 = x1;
        F0 = F1;

        x1 = resta(x0, mult(J_inv, F0));
        switch (opc){
        case 1:
            F1[0][0] = broyden1_funcion1(x1[0][0], x1[1][0]);
            F1[1][0] = broyden1_funcion2(x1[0][0], x1[1][0]);
            break;

        case 2:
            F1[0][0] = broyden2_funcion1(x1[0][0], x1[1][0]);
            F1[1][0] = broyden2_funcion2(x1[0][0], x1[1][0]);
            break;

        case 3:
            F1[0][0] = broyden3_funcion1(x1[0][0], x1[1][0], x1[2][0]);
            F1[1][0] = broyden3_funcion2(x1[0][0], x1[1][0], x1[2][0]);
            F1[2][0] = broyden3_funcion3(x1[0][0], x1[1][0], x1[2][0]);
            break;

        case 4:
            F1[0][0] = broyden4_funcion1(x1[0][0], x1[1][0], x1[2][0]);
            F1[1][0] = broyden4_funcion2(x1[0][0], x1[1][0], x1[2][0]);
            F1[2][0] = broyden4_funcion3(x1[0][0], x1[1][0], x1[2][0]);
            break;
        }
        
        iter++;
        if(n == 2)
        cout << " x" << iter << " = " << x1[0][0] 
             << ", y" << iter << " = " << x1[1][0] << "\n";
        else if(n == 3)
        cout << " x" << iter << " = " << x1[0][0] 
             << ", y" << iter << " = " << x1[1][0] 
             << ", z" << iter << " = " << x1[2][0] << "\n";
    }

    if(n == 2){
        cout << "Solucion encontrada: x = " << x1[0][0] << ", y = " << x1[1][0] << "\n";
    }
    else if(n == 3){
        cout << "Solucion encontrada: x = " << x1[0][0] << ", y = " << x1[1][0] << ", z = " << x1[2][0] << "\n";
    }
    cout << "Numero de iteraciones: " << iter << endl;
}


void Newton(int opc){
    int n;
    if(opc == 1 || opc == 2) n = 2;
    else if(opc == 3 || opc == 4) n = 3;
    int max_iter;

    Mat J(n, vector<db>(n));
    Mat J_inv(n, vector<db>(n));
    Mat x0(n, vector<db>(1));
    Mat x1(n, vector<db>(1));
    Mat F0(n, vector<db>(1));
    Mat F1(n, vector<db>(1));      
    db error;
    switch(opc){
    case 1:
        cout << "f1(x,y)=x^2+xy+2y^2-5\nf2(x,y)=5y-2xy^2+3\n";
        cout << "Ingrese x0: ";  cin >> x0[0][0];
        cout << "Ingrese y0: ";  cin >> x0[1][0];
        cout << "Ingrese tolerancia: "; cin >> error;
        cout << "Ingerese numero maximo de iteraciones: "; cin >> max_iter;
        J[0][0] = 2*x0[0][0] + x0[1][0];
        J[0][1] = x0[0][0] + 4*x0[1][0];    
        J[1][0] = -2*x0[1][0]*x0[1][0];
        J[1][1] = 5 - 4*x0[0][0]*x0[1][0];

        J_inv = inversa(J);
        if(J_inv.empty()){
            cout << "Jacobiano singular inicial.\n";
            return;
        }

        F0[0][0] = broyden1_funcion1(x0[0][0], x0[1][0]);
        F0[1][0] = broyden1_funcion2(x0[0][0], x0[1][0]);
        x1 = resta(x0, mult(J_inv, F0));
        F1[0][0] = broyden1_funcion1(x1[0][0], x1[1][0]);
        F1[1][0] = broyden1_funcion2(x1[0][0], x1[1][0]);

        break;
    case 2:
        cout << "f1(x,y)=x^2-3y^2-10\nf2(x,y)=2y^2-3xy+1\n";
        cout << "Ingrese x0: ";  cin >> x0[0][0];
        cout << "Ingrese y0: ";  cin >> x0[1][0];
        cout << "Ingrese tolerancia: "; cin >> error;
        cout << "Ingerese numero maximo de iteraciones: "; cin >> max_iter;

        J[0][0] = 2*x0[0][0];
        J[0][1] = -6*x0[1][0];    
        J[1][0] = -3*x0[1][0];
        J[1][1] = 4*x0[1][0] - 3*x0[0][0];

        J_inv = inversa(J);
        if(J_inv.empty()){
            cout << "Jacobiano singular inicial.\n";
            return;
        }
        F0[0][0] = broyden2_funcion1(x0[0][0], x0[1][0]);
        F0[1][0] = broyden2_funcion2(x0[0][0], x0[1][0]);
        x1 = resta(x0, mult(J_inv, F0));
        F1[0][0] = broyden2_funcion1(x1[0][0], x1[1][0]);
        F1[1][0] = broyden2_funcion2(x1[0][0], x1[1][0]);

        break;
    case 3:
        cout << "f1(x,y,z)=3x^2+y^2-8y+2z^2-5\nf2(x,y,z)=-2x^2-12x+y^2-3z^2+10\nf3(x,y,z)=-x+2y+5z\n";
        cout << "Ingrese x0: ";  cin >> x0[0][0];
        cout << "Ingrese y0: ";  cin >> x0[1][0];
        cout << "Ingrese z0: ";  cin >> x0[2][0];
        cout << "Ingrese tolerancia: "; cin >> error;
        cout << "Ingerese numero maximo de iteraciones: "; cin >> max_iter;
        J[0][0] = 6*x0[0][0];
        J[0][1] = 2*x0[1][0] - 8;    
        J[0][2] = 4*x0[2][0];
        J[1][0] = -4*x0[0][0] - 12;
        J[1][1] = 2*x0[1][0];
        J[1][2] = -6*x0[2][0];
        J[2][0] = -1;
        J[2][1] = 2;
        J[2][2] = 5;

        J_inv = inversa(J);
        if(J_inv.empty()){
            cout << "Jacobiano singular inicial.\n";
            return;
        }

        F0[0][0] = broyden3_funcion1(x0[0][0], x0[1][0], x0[2][0]);
        F0[1][0] = broyden3_funcion2(x0[0][0], x0[1][0], x0[2][0]);
        F0[2][0] = broyden3_funcion3(x0[0][0], x0[1][0], x0[2][0]);
        x1 = resta(x0, mult(J_inv, F0));
        F1[0][0] = broyden3_funcion1(x1[0][0], x1[1][0], x1[2][0]);
        F1[1][0] = broyden3_funcion2(x1[0][0], x1[1][0], x1[2][0]);
        F1[2][0] = broyden3_funcion3(x1[0][0], x1[1][0], x1[2][0]);

        break;
    case 4:
        cout << "f1(x,y,z)=x^2+y^2-2z+3\nf2(x,y,z)=x+y+z-5\nf3(x,y,z)=x^2-y^2+z^2-9\n";
        cout << "Ingrese x0: ";  cin >> x0[0][0];
        cout << "Ingrese y0: ";  cin >> x0[1][0];
        cout << "Ingrese z0: ";  cin >> x0[2][0];
        cout << "Ingrese tolerancia: "; cin >> error;
        cout << "Ingerese numero maximo de iteraciones: "; cin >> max_iter;
        J[0][0] = 2*x0[0][0];
        J[0][1] = 2*x0[1][0];
        J[0][2] = -2;
        J[1][0] = 1;
        J[1][1] = 1;
        J[1][2] = 1;
        J[2][0] = 2*x0[0][0];
        J[2][1] = -2*x0[1][0];
        J[2][2] = 2*x0[2][0];

        J_inv = inversa(J);
        if(J_inv.empty()){
            cout << "Jacobiano singular inicial.\n";
            return;
        }

        F0[0][0] = broyden4_funcion1(x0[0][0], x0[1][0], x0[2][0]);
        F0[1][0] = broyden4_funcion2(x0[0][0], x0[1][0], x0[2][0]);
        F0[2][0] = broyden4_funcion3(x0[0][0], x0[1][0], x0[2][0]);
        x1 = resta(x0, mult(J_inv, F0));
        F1[0][0] = broyden4_funcion1(x1[0][0], x1[1][0], x1[2][0]);
        F1[1][0] = broyden4_funcion2(x1[0][0], x1[1][0], x1[2][0]);
        F1[2][0] = broyden4_funcion3(x1[0][0], x1[1][0], x1[2][0]);

        break;

    }
    

    if(n == 2)
        cout << " x1 = " << x1[0][0] << ", y1 = " << x1[1][0] << "\n";
    else if(n == 3)
        cout << " x1 = " << x1[0][0] << ", y1 = " << x1[1][0] << ", z1 = " << x1[2][0] << "\n";

    broyden(x1, F1, J_inv, error, x0, F0, opc, max_iter);
    cout<<"Desea probar otros datos iniciales? (s/n): ";
    char cont; cin>>cont;
    if(cont == 's' || cont == 'S'){
        Newton(opc);
    }
}


int main(){
    while(true){
        cout<<"\n\nMetodo de Broyden para sistemas no lineales\n\n";
        cout<<"Desarrollado por: \nGomez Perez Vader Ali\nRamos Renteria Emiliano\nAlmaraz Remigio Luis\n\n";
        cout<<"Seleccione el sistema de ecuaciones a resolver:\n\n";    
        cout<<"1.f1(x,y)=x^2+xy+2y^2-5\n  f2(x,y)=5y-2xy^2+3\n\n";
        cout<<"2.f1(x,y)=x^2-3y^2-10\n  f2(x,y)=2y^2-3xy+1\n\n";
        cout<<"3.f1(x,y,z)=3x^2+y^2-8y+2z^2-5\n  f2(x,y,z)=-2x^2-12x+y^2-3z^2+10\n  f3(x,y,z)=-x+2y+5z\n\n";
        cout<<"4.f1(x,y,z)=x^2+y^2-2z+3\n  f2(x,y,z)=x+y+z-5\n  f3(x,y,z)=x^2-y^2+z^2-9\n\n";
        cout<<"5.Salir\n\nOpcion: ";
        int opc; 
        cin>>opc;
        if(opc < 1 || opc > 4){
            cout << "Saliendo.\n";
            return 0;
        }
        Newton(opc);
        cout<<"Desea resolver otro sistema? (s/n): ";
        char cont; cin>>cont;
        if(cont == 'n' || cont == 'N'){
            cout << "Saliendo.\n";  
            break;
        }
    }   

    return 0;
}
