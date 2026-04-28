#include <bits/stdc++.h>
using namespace std;
void splineCubico() {
    char resp = 's';
    while (resp=='s' || resp=='S')
    {
    char resp2 = 's';
    int n;
    cout << "Ingrese la cantidad de puntos: ";
    cin >> n;
    if(cin.fail()) {
        cout << "Error: Ingrese solo numeros.\n";
        cin.clear();
        cin.ignore(10000, '\n');
        return;
    }
    cout << "Ingrese los puntos (x y) separados por espacio:\n";
    vector<vector<float>> mat(n, vector<float>(2));
    for(int i = 0; i < n; i++){
        cin >> mat[i][0] >> mat[i][1];
        if(cin.fail()) {
            cout << "Error: Ingrese solo numeros.\n";
            cin.clear();
            cin.ignore(10000, '\n');
            return;
        }
    }
    if(n < 2) {
        cout << "Se necesitan al menos 2 puntos para construir un spline cubico.\n";
        return;
    }
     //verificar que si ingrese n puntos
     if(mat.size() !=n){
        cout << "Ingresaste un numero diferente de puntos al que dijiste que ibas a ingresar.\n";
        cout<<"Deseas ingresar otro grupo de puntos? (s/n): ";        cin>>resp;
        if(resp == 'n' || resp == 'N') return;
        else continue;  
     }
    cout<<"Los puntos x_i , y_i, son correctos (s/n)?"<<"\n";
    for(int i = 0; i<n; i++){
        cout<<i<<": "<< mat[i][0]<<", "<< mat[i][1];
        cout<<"\n";
    }
    cin>>resp2;
    if(resp2== 'n'|| resp2== 'N'){
        cout<<"Ingresa el indice a corregir: ";
        int idx;
        cin>>idx;
        cout<<"Ingresa el nuevo valor de x: ";
        cin>>mat[idx][0];
        cout<<"Ingresa el nuevo valor de y: ";
        cin>>mat[idx][1];
    }
    vector<float> h(n-1), f(n-1);
// se calculan las h y las diferencias 
    for(int i = 0; i < n-1; i++){
        h[i] = mat[i+1][0] - mat[i][0];
        f[i] = (mat[i+1][1] - mat[i][1]) / h[i];
    }
   int m= n-2;
   vector<vector<double>> A(m, vector<double>(m, 0)); //matriz tridiagonal
   vector<double> B(m); // matriz de resultados
   for(int i=0; i<m; i++) {
    A[i][i] = 2*(h[i]+h[i+1]); // diagonal principal
    if(i > 0) A[i][i-1] = h[i]; // diagonal inferior
    if(i < m-1) A[i][i+1] = h[i+1]; // diagonal superior
    B[i]= 6*(f[i+1]-f[i]); 
   }
   cout << "\nMatriz A:\n";
   for(int i=0; i<m; i++) {
    for(int j=0; j<m; j++) cout << A[i][j] << " ";
    cout << "\n";
   }
    cout << "\nVector B:\n";
    for(int i=0; i<m; i++) cout << B[i] << "\n";
  
   // B a matriz columna
   Mat B_mat(m, vector<double>(1, 0));
   for(int i = 0; i < m; i++) {
       B_mat[i][0] = B[i];
   }
   Mat A_inv = inversa(A);
   
   if(A_inv.empty()) {
       cout << "Error: La matriz A no es invertible.\n";
       return;
   }
   
   cout  << "Resolviendo el sistema, tenemos que S):\n";
   Mat X = mult(A_inv, B_mat);
   vector<double> S(m);
    for(int i = 0; i < m; i++) {
         S[i] = X[i][0];
    }
    for(int i = 0; i < m; i++) {
        cout << "\nS[" << i+1 << "] = " << S[i] << "\n";
    }
    vector<double> a(n-1), b(n-1), c(n-1), d(n-1);

    vector<double> S_b(n);
    S_b[0] = 0;
    for(int i = 0; i < m; i++) {
        S_b[i+1] = S[i];
    }
    S_b[n-1] = 0;
    
    for ( int i = 0; i < n-1; i++)
    {
        a[i] = (S_b[i+1] -S_b[i]) / (6 * h[i]); cout << h[i];;
        b[i] = S_b[i] / 2;
        c[i] = f[i] - ((S_b[i+1] + 2*S_b[i]) * h[i] / 6);
        d[i] = mat[i][1];

    }
    cout << "\nCoeficientes:\n";
    for(int i = 0; i < n-1; i++) {
        cout << "a[" << i+1 << "] = " << a[i] << "\n";
        cout << "b[" << i+1 << "] = " << b[i] << "\n";
        cout << "c[" << i+1 << "] = " << c[i] << "\n";
        cout << "d[" << i+1 << "] = " << d[i] << "\n";
    }
    cout << "\nLas funciones de spline cubicas son:\n";
    for(int i = 0; i < n-1; i++) {
        cout << "S" << i+1 << "(x) = " << a[i] << "*(x - " << mat[i][0] << ")^3 + " 
             << b[i] << "*(x - " << mat[i][0] << ")^2 + "
             << c[i] << "*(x - " << mat[i][0] << ") + "
             << d[i] << "\n";
    }
    cout << "\nDesea ingresar otro conjunto de puntos? (s/n): ";
    cin >> resp;

}
}