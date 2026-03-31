#pragma once
#include <bits/stdc++.h>
using namespace std;

#ifndef DB_DEFINED
#define DB_DEFINED
#define db double
using Mat = vector<vector<db>>;
#endif

void newton(){
    int n; 
    cout << "Ingrese la cantidad de puntos en la tabla: "; cin >> n;
    vector<db> x(n), y(n);
    cout << "Ingrese los puntos (x y) separados por espacio:\n";
    for(int i = 0; i < n; i++) cin >> x[i] >> y[i];

    cout << "x: "; for(int i = 0; i < n; i++) cout << x[i] << " ";
    cout << "\ny: "; for(int i = 0; i < n; i++) cout << y[i] << " ";
    cout << "\nConfirme que los puntos son correctos (s/n): ";  
    char confirm; cin >> confirm;
    if(confirm != 's') {
        int ind;
        cout << "Cual es la fila del valor erroneo? "; cin >> ind;
        cout << "Ingrese el nuevo valor (x y): "; cin >> x[ind] >> y[ind];
        return;
    }

    for(int i = 1; i < n; i++)
        if(x[i] <= x[i-1]) {
            swap(x[i], x[i-1]);
            swap(y[i], y[i-1]);
        }

    db h = x[1] - x[0];
    for(int i = 2; i < n; i++){
        if(fabs((x[i] - x[i-1]) - h) > 1e-9) {
            cout << "Los puntos no son equidistantes, no se puede aplicar el metodo de diferencias de Newton.\n";
            return;
        }
    }

    vector<vector<db>> tabla(n, vector<db>(n, 0.0));
    for(int i = 0; i < n; i++) tabla[i][0] = y[i];
    for(int j = 1; j < n; j++)
        for(int i = 0; i < n - j; i++)
            tabla[i][j] = tabla[i+1][j-1] - tabla[i][j-1];

    cout << "\nTabla de diferencias finitas:\n";
    cout << fixed << setprecision(6);
    for(int i = 0; i < n; i++){
        cout << x[i] << "\t" << tabla[i][0];
        for(int j = 1; j < n - i; j++) cout << "\t" << tabla[i][j];
        cout << "\n";
    }

    while(true){
        cout << "\nSeleccione el metodo:\n";
        cout << "1. Progresivo\n";
        cout << "2. Regresivo\n";
        cout << "Opcion: ";
        int metodo; cin >> metodo;
        if(metodo != 1 && metodo != 2){
            cout << "Opcion no valida.\n";
            continue;
        }

        db valor;
        cout << "Ingrese el valor a interpolar: "; cin >> valor;
        if(valor < x[0] || valor > x[n-1]) {
            cout << "El valor a interpolar esta fuera del rango de los puntos.\n";
            continue;
        }

        int grado;
        cout << "Ingrese el grado del polinomio interpolante: "; cin >> grado;
        if(grado >= n) {
            cout << "El grado debe ser menor que la cantidad de puntos (" << n << ").\n";
            continue;
        }

        db resultado;

        if(metodo == 1){
            db s = (valor - x[0]) / h;
            resultado = tabla[0][0];
            db coef = 1.0;
            db factorial = 1.0;

            for(int k = 1; k <= grado; k++){
                coef *= (s - (k - 1));
                factorial *= k;
                resultado += (coef / factorial) * tabla[0][k];
            }
            cout << "Metodo: Newton Progresivo\n";
        }
        else{
            int base = n - 1; 
            db s = (valor - x[base]) / h;
            resultado = tabla[base][0];
            db coef = 1.0;
            db factorial = 1.0;

            for(int k = 1; k <= grado; k++){
                coef *= (s + (k - 1));
                factorial *= k;
               
                resultado += (coef / factorial) * tabla[base - k][k];
            }
            cout << "Metodo: Newton Regresivo\n";
        }

        cout << "El valor interpolado en x = " << valor << " es: " << resultado << "\n";

        cout << "Desea interpolar otro valor? (s/n): ";
        char cont; cin >> cont;
        if(cont == 'n' || cont == 'N') {
            cout << "Desea cambiar los puntos de la tabla? (s/n): ";
            char cambiar; cin >> cambiar;
            if(cambiar == 's' || cambiar == 'S') newton();
            else cout << "Saliendo.\n";
            break;
        }
    }
}

