#include <bits/stdc++.h>
#include "Broyden.h"
#include "Diferencias_Newton.h"
#include "Spline_Cubico.h"
using namespace std;

int main() {
    bool salir = false;
    
    cout << "Menu Principal - Metodos Numericos 2" << endl;
    cout<<"Desarrollado por: \nGomez Perez Vader Ali\nRamos Renteria Emiliano\nAlmaraz Remigio Luis\n\n";

    while(!salir) {
        cout << "1. Sistemas de ecuaciones no lineales (Broyden)" << endl;
        cout << "2. Interpolacion polinomial (Diferencias de Newton)" << endl;
        cout << "3. Spline Cubico" << endl;
        cout << "4. Salir" << endl;
        cout << "Seleccione una opcion: ";
        int opcion; cin >> opcion;
        switch(opcion) {
            case 1:
                cout << "Sistemas de ecuaciones no lineales seleccionados." << endl;
                broyden();
                break;

            case 2:
                cout << "Interpolacion polinomial seleccionada." << endl;
                newton();
                break;
            case 3:
                cout << "Spline Cubico seleccionada." << endl;
                splineCubico();
                break;
            case 4:

                cout << "Saliendo." << endl;
                salir = true;

                break;
            default:
                cout << "Opcion no valida." << endl;
        }
    }

    return 0;
}