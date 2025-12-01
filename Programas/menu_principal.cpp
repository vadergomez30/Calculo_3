#include <bits/stdc++.h>
#include "sistemalineales_primerosmetodos.h"
#include "M_Gauss.h"
#include "Particionado.h"
#include "GaussP.h"
#include "intercambioCompleto.h"
#include "Jacobi.h"
#include "Gauss-Seidel.h"
#include "Cholesky.h"
using namespace std;

int main() {
    int opc;

    do {
        cout << "\n===== MENU PRINCIPAL =====\n";
        cout << "1. Solucion numerica de ecuaciones de una sola variable\n";
        cout << "2. Soluciones de sistemas de ecuaciones lineales\n";
        cout << "3. Metodos iterativos para sistemas de ecuaciones lineales\n";
        cout << "4. Factorizacion LU\n";
        cout << "5. Salir\n";
        cout << "Seleccione una opcion (1-5): ";
        cin >> opc;

        switch (opc) {
            case 1: {
                int opc2;
                int repetir; // 1 = seguir en este submenu, 0 = regresar al menu principal
                do {
                    cout << "\n--- Opcion 1: Ecuaciones de una sola variable ---\n";
                    cout << "Seleccione el metodo a utilizar: \n";
                    cout << "1. Metodo de la Biseccion\n";
                    cout << "2. Metodo de la Falsa Posicion\n";
                    cout << "3. Metodo de Newton Raphson\n";
                    cout << "4. Metodo de la Secante\n";
                    cout << "5. Volver al menu principal\n";
                    cout << "Seleccione una opcion (1-5): ";
                    cin >> opc2;

                    switch (opc2) {
                        case 1:
                            ejecutarMetodoBiseccion();
                            break;
                        case 2:
                            ejecutarFalsaPosicion();
                            break;
                        case 3:
                            ejecutarNewtonRaphson();
                            break;
                        case 4:
                            ejecutarMetodoSecante();
                            break;
                        case 5:
                            cout << "Regresando al menu principal...\n";
                            break;
                        default:
                            cout << "Opcion no valida\n";
                    }

                    if (opc2 == 5) {
                        repetir = 0; // salir del do de esta seccion
                    } else {
                        cout << "\n¿Desea escoger otro metodo de esta seccion? (1=Si, 0=No): ";
                        cin >> repetir;
                    }

                } while (repetir == 1);
                break;
            }

            case 2: {
                int opc2;
                int repetir;
                do {
                    cout << "\n--- Opcion 2: Sistemas de ecuaciones lineales ---\n";
                    cout << "Selecciona el metodo a utilizar: \n";
                    cout << "1. Metodo de Gauss\n";
                    cout << "2. Metodo de Gauss-Jordan\n";
                    cout << "3. Metodo de la Inversa\n";
                    cout << "4. Inversa particionado\n";
                    cout << "5. Gauss Particionado\n";
                    cout << "6. Intercambio\n";
                    cout << "7. Volver al menu principal\n";
                    cout << "Seleccione una opcion (1-7): ";
                    cin >> opc2;

                    switch (opc2) {
                        case 1:
                            Gauss();
                            break;
                        case 2:
                            GaussJordan();
                            break;
                        case 3:
                            Inversa();
                            break;
                        case 4:
                            Particionado();
                            break;
                        case 5:
                            GaussParticionado();
                            break;
                        case 6:
                            Intercambio();
                            break;
                        case 7:
                            cout << "Regresando al menu principal...\n";
                            break;
                        default:
                            cout << "Opcion no valida\n";
                    }

                    if (opc2 == 7) {
                        repetir = 0;
                    } else {
                        cout << "\n¿Desea escoger otro metodo de esta seccion? (1=Si, 0=No): ";
                        cin >> repetir;
                    }

                } while (repetir == 1);
                break;
            }

            case 3: {
                int opc2;
                int repetir;
                do {
                    cout << "\n--- Opcion 3: Metodos iterativos ---\n";
                    cout << "Selecciona el metodo a utilizar: \n";
                    cout << "1. Metodo de Jacobi\n";
                    cout << "2. Metodo de Gauss-Seidel\n";
                    cout << "3. Metodo de relajacion\n";
                    cout << "4. Volver al menu principal\n";
                    cout << "Seleccione una opcion (1-4): ";
                    cin >> opc2;

                    switch (opc2) {
                        case 1:
                            Jacobi();
                            break;
                        case 2:
                            GaussSeidel();
                            break;
                        case 3:
                            M_relajacion();
                            break;
                        case 4:
                            cout << "Regresando al menu principal...\n";
                            break;
                        default:
                            cout << "Opcion no valida\n";
                    }

                    if (opc2 == 4) {
                        repetir = 0;
                    } else {
                        cout << "\n¿Desea escoger otro metodo de esta seccion? (1=Si, 0=No): ";
                        cin >> repetir;
                    }

                } while (repetir == 1);
                break;
            }

            case 4: {
                int repetir;
                do {
                    cout << "\n--- Opcion 4: Factorizacion LU ---\n";
                    menuLU();
                    cout << "\n¿Desea realizar otra factorizacion LU? (1=Si, 0=No): ";
                    cin >> repetir;
                } while (repetir == 1);
                break;
            }

            case 5:
                cout << "Saliendo del programa...\n";
                break;

            default:
                cout << "Opcion no valida\n";
        }

    } while (opc != 5);

    return 0;
}
