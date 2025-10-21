#include "sistemalineales_primerosmetodos.h"
#include <iostream>
#include <cmath>
#include "M_Gauss.h"
using namespace std;
int main(){
    int opc;
    cout<<"Menu Principal\n";
    cout<<"1. Solucion numerica de ecuaciones de una sola variable"<<endl;
    cout<<"2.Soluciones de sistemas ecuaciones lineales"<<endl;
    cout<<"Seleccione una opcion (1-2): ";
    cin>>opc;
    switch(opc){
        case 1: {
            int opc2;
            cout<<"Has seleccionado la opcion 1"<<endl;
            cout<<"Seleccione el metodo a utilizar: "<<endl;
            cout<<"1. Metodo de la Biseccion"<<endl;
            cout<<"2. Metodo de la Falsa Posicion"<<endl;
            cout<<"3. Metodo de Newton Raphson"<<endl;
            cout<<"4. Metodo de la Secante"<<endl;
            cout<<"Seleccione una opcion (1-4): ";
            cin>>opc2;
            switch(opc2){
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
                default:
                    cout<<"Opcion no valida"<<endl;
            }
            break;
        }
        case 2: {
            int opc3;
            cout<<"Has seleccionado la opcion 2"<<endl;
            cout<<"Selecciona el metodo a utilizar: "<<endl;
            cout<<"1. Metodo de Gauss"<<endl;
            cout<<"2. Metodo de Gauss-Jordan"<<endl;
            cout<<"3. Metodo de la Inversa"<<endl;
            cout<<"Seleccione una opcion (1-3): ";
            cin>>opc3;
            switch(opc3){
                case 1:
                    Gauss();
                    break;
                case 2:
                    GaussJordan();
                    break;
                case 3:
                    Inversa();
                    break;
                default:
                    cout<<"Opcion no valida"<<endl;
            }
            break;
        }
        default:
            cout<<"Opcion no valida"<<endl;
    }
}