#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define db double
#define MAX_ITER 2000
#define kn 1000
db f1(db x) { return sin(x) + cos(x); }
db df1(db x) { return cos(x) - sin(x); }

db f2(db x) { return (x - 5) * (x + 3) * pow(cos(x), 3); }
db df2(db x) { return cos(x) * cos(x) * (2 * (x - 1) * cos(x) - 3 * (x * x - 2 * x - 15) * sin(x)); }

db f3(db x) { return x * x * sin(x) * sin(x); }
db df3(db x) { return 2 * x * sin(x) * (sin(x) + x * cos(x)); }

db f4(db x) { return exp(x) * sin(x); }
db df4(db x) { return exp(x) * (sin(x) + cos(x)); }

void ejecutarNewtonRaphson() { 
    int opc;
    
    do {
        
        cout << "Metodo de Newton-Raphson\n\nSelecciona una funcion (1-4):\n\n";
        cout << "1. f(x)=sen(x)+cos(x)\n";
        cout << "2. f(x)=(x-5)(x+3)cos^3(x)\n";
        cout << "3. f(x)=x^2 sen^2(x)\n";
        cout << "4. f(x)=e^x sen(x)\n";
        cout << "5. Salir\n";
        cout << "Opcion: ";
        cin >> opc;
        
        if (opc >= 1 && opc <= 4) {
            db x, error;
            
           
            switch(opc) {
                case 1: cout << "\nf(x)=sen(x)+cos(x)\n"; break;
                case 2: cout << "\nf(x)=(x-5)(x+3)cos^3(x)\n"; break;
                case 3: cout << "\nf(x)=x^2 sen^2(x)\n"; break;
                case 4: cout << "\nf(x)=e^x sen(x)\n"; break;
            }
            
            cout << "Ingrese un valor aproximado a la raiz: ";
            cin >> x;
            cout << "Ingrese el error de tolerancia: ";
            cin >> error;
            
            
            bool exito = true;
            int iteraciones = 0;
            db (*f)(db) = nullptr;
            db (*df)(db) = nullptr;
            
           
            switch(opc) {
                case 1: f = f1; df = df1; break;
                case 2: f = f2; df = df2; break;
                case 3: f = f3; df = df3; break;
                case 4: f = f4; df = df4; break;
            }
            
           
            while (abs(f(x)) > error && iteraciones < 1000) {
                if (df(x) == 0) {
                    cout << "\nError: Derivada cero. Da otro valor inicial.\n";
                    exito = false;
                    break;
                } else {
                    x = x - (f(x) / df(x));
                }
                iteraciones++;
            }
            
           
            if (exito && iteraciones < 1000) {
                cout << "\nLa raiz es: " << x << '\n';
                cout << "Iteraciones realizadas: " << iteraciones << '\n';
            } else if (iteraciones >= 1000) {
                cout << "\nError: No converge en 1000 iteraciones.\n";
            }
            cout << "----------------------------------------\n";
            
        } else if (opc != 5) {
            cout << "Opcion no valida. Intente de nuevo.\n";
        }
        
    } while (opc != 5);
    
    cout << "Programa terminado.\n";
}

double func1(float x) { return sin(x) + cos(x); }
double func2(float x) { return (x - 5) * (x + 3) * pow(cos(x), 3); }
double func3(float x) { return x * x * pow(sin(x), 2); }
double func4(float x) { return exp(x) * sin(x); }

void ejecutarMetodoSecante() {
    int opc, iteraciones;
    float x0, x1, error;
    
    cout << "Metodo de la secante\n";
    cout << "Elige la funcion\n";
    cout << "1. f(x)=sin(x)+cos(x)\n";
    cout << "2. f(x)=(x-5)(x+3)cos^3(x)\n";
    cout << "3. f(x)=x^2*sin^2(x)\n";
    cout << "4. f(x)=e^x*sin(x)\n";
    cin >> opc;
    
    cout << "Ingresa el valor de x0\n";
    cin >> x0;
    cout << "Ingresa el valor de x1\n";
    cin >> x1;
    cout << "Ingresa el error\n";
    cin >> error;
    cout << "Ingresa el numero maximo de iteraciones\n";
    cin >> iteraciones;
    
    int i = 0;
    bool raizEncontrada = false;
    
    double (*func)(float) = nullptr;
    
    switch (opc) {
        case 1: func = func1; break;
        case 2: func = func2; break;
        case 3: func = func3; break;
        case 4: func = func4; break;
        default: 
            cout << "Opcion no valida" << endl;
            return;
    }
    
   
    while (i < iteraciones) {
        float f0 = func(x0);
        float f1 = func(x1);
        
       
        if (fabs(f0 - f1) < 1e-10) {
            cout << "Error: Division por cero. Los puntos tienen valores de funcion muy cercanos.\n";
            break;
        }
        
        float x2 = x1 - ((f1 * (x0 - x1)) / (f0 - f1));
        
        if (fabs(x2 - x1) < error) {
            cout << "La raiz es: " << x2 << "\n";
            cout << "Numero de iteraciones: " << i + 1 << endl;
            raizEncontrada = true;
            break;
        }
        
        x0 = x1;
        x1 = x2;
        i++;
    }
    
    if (!raizEncontrada && i >= iteraciones) {
        cout << "No se encontro la raiz en el numero maximo de iteraciones (" << iteraciones << ")\n";
        cout << "Ultima aproximacion: " << x1 << endl;
    }
}

float funcion1(float x) { return (sin(x) + cos(x)); }
float funcion2(float x) { return ((x - 5) * (x + 3) * pow(cos(x), 3)); }
float funcion3(float x) { return (2 * x * pow(sin(x), 2) + 2 * pow(x, 2) * cos(x) * sin(x)); }
float funcion4(float x) { return (exp(x) * sin(x)); }


void ejecutarMetodoBiseccion() {
    
    float (*funciones[4])(float) = { funcion1, funcion2, funcion3, funcion4 };

    float a, b, x, fa, fb, fx;
    float error;
    int op, para = 0;
    bool raiz = false;

    cout << "\t--Metodo de biseccion--\n";
    cout << "Funciones disponibles: \n1.f(x)=Sen(x)+Cos(x)\n2.f(x)=(x-5)(x+3)Cos3(x)\n";
    cout << "3.f(x)=x^2*Sen^2(x)\n4.f(x)e^x*Sen(x)\n";

    cout << "Ingrese la funcion: ";
    cin >> op;

    if (op < 1 || op > 4) {
        cout << "La funcion no existe...\n";
        return;
    }

    
    cout << "Ingrese el intervalo [a,b]\n";
    cout << "Ingrese a: ";
    cin >> a;
    cout << "Ingresa b: ";
    cin >> b;
    cout << "Ingresa el error: ";
    cin >> error;

   
    fa = funciones[op - 1](a);
    fb = funciones[op - 1](b);
    
    if (fa * fb >= 0) {
        cout << "\nNo existe cambio de signo. No se puede aplicar el metodo...\n";
        return;
    }
    
   
    if (op == 3) {
        if (funciones[2](b) < 0) {
            cout << "\nNo existe cambio de signo. No se puede aplicar el metodo...\n";
            return;
        }
    }

    
    do {
        fa = funciones[op - 1](a);
        fb = funciones[op - 1](b);        

        if (abs(fa) < error) {
            cout << "La raiz es: " << fixed << setprecision(8) << a << "\n";
            raiz = true;
        }
        if (abs(fb) < error) {
            cout << "La raiz es: " << fixed << setprecision(8) << b << "\n";
            raiz = true;
        }

        
        x = (a + b) / 2;
        fx = funciones[op - 1](x);

        if (abs(fx) < error) {
            cout << "La raiz es: " << fixed << setprecision(8) << x << "\n";
            raiz = true;
        }

        if (!raiz) {
            cout << "En el intervalo: [" << a << "," << b << "] aun no se encuentra la raiz.\n";
            if (fa * fx < 0) {
                b = x;
                cout << "El nuevo intervalo es: [" << a << "," << b << "]\n";
            } else {
                a = x;
                cout << "El nuevo intervalo es: [" << a << "," << b << "]\n";
            }
            cout << "\n";
        }

        para++;
        if (para > 49) {
            cout << "Se alcanzo el maximo de iteraciones...\n";
            break;
        }
    } while (!raiz && para < 50);
}
double f1(float x){
    return sin(x)+cos(x);
}
double f1_1(float x){
    return cos(x)-sin(x);
}
double f2_a(float x){  
    return (x-5)*(x+3)*pow(cos(x),3);
}
double f2_1d(float x){  
    return pow(cos(x),2)*( 2*(x-1)*cos(x)-3*(x*x-2*x-15)*sin(x) );
}
double f3(float x){
    return x*x*pow(sin(x),2);
}
double f3_1(float x){
    return 2*x*sin(x)*(sin(x)+x*cos(x));
}
double f4_a(float x){
    return exp(x)*sin(x);
}
double f4_1d(float x){
    return exp(x)*(sin(x)+cos(x));
}

void ejecutarFalsaPosicion() {
    double a,b,error;
    int opc;
    cout<<"Metodo de la Falsa Posicion"<<endl;
    cout<<"Funcion 1. f(x)=sin(x)+cos(x)"<<endl;
    cout<<"Funcion 2. f(x)=(x-5)(x+3)Cos^3(x)"<<endl;
    cout<<"Funcion 3. f(x)=x^2*sen^2(x)"<<endl;
    cout<<"Funcion 4. f(x)=e^x*sen(x)"<<endl;
    cout<<"Seleccione una funcion (1-4): ";
    cin>>opc;
    
    switch(opc){
        case 1:
            cout<<"Has seleccionado la funcion f(x)=sin(x)+cos(x)"<<endl;
            cout<<"Ingrese el valor de a: ";cin>>a;
            cout<<"Ingrese el valor de b: ";cin>>b;
            cout<<"Ingrese el valor del error: ";cin>>error;
            if(f1(a)*f1(b)<0 && f1_1(a)*f1_1(b)>0){
                cout<<"Se puede aplicar el metodo"<<endl;
            }
            else{
                cout<<"No se puede aplicar el metodo"<<endl;
                return;
            }
            for(int j=0;j<kn;j++){
                double xr=a-f1(a)*(b-a)/(f1(b)-f1(a));
                if(f1(a)*f1(xr)<0){
                    b=xr;
                }
                else{
                    a=xr;
                }
                if(fabs(f1(xr))<error){
                    cout<<"La raiz es: "<<xr<<endl;
                    return;
                }
            }
            cout<<"No se encontro la raiz en "<<kn<<" iteraciones"<<endl;
            return;
            
        case 2:
            cout<<"Has seleccionado la funcion f(x)=(x-5)(x+3)Cos^3(x)"<<endl;
            cout<<"Ingrese el valor de a: ";cin>>a;
            cout<<"Ingrese el valor de b: ";cin>>b;
            cout<<"Ingrese el valor del error: ";cin>>error;
            if(f2(a)*f2(b) < 0){  
                cout<<"Se puede aplicar el metodo"<<endl;
            }
            else{
                cout<<"No se puede aplicar el metodo"<<endl;
                return;
            }
            for(int j=0;j<kn;j++){
                long double xr=a-f2_a(a)*(b-a)/(f2_a(b)-f2_a(a));  
                if(f2_a(a)*f2_a(xr)<0){  
                    b=xr;
                    cout<<b<<endl;
                }
                else{
                    a=xr;
                    cout<<a<<endl;
                }
                if(fabs(f2_a(xr))<error){  
                    cout<<"La raiz es: "<<xr<<endl;
                    return;
                }
            }
            cout<<"No se encontro la raiz en "<<kn<<" iteraciones"<<endl;
            return;
            
        case 3:
            cout<<"Has seleccionado la funcion f(x)=x^2*sen^2(x)"<<endl;
            cout<<"Ingrese el valor de a: ";cin>>a; 
            cout<<"Ingrese el valor de b: ";cin>>b;
            cout<<"Ingrese el valor del error: ";cin>>error;
            if(f3(a)*f3(b)<0 && f3_1(a)*f3_1(b)>0){
                cout<<"Se puede aplicar el metodo"<<endl;
            }
            else{
                cout<<"No se puede aplicar el metodo"<<endl;
                return;
            }
            for(int j=0;j<kn;j++){
                double xr=a-f3(a)*(b-a)/(f3(b)-f3(a));
                if(f3(a)*f3(xr)<0){
                    b=xr;
                }
                else{
                    a=xr;
                }
                if(fabs(f3(xr))<error){
                    cout<<"La raiz es: "<<xr<<endl;
                    return;
                }
            }
            cout<<"No se encontro la raiz en "<<kn<<" iteraciones"<<endl;
            return;
            
        case 4:
            cout<<"Has seleccionado la funcion f(x)=e^x*sen(x)"<<endl;
            cout<<"Ingrese el valor de a: ";cin>>a;
            cout<<"Ingrese el valor de b: ";cin>>b;
            cout<<"Ingrese el valor del error: ";cin>>error;
            if(f4_a(a)*f4_a(b)<0 && f4_1d(a)*f4_1d(b)>0){
                cout<<"Se puede aplicar el metodo"<<endl;
            }
            else{
                cout<<"No se puede aplicar el metodo: f(a) y f(b) deben tener signos opuestos"<<endl;
                return;
            }
            for(int j=0;j<kn;j++){
                long double xr=a-(f4_a(a)*(b-a)/(f4_a(b)-f4_a(a)));
                if(f4(a)*f4(xr)<0){
                    b=xr;
                }
                else{
                    a=xr;
                }
                if(fabs(f4_a(xr))<error){
                    cout<<"La raiz es: "<<xr<<endl;
                    return;
                }
            }
            cout<<"No se encontro la raiz en "<<kn<<" iteraciones"<<endl;
            return;
            
        default:
            cout<<"Opcion no valida"<<endl;
            return;
    }
}