//g++ -std=c++20 Proyecto_Final.cpp -o Proyecto_Final.exe
#include <bits/stdc++.h>
using namespace std;

#define ll long long 
using Mat = vector<vector<ll>>;
static const double EPS = 1e-12;

Mat multiplicar(const Mat& A, const Mat& B, int mod) {
    int n = A.size();
    Mat C(n, vector<ll>(n, 0));
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
    return C;
}

Mat suma(const Mat& A, const Mat& B){
    int n = A.size(), m = A[0].size();
    Mat C(n, vector<ll>(m,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

int main() {
	int ver, lin, tipoGrafica;
    unordered_set<char> vertices;
    unordered_map<char,int> indice;
    int mod = 1e9 + 7; 
    char masmenos = 241;
    char aux;
    aux = 'A';
    bool simple=true, conectada=true, regular=true;
    cout<<"Ingresa la cantidad de vertices: "; cin>>ver; cout<<'\n';
    cout<<"Ingresa la cantidad de lineas: "; cin>>lin; cout<<'\n';
    cout<<"La Grafica es no dirigida(1) o dirigida(2)?: "; cin>>tipoGrafica; cout<<'\n';
    Mat grados2 (ver,vector<ll>(2,0));
    vector<ll> grados1 (ver,0);
    Mat incidencia(ver, vector<ll>(lin,0));
    Mat adyacencia(ver,vector<ll>(ver,0));
    Mat accesibilidad(ver, vector<ll>(ver,0));
    switch(tipoGrafica){
        
        case 1:

            cout<<"Ingrese la relacion de cada linea (ambos vertices separados por un espacio) :"<<'\n';        
            for(int i=0; i<lin; i++){
                char aux1, aux2;
                cout<<"Linea "<<i+1<<" :"<<'\n';
                cin>>aux1>>aux2;
                if(!vertices.contains(aux1)){
                    if(vertices.size()<ver){
                        indice[aux1] = vertices.size();
                        vertices.insert(aux1);
                    }
                    else {
                        cout<<"Numero de vertices excedido("<<ver<<")"<<'\n';
                        return 0;
                    }
                }
            if(!vertices.contains(aux2)){
                    if(vertices.size()<ver){
                        indice[aux2] = vertices.size();
                        vertices.insert(aux2);
                    }
                    else {
                        cout<<"Numero de vertices excedido("<<ver<<")"<<'\n';
                        return 0;
                    }
                }
                incidencia[indice[aux1]][i] = 1;
                incidencia[indice[aux2]][i] = 1;
                adyacencia[indice[aux1]][indice[aux2]]=1;
                adyacencia[indice[aux2]][indice[aux1]]=1;

            } 

            while (vertices.size() != ver) {
                if (!vertices.contains(aux)) {
                    indice[aux] = vertices.size();
                    vertices.insert(aux);
                }
                aux++;
            }
            cout<<'\n';
            

        break;

        case 2:
            cout<<"Ingrese la relacion de cada linea (ambos vertices separados por un espacio) :"<<'\n';
            cout<<"*Primero el vertice del que sale la linea y depues en el que incide*"<<'\n';
            for(int i=0; i<lin; i++){
                char aux1, aux2;
                cout<<"Linea "<<i+1<<" :"<<'\n';
                cin>>aux1>>aux2;
                if(!vertices.contains(aux1)){
                    if(vertices.size()<ver){
                        indice[aux1] = vertices.size();
                        vertices.insert(aux1);
                    }
                    else {
                        cout<<"Numero de vertices excedido("<<ver<<")"<<'\n';
                        return 0;
                    }
                }
            if(!vertices.contains(aux2)){
                    if(vertices.size()<ver){
                        indice[aux2] = vertices.size();
                        vertices.insert(aux2);
                    }
                    else {
                        cout<<"Numero de vertices excedido("<<ver<<")"<<'\n';
                        return 0;
                    }
                }
                if(aux1==aux2)incidencia[indice[aux1]][i] = 2;
                else {
                    incidencia[indice[aux1]][i] = 1;
                    incidencia[indice[aux2]][i] = -1;
                }
                adyacencia[indice[aux1]][indice[aux2]] = 1;
            }
            while (vertices.size() != ver) {
                if (!vertices.contains(aux)) {
                    indice[aux] = vertices.size();
                    vertices.insert(aux);
                }
                aux++;
            }
            cout<<'\n';
            
            break;


        default:

            cout<< "Opcion no valida";
            return 0;
    }
    long long exp = (ver * ver + ver)/2;
    for(int i=0; i<exp; i++){
        accesibilidad = adyacencia;
        accesibilidad = suma(accesibilidad,multiplicar(accesibilidad, adyacencia, mod));
    }

    vector<char> v(vertices.begin(), vertices.end());
    reverse(v.begin(),v.end());
    cout<<"Las matrices correspondientes son las siguientes :\n"<<'\n';
    cout<<"Mat incidencia"<<'\n'<<"* ";
    if(lin==0)cout<<"No hay lineas por lo que no hay matriz de adyacencia";
    for(int i=0; i<lin; i++)cout<<i+1<<" ";
    cout<<'\n';
    for(int i=0; i<ver; i++){
        for(int j=0; j<lin; j++){
            if(j==0)cout<<v[i]<<" ";
            if(incidencia[i][j]==2)cout<<masmenos<<" ";
            else cout<<incidencia[i][j]<<" ";
            if(tipoGrafica==2){
                if(incidencia[i][j]==2){
                    grados2[i][0]+=1;
                    grados2[i][1]+=1;
                }
                else{
                    if(incidencia[i][j]==-1) grados2[i][0]+=1;
                    else if(incidencia[i][j]==1)grados2[i][1]+=1;
                }
            }
            else grados1[i]+= incidencia[i][j];
        }
        cout <<'\n'; 
    }
    cout<<'\n';

    cout<<"Mat adyacencia"<<'\n'<<"* ";
    for(int i=0; i<ver; i++)cout<<v[i]<<" ";
    cout<<'\n';
    for(int i=0; i<ver; i++){
        for(int j=0; j<ver; j++){
            if(j==0)cout<<v[i]<<" ";
            cout<<adyacencia[i][j]<<" ";
        }
        cout <<'\n'; 
    }
    cout<<'\n';

    cout<<"Mat accesibilidad"<<'\n'<<"* ";
    for(int i=0; i<ver; i++)cout<<v[i]<<" ";
    cout<<'\n';
    for(int i=0; i<ver; i++){
        for(int j=0; j<ver; j++){
            if(j==0)cout<<v[i]<<" ";
            if(accesibilidad[i][j] != 0) cout<<"+"<<" ";
            else cout<<accesibilidad[i][j]<<" ";
        }
        cout <<'\n'; 
    }
    cout<<'\n';

    cout<<"Informacion de los vertices: "<<'\n';
    if(tipoGrafica==1){ 
        cout<<"Grado: \n";
        for(int i=0; i<ver; i++){
            cout<<v[i]<<": "<<grados1[i]<<'\n';
        }
        cout<<'\n';
        
        cout<<"Caracteristicas: \n";
        for(int i=0; i<ver; i++){
            if(grados1[i]==0)cout<<v[i]<<" Es aislado\n";
            else if(grados1[i]==1)cout<<v[i]<<" Es colgante\n"; 
        }
    }
    else{
        cout<<"Grado: \n";
        cout<<"  Int Ext\n";
        for(int i=0; i<ver; i++){
            cout<<v[i]<<": ";
            for(int j=0; j<2; j++){
                cout<<grados2[i][j]<<"    ";
            }
            cout<<'\n';
        }
        cout<<"Caracteristicas: \n";
        for(int i=0; i<ver; i++){
            if(grados2[i][0]==0 && grados2[i][1]==0)cout<<v[i]<<" Es aislado\n";
            else if(grados2[i][1]==0)cout<<v[i]<<" Es final\n";
            else if(grados2[i][0]==0)cout<<v[i]<<" Es inicial\n";
        }
    }
    cout<<'\n';

    cout<<"Informacion de lineas: \n";
    map<vector<ll>, vector<int>> paralelas;

    for(int j = 0; j < lin; j++){
        vector<ll> columna(ver);
        for(int i = 0; i < ver; i++){
            columna[i] = incidencia[i][j];
        }        
        paralelas[columna].push_back(j);
    }

    for(const auto& [col, indices] : paralelas){
        if(indices.size() > 1){
            simple=false;
            cout<<"Las lineas ";
            for(int ind : indices)
                cout << ind+1 << " ";
            cout<<"son paralelas";
            cout << '\n';
        }
    }    
    if(tipoGrafica==2){
        for(int i=0; i<ver; i++){
            for(int j=0; j<lin; j++){
                if(incidencia[i][j]==2){
                    simple=false;
                    cout<<"La linea "<<j+1<<" es un bucle\n";
                }
            }
        }
    }
    else{
        for(int j=0; j<lin; j++){
            int sum=0;
            for(int i=0; i<ver; i++){
                sum+=incidencia[i][j];
            }
            if(sum==1){
                simple=false;
                cout<<"La linea "<<j+1<<" es un bucle\n";
            }
        }
    }
    if(tipoGrafica==1)
    for(int i=0; i<ver; i++){
        if(grados1[i]==2){
            for(int j=0; j<lin; j++){
                if(incidencia[i][j]==1){
                    vector<ll> colum(ver);
                    for(int h=0; h<ver; h++){
                        colum[h]=incidencia[h][j];
                    }
                    if(paralelas[colum].size()==1){
                        for(int z=j+1; z<ver; z++){
                            if(incidencia[i][z]==1){
                                cout<<"Las lineas "<<j+1<<" y "<<z+1<<" son en serie\n";
                                z=ver;
                            }
                        }
                    }
                }
            }
        }
    }
    
    cout<<"\nClasificacion de la Grafica: \n";
    if(!simple)cout<<"La Grafica es General\n";
    else cout<<"La Grafica es simple\n";
    if(lin==0)cout<<"La Grafica es nula\n";
    for(int i=0; i<ver; i++){
        for(int j=0; j<ver; j++){
            if(accesibilidad[i][j]==0){
                conectada=false;
                goto salir;
            }
        }
    }
    salir:
    if(conectada)cout<<"La Grafica es conectada\n";
    else cout<<"La Grafica es desconectada\n";

    if(tipoGrafica==1){
        for(int i=1; i<ver; i++){
            if(grados1[i]!=grados1[i-1]){
                regular=false;
                goto salir1;
            } 
        }
    }
    else{
        ll grad=grados2[0][0];
        for(int i=0; i<ver; i++){
            if(grad!=grados2[i][0] || grad!=grados2[i][1]){
                regular=false;
                goto salir1;
            }   
        }
    }

    salir1:
    if(regular)cout<<"La Grafica es regular\n";
    
    if(tipoGrafica==1){
        if(simple && lin==(ver*(ver-1)/2)){
            cout<<"La Grafica es completa\n";
        }
    }
    else{
        if(simple && lin==(ver*(ver-1))){
            cout<<"La Grafica es completa\n";
        }
    }



    return 0;
}
