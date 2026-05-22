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
// Agregar después de las funciones suma() y multiplicar()

bool esPuente(Mat& adj, int u, int v, int ver) {
    // Contar vértices alcanzables desde u
    vector<bool> visitado(ver, false);
    queue<int> q;
    q.push(u);
    visitado[u] = true;
    int count1 = 0;
    while (!q.empty()) {
        int node = q.front(); q.pop();
        count1++;
        for (int i = 0; i < ver; i++)
            if (adj[node][i] > 0 && !visitado[i]) {
                visitado[i] = true;
                q.push(i);
            }
    }

    // Quitar la arista y contar de nuevo
    adj[u][v]--;
    adj[v][u]--;
    fill(visitado.begin(), visitado.end(), false);
    q.push(u);
    visitado[u] = true;
    int count2 = 0;
    while (!q.empty()) {
        int node = q.front(); q.pop();
        count2++;
        for (int i = 0; i < ver; i++)
            if (adj[node][i] > 0 && !visitado[i]) {
                visitado[i] = true;
                q.push(i);
            }
    }

    // Restaurar la arista
    adj[u][v]++;
    adj[v][u]++;

    // Es puente si al quitarla hay menos vértices alcanzables
    return count2 < count1;
}
void fleury(Mat adyacencia, vector<char>& v, int ver, int tipoGrafica, vector<ll>& grados1) {
    
    // Usar grados1 ya calculado en lugar de recalcular desde adyacencia
    int sumaGrados = 0;
    int primerImpar = -1, segundoImpar = -1;
    int totalImpares = 0;

    for (int i = 0; i < ver; i++) {
        sumaGrados += grados1[i];
        if (grados1[i] % 2 != 0) {
            totalImpares++;
            if (primerImpar == -1) primerImpar = i;
            else if (segundoImpar == -1) segundoImpar = i;
        }
    }
    cout<<totalImpares<<" vertices de grado impar\n";
    // Determinar tipo y vertice de inicio
    int inicio;
    if (totalImpares == 0) {
        // 0 vertices impares = circuito euleriano
        inicio = -1;
        for (int i = 0; i < ver; i++)
            if (grados1[i] > 0) { inicio = i; break; }
        if (inicio == -1) { cout << "No hay aristas\n"; return; }
        cout << "Circuito Euleriano: ";
    }
    else if (totalImpares == 2) {
        // exactamente 2 vertices impares = camino euleriano
        inicio = primerImpar;
        cout << "Camino Euleriano (Unicursal): ";
    }
    else {
        // 4, 6, 8... vertices impares = no existe camino euleriano
        cout << "No existe camino ni circuito euleriano ("<<totalImpares<<" vertices de grado impar)\n";
        return;
    }
// Algoritmo de Fleury
    Mat adjCopia = adyacencia;
    vector<int> camino;
    int actual = inicio;
    camino.push_back(actual);

    while (true) {
        int siguiente = -1;
        for (int i = 0; i < ver; i++) {
            if (adjCopia[actual][i] > 0) {
                if (!esPuente(adjCopia, actual, i, ver)) {
                    siguiente = i;
                    break;
                }
                if (siguiente == -1) siguiente = i;
            }
        }

        if (siguiente == -1) break;

        camino.push_back(siguiente);
        adjCopia[actual][siguiente]--;
        adjCopia[siguiente][actual]--;
        actual = siguiente;
    }

    // ✅ Si es circuito euleriano, cerrar el ciclo regresando al inicio
    if (totalImpares == 0 && camino.front() != camino.back())
        camino.push_back(inicio);

    // Imprimir camino
    for (int i = 0; i < (int)camino.size(); i++) {
        cout << v[camino[i]];
        if (i != (int)camino.size() - 1) cout << " -> ";
    }
    cout << '\n';
}
// BFS para encontrar camino más corto entre dos vertices
vector<int> bfs(Mat& adj, int inicio, int fin, int ver) {
    vector<int> padre(ver, -1);
    vector<bool> visitado(ver, false);
    queue<int> q;
    q.push(inicio);
    visitado[inicio] = true;

    while (!q.empty()) {
        int actual = q.front(); q.pop();
        if (actual == fin) break;
        for (int i = 0; i < ver; i++) {
            if (adj[actual][i] > 0 && !visitado[i]) {
                visitado[i] = true;
                padre[i] = actual;
                q.push(i);
            }
        }
    }

    // Reconstruir camino
    vector<int> camino;
    for (int actual = fin; actual != -1; actual = padre[actual])
        camino.push_back(actual);
    reverse(camino.begin(), camino.end());

    // Si no hay camino
    if (camino[0] != inicio) return {};
    return camino;
}

void eulerizar(Mat& adyacencia, vector<char>& v, int ver, vector<ll>& grados1, vector<pair<char,char>>& aristas) {

    vector<int> impares;
    for (int i = 0; i < ver; i++)
        if (grados1[i] % 2 != 0)
            impares.push_back(i);

    if (impares.empty()) {
        cout << "La grafica ya es euleriana, no necesita eulerización\n";
        return;
    }

    cout << "\nEulerizacion:\n";
    cout << "Vertices de grado impar: ";
    for (int i : impares) cout << v[i] << " ";
    cout << '\n';

    cout << "Aristas a duplicar:\n";

    for (int i = 0; i < (int)impares.size(); i += 2) {
        int u = impares[i];
        int w = impares[i + 1];

        vector<int> camino = bfs(adyacencia, u, w, ver);

        if (camino.empty()) {
            cout << "No hay camino entre " << v[u] << " y " << v[w] << '\n';
            continue;
        }

        for (int j = 0; j < (int)camino.size() - 1; j++) {
            int a = camino[j], b = camino[j + 1];
            adyacencia[a][b]++;
            adyacencia[b][a]++;
            grados1[a]++;
            grados1[b]++;
            cout << "  " << v[a] << " -- " << v[b] << '\n';

            // ✅ Agregar la arista duplicada al vector
            aristas.push_back({v[a], v[b]});
        }
    }

    cout << "\nGrados despues de eulerizar:\n";
    for (int i = 0; i < ver; i++)
        cout << v[i] << ": " << grados1[i] << '\n';

    cout << "\nRecorrido Euleriano tras eulerizar:\n";
    fleury(adyacencia, v, ver, 1, grados1);
}
void generarDOT(const vector<char>& vertices,const vector<pair<char,char>>& aristas,int tipoGrafica){
    ofstream archivo("grafo.dot");

    if(tipoGrafica == 1)
        archivo << "graph G {\n";
    else
        archivo << "digraph G {\n";

    archivo << "    node [shape=circle];\n";

    for(char v : vertices){
        archivo << "    " << v << ";\n";
    }

    for(auto arista : aristas){

        char a = arista.first;
        char b = arista.second;

        if(tipoGrafica == 1)
            archivo << "    " << a << " -- " << b << ";\n";
        else
            archivo << "    " << a << " -> " << b << ";\n";
    }

    archivo << "}\n";

    archivo.close();
}

int main() {
    int ver, lin, tipoGrafica;
    unordered_set<char> vertices;
    unordered_map<char,int> indice;
    vector<char> v;
    int mod = 1e9 + 7; 
    char masmenos = 241;
    char aux;
    aux = 'A';
    bool simple=true, conectada=true, regular=true, simetrica=true, balanceada=true;
    cout<<"Ingresa la cantidad de vertices: "; cin>>ver; cout<<'\n';
    if(ver <= 0){
        cout<<"Numero de vertices debe ser mayor que 0"<<'\n';
        return 0;
    }
    size_t verCount = static_cast<size_t>(ver);
    cout<<"Ingresa la cantidad de lineas: "; cin>>lin; cout<<'\n';
    cout<<"La Grafica es no dirigida(1) o dirigida(2)?: "; cin>>tipoGrafica; cout<<'\n';
    Mat grados2 (ver,vector<ll>(2,0));
    vector<ll> grados1 (ver,0);
    vector<pair<char,char>> aristas;
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
                if(!vertices.count(aux1)){
                    if(vertices.size()<verCount){
                        indice[aux1] = vertices.size();
                        vertices.insert(aux1);
                        v.push_back(aux1);
                    }
                    else {
                        cout<<"Numero de vertices excedido("<<ver<<")"<<'\n';
                        return 0;
                    }
                }
            if(!vertices.count(aux2)){
                    if(vertices.size()<verCount){
                        indice[aux2] = vertices.size();
                        vertices.insert(aux2);
                        v.push_back(aux2);
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
                aristas.push_back({aux1, aux2});
            } 

            while (vertices.size() != verCount) {
                if (!vertices.count(aux)) {
                    indice[aux] = vertices.size();
                    vertices.insert(aux);
                    v.push_back(aux);
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
                if(!vertices.count(aux1)){
                    if(vertices.size()<verCount){
                        indice[aux1] = vertices.size();
                        vertices.insert(aux1);
                        v.push_back(aux1);
                    }
                    else {
                        cout<<"Numero de vertices excedido("<<ver<<")"<<'\n';
                        return 0;
                    }
                }
            if(!vertices.count(aux2)){
                    if(vertices.size()<verCount){
                        indice[aux2] = vertices.size();
                        vertices.insert(aux2);
                        v.push_back(aux2);
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
                aristas.push_back({aux1, aux2});
            }
            while (vertices.size() != verCount) {
                if (!vertices.count(aux)) {
                    indice[aux] = vertices.size();
                    vertices.insert(aux);
                    v.push_back(aux);
                }
                aux++;
            }
            cout<<'\n';

            break;


        default:

            cout<< "Opcion no valida";
            return 0;
    }
    long long exp = ((ver * ver + ver)/2)+10;
    accesibilidad = adyacencia;
    for(int i=0; i<exp; i++){
        accesibilidad = suma(accesibilidad,multiplicar(accesibilidad, adyacencia, mod));
    }

    cout<<'\n';
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
                        for(int z=j+1; z<lin; z++){
                            if(incidencia[i][z]==1){
                                cout<<"Las lineas "<<j+1<<" y "<<z+1<<" son en serie\n";
                                z=lin;
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
            if(grados2[i][0]!=grados2[i][1]){
                balanceada=false;
                regular=false;
                goto salir1;
            }
            if(grad!=grados2[i][0] || grad!=grados2[i][1]){
                regular=false;
            }   
        }
    }

    salir1:
    if(regular)cout<<"La Grafica es regular\n";

    if(tipoGrafica==1){
        if(ver>1 && simple && lin==(ver*(ver-1)/2)){
            cout<<"La Grafica es completa\n";
        }
    }
    else{
        if(ver>1 && simple && lin==(ver*(ver-1))){
            cout<<"La Grafica es completa\n";
        }
    }

    if(tipoGrafica==1 && ver>1 && conectada && lin==ver-1)cout<<"La Grafica es un Arbol\n";

    if(tipoGrafica==2)
    for(int i=0; i<ver; i++){
        for(int j=0; j<ver; j++){
            if(adyacencia[i][j]!=adyacencia[j][i]){
                simetrica=false;
                goto salir2;
            }
        }
    }
    
    salir2: 
    if(lin>0 && tipoGrafica==2 && simetrica)cout<<"La Grafica es simetrica\n";

    if(tipoGrafica==2 && balanceada)cout<<"La Grafica es balanceda\n";

    if(tipoGrafica==1 && conectada){
        int gradpar=0;
        for(int i=0; i<ver; i++){
            if(grados1[i]%2==0){
                gradpar++;
            }
        }
        if(gradpar==ver)cout<<"La Grafica es euleriana"<<'\n';
        else if(gradpar==ver-2)cout<<"La Grafica es unicursal"<<'\n';
    }
    if(tipoGrafica==2 && balanceada)cout<<"La Grafica es euleriana"<<'\n';
    else{
        int gradigual=0;
        for(int i=0; i<ver; i++){
            if(grados2[i][0]==grados2[i][1])gradigual++;
        }
        if(gradigual==ver-2)cout<<"La Grafica es unicursal"<<'\n';
    }
    
    generarDOT(v, aristas, tipoGrafica);
    system("dot -Tpng grafo.dot -o grafo.png");
    cout << "\nSe genero la imagen: grafo.png\n";
    system("pause");
    
    cout << "\nRecorrido Euleriano:\n";
    if (tipoGrafica == 1 && conectada) {
        int totalImpares = 0;
        for (int i = 0; i < ver; i++)
            if (grados1[i] % 2 != 0) totalImpares++;

        if (totalImpares == 0 || totalImpares == 2) {
            fleury(adyacencia, v, ver, tipoGrafica, grados1);
        } else {
            cout << "La grafica no es euleriana ni unicursal\n";
            eulerizar(adyacencia, v, ver, grados1, aristas);
        }
    }

    generarDOT(v, aristas, tipoGrafica);

    system("dot -Tpng grafo.dot -o grafo.png");

    cout << "\nSe genero la imagen: grafo.png\n";

    return 0;
}
