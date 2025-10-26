#include <iostream>
#include <vector>
using namespace std;

vector<vector<double>> getSubmatrix(const vector<vector<double>>& mat, int row, int col) {
    int n = mat.size();
    vector<vector<double>> submat(n - 1, vector<double>(n - 1));
    int subi = 0;
    for (int i = 0; i < n; i++) {
        if (i == row) continue;
        int subj = 0;
        for (int j = 0; j < n; j++) {
            if (j == col) continue;
            submat[subi][subj] = mat[i][j];
            subj++;
        }
        subi++;
    }
    return submat;
}

double determinant(const vector<vector<double>>& mat) {
    int n = mat.size();
    
    if (n == 1) return mat[0][0];
    if (n == 2) return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]; // caso base 2x2

    double det = 0;
    for (int j = 0; j < n; j++) {
        vector<vector<double>> submat = getSubmatrix(mat, 0, j);
        det += ( (j % 2 == 0 ? 1 : -1) * mat[0][j] * determinant(submat) );
    }
    return det;
}

vector<vector<double>> identidad(int n) {
    vector<vector<double>>pos(n,vector<double>(n));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j)pos[i][j]=1;
            else pos[i][j]=0;
        }
    }
    return pos;
}

vector<vector<double>> inversa(vector<vector<double>>& mat, vector<vector<double>>& iden, int k, int n) {
    if (k >= n) return iden;
    double factor = mat[k][k];
    for (int j = 0; j < n; j++) {
        mat[k][j] /= factor;
        iden[k][j] /= factor;
    }

    for (int i = 0; i < n; i++) {
        if (i == k) continue;
        factor = mat[i][k];
        for (int j = 0; j < n; j++) {
            mat[i][j] -= factor * mat[k][j];
            iden[i][j] -= factor * iden[k][j];
        }
    }

    return inversa(mat, iden, k + 1, n);
}

vector<vector<double>> mMUlt(const vector<vector<double>>& A, const vector<vector<double>>& B){
    int n=A.size();
    int m=B[0].size();
    vector<vector<double>> C(n,vector<double>(m,0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            for (int k = 0; k < B.size(); k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

