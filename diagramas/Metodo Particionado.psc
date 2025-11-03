Proceso ParticionadoDeMatrices
	
    Particionado
	
FinProceso


SubProceso Particionado
    Definir n, P, i, j Como Entero
    Definir detA, detA22 Como Real
	
    Escribir "Metodo de Particionado de Matrices"
    Escribir ""
    Escribir "Ingrese el tamano de la matriz (n >= 4): "
    Leer n
	
    P <- n - 3
    Si n < 4 Entonces
        Escribir ""
        Escribir "La cantidad de variables no es suficiente para aplicar este metodo."
    FinSi

    Dimension mat[n,n]
    Dimension ind[n]
	
    Escribir ""
    Escribir "Ingrese los elementos de su matriz de coeficientes (", n, " x ", n, "):"
    Para i <- 0 Hasta n-1
        Escribir "Fila ", i+1, ":"
        Para j <- 0 Hasta n-1
            Leer mat[i,j]
        FinPara
    FinPara
	
    detA <- Determinante(mat, n)
    Si Abs(detA) < 1*e-12 Entonces
        Escribir ""
        Escribir "La matriz no tiene solucion unica (det=0). No se puede aplicar el metodo."
    FinSi
	
    Escribir ""
    Escribir "Ingrese los elementos del vector independiente (", n, "):"
    Para i <- 0 Hasta n-1
        Escribir "Independiente ", i+1, ": "
        Leer ind[i]
    FinPara
	
    Para i <- 0 Hasta n-1
        Si Abs(mat[i,i]) < 1*e-12 Entonces
            Escribir ""
            Escribir "La matriz contiene un 0 en la diagonal principal (posicion ", i, ",", i, ")."
        FinSi
    FinPara
	
    // A22 (3x3) y su inversa
    Dimension A22[3,3]
    Para i <- 0 Hasta 2
        Para j <- 0 Hasta 2
            A22[i,j] <- mat[i+P,j+P]
        FinPara
    FinPara
	
    detA22 <- Determinante(A22, 3)
    Si Abs(detA22) < 1*e-12 Entonces
        Escribir ""
        Escribir "El bloque A22 no es invertible. No se puede aplicar el metodo."
    FinSi
	
    // Identidad 3x3
    Definir dummyK Como Entero
    dummyK <- 0
    Dimension I3[3,3]
    I3 <- Identidad(3)
	
    Definir A22inv Como Real
    Dimension A22inv[3,3]
    A22inv <- Inversa(A22, I3, dummyK, 3)
	
    // Construccion de D, E, C, F
    Definir D, E, C, F Como Real
	
    D <- getD(mat, n, n, P, A22inv, 3, 3)          // D es P x 3
    E <- getE(mat, n, n, P, A22inv, 3, 3)          // E es 3 x P
    C <- getC(mat, n, n, D, P, 3)                  // C es P x P
    F <- getF(D, P, 3, C, P, P, E, 3, P, A22inv, 3, 3)
	
    // Mostrar inversa por bloques y resolver el sistema
    resultado(D, P, 3, E, 3, P, C, P, P, F, 3, 3, ind, n, P)
	
FinSubProceso


// -------------------- UTILIDADES DE MATRICES --------------------

Funcion M <- Identidad(n)
    Definir i, j Como Entero
    Dimension M[n,n]
    Para i <- 0 Hasta n-1
        Para j <- 0 Hasta n-1
            Si i = j Entonces
                M[i,j] <- 1
            SiNo
                M[i,j] <- 0
            FinSi
        FinPara
    FinPara
FinFuncion


// Inversa por Gauss-Jordan (firma compatible con la pedida)
// mat: n x n, iden: n x n (se ignora su contenido inicial; se recalcula identidad)
// k no se usa (solo para compatibilidad de firma)
Funcion Inv <- Inversa(mat, iden, k, n)
    Definir i, j, r, c Como Entero
    Definir piv, tmp Como Real
	
    Dimension A[n,n]
    Dimension Inv[n,n]
	
    // Copiar mat en A y construir identidad en Inv
    Para i <- 0 Hasta n-1
        Para j <- 0 Hasta n-1
            A[i,j] <- mat[i,j]
        FinPara
    FinPara
    Inv <- Identidad(n)
	
    Para c <- 0 Hasta n-1
		
        // Pivoteo parcial: buscar mejor fila desde c
        Definir filaPiv Como Entero
        filaPiv <- c
        Definir maxVal Como Real
        maxVal <- Abs(A[c,c])
        Para r <- c+1 Hasta n-1
            Si Abs(A[r,c]) > maxVal Entonces
                maxVal <- Abs(A[r,c])
                filaPiv <- r
            FinSi
        FinPara
		
        Si maxVal < 1*e-15 Entonces
            // Matriz no invertible
            // Retornar matriz de ceros para señalar fallo
            Para i <- 0 Hasta n-1
                Para j <- 0 Hasta n-1
                    Inv[i,j] <- 0
                FinPara
            FinPara
        FinSi
		
        // Intercambiar filas si es necesario
        Si filaPiv <> c Entonces
            Para j <- 0 Hasta n-1
                tmp <- A[c,j]
                A[c,j] <- A[filaPiv,j]
                A[filaPiv,j] <- tmp
				
                tmp <- Inv[c,j]
                Inv[c,j] <- Inv[filaPiv,j]
                Inv[filaPiv,j] <- tmp
            FinPara
        FinSi
		
        // Normalizar fila pivote
        piv <- A[c,c]
        Para j <- 0 Hasta n-1
            A[c,j] <- A[c,j] / piv
            Inv[c,j] <- Inv[c,j] / piv
        FinPara
		
        // Eliminar en las demas filas
        Para r <- 0 Hasta n-1
            Si r <> c Entonces
                tmp <- A[r,c]
                Para j <- 0 Hasta n-1
                    A[r,j] <- A[r,j] - tmp * A[c,j]
                    Inv[r,j] <- Inv[r,j] - tmp * Inv[c,j]
                FinPara
            FinSi
        FinPara
    FinPara
FinFuncion


// Multiplicacion de matrices: C = A(fA x cA) * B(fB x cB)
Funcion C <- mMUlt(A, fA, cA, B, fB, cB)
    Definir i, j, k Como Entero
    Dimension C[fA,cB]
    // Inicializar
    Para i <- 0 Hasta fA-1
        Para j <- 0 Hasta cB-1
            C[i,j] <- 0
        FinPara
    FinPara
	
    // Producto
    Para i <- 0 Hasta fA-1
        Para j <- 0 Hasta cB-1
            Para k <- 0 Hasta cA-1
                C[i,j] <- C[i,j] + A[i,k] * B[k,j]
            FinPara
        FinPara
    FinPara
FinFuncion


// Determinante por eliminacion gaussiana con pivoteo parcial
Funcion det <- Determinante(A, n)
    Definir i, j, r, c Como Entero
    Definir piv, factor Como Real
    Definir signo Como Entero
    signo <- 1
	
    Dimension M[n,n]
    Para i <- 0 Hasta n-1
        Para j <- 0 Hasta n-1
            M[i,j] <- A[i,j]
        FinPara
    FinPara
	
    Para c <- 0 Hasta n-1
        // Buscar pivote maximo por debajo/igual a c
        Definir filaPiv Como Entero
        Definir maxVal Como Real
        filaPiv <- c
        maxVal <- Abs(M[c,c])
        Para r <- c+1 Hasta n-1
            Si Abs(M[r,c]) > maxVal Entonces
                maxVal <- Abs(M[r,c])
                filaPiv <- r
            FinSi
        FinPara
		
        Si maxVal < 1*e-15 Entonces
            det <- 0
        FinSi
		
        // Intercambio de filas si corresponde
        Si filaPiv <> c Entonces
            Para j <- 0 Hasta n-1
                piv <- M[c,j]
                M[c,j] <- M[filaPiv,j]
                M[filaPiv,j] <- piv
            FinPara
            signo <- -signo
        FinSi
		
        // Eliminar debajo
        Para r <- c+1 Hasta n-1
            factor <- M[r,c] / M[c,c]
            Para j <- c Hasta n-1
                M[r,j] <- M[r,j] - factor * M[c,j]
            FinPara
        FinPara
    FinPara
	
    // Producto de la diagonal (con signo por intercambios)
    det <- signo
    Para i <- 0 Hasta n-1
        det <- det * M[i,i]
    FinPara
FinFuncion


// -------------------- CONSTRUCCION DE BLOQUES --------------------

// D = A12 * A22inv
// A12 es P x 3 (toma columnas P..P+2 de filas 0..P-1)
Funcion D <- getD(mat, fM, cM, P, A22inv, fI, cI)
    Definir i, j Como Entero
    Dimension A12[P,3]
    Para i <- 0 Hasta P-1
        Para j <- 0 Hasta 2
            A12[i,j] <- mat[i,j+P]
        FinPara
    FinPara
    D <- mMUlt(A12, P, 3, A22inv, 3, 3)  // D es P x 3
FinFuncion


// E = A22inv * A21
// A21 es 3 x P (toma filas P..P+2 y columnas 0..P-1)
Funcion E <- getE(mat, fM, cM, P, A22inv, fI, cI)
    Definir i, j Como Entero
    Dimension A21[3,P]
    Para i <- 0 Hasta 2
        Para j <- 0 Hasta P-1
            A21[i,j] <- mat[i+P,j]
        FinPara
    FinPara
    E <- mMUlt(A22inv, 3, 3, A21, 3, P)  // E es 3 x P
FinFuncion


// C = A11 - (D * A21) , donde A11 es P x P y A21 es 3 x P
Funcion C <- getC(mat, fM, cM, D, fD, cD)
    Definir i, j Como Entero
    Definir P Como Entero
    P <- fD  // D es P x 3
	
    Dimension A21[3,P]
    Para i <- 0 Hasta 2
        Para j <- 0 Hasta P-1
            A21[i,j] <- mat[i+P,j]
        FinPara
    FinPara
	
    Dimension A11[P,P]
    Para i <- 0 Hasta P-1
        Para j <- 0 Hasta P-1
            A11[i,j] <- mat[i,j]
        FinPara
    FinPara
	
    Definir M Como Real
    M <- mMUlt(D, P, 3, A21, 3, P)  // P x P
	
    Dimension C[P,P]
    Para i <- 0 Hasta P-1
        Para j <- 0 Hasta P-1
            C[i,j] <- A11[i,j] - M
        FinPara
    FinPara
FinFuncion


// F = A22inv + E * C^{-1} * D
Funcion F <- getF(D, fD, cD, C, fC, cC, E, fE, cE, A22inv, fA22, cA22)
    Definir i, j, dummyK Como Entero
    Definir P Como Entero
    P <- fC    // C es P x P
	
    // C^{-1}
    Dimension IdenC[P,P]
    IdenC <- Identidad(P)
    dummyK <- 0
    Definir Cinv Como Real
    Cinv <- Inversa(C, IdenC, dummyK, P)
	
    // E * C^{-1}
    Definir ECinv Como Real
    ECinv <- mMUlt(E, fE, cE, Cinv, P, P)  // (3 x P) * (P x P) = 3 x P
	
    // (E C^{-1}) * D
    Definir ECinvD Como Real
    ECinvD <- mMUlt(ECinv, 3, P, D, P, 3)  // (3 x P) * (P x 3) = 3 x 3
	
    // F = A22inv + ECinvD
    Dimension F[3,3]
    Para i <- 0 Hasta 2
        Para j <- 0 Hasta 2
            F[i,j] <- A22inv[i,j] + ECinvD
        FinPara
    FinPara
FinFuncion


// -------------------- SALIDA (INVERSA POR BLOQUES) Y SOLUCION --------------------

// Construye la inversa completa por bloques y muestra la solucion del sistema
SubProceso resultado(D, fD, cD, E, fE, cE, C, fC, cC, F, fF, cF, ind, n, P)
    Definir i, j, dummyK Como Entero
	
    // C^{-1}, C^{-1}D, E C^{-1}
    Dimension IdenC[fC,cC]
    IdenC <- Identidad(fC)
    dummyK <- 0
	
    Definir Cinv, CinvD, ECinv Como Real
    Cinv <- Inversa(C, IdenC, dummyK, fC)                  // P x P
    CinvD <- mMUlt(Cinv, fC, cC, D, fD, cD)                // (P x P) * (P x 3) = P x 3
    ECinv <- mMUlt(E, fE, cE, Cinv, fC, cC)                // (3 x P) * (P x P) = 3 x P
	
    // Matriz inversa completa (n x n) por bloques
    Dimension res[n,n]
	
    // Res[0..P-1, 0..P-1] = Cinv
    Para i <- 0 Hasta P-1
        Para j <- 0 Hasta P-1
            res[i,j] <- Cinv
        FinPara
    FinPara
	
    // Res[0..P-1, P..P+2] = - CinvD
    Para i <- 0 Hasta P-1
        Para j <- 0 Hasta 2
            res[i,j+P] <- - CinvD
        FinPara
    FinPara
	
    // Res[P..P+2, 0..P-1] = - ECinv
    Para i <- 0 Hasta 2
        Para j <- 0 Hasta P-1
            res[i+P,j] <- - ECinv
        FinPara
    FinPara
	
    // Res[P..P+2, P..P+2] = F
    Para i <- 0 Hasta 2
        Para j <- 0 Hasta 2
            res[i+P,j+P] <- F[i,j]
        FinPara
    FinPara
	
    // Mostrar inversa
    Escribir ""
    Escribir "Matriz inversa:"
    Para i <- 0 Hasta n-1
        Para j <- 0 Hasta n-1
            Si Abs(res[i,j]) < 1*e-10 Entonces
                Escribir Sin Saltar 0, " "
            SiNo
                Escribir Sin Saltar res[i,j], " "
            FinSi
        FinPara
        Escribir ""
    FinPara
	
    // Solucion x = res * ind
    Definir solucion Como Real
    Dimension solucion[n]
    Para i <- 0 Hasta n-1
        solucion[i] <- 0
        Para j <- 0 Hasta n-1
            solucion[i] <- solucion[i] + res[i,j] * ind[j]
        FinPara
    FinPara
	
    Escribir ""
    Escribir "La solucion del sistema de ecuaciones es:"
    Para i <- 0 Hasta n-1
        Escribir "x", i+1, " = ", solucion[i]
    FinPara
    Escribir ""
FinSubProceso
