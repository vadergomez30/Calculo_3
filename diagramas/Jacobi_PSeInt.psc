    // Valor absoluto para reales
Funcion rabs <- AbsR(x Por Valor)
        Definir rabs Como Real
        Si x >= 0 Entonces
            rabs <- x
        SiNo
            rabs <- -x
        FinSi
FinFuncion

// Copia de matriz A -> B (n x n)
SubProceso CopiarMatriz(A, n, B)
	Definir i, j Como Entero
	Para i <- 1 Hasta n
		Para j <- 1 Hasta n
			B[i,j] <- A[i,j]
		FinPara
	FinPara
FinSubProceso

// Intercambio de filas en matriz M (1..n)
SubProceso SwapFilas(M, n, f1, f2)
	Definir j Como Entero
	Definir tmp Como Real
	Para j <- 1 Hasta n
		tmp <- M[f1,j]
		M[f1,j] <- M[f2,j]
		M[f2,j] <- tmp
	FinPara
FinSubProceso

// Determinante por eliminación gaussiana con pivoteo parcial
Funcion det <- Determinante(A, n)
	Definir det, maxv, factor Como Real
	Definir i, j, k, p Como Entero
	Dimension T[n,n]
	CopiarMatriz(A, n, T)
	
	Definir swaps Como Entero
	swaps <- 0
	
	Para k <- 1 Hasta n
		// Buscar pivote máximo en columna k, filas k..n
		p <- k
		maxv <- absr (T[k,k]) 
		Para i <- k+1 Hasta n
			Si absr(T[i,k]) > maxv Entonces
				maxv <- absr(T[i,k])
				p <- i
			FinSi
		FinPara
		
		// Si pivote ~ 0, determinante = 0
		Si maxv < EPS Entonces
			det <- 0
		FinSi
		
		// Intercambiar fila k con p si es necesario
		Si p <> k Entonces
			SwapFilas(T, n, k, p)
			swaps <- swaps + 1
		FinSi
		
		// Eliminar entradas debajo del pivote
		Para i <- k+1 Hasta n
			factor <- T[i,k] / T[k,k]
			Para j <- k Hasta n
				T[i,j] <- T[i,j] - factor * T[k,j]
			FinPara
		FinPara
	FinPara
	
	// Producto de la diagonal
	det <- 1
	Para i <- 1 Hasta n
		det <- det * T[i,i]
	FinPara
	
	// Ajuste por número de swaps (cada swap cambia el signo)
	Si (swaps % 2) = 1 Entonces
		det <- -det
	FinSi
FinFuncion



// Verifica dominancia diagonal por filas y no singularidad (det != 0)
Funcion ok <- Supuestos(A, n)
	Definir ok Como Logico
	Definir i, j Como Entero
	Definir diag, suma, d Como Real
	
	Para i <- 1 Hasta n
		diag <- absr(A[i,i])
		suma <- 0
		Para j <- 1 Hasta n
			Si j <> i Entonces
				suma <- suma + absr(A[i,j])
			FinSi
		FinPara
		Si diag < suma Entonces
			ok <- Falso
		FinSi
	FinPara
	
	d <- Determinante(A, n)
	Si absr(d) < EPS Entonces
		ok <- Falso
	SiNo
		ok <- Verdadero
	FinSi
FinFuncion

// Criterio de paro: ||x_nue - x_ant||_inf <= error (componente a componente)
Funcion listo <- criterioParo(x_ant, x_nue, n, error)
	Definir listo Como Logico
	Definir i Como Entero
	listo <- Verdadero
	Para i <- 1 Hasta n
		Si absr(x_nue[i] - x_ant[i]) > error Entonces
			listo <- Falso
		FinSi
	FinPara
FinFuncion

// Una iteración (versión recursiva como en el C++)
SubProceso IterativoJacobi(A, n, b, x_ant, x_nue, error, count Por Referencia)
	Definir i, j Como Entero
	Definir res Como Real
	
	Para i <- 1 Hasta n
		res <- b[i]
		Para j <- 1 Hasta n
			Si j <> i Entonces
				res <- res - A[i,j] * x_ant[j]
			FinSi
		FinPara
		x_nue[i] <- res / A[i,i]
	FinPara
	
	count <- count + 1
	
	Escribir ""
	Escribir "Iteracion ", count
	Escribir "Solucion X: " Sin Saltar
	Para i <- 1 Hasta n
		Escribir x_nue[i], " " Sin Saltar
	FinPara
	Escribir ""
	
	Si criterioParo(x_ant, x_nue, n, error) Entonces
		Escribir "Convergencia alcanzada despues de ", count, " iteraciones."
	SiNo
		// x_ant <- x_nue
		Para i <- 1 Hasta n
			x_ant[i] <- x_nue[i]
		FinPara
		IterativoJacobi(A, n, b, x_ant, x_nue, error, count)
	FinSi
FinSubProceso

// ================== Programa principal (equivale a Jacobi()) ==================
Algoritmo Jacobi_PSeInt
	Definir EPS Como Real
    EPS <- 0.000000000001	
Definir n, i, j, count Como Entero
Definir error Como Real
Escribir "Ingrese el tamano n de la matriz (n x n):"
Leer n
Escribir "Ingrese el error tolerable:"
Leer error

Dimension A[n,n]
Dimension b[n]
Dimension x_ant[n]
Dimension x_nue[n]

Escribir "Ingrese la matriz A (", n, "x", n, "), por filas:"
Para i <- 1 Hasta n
	Para j <- 1 Hasta n
		Leer A[i,j]
	FinPara
FinPara

Escribir "Ingrese el vector b (", n, " valores):"
Para i <- 1 Hasta n
	Leer b[i]
FinPara

Si No Supuestos(A, n) Entonces
	Escribir "La matriz no cumple los supuestos necesarios (dominancia diagonal y no singularidad)."
FinSi

Para i <- 1 Hasta n
	x_ant[i] <- 0
	x_nue[i] <- 0
FinPara

count <- 0
IterativoJacobi(A, n, b, x_ant, x_nue, error, count)

FinAlgoritmo
