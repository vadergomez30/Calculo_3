Algoritmo MetodoGauss
	// --- Declaraci�n de variables ---
	Definir n, i, j, k Como Entero
	Definir factor, suma, det Como Real
	Dimensionar mat(10,10)
	Dimensionar b(10)
	Dimensionar Res(10)
	// --- Entrada de datos ---
	Escribir 'Ingrese el tama�o de la matriz cuadrada (n x n): '
	Leer n
	Escribir 'Ingrese los elementos de la matriz fila por fila:'
	Para i<-1 Hasta n Hacer
		Para j<-1 Hasta n Hacer
			Escribir 'Elemento [', i, '][', j, ']: '
			Leer mat[i,j]
		FinPara
	FinPara
	Escribir ''
	Escribir 'La matriz ingresada es:'
	Para i<-1 Hasta n Hacer
		Para j<-1 Hasta n Hacer
			Escribir mat[i,j], '   'Sin Saltar
		FinPara
		Escribir ''
	FinPara
	// --- C�lculo del determinante (hasta 3x3 para ejemplo pr�ctico) ---
	Si n=1 Entonces
		det <- mat[1,1]
	SiNo
		Si n=2 Entonces
			det <- mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1]
		SiNo
			// Determinante para matriz 3x3
			det <- mat[1,1]*(mat[2,2]*mat[3,3]-mat[2,3]*mat[3,2])-mat[1,2]*(mat[2,1]*mat[3,3]-mat[2,3]*mat[3,1])+mat[1,3]*(mat[2,1]*mat[3,2]-mat[2,2]*mat[3,1])
		FinSi
	FinSi
	Escribir ''
	Escribir 'El determinante es: ', det
	Si det=0 Entonces
		Escribir 'No es posible aplicar el m�todo de Gauss (determinante = 0).'
	SiNo
		// --- Lectura del vector b ---
		Escribir ''
		Escribir 'Ingrese el vector de t�rminos independientes b (', n, 'x1):'
		Para i<-1 Hasta n Hacer
			Escribir 'b[', i, ']: '
			Leer b[i]
		FinPara
		// --- Eliminaci�n Gaussiana ---
		Para k<-1 Hasta n-1 Hacer
			Para i<-k+1 Hasta n Hacer
				factor <- mat[i,k]/mat[k,k]
				Para j<-k Hasta n Hacer
					mat[i,j]<-mat[i,j]-factor*mat[k,j]
				FinPara
				b[i] <- b[i]-factor*b[k]
			FinPara
		FinPara
		Escribir ''
		Escribir '--- Matriz transformada (forma triangular superior) ---'
		Para i<-1 Hasta n Hacer
			Para j<-1 Hasta n Hacer
				Si ABS(mat[i,j])<0.0000001 Entonces
					Escribir '0   'Sin Saltar
				SiNo
					Escribir mat[i,j], '   'Sin Saltar
				FinSi
			FinPara
			Escribir ''
		FinPara
		// --- Sustituci�n regresiva ---
		Para i<-n Hasta 1 Con Paso -1 Hacer
			suma <- 0
			Para j<-i+1 Hasta n Hacer
				suma <- suma+mat[i,j]*Res[j]
			FinPara
			Res[i] <- (b[i]-suma)/mat[i,i]
		FinPara
		// --- Resultados ---
		Escribir ''
		Escribir 'Soluciones del sistema:'
		Para i<-1 Hasta n Hacer
			Escribir 'x[', i, '] = ', Res[i]
		FinPara
	FinSi
FinAlgoritmo
