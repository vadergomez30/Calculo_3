Algoritmo MetodoIntercambio_4x4
    Definir i, j, it, l, k Como Entero
    Definir best, alk, EPS Como Real
    EPS <- 0.000000000001
    
    Dimension A[4,4]
    Dimension b[4]
    Dimension x[4]
    // ---- etiquetas para metodo de intercambio ----
    Dimension rowLab[4]
    Dimension colLab[4]
    Para i <- 1 Hasta 4
        rowLab[i] <- i
        colLab[i] <- i
    FinPara
    
    Escribir "Captura A (4x4) por filas:"
    Para i <- 1 Hasta 4
        Para j <- 1 Hasta 4
            Escribir "A[", i, ",", j, "] = "
            Leer A[i,j]
        FinPara
    FinPara
    
    Escribir "Captura b (4 componentes):"
    Para i <- 1 Hasta 4
        Escribir "b[", i, "] = "
        Leer b[i]
    FinPara
    
    Para it <- 1 Hasta 4
        best <- 0
        l <- 1
        k <- 1
        
        Para i <- 1 Hasta 4
            Para j <- 1 Hasta 4
                Si Abs(A[i,j]) > best Entonces
                    best <- Abs(A[i,j])
                    l <- i
                    k <- j
                FinSi
            FinPara
        FinPara
        
        alk <- A[l,k]
        Si Abs(alk) < EPS Entonces
            Escribir "Fallo: pivote ~ 0 (matriz singular)."
			//FinAlgoritmo
		FinSi
		
		Escribir ""
		Escribir "Iteracion ", it, ": pivote A[", l, ",", k, "] = ", alk
		
		Para i <- 1 Hasta 4
			Si i <> l Entonces
				Para j <- 1 Hasta 4
					Si j <> k Entonces
						A[i,j] <- A[i,j] - (A[i,k] * A[l,j]) / alk
					FinSi
				FinPara
			FinSi
		FinPara
		
		Para i <- 1 Hasta 4
			Si i <> l Entonces
				A[i,k] <- A[i,k] / alk
			FinSi
		FinPara
		
		Para j <- 1 Hasta 4
			Si j <> k Entonces
				A[l,j] <- -A[l,j] / alk
			FinSi
		FinPara
		
		A[l,k] <- 1 / alk
		
		Escribir "Matriz despues de iteracion ", it, ":"
		Para i <- 1 Hasta 4
			Para j <- 1 Hasta 4
				Escribir Sin Saltar A[i,j], " "
			FinPara
			Escribir ""
		FinPara
		
		// ---- intercambio de etiquetas (NO cambia la logica numerica) ----
		//Definir tmp Como Entero
		tmp <- rowLab[l]
		rowLab[l] <- colLab[k]
		colLab[k] <- tmp
		
	FinPara  // fin del bucle it=1..4
	
	Escribir ""
	Escribir "Solucion x = A^{-1}b (respetando etiquetas):"
	// A contiene la "inversa" con filas etiquetadas por rowLab y columnas por colLab.
	// x[rowLab[i]] = sum_j A[i,j] * b[colLab[j]]
	Para i <- 1 Hasta 4
		best <- 0
		Para j <- 1 Hasta 4
			best <- best + A[i,j] * b[colLab[j]]
		FinPara
		x[rowLab[i]] <- best
	FinPara
	
	Para i <- 1 Hasta 4
		Escribir "x", i, " = ", x[i]
	FinPara
FinAlgoritmo