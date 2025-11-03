Algoritmo MetodoIntercambio_4x4
    Definir i, j, it, l, k Como Entero
    Definir tmp Como Entero
    Definir best, alk, EPS Como Real
    EPS <- 0.000000000001
    Dimension A[5,5]
    Dimension b[5]
    Dimension x[5]
    Dimension rowLab[5]
    Dimension colLab[5]
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
        tmp <- rowLab[l]
        rowLab[l] <- colLab[k]
        colLab[k] <- tmp
		
    FinPara 
	
    Escribir ""
    Escribir "Solucion x = A^{-1}b (respetando etiquetas):"
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
