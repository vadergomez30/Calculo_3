Proceso MenuPrincipal
    Definir opc, opc2 Como Entero
	
    Escribir "Menu Principal"
    Escribir "1. Solucion numerica de ecuaciones de una sola variable"
    Escribir "2. Soluciones de sistemas de ecuaciones lineales"
    Escribir "Seleccione una opcion (1-2): " Sin Saltar
    Leer opc
	
    Segun opc Hacer
        1:
            Escribir "Has seleccionado la opcion 1"
            Escribir "Seleccione el metodo a utilizar:"
            Escribir "1. Metodo de la Biseccion"
            Escribir "2. Metodo de la Falsa Posicion"
            Escribir "3. Metodo de Newton Raphson"
            Escribir "4. Metodo de la Secante"
            Escribir "Seleccione una opcion (1-4): " Sin Saltar
            Leer opc2
			
            Segun opc2 Hacer
                1:
                    EjecutarMetodoBiseccion
                2:
                    EjecutarFalsaPosicion
                3:
                    EjecutarNewtonRaphson
                4:
                    EjecutarMetodoSecante
                De Otro Modo:
                    Escribir "Opcion no valida"
            FinSegun
			
        2:
            Escribir "En proceso"
			
        De Otro Modo:
            Escribir "Opcion no valida"
    FinSegun
FinProceso
//------------------------------------------------
// FUNCIONES DE APOYO
//------------------------------------------------

Funcion r <- f1(x)
    r <- sen(x) + cos(x)
FinFuncion

Funcion r <- f2(x)
    r <- (x - 5) * (x + 3) * (cos(x))^3
FinFuncion

Funcion r <- f3(x)
    r <- (x^2) * (sen(x))^2
FinFuncion

Funcion r <- f4(x)
    r <- exp(x) * sen(x)
FinFuncion

//------------------------------------------------
// METODO DE BISECCION
//------------------------------------------------

SubProceso ejecutarMetodoBiseccion
    Definir a,b,x,error,fa,fb,fx Como Real
    Definir opcion, iteraciones Como Entero
    Definir raizEncontrada Como Logico
	
    Escribir "== Metodo de la Biseccion =="
    Escribir "1. f(x)=sen(x)+cos(x)"
    Escribir "2. f(x)=(x-5)(x+3)cos^3(x)"
    Escribir "3. f(x)=x^2*sen^2(x)"
    Escribir "4. f(x)=e^x*sen(x)"
    Escribir "Seleccione una funcion: " Sin Saltar
    Leer opcion
	
    Escribir "Ingrese el intervalo [a,b]:"
    Escribir "a = " Sin Saltar
    Leer a
    Escribir "b = " Sin Saltar
    Leer b
    Escribir "Ingrese el error permitido: " Sin Saltar
    Leer error
	
    Segun opcion Hacer
        1: fa <- f1(a); fb <- f1(b)
        2: fa <- f2(a); fb <- f2(b)
        3: fa <- f3(a); fb <- f3(b)
        4: fa <- f4(a); fb <- f4(b)
    FinSegun
	
    Si fa * fb > 0 Entonces
        Escribir "No hay cambio de signo. No se puede aplicar el metodo."
    FinSi
	
    iteraciones <- 0
    raizEncontrada <- Falso
	
    Mientras raizEncontrada = Falso Y iteraciones < 50 Hacer
        x <- (a + b) / 2
		
        Segun opcion Hacer
            1: fx <- f1(x)
            2: fx <- f2(x)
            3: fx <- f3(x)
            4: fx <- f4(x)
        FinSegun
		
        Si Abs(fx) < error Entonces
            Escribir "La raiz aproximada es: ", x
            raizEncontrada <- Verdadero
        Sino
            Si fa * fx < 0 Entonces
                b <- x
            Sino
                a <- x
            FinSi
        FinSi
        iteraciones <- iteraciones + 1
    FinMientras
	
    Si raizEncontrada = Falso Entonces
        Escribir "No se encontro la raiz en 50 iteraciones."
    FinSi
FinSubproceso

//------------------------------------------------
// METODO DE FALSA POSICION
//------------------------------------------------

SubProceso ejecutarFalsaPosicion
    Definir a,b,xr,error,fa,fb,fx Como Real
    Definir opcion,j Como Entero
	
    Escribir "== Metodo de la Falsa Posicion =="
    Escribir "1. f(x)=sen(x)+cos(x)"
    Escribir "2. f(x)=(x-5)(x+3)cos^3(x)"
    Escribir "3. f(x)=x^2*sen^2(x)"
    Escribir "4. f(x)=e^x*sen(x)"
    Escribir "Seleccione una funcion: " Sin Saltar
    Leer opcion
	
    Escribir "Ingrese a: " Sin Saltar
    Leer a
    Escribir "Ingrese b: " Sin Saltar
    Leer b
    Escribir "Ingrese error de tolerancia: " Sin Saltar
    Leer error
	
    j <- 0
    Repetir
        Segun opcion Hacer
            1: fa <- f1(a); fb <- f1(b)
            2: fa <- f2(a); fb <- f2(b)
            3: fa <- f3(a); fb <- f3(b)
            4: fa <- f4(a); fb <- f4(b)
        FinSegun
		
        xr <- a - fa * (b - a) / (fb - fa)
		
        Segun opcion Hacer
            1: fx <- f1(xr)
            2: fx <- f2(xr)
            3: fx <- f3(xr)
            4: fx <- f4(xr)
        FinSegun
		
        Si Abs(fx) < error Entonces
            Escribir "Raiz aproximada: ", xr
        FinSi
		
        Si fa * fx < 0 Entonces
            b <- xr
        Sino
            a <- xr
        FinSi
		
        j <- j + 1
    Hasta Que j >= 1000
    Escribir "No se encontro la raiz en 1000 iteraciones."
FinSubproceso

//------------------------------------------------
// METODO DE NEWTON-RAPHSON
//------------------------------------------------

SubProceso ejecutarNewtonRaphson
    Definir x,error,fx,dfx Como Real
    Definir opcion,iteraciones Como Entero
	
    Escribir "== Metodo de Newton-Raphson =="
    Escribir "1. f(x)=sen(x)+cos(x)"
    Escribir "2. f(x)=(x-5)(x+3)cos^3(x)"
    Escribir "3. f(x)=x^2*sen^2(x)"
    Escribir "4. f(x)=e^x*sen(x)"
    Escribir "Seleccione una funcion: " Sin Saltar
    Leer opcion
	
    Escribir "Ingrese el valor inicial (x0): " Sin Saltar
    Leer x
    Escribir "Ingrese el error de tolerancia: " Sin Saltar
    Leer error
	
    iteraciones <- 0
	
    Repetir
        Segun opcion Hacer
            1: fx <- f1(x); dfx <- cos(x) - sen(x)
            2: fx <- f2(x); dfx <- cos(x)^2*(2*(x-1)*cos(x)-3*(x^2-2*x-15)*sen(x))
            3: fx <- f3(x); dfx <- 2*x*sen(x)*(sen(x)+x*cos(x))
            4: fx <- f4(x); dfx <- exp(x)*(sen(x)+cos(x))
        FinSegun
		
        Si dfx = 0 Entonces
            Escribir "Error: derivada cero."
        FinSi
		
        x <- x - (fx / dfx)
        iteraciones <- iteraciones + 1
    Hasta Que Abs(fx) < error O iteraciones >= 1000
	
    Si iteraciones >= 1000 Entonces
        Escribir "No converge en 1000 iteraciones."
    Sino
        Escribir "Raiz aproximada: ", x
        Escribir "Iteraciones realizadas: ", iteraciones
    FinSi
FinSubproceso

//------------------------------------------------
// METODO DE LA SECANTE
//------------------------------------------------

SubProceso ejecutarMetodoSecante
    Definir x0,x1,x2,fx0,fx1,error Como Real
    Definir opcion,iteraciones,i Como Entero
	
    Escribir "== Metodo de la Secante =="
    Escribir "1. f(x)=sen(x)+cos(x)"
    Escribir "2. f(x)=(x-5)(x+3)cos^3(x)"
    Escribir "3. f(x)=x^2*sen^2(x)"
    Escribir "4. f(x)=e^x*sen(x)"
    Escribir "Seleccione una funcion: " Sin Saltar
    Leer opcion
	
    Escribir "Ingrese x0: " Sin Saltar
    Leer x0
    Escribir "Ingrese x1: " Sin Saltar
    Leer x1
    Escribir "Ingrese error: " Sin Saltar
    Escribir "Numero maximo de iteraciones: " Sin Saltar
    Leer iteraciones
	
    i <- 0
	
    Mientras i < iteraciones Hacer
        Segun opcion Hacer
            1: fx0 <- f1(x0); fx1 <- f1(x1)
            2: fx0 <- f2(x0); fx1 <- f2(x1)
            3: fx0 <- f3(x0); fx1 <- f3(x1)
            4: fx0 <- f4(x0); fx1 <- f4(x1)
        FinSegun
		
        Si Abs(fx0 - fx1) < 0.0000000001 Entonces
            Escribir "Error: division por cero."
        FinSi
		
        x2 <- x1 - (fx1 * (x0 - x1)) / (fx0 - fx1)
		
        Si Abs(x2 - x1) < error Entonces
            Escribir "Raiz aproximada: ", x2
            Escribir "Iteraciones: ", i + 1
        FinSi
		
        x0 <- x1
        x1 <- x2
        i <- i + 1
    FinMientras
	
    Escribir "No se encontro la raiz en ", iteraciones, " iteraciones."
    Escribir "Ultima aproximacion: ", x1
FinSubproceso
