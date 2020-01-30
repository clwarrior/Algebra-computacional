︠7da57ba0-89d1-4a98-af09-853c59b6e5d2s︠
#SECCIÓN CON CÓDIGO DE OTROS APARTADOS, QUE UTILIZAREMOS PARA EL ALGORITMO
#########################################################################################
#Definimos la división entre dos elementos
def div(a, b): return a.quo_rem(b)

#Definimos tres métodos, uno para calcular el algoritmo de euclides para enteros, otro para polinomios sobre Q y otro para polinomios sobre un anillo de restos
def ea_ext_int(a, b, R): return Euclidean_Algorithm_Ext(a, b, div, R)
def ea_ext_fp(a, b, R): return Euclidean_Algorithm_Ext(a, b, div, R)
def ea_ext_qx(a, b, R): return Euclidean_Algorithm_Ext(a, b, div, R)

#Recibe dos elementos "a" y "b", la función de división "f", el dominio en el que trabajamos "K" y el criterio de parada "p"
def Euclidean_Algorithm_Ext(a,b,f,K):
    r0=a
    r1=b
    s0=1
    s1=0
    t0=0
    t1=1
    while r1 != 0:
        q,_ = f(r0, r1)
        auxr1 = r1
        auxs1 = s1
        auxt1 = t1
        r1=K(r0-q*r1)
        s1=s0-q*s1
        t1=t0-q*t1
        r0 = auxr1
        s0 = auxs1
        t0 = auxt1
    
    return r0, s0, t0
#########################################################################################

#Recibe dos listas de elementos: "x" e "y". Los elementos de "x" son coprimos dos a dos.
#Devuelve un elemento tal que al dividirlo entre cada elemento de x tenga como resto el respectivo y.
#NO TIENE SENTIDO HACERLO SOBRE CUERPOS FINITOS NI SOBRE LOS ENTEROS DE GAUSS PORQUE TODOS SON UNIDADES. NO EXISTE EL CONCEPTO DE SER COPRIMOS

#Definimos tres métodos, uno para enteros, otro para polinomios sobre el anillo de restos y otro para polinomios sobre Q
def Chinese_int(x, y): return ChineseRemainderAlgorithm(x, y, ea_ext_int, div, Integer)
def Chinese_fp(x, y, R): return ChineseRemainderAlgorithm(x, y, ea_ext_fp, div, R)
def Chinese_qx(x, y):
    return ChineseRemainderAlgorithm(x, y, ea_ext_qx, div, QQ['x'])
   



#El método genérico recibe las listas "x" e "y", el algoritmo de euclides extendido "f", el algoritmo de división "d" y el dominio sobre el que trabajamos "G"
def ChineseRemainderAlgorithm(x, y, f, d, G):
    m=1
    #Calculamos el producto de los elementos de x
    for i in x :
        m *= i
    c = 0
    for i in range(len(x)) :
        #a = m // x[i]
        a, _ = d(m, x[i])
        #Calculamos el algoritmo de euclides extendido para a y x[i]
        g, s, t = f( a, x[i], G)
        #Si en vez de g == 1, g está multiplicado por una unidad, divido la solución entre esa unidad
        s = s // g
        t = t // g
        #aux1 = y[i]*s % x[i]
        #aux2 = m // x[i]
        _, aux1 = d(y[i]*s, x[i])
        aux2, _ = d(m, x[i])
        #acumulo el valor en c
        c += aux1*aux2
    return c

#EJEMPLOS
########################################################
#Ejemplo para Enteros
x = [4, 7, 9]
y = [2, 5, 3]
res = Chinese_int(x, y)
print("El resultado es {}".format(res))
for i in range(len(x)):
    print("{} mod {} = {}".format(res, x[i], res % x[i]))
︡566d3077-04b9-4fdb-89bf-4c22f1ae978c︡{"stdout":"El resultado es 390\n"}︡{"stdout":"390 mod 4 = 2\n390 mod 7 = 5\n390 mod 9 = 3\n"}︡{"done":true}
︠7a2457aa-bb3b-479d-93c0-95d91d8554bes︠
#Ejemplo para el anillo de polinomios sobre Q
R.<x> = QQ[] #Definimos el anillo de polinomios sobre Q
X = [x**2+1, x**2-1]
res = Chinese_qx(X, [x, x+1]) #Funciona pero salvo multiplicacion por unidades

print("El resultado es {}".format(res))
for i in range(len(X)):
    print("{} mod {} = {}".format(res, X[i], res % X[i]))
︡6a307cfa-11d8-4a65-9422-30d75179f8c1︡{"stdout":"El resultado es 1/2*x^2 + x + 1/2\n"}︡{"stdout":"1/2*x^2 + x + 1/2 mod x^2 + 1 = x\n1/2*x^2 + x + 1/2 mod x^2 - 1 = x + 1\n"}︡{"done":true}
︠69a5f5e6-00ff-467e-85f0-2352cce935c3s︠
R.<x> = QQ[] #Definimos el anillo de polinomios sobre Q
X = [x**4+1, x**6+1]
res = Chinese_qx(X, [x**2, x+1]) #Funciona pero salvo multiplicacion por unidades
print("\nEl resultado es {}".format(res))
for i in range(len(X)):
    print("{} mod {} = {}".format(res, X[i], res % X[i]))
︡779913ee-0e95-4eb6-9d93-a83e2ac5b8bf︡{"stdout":"\nEl resultado es -1/2*x^9 - 1/2*x^7 - x^6 - 1/2*x^3 + 1/2*x\n"}︡{"stdout":"-1/2*x^9 - 1/2*x^7 - x^6 - 1/2*x^3 + 1/2*x mod x^4 + 1 = x^2\n-1/2*x^9 - 1/2*x^7 - x^6 - 1/2*x^3 + 1/2*x mod x^6 + 1 = x + 1\n"}︡{"done":true}
︠9656fe88-dd65-40cb-8711-e32d33628fe8s︠
    
#ANILLO DE POLINOMIOS SOBRE UN CUERPO FINITO
R.<x> = PolynomialRing(FiniteField(7))
X = [6*x**5 + 4*x**4 + 2*x**3 + 3*x**2 + 3, x**3 + 4*x**2 + 1]
Y = [x**4+5*x**3+3*x, x**2 + x + 1]
res = Chinese_fp(X, Y, R)
print("\nEl resultado es {}".format(res))
for i in range(len(X)):
    print("{} mod {} = {}".format(res, X[i], res % X[i]))
︡2e3060fb-533b-4a7c-9ad7-42c3b2c7e13b︡{"stdout":"\nEl resultado es 2*x^7 + 2*x^6 + 2*x^5 + x^4 + 2*x^3 + 3*x^2 + x + 2\n"}︡{"stdout":"2*x^7 + 2*x^6 + 2*x^5 + x^4 + 2*x^3 + 3*x^2 + x + 2 mod 6*x^5 + 4*x^4 + 2*x^3 + 3*x^2 + 3 = x^4 + 5*x^3 + 3*x\n2*x^7 + 2*x^6 + 2*x^5 + x^4 + 2*x^3 + 3*x^2 + x + 2 mod x^3 + 4*x^2 + 1 = x^2 + x + 1\n"}︡{"done":true}









