︠d26c601b-53f5-4805-82a4-39f7b4d9200es︠
#Sección 14 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard

import random
#El cuerpo finito Fp[x] / (f(x)) tiene orden p^n siendo p un primo (Fp) y n el grado del polinomio f, irreducible en Fp

##ALGORITMO DE CANTOR ZASSENHAUS: Factorización en GF(p)[x]

#PARTE 1: DESCOMPOSICIÓN LIBRE DE CUADRADOS
#Entrada: polinomio f en GF(p)[x]
#Salida: una lista L de tuplas, tal que el primer elemento de la tupla es un factor de f libre de cuadrados y el segundo su multiplicidad en f. Es decir si L = [(f1,e1),...m(fn,en)], tenemos que f1,...,fn están libres de cuadros y f = f1^e1·...·fn^en

def cuerpo_bucle(f, g, s, L, p):
    R.<x> = PolynomialRing(FiniteField(p))
    j=1
    f = R(f)
    #Calculamos g el polinomio libre de cuadrados asociado a f (tiene las mismas raíces complejas que f , pero todas de multiplicidad 1)
    g, _ = f.quo_rem(f.gcd(f.derivative()))
    g = R(g)
    while (g != R(1)):
        #Calculamos un posible factor libre de cuadrados de f
        f, _ = f.quo_rem(g)
        f = R(f)
        h = gcd(f, g)
        m, _ = g.quo_rem(h)
        m = R(m)
        #Si el término calculado es válido (1 no lo es por ser el factor trivial), añadimos el nuevo factor libre de cuadrados y su multiplicidad a la lista L
        if (m != R(1)): L.append((m, j*s))
        g = h
        j = j+1
    if (f != R(1)): #f es una potencia p-ésima
        # Computamos la raíz p-ésima de f
        aux=0
        exponentes = f.exponents()
        coeficientes = f.coefficients()
        l=len(exponentes)
        for i in range(l):
            aux = aux + coeficientes[i]*x^(exponentes[i] // p)
        f = R(aux)
        #Actualizamos s
        s=p*s
    return f, g, s, L, p

def SFD(f, p):
    L=[]
    R.<x> = PolynomialRing(FiniteField(p))
    f = R(f.list())
    #Inicializamos las variables y ejecutamos el cuerpo del bucle en estructura do-while
    s=1
    g = 0
    f, g, s, L, p = cuerpo_bucle(f, g, s, L, p)
    while (f!=R(1)):
        f, g, s, L, p = cuerpo_bucle(f, g, s, L, p)
    return L



#PARTE 2: DESCOMPOSICIÓN DE DISTINTO GRADO

#Método auxiliar para calcular b^e mod m de forma más eficiente que con fuerza bruta

#Algoritmo para calcular rápidamente b^e mod m
def fast_exp_mod(b,e,m):
    r = 1
    _, b = b.quo_rem(m)
    while e>0:
        if e%2 == 1:
            _, r = (r*b).quo_rem(m)
        e = e//2
        _, b  = (b*b).quo_rem(m)
    return r

#Entrada: Un polinomio mónico libre de cuadrados f en GF(q)[x] de grado n>0
#Salida: Una lista [g_1,...,g_t] con los términos resultantes de la descomposición de distinto grado de f. Es decir cada término es el producto de polinomios mónicos irreducibles de grado i y divide a f.

#Este algoritmo se basa en que, siendo (G_1,...,G_t) la descomposición de distinto grado de f con g_t!=1, se tiene para todo i>=0 que:
#    hi = x^q^i mod f
#    fi = G_(i+1)·...·G_t
#    g_i = G_i si i>=1


def dd_fac(f, p):
    R.<x> = PolynomialRing(FiniteField(p))
    f = R(f)
    #Comprobamos que el polinomio de entrada sea libre de cuadrados
    if gcd(f,derivative(f, x)) != 1:
        print("El polinomio no es libre de cuadrados")
        return -1
    #Inicializamos las variables
    h0 = x
    f0 = R(f)
    i = 1
    fi = 0
    G = []
    #Damos la primera vuelta del bucle para conseguir una estructura de do-while
    #Calculamos hi
    hi = fast_exp_mod(h0, p, f)
    hi = R(hi)
    #Calculamos gi y lo añadimos a la lista
    aux = hi - x
    aux = R(aux)
    gi = gcd(aux, f0)
    gi = R(gi)
    G.append(gi)
    #Calculamos fi
    fi, _ = f0.quo_rem(gi)
    fi = R(fi)
    #Actualizamos los valores para la siguiente vuelta
    f0=fi
    h0=hi
    while (fi!=R(1)):
        i = i + 1
        #Calculamos hi
        hi = fast_exp_mod(h0, p, f)
        hi = R(hi)
        #Calculamos gi y lo añadimos a la lista
        aux = hi - x
        aux = R(aux)
        gi = gcd(aux, f0)
        gi = R(gi)
        G.append(gi)
        #Calculamos fi
        fi, _ = f0.quo_rem(gi)
        fi = R(fi)
        #Actualizamos los valores para la siguiente vuelta
        f0=fi
        h0=hi
    return G

#PARTE 3: DESCOMPOSICIÓN DE IGUAL GRADO

#EQUAL-DEGREE SPLITTING: Algoritmo auxiliar que nos da una factorización en dos factores de un polinomio

#Entrada:
#    f = polinomio mónico libre de cuadrados en GF(p)[x] de grado n>0
#    p = potencia impar de un primo
#    d = divisor de n tal que todos los factores irreducibles de f tienen grado d
#Salida: Un factor mónico apropiado en GF(p)[x] o 'Fallo'

def eq_deg_sp(f, p, d):
    R.<x> = PolynomialRing(FiniteField(p))
    #Escogemos a en GF(p)[x] con grado menor estricto que n de forma aleatoria
    grado = random.randint(1, f.degree())
    a = R.random_element(grado)
    #Comprobamos que el elemento aleatorio sea no constante
    if (a.degree()==0): return "Fallo"
    #Si encontramos un factor válido, lo devolvemos
    g1 = gcd(a, f)
    if (g1 != R(1)): return g1
    #Calculamos b = a^((q^d-1)/2) mod f
    aux = (p^d-1)//2
    b = fast_exp_mod(a, aux, f)
    #Calculamos g2 = mcd(b-1,f)
    aux2 = b - 1
    aux2 = R(aux2)
    g2 = gcd(aux2, f)
    g2 = R(g2)
    #Si encontramos un factor válido, lo devolvemos
    if (g2 != 1 and g2 != f): return g2
    #Si no hemos encontrado ningún factor, el algoritmo falla
    else: return "Fallo"

#EQUAL-DEGREE FACTORIZATION

#Haciendo uso del algoritmo anterior recursivamente obtenemos los factores de igual grado del polinomio

#Entrada:
#    f = polinomio mónico libre de cuadrados en GF(p)[x] de grado n>0
#    p = potencia impar de un primo
#    d = divisor de n tal que todos los factores irreducibles de f tienen grado d
#    L = lista de factores a devolver
#Salida: Los factores mónicos irreducibles de f en GF(p)[x]

def eq_deg_fac(f, p, d, L):
    #Si f tiene el grado de los factores que buscamos, lo añadimos a la lista
    n = f.degree()
    R.<x> = PolynomialRing(FiniteField(p))
    f = R(f)
    if (n==d):
        L.append(f)
        return L
    #Llamamos a equal_degree_splitting hasta que devuelva un factor de f
    factor = eq_deg_sp(f, p, d)
    while (factor == "Fallo"):
        factor = eq_deg_sp(f, p, d)
    factor = R(factor)
    #Llamamos al algoritmo equal_degree_factorization recursivamente con entrada g y f/factor añadiendo los resultados de ambas llamadas a la lista L (los parámetros se pasan por referencia)
    f = R(f)
    aux, _ = f.quo_rem(factor)
    aux = R(aux)
    eq_deg_fac(factor, p, d, L)
    eq_deg_fac(aux, p, d, L)
    #Devolvemos la lista de factores
    return L

'''Recordamos que para las lista de tuplas el primer termino es el factor y el segundo su multiplicidad'''


#EJEMPLOS
########################################################################################
p = 11
R.<x> = GF(p)[x]

#EJEMPLO SQUARE-FREE DECOMPOSITION
g = (x^2+1)*(x^2 +1)*(-x^2 + x -1)
L1 = SFD(g, p)
print("El polinomio de entrada es: {}".format(g))
print("Resultado de SFD: {}".format(L1))

#EJEMPLO DISTINCT-DEGREE FACTORIZATION
#Llamamos al algoritmo para cada término dado por SFD
L2 = [(dd_fac(g, p), a) for [g, a] in L1]
print("Resultado de DDF: {}".format(L2))

#Eliminamos los factores triviales
aux = []
for a in range(len(L2)):
    touple = L2[a]
    List = touple[0]
    V = touple[1]
    for b in range(len(List)):
        if List[b]!=1: aux.append((List[b], V))
print("Resultado de DDF sin factores triviales: {}".format(L1))
︡37300d1d-1c5f-4072-9cc0-4775acb1be87︡{"stdout":"'Recordamos que para las lista de tuplas el primer termino es el factor y el segundo su multiplicidad'\n"}︡{"stdout":"El polinomio de entrada es: 10*x^6 + x^5 + 8*x^4 + 2*x^3 + 8*x^2 + x + 10\n"}︡{"stdout":"Resultado de SFD: [(10*x^2 + x + 10, 1), (x^2 + 1, 2)]\n"}︡{"stdout":"Resultado de DDF: [([1, 10*x^2 + x + 10], 1), ([1, x^2 + 1], 2)]\n"}︡{"stdout":"Resultado de DDF sin factores triviales: [(10*x^2 + x + 10, 1), (x^2 + 1, 2)]\n"}︡{"done":true}
︠be0bcd2a-63f6-403a-b48c-73aa126ca814s︠

#EJEMPLO EQUAL-DEGREE FACTORIZATION
f = (x^2+1)*(x^2 +x +1)*(-x^2 + x -1)
L=[]
print("El polinomio de entrada es: {}".format(f))
print("El resultado de  EDF es: {}".format(eq_deg_fac(f, p, 2, L)))
︡f12052f2-6b45-4a98-8823-aec87f291e2b︡{"stdout":"El polinomio de entrada es: 10*x^6 + 9*x^4 + 9*x^2 + 10\n"}︡{"stdout":"El resultado de  EDF es: [x^2 + 1, x^2 + 10*x + 1, 10*x^2 + 10*x + 10]\n"}︡{"done":true}
︠8fa8b451-7b31-4bac-95ae-2277337aab38s︠


################################
#Ejemplo para GF(9)
p = 9
R.<x> = GF(p)[x]

#EJEMPLO SQUARE-FREE DECOMPOSITION
cuadrado = R.random_element(2, 2)
g = cuadrado*cuadrado*R.random_element(2, 2)
L1 = SFD(g, p)
print("El polinomio de entrada es: {}".format(g))
print("Resultado de SFD: {}".format(L1))

#EJEMPLO DISTINCT-DEGREE FACTORIZATION
#Llamamos al algoritmo para cada término dado por SFD
L2 = [(dd_fac(g, p), a) for [g, a] in L1]
print("Resultado de DDF: {}".format(L2))

#Eliminamos los factores triviales
aux = []
for a in range(len(L2)):
    touple = L2[a]
    List = touple[0]
    V = touple[1]
    for b in range(len(List)):
        if List[b]!=1: aux.append((List[b], V))
print("Resultado de DDF sin factores triviales: {}".format(L1))
︡6395ac33-780d-447e-80f4-bb566dd5cabe︡{"stdout":"El polinomio de entrada es: 2*z2*x^6 + (2*z2 + 1)*x^5 + 2*z2*x^4 + x^3 + (z2 + 2)*x^2 + (z2 + 2)*x + 2*z2 + 1\n"}︡{"stdout":"Resultado de SFD: [(2*z2*x^2 + (z2 + 2)*x + z2 + 2, 1), (x^2 + (2*z2 + 2)*x + z2 + 1, 2)]\n"}︡{"stdout":"Resultado de DDF: [([1, 2*z2*x^2 + (z2 + 2)*x + z2 + 2], 1), ([1, x^2 + (2*z2 + 2)*x + z2 + 1], 2)]\n"}︡{"stdout":"Resultado de DDF sin factores triviales: [(2*z2*x^2 + (z2 + 2)*x + z2 + 2, 1), (x^2 + (2*z2 + 2)*x + z2 + 1, 2)]\n"}︡{"done":true}
︠8c73819c-ba3f-411c-b325-4af35adfe34bs︠
p = 27
R.<x> = GF(p)[x]
#EJEMPLO EQUAL-DEGREE FACTORIZATION
f = R.random_element(2)*R.random_element(2)*R.random_element(2)
L=[]
print("El polinomio de entrada es: {}".format(f))
print("El resultado de  EDF es: {}".format(eq_deg_fac(f, p, 2, L)))
︡0c0e61c3-dfaf-44cb-8f3b-b2208942918f︡{"stdout":"El polinomio de entrada es: 2*z3*x^6 + (2*z3^2 + 2)*x^5 + z3^2*x^4 + (z3^2 + 2)*x^3 + (2*z3 + 1)*x^2 + (2*z3^2 + 2*z3 + 1)*x + 2*z3^2 + z3 + 2\n"}︡{"stdout":"El resultado de  EDF es: [x^2 + (2*z3^2 + z3 + 1)*x + 2, x^2 + (z3^2 + 2)*x + 2, 2*z3*x^2 + 2*x + 2*z3^2 + z3 + 2]\n"}︡{"done":true}
︠8823a4c3-f309-4535-959e-8678f476cd9a︠









