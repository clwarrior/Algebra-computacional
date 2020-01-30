︠1dfbba85-2258-4d63-bcb0-787fd368ed34︠
#Algoritmo auxiliar para comprobar si un polinomio es primitivo
#Polinomio primitivo : mcd(coeficientes) = 1
def primitivo(f):
    c = f.coefficients()
    if len(c)==1: return True
    i = 2
    g = gcd(c[0],c[1])
    while i<(len(c)-1) and g==1:
        g = gcd(g, c[i])
        i = i + 1
    return g==1

#Algoritmo auxiliar que calcula la norma infinito de un polinomio
def norma_infinito(f):
    x = ZZ['x'].0
    u = vector(QQ, f.coefficients())
    return u.norm(Infinity)

#Algoritmo auxiliar que calcula la norma uno de un polinomio
def norma_uno(f):
    x = ZZ['x'].0
    u = vector(QQ, f.coefficients())
    return u.norm(1)

#Algorimo auxiliar que transforma un polinomio de Z(x) a Z_p(x)
def modulo(f,p):
    c = f.coefficients()
    e = f.exponents()
    return sum((c[i]%p)*x^e[i] for i in range(len(c)))

def min_repr(n,p):
    if abs(ZZ(n%p)) < p / 2:
        return ZZ(n%p)
    else:
        return ZZ(n%p) - p

def min_norm(g,p):
    return ZZ['x']([min_repr(i, p) for i in g])

#Entrada: un polinomio primitivo libre de cuadrados f en Z[x] de grado n>=1 con lc(f)>0 y norma infinito igual a A
#Salida: los factores irreducibles de f
#Basado en el algoritmo 15.2 del libro Modern Computer Algebra
def factZx(f):
    x = ZZ['x'].0
    f = ZZ['x'](f)
    n = f.degree()
    #Si el polinomio es de grado 1, ya está factorizado
    if n == 1: return f

    #Si comparte un factor con su derivada, no es square free
    if gcd(f, derivative(f))!=1:
        print("El polinomio no es square-free")
        return "Fallo"

    #Inicializamos las variables
    b = f.lc()    #leading coefficient
    A = norma_infinito(f)
    B = sqrt(n+1)*(2^n)*A*b

    #En estructura repeat-until ponemos primero una vuelta del cuerpo del bucle
    #Queremos encontrar un primo p tal que f módulo p, y la deerivada de f módulo p sean coprimos

    #Escogemos un primo impar aleatorio p tal que 2B < p < 4B
    p = Primes().next(ceil(2*B))
    if p==2:
        if 3 >= 4*B: return "Fallo"
        else: p=3

    #Calculamos f módulo dicho primo
    R.<x> = PolynomialRing(FiniteField(p))
    ff = modulo(f,p)
    ff = R(ff)

    #Mientras no se cumpla la condición, seguimos buscando primos
    while gcd(ff, derivative(ff))!=1:
        p = Primes().next(p)
        if p >= 4*B: return "Fallo"
        R.<x> = PolynomialRing(FiniteField(p))
        ff = modulo(f,p)
        ff = R(ff)

    #Factorización modular
    f = GF(p)[x](f)
    factors = [a[0] for a in list(f.factor())]

    #Inicializar el conjunto de índices T  de factores modulares aún sin tratar, el conjunto G de factores encontrados y el polinomio f_star aún sin factorizar
    T = [0..len(factors)-1]
    s = 1
    G = Set()
    f_star = f
    #Combinación de factores
    while 2*s<=len(T):
        for S in Subsets(T):
            if len(S)==s:
                g_star = b*prod(modulo(factors[i],p) for i in S)
                h_star = b*prod(modulo(factors[i],p) for i in list(Set(T)-Set(S)))
                #Nos aseguramos de que los coeficientes de g_star y h_star están entre -p/2 y p/2
                g_star = min_norm(g_star, p)
                h_star = min_norm(h_star, p)
                if norma_infinito(g_star) >= p/2 or norma_infinito(h_star) >= p/2: return "Fallo con g_star y h_star"

                if norma_uno(g_star)*norma_uno(h_star) <= B:
                    T = list(Set(T)-S)
                    G = G + Set({g_star})
                    f_star = h_star
                    b = f_star.lc()
                    break;
        s = s + 1
    return G + Set({f_star})



#EJEMPLOS
#######################################################

#f = ZZ['x'].random_element(8)
#if f.lc()<0: f=-f
#while primitivo(f)==False or f.is_irreducible()==True:
#    f = ZZ['x'].random_element(8)
#    if f.lc()<0: f=-f

f = x^8+x^6-2*x^5-2*x^4-x^3+x-1
x = ZZ['x'].0
print("El polinomio es {}".format(f))
print("Su factorización es {}".format(f.factor()))
print("Según nuestro algoritmo es {}".format(factZx(f)))
︡1f07767f-dc54-4ff8-ba63-4b3457aa1645︡{"stdout":"El polinomio es x^8 + x^6 - 2*x^5 - 2*x^4 - x^3 + x - 1\n"}︡{"stdout":"Su factorización es (x^6 - x^5 + x^4 - 2*x^3 - x^2 + 2*x - 1)*(x^2 + x + 1)\n"}︡{"stdout":"Según nuestro algoritmo es {x^2 + x + 1, x^6 - x^5 + x^4 - 2*x^3 - x^2 + 2*x - 1}"}︡{"stdout":"\n"}︡{"done":true}
︠83a4757c-5430-4a07-8dd9-0197af7043f3︠









