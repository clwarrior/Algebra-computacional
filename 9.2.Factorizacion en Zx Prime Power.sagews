︠fbdb4367-9950-471c-bcee-6b953f80598d︠
def min_repr(n,p):
    if abs(ZZ(n%p)) < p / 2:
        return ZZ(n%p)
    else:
        return ZZ(n%p) - p

def min_norm(g,p):
    return ZZ['x']([min_repr(i, p) for i in g])

#Recibe m, elemento de un anillo R, y polinomios f, g, h, s, t, con coeficientes en R.
#Tales que f = g*h mod m y s*g + t*h = 1 mod m.
#lc(f) no es diviso de cero módulo m, h es mónico, deg(f) = n = deg(g) + deg(h), deg(s) < deg(h) y deg(t) < deg(g)

#Devuelve g1, h1, s1, t* con coeficientes en R, tales que
#f = g1*h1 mod m^2 y s1*g1 + t1*h = 1 mod m^2.
#h1 es mónico, todos los polinomios de salida son congruentes con su respectivo polinomio de entrada módulo m
#deg(g1) = deg(g), deg(h1) = deg(h), deg(s1) < deg(h1), deg(t1) < deg(g1)

#BASADO EN EL ALGORITMO 15.10 de Modern Computer Algebra
def hensel_step(f,m,g,h,s,t):
    #Casteamos los tipos inicialmente, por si acaso
    m = ZZ(m)
    GFm2X = GF(m^2)[x]
    f = GFm2X(f)
    g = GFm2X(g)
    s = GFm2X(s)
    t = GFm2X(t)

    # calculamos valores auxiliares
    e = f - g*h
    q = (s*e) // h
    r = (s*e) % h
    g1 = g + t*e + q*g
    h1 = h + r

    # Segundo paso de hensel-step
    b = s*g1 + t*h1 - 1
    c = (s*b) // h1
    d = (s*b) % h1

    s1 = s - d
    t1 = t - t*b - c* g1

    # Casteamos de nuevo
    g1 = ZZX(g1)
    h1 = ZZX(h1)
    s1 = ZZX(s1)
    t1 = ZZX(t1)

    return (g1, h1, s1, t1)

#Recibe un elemento p de un anillo. f, polinomio con coeficientes en ese anillo, tal que
#lc(f) es una unidad módulo p. Recibe también una lista de polinomios mónicos no constatntes
#coprimos dos a dos módulo p, tales que f = lc(f)*f1*....*fr mod p. Recibe un natural l

#Devuelve una lista de poinomios mónicos f*1.... f*r tales que f = lc(f)*f*1*...*f*r mod p^l
#Y f*i = fi mod p para todo i

#Basado en el algoritmo 15.17 de Modern Computer Algebra
def multifactor_hensel_lifting(f, factors_mod_p, p, l):
    ZZX = ZZ[x]
    GFpX = GF(p)[x]
    f = ZZX(f)

    if len(factors_mod_p) == 1: return [ZZX(f / f.leading_coefficient())], p
    k = len(factors_mod_p) // 2
    aux1 = 1

    #Calculamos g0 = lc(f)*f1*...*fk mod p
    #h0 = f(k+1)*...*fr mod p
    for i in range(k):
        aux1 = aux1*factors_mod_p[i]
    aux2 = 1
    for i in range(k, len(factors_mod_p)):
        aux2 = aux2*factors_mod_p[i]
    g = aux1*f.leading_coefficient()
    h = aux2

    #Casteamos
    g = GFpX(g)
    h = GFpX(h)

    #Buscamos s, t tales que s*g +t*h = 1 mod p
    r,s,t = g.xgcd(h)
    s = ZZX(s); t = ZZX(t)

    #Con j = 1....d
    #Llamamos a Hensel Step con m = p^(2^(j-1))
    #Para ir elevando las congruencias de la iteración anterior
    m = p
    while m < p^l:
        g,h,s,t = hensel_step(f,m,g,h,s,t)
        m = m^2

    #Llamamos al algoritmo de forma recursiva con g y los primeros k factores
    g_factors_mod_p_l, _ = multifactor_hensel_lifting(g, factors_mod_p[:k], p, l)
    #Llamamos también con h y los factores de k+1 al último
    h_factors_mod_p_l, _ = multifactor_hensel_lifting(h, factors_mod_p[k:], p, l)

    #Devolvemos la concatenación de ambas llamadas
    return g_factors_mod_p_l + h_factors_mod_p_l, m


# Modular Factorization Algorithm

import itertools

#Recibe un polinomio primitivo libre de cuadrados de Z[x] de grado n >= 1, con lc(f) > 0 y norma-infinito = A

#Devuelve los factores irreducibles f1, ... , fk, en Z[x] de f
def factZx(f):
    if f.degree() == 1: return [f]
    n = f.degree()
    A = f.norm(infinity)
    b = f.leading_coefficient()
    B = sqrt(n+1)*2^n*A*b
    C = (n+1)^(2*n)*A^(2*n-1)
    gamma = ceil(2*log(C,2))
    p_bound = ceil(2*gamma*ln(gamma))

    # Elegimos un primo p <= 2*gamma*ln(gamma)
    #Tal que p no divida a b y gcd(f mod p, (f mod p)') = 1 en Fp
    while(True):
        p = random_prime(p_bound) # pick prime less than p_bound
        GFpX = GF(p)[x]
        f_mod_p = GFpX(f)
        df_mod_p = f_mod_p.derivative()
        # Si se cumple la condición, lo hemos encontrado
        if(b % p != 0 and gcd(f_mod_p, df_mod_p)==1): break



    # Factorizamos módulo en Zp[x] h1 ... hr
    factors_mod_p = [a[0] for a in list(f_mod_p.factor())]
    #Casteamos
    factors_mod_p = [ZZX(g) for g in factors_mod_p]


    # Llamamos a hensel lifting para buscar una factorización en p^l
    # Y tales que cada factor gi = hi mod p
    l = ceil(log(2*B+1, p))
    factors_mod_m, m = multifactor_hensel_lifting(f, factors_mod_p, p, l)


    # Paso de combinación de los factores
    factors_mod_m = set(factors_mod_m)
    true_factors = []
    s = 1; f_star = f
    while(2*s <= len(factors_mod_m)):
        #Para cada subconjunto S de factores de cardinalidad s
        for factors_mod_m_subset in itertools.combinations(factors_mod_m, s):

            aux1 = 1
            for i in range(len(factors_mod_m_subset)):
                aux1 = aux1 * factors_mod_m_subset[i]
            g_star = b * aux1
            complement = tuple(factors_mod_m - set(factors_mod_m_subset))

            aux1 = 1
            for i in range(len(complement)):
                aux1 = aux1 * complement[i]
            h_star = b*aux1

            # Después de generar g* y h*, nos aseguramos de que sus normas son como máximo p^l/2
            g_star = min_norm(g_star, m)
            assert g_star.norm(infinity) < m/2, "La norma de g* después de minimizar es {} >= {}/2!".format(g_star.norm(1), m)
            h_star = min_norm(h_star, m)
            assert h_star.norm(infinity) < m/2, "La norma de h* después de minimizar es {} >= {}/2!".format(h_star.norm(1), m)

            #Si la norma 1 de g* * norma 1 de h* es <= B
            if (g_star.norm(1) * h_star.norm(1) <= B):
                factors_mod_m = factors_mod_m.difference(factors_mod_m_subset)
                true_factors.append(g_star)
                f_star = h_star
                b = f_star.leading_coefficient()
                break
            else: s += 1

    # Si f* es un polinomio no constante, lo añadimo s la lista a devolver
    if(f_star.degree() > 0): true_factors.append(f_star)
    return true_factors


# EJEMPLO
#######################################################################
ZZX = ZZ[x]
f = (x+2)*(x-2)*(x**2+2*x+2)
f = ZZX(f)
assert f.is_squarefree()
factors = []
while(len(factors) <= 1):
    factors = factZx(f)
print("El polinomio es {}".format(f))
print("Su factorización es {}".format(f.factor()))
print("Según nuestro algoritmo es {}".format(factors))
︡db0404a5-9215-4cc0-87fd-60a4b6df7671︡{"stdout":"El polinomio es x^4 + 2*x^3 - 2*x^2 - 8*x - 8\n"}︡{"stdout":"Su factorización es (x - 2) * (x + 2) * (x^2 + 2*x + 2)\n"}︡{"stdout":"Según nuestro algoritmo es [x - 2, x + 2, x^2 + 2*x + 2]\n"}︡{"done":true}
︠a83aee52-fd02-4033-9f8e-9c093ee587d1︠









