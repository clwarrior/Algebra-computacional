︠8693f57e-49bf-4a14-8f94-95fa8b511424︠
#Algoritmo 21.33 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard
from sage.rings.polynomial.toy_buchberger import spol

#Método auxiliar que implementa la definición de polinomio-S de g y h, donde g y h en R=F[x1,...,xn]
def pol_S(g,h,F):
    a = g.lt().degrees() # degrees() devuelve una n-pla  con el grado en cada variable
    b = h.lt().degrees()
    c = [max(a[i],b[i]) for i in range(len(a))]
    aux = prod(F[i]^c[i] for i in range(len(c)))
    return aux*g//g.lt()-aux*h//h.lt()

#Algoritmo de división generalizada de un polinomio multivariable f entre varios polinomios (dados en la lista pol=[p1,...,pk]).
#Devuelve una lista l=[l1,...,lk] y r tales que f = p1*l1 + ... + pk*lk + r

def divide(f, pols):
    r = 0
    p = R(f)
    l = []
    for i in range(len(pols)): l.append(0)
    while p != 0:
        j = 0
        while (j < len(pols) and p.lt()%pols[j].lt() != 0): j = j +1
        if (j < len(pols)):
            l[j] = l[j] + p.lt()//pols[j].lt()
            p = p - (p.lt()//pols[j].lt())*pols[j]
        else:
            r = r + p.lt()
            p = p - p.lt()
    return l, r

##ALGORITMO DE BUCHBERGER: CÁLCULO DE UNA BASE DE GRÖBNER
#Entrada: Una lista F=[f1,...,fs] con elementos en R=F[x1,...,xn] (anillo de polinomios en n variables sobre un cuerpo F) y un orden monomial dado en forma de función con dos argumentos
#Salida: Una base de Gröbner G subconjunto de R del ideal I = <f1,...,fs> con respecto al orden dado, con f1,...,fs pertenecientes a G

def Buchberger(F):
    #Iniclizamos G a los elementos de entrada y P a ciertos pares de dichos elementos
    G = Set(F)
    P = Set({})
    for i in range(len(F)):
        for j in range(i+1, len(F)):
            P = P + Set({(G[i], G[j])})
    #Para cada par de P calculamos el polinomio-S módulo los polinomios de entrada
    while not(P.is_empty()):
        aux = P[0]
        P = P - {aux}
        _, h = divide(pol_S(aux[0], aux[1], [x, y]), G)
        h = R(h)
        #Si este resto no es 0, añadimos a P los pares con primer término dicho resto mónico, y como segundo término cada polinomio de entrada
        #Añadimos h al conjunto de los polinomios de entrada
        if h != 0:
            h = h // h.lc()
            for i in range(len(G)):
                P = P + Set({(h, G[i])})
            G = G + Set({h})
    return G

#EJEMPLOS
##########################################################################

#Le damos al anillo orden lexicográfico
R.<x,y> = PolynomialRing(QQ, 2, order='lex')

#EJEMPLO 1
f = R(x^2+3*x+5+y^3+x*y)
g = R(x^2+y^2)
print("Polinomios de entrada: \n{}\n{}".format(f,g))
I = Ideal([f, g])
S = Buchberger([f,g])
#Como los elementos de la base de Groebner tienen la propiedad de que divididos cada uno entre el resto (con la división multivariable) dan resto cero. Por tanto para hacer la base irreducible, dividimos cada elemento de la base dada por nuestro algoritmo entre los demás y nos quedamos únicamente con los que dan un resto no nulo.
Basis = []
for i in range(len(S)):
    _, r = divide(S[i], S - Set({S[i]}))
    if r != 0:
        Basis.append(r)
print("Base de Grobner esperada es: {}".format(I.groebner_basis()))
print("Base de Grobner calculada es: {}".format(Basis))
print("Algoritmo correcto: {}".format(Set(Basis) == Set(I.groebner_basis())))
︡afd787d6-5506-4f2d-8200-4ccc6490b43c︡{"stdout":"Polinomios de entrada: \nx^2 + x*y + 3*x + y^3 + 5\nx^2 + y^2\n"}︡{"stdout":"Base de Grobner esperada es: [x + 1/31*y^5 - 5/31*y^4 + 17/31*y^3 - 4/31*y^2 - 20/31*y + 60/31, y^6 - 2*y^5 + 2*y^4 + 16*y^3 - y^2 + 25]\n"}︡{"stdout":"Base de Grobner calculada es: [y^6 - 2*y^5 + 2*y^4 + 16*y^3 - y^2 + 25, x + 1/31*y^5 - 5/31*y^4 + 17/31*y^3 - 4/31*y^2 - 20/31*y + 60/31]\n"}︡{"stdout":"Algoritmo correcto: True\n"}︡{"done":true}
︠2dc68fb9-ae1b-4531-96f9-2e3ffc016649︠
#EJEMPLO 2
f = R(x^3+y^3+5*x*y)
g = R(3*x^2-4)
print("Polinomios de entrada: \n{}\n{}".format(f,g))
I = Ideal([f, g])
S = Buchberger([f,g])
#Igual que en el ejemplo anterior
Basis = []
for i in range(len(S)):
    _, r = divide(S[i], S - Set({S[i]}))
    if r != 0:
        Basis.append(r)
print("Base de Grobner esperada es: {}".format(I.groebner_basis()))
print("Base de Grobner calculada es: {}".format(Basis))
print("Algoritmo correcto: {}".format(Set(Basis) == Set(I.groebner_basis())))
︡30cc3c15-58d2-4a70-b1f6-3d427338c7ba︡{"stdout":"Polinomios de entrada: \nx^3 + 5*x*y + y^3\n3*x^2 - 4\n"}︡{"stdout":"Base de Grobner esperada es: [x + 675/64*y^5 - 45/16*y^4 + 3/4*y^3 - 5625/16*y - 375/4, y^6 - 100/3*y^2 - 160/9*y - 64/27]\n"}︡{"stdout":"Base de Grobner calculada es: [x + 675/64*y^5 - 45/16*y^4 + 3/4*y^3 - 5625/16*y - 375/4, y^6 - 100/3*y^2 - 160/9*y - 64/27]\n"}︡{"stdout":"Algoritmo correcto: True\n"}︡{"done":true}
︠040b7ef9-5e2a-42b0-a328-5b5aad441ba0︠









