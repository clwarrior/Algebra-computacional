︠7fd2bad2-3491-4ce2-977f-84c4fb5e3895s︠
#SECCIÓN CON CÓDIGO DE OTROS APARTADOS, QUE UTILIZAREMOS PARA EL ALGORITMO
#########################################################################################
#Definimos la división entre dos elementos
def div(a, b): return a.quo_rem(b)

#Definimos tres métodos, uno para calcular el algoritmo de euclides para enteros, otro para polinomios sobre Q y otro para polinomios sobre un anillo de restos
def ea_ext_int(a, b, R): return Euclidean_Algorithm_Ext(a, b, div, R)
def ea_ext_fp(a, b, R): return Euclidean_Algorithm_Ext(a, b, div, R)

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


# Si K = R/<m> es un cuerpo, con R un domineo euclideo y <m> ideal maximal, se puede calcular el inverso aplicando el
# algoritmo de euclides extendido: mcd(k,m)=1, ya que m es primo, y si 1 = sk+tm, s es el inverso de k, ya que sk=sk+tm=1 mod m

def inverso_zp(a, p):
    _, s, _ = ea_ext_int(a, p, Integer)
    return s % p

def inverso_f(a, f, p):
    R.<x> = PolynomialRing(FiniteField(p))
    a = R(a.list())
    f = R(f.list())
    aux, s, _ = ea_ext_fp(a, f, R)
    s = s // aux
    return s



#EJEMPLOS
########################################################
p = 11 #Sea p cualquier primo
Zp = Integers(p)
aux1 = ZZ.random_element(0, 5*p) % p
res = inverso_zp(aux1, p)
print("El inverso de {} es {}: {}*{} = {}".format(aux1, res, aux1, res, aux1*res % p))
︡8980fbdd-0f9c-4d4c-bcbe-1cb6f6273fee︡{"stdout":"El inverso de 2 es 6: 2*6 = 1\n"}︡{"done":true}
︠51124739-bb39-47ac-9489-76b96a4df40c︠
#Segundo ejemplo
p = 11
R.<x> = PolynomialRing(FiniteField(p))
f = R.irreducible_element(5)
I = f*R
G.<x> = R.quotient_by_principal_ideal(I)
g = G.random_element()
print("Trabajamos en polinomios sobre F_{} / {}. Nuestro elemento será g = {}".format(p, f, g))
out = inverso_f(g, f, p)
g = R(g.list())
out = R(out.list())
print("El resultado es: {}. g*{} = {}".format(out, out, (g*out) % f))
︡c4ef0c17-62e6-4d31-ac7a-bca9849801fa︡{"stdout":"Trabajamos en polinomios sobre F_11 / x^5 + 10*x^2 + 9. Nuestro elemento será g = 9*x^4 + x^3 + 2*x + 5\n"}︡{"stdout":"El resultado es: 7*x^4 + 6*x^3 + 9*x^2 + 8*x + 3. g*7*x^4 + 6*x^3 + 9*x^2 + 8*x + 3 = 1\n"}︡{"done":true}
︠92978bc6-c516-4b76-b7b0-2d746e9b3a92︠

︠f9cf021f-fb67-4def-be3b-70c099373a6b︠









