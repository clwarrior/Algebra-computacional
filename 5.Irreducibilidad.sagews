︠1fd434ed-f090-4d1a-8f21-2c59f383ec24s︠
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
#Algoritmo 14.36 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard

#TEST DE IRREDUCIBILIDAD SOBRE Fq[x]
#Se basa en resultado que nos dice que un polinomio en GF(q)[x] de grado n>=1 es  irreducible sii se cumplen ambas:
#    (i) f divide a x^(q^n)-x
#    (ii) mcd(x^q^(n/t)-x, f) = 1 para todos los divisores primos t de n

#A su vez este resultado es la consecuencia de que x^q^n-x es el producto de todos los polinomios mónicos e irreducibles en GF(q) cuyo grado divide a n


#Entrada: f polinomio en la variable x con coeficientes en el cuerpo finito GF(q) con q potencia de un primo p y dicho q
#Salida: 'Reducible' o 'Irreducible'

def irreducible(f, q):
    n = f.degree()
    a = fast_exp_mod(R(x), Integer(q^n), R(f))
    if a!=x: return "Reducible"

    i = 2
    #Para los divisores primos de n
    while i <= n:
        if i % n == 0:
            _ , b = (x^(q^(n/i))).quo_rem(f)
            if gcd(b-x, f) != 1:
                return "Reducible"
                break
        i = next_prime(i)
    if i > n:
        return "Irreducible"


#Para contrastar el veredicto del algoritmo comparamos con el resultado del método is_irreducible() soportado por sage
def rdo(f):
    if f.is_irreducible():
        return "Irreducible"
    else:
        return "Reducible"

#EJEMPLOS
#####################################################################################
#Ejemplos de polinomios de Conway irreducibles
p = 5
R.<x> = GF(p)[x]
f = conway_polynomial(5, 4)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,p)))

︡8b5d7f56-19c6-45f7-963a-1b1d3052f87c︡{"stdout":"Datos de entrada: x^4 + 4*x^2 + 4*x + 2\n"}︡{"stdout":"Resultado esperado: Irreducible\n"}︡{"stdout":"Resultado obtenido: Irreducible\n"}︡{"done":true}
︠0e8537f9-baa0-4a1f-899e-a0a929296ee2s︠

p = 5
R.<x> = GF(p)[x]
f = conway_polynomial(5, 2)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,p)))
︡d84f15b2-c865-4a43-9272-857bd3963421︡{"stdout":"Datos de entrada: x^2 + 4*x + 2\n"}︡{"stdout":"Resultado esperado: Irreducible\n"}︡{"stdout":"Resultado obtenido: Irreducible\n"}︡{"done":true}
︠db7b0460-dcda-4d63-a455-5c3cda85fb2as︠
p = 7
R.<x> = GF(p)[x]
f = conway_polynomial(7, 6)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,p)))
︡5c1b082e-28a0-4916-96a0-41f3ed46d7c5︡{"stdout":"Datos de entrada: x^6 + x^4 + 5*x^3 + 4*x^2 + 6*x + 3\n"}︡{"stdout":"Resultado esperado: Irreducible\n"}︡{"stdout":"Resultado obtenido: Irreducible\n"}︡{"done":true}
︠9d9c3d20-67f6-4d84-b501-73c0cd8b0007s︠

#Ejemplos de polinomios reducibles
p = 7
R.<x> = GF(p)[x]
f = x^6+x^5+x^4+x^3+x^2+x+1
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,p)))
︡b23b2ff3-274a-4bf0-be46-0dd486056552︡{"stdout":"Datos de entrada: x^6 + x^5 + x^4 + x^3 + x^2 + x + 1\n"}︡{"stdout":"Resultado esperado: Reducible\n"}︡{"stdout":"Resultado obtenido: Reducible\n"}︡{"done":true}
︠3e434a3e-de41-4ace-a799-f083ac60e9b4s︠
p = 5
R.<x> = GF(p)[x]
f = x^4+x^3+x^2+x+1
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,p)))
︡a7b735fc-f566-4ce6-b246-88dd4b5c7020︡{"stdout":"Datos de entrada: x^4 + x^3 + x^2 + x + 1\n"}︡{"stdout":"Resultado esperado: Reducible\n"}︡{"stdout":"Resultado obtenido: Reducible\n"}︡{"done":true}
︠ffdbf15f-e00e-4157-9bba-cdf120897b81s︠


#Ejemplos sobre cuerpos finitos
q = 25
R.<x> = GF(q)[x]
f = R.irreducible_element(2)*R.irreducible_element(3)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,q)))
︡4c2a8937-02ae-4dda-88ac-e661c52afac9︡{"stdout":"Datos de entrada: x^5 + (z2 + 3)*x^4 + (2*z2 + 3)*x^3 + 4*z2*x^2 + (z2 + 2)*x + 4\n"}︡{"stdout":"Resultado esperado: Reducible\n"}︡{"stdout":"Resultado obtenido: Reducible\n"}︡{"done":true}
︠40d17b8b-25c5-49d2-96a2-1d64262cdf5c︠


q = 9
R.<x> = GF(q)[x]
f = R.irreducible_element(5)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,q)))
︡1be26c5b-a4fa-4628-a1a7-c471fab89ec3︡{"stdout":"Datos de entrada: x^5 + (2*z2 + 1)*x^2 + z2*x + 2\n"}︡{"stdout":"Resultado esperado: Irreducible\n"}︡{"stdout":"Resultado obtenido: Irreducible\n"}︡{"done":true}
︠91fc4711-44b3-43f9-b36d-5ec9dcbdf6fb︠


q = 27
R.<x> = GF(q)[x]
f = R.irreducible_element(6)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,q)))
︡060e469e-ad55-4730-a677-79eb04c654d6︡{"stdout":"Datos de entrada: x^6 + (z3^2 + 2*z3)*x^4 + (z3^2 + z3 + 1)*x^2 + (z3^2 + z3 + 2)*x + z3^2\n"}︡{"stdout":"Resultado esperado: Irreducible\n"}︡{"stdout":"Resultado obtenido: Irreducible\n"}︡{"done":true}
︠27c58883-f918-446a-9455-423508485982s︠
q = 27
R.<x> = GF(q)[x]
f = R.irreducible_element(3)*R.irreducible_element(5)
print("Datos de entrada: {}".format(f))
print("Resultado esperado: {}".format(rdo(f)))
print("Resultado obtenido: {}".format(irreducible(f,q)))
︡d2c358af-7774-4e0e-b7f6-18a94c773819︡{"stdout":"Datos de entrada: x^8 + (z3 + 2)*x^7 + (2*z3^2 + 1)*x^6 + (z3 + 2)*x^5 + (z3 + 1)*x^4 + (2*z3^2 + z3)*x^3 + 2*z3*x^2 + (z3^2 + 1)*x + 2*z3^2\n"}︡{"stdout":"Resultado esperado: Reducible\n"}︡{"stdout":"Resultado obtenido: Reducible\n"}︡{"done":true}
︠95b56210-b492-4273-b7f3-65e4e8f5c4b4︠









