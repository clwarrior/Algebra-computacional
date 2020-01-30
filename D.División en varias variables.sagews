︠70a8a22c-4c7b-49fc-b4bd-7727cb75ee7cs︠
#Demostración de corrección en página 598 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard

#ALGORITMO DE DIVISIÓN MULTIVARIABLE
#Algoritmo de división generalizada de un polinomio multivariable f entre varios polinomios.

#Entrada: El polinomio dividendo f y la lista pol=[p1,...,pk] de polinomios divisores. Todos ellos no nulos
#Salida: Una lista l=[l1,...,lk] y r tales que f = p1*l1 + ... + pk*lk + r

def divide(f, pols):
    #Inicialización de variables, p indica la parte que queda por dividir
    r = 0
    p = f
    #Inicializamos la lista de coeficientes a 0
    l = []
    for i in range(len(pols)): l.append(0)
    #Mientras quede polinomio por dividir
    while p != 0:
        #Buscamos j tal que pols[j].lt() divida a p.lt()
        j = 0
        while (j < len(pols) and p.lt()%pols[j].lt() != 0): j = j +1
        #Si lo encontramos, lo usamos para actualizar el cociente correspondiente a dicho polinomio y la parte de polinomio que queda por dividir
        if (j < len(pols)):
            l[j] = l[j] + p.lt()//pols[j].lt()
            p = p - (p.lt()//pols[j].lt())*pols[j]
        #Si no lo encontramos, actualizamos el resto y la parte de polinomio que queda por dividir
        else:
            r = r + p.lt()
            p = p - p.lt()
    return l, r
︡b02cb1ec-918d-40f5-8b97-3b53c01745f1︡{"done":true}
︠0bd70ae8-6ea8-4182-a16f-53c840c4a3fc︠
#EJEMPLOS
#######################################################################################3
#Damos orden lexicográfico al anillo
P.<x,y> = PolynomialRing(QQ, 2, order='lex')
f = x*y**2 +1
f1 = x*y+1
f2 = y+1
F=[f1,f2]
print("Dividendo:{}\nDivisores: {}, {}".format(f, f1,f2))
l, r = divide(f, F)
print("Cocientes: {}, {}\nResto: {}".format(l[0],l[1], r))
#Comprobemos que es correcto
aux = r
for i in range(len(l)):
    aux = aux + l[i]*F[i]
print("Correcto: {}".format(aux==f))
︡a01baa2e-9d85-44dc-9055-8f4b9f64df9b︡{"stdout":"Dividendo:x*y^2 + 1\nDivisores: x*y + 1, y + 1\n"}︡{"stdout":"Cocientes: y, -1\nResto: 2\n"}︡{"stdout":"Correcto: True\n"}︡{"done":true}
︠101a9166-5b10-4d4b-a2ba-afc3dea04f73︠
#SEGUNDO EJEMPLO

f = x^3 + 4*y^2 + 3*x^2*y + 5*x*x^2 + 3*x + 27
f1 = x^2 +3
f2 = x^3 + x*y
F=[f1,f2]
print("Dividendo:{}\nDivisores: {},{}".format(f, f1,f2))
l, r = divide(f, F)
print("Cocientes: {}\nResto: {}".format(l, r))
#Comprobemos que es correcto
aux = r
for i in range(len(l)):
    aux = aux + l[i]*F[i]
print("Correcto: {}".format(aux==f))
︡c1699695-ab6c-49fa-ab46-6b2cd843f672︡{"stdout":"Dividendo:6*x^3 + 3*x^2*y + 3*x + 4*y^2 + 27\nDivisores: x^2 + 3,x^3 + x*y\n"}︡{"stdout":"Cocientes: [6*x + 3*y, 0]\nResto: -15*x + 4*y^2 - 9*y + 27\n"}︡{"stdout":"Correcto: True\n"}︡{"done":true}
︠c775f457-e726-423e-91ce-cc59d3f90f59︠









