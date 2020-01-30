︠fcd50249-c54b-4d9f-bb6f-435acac2bb9d︠
#Sección 2.1 de Introduction to Computer Algebra de Anne Frühbis-Krüger y Christoph Lossen

#División en el anillo de los enteros de Gauss
def div_gauss(a, b):
    a, b = ZZ[i](a), ZZ[i](b)
    q, _ = a.quo_rem(b)
    Q = q.imag().round()*I + q.real().round()
    r = a - Q*b
    return Q, r

#División genérica (soportada por sage)
def div(a, b): return a.quo_rem(b)

#Definimos distintas funciones para poder probar el algoritmo con distintos D.E., en cada una llamamos a Euclidean Algorithm con la operación de división correspondiente como tercer parámetro
def ea_int(a, b): return Euclidean_Algorithm(a, b, div)
def ea_fp(a, b): return Euclidean_Algorithm(a, b, div)
def ea_qx(a, b): return Euclidean_Algorithm(a, b, div)
def ea_g(a, b): return Euclidean_Algorithm(a, b, div_gauss)
def ea_zp(a, b, p):
    Zp = Integers(p)
    return Euclidean_Algorithm(Integer(Zp(a)), Integer(Zp(b)), div)


#El algoritmo de Euclides se basa en la la igualdad mcd(a,b) = mcd(a,b%a). Hacemos uso de esta igualdad sucesivamente hasta llegar al mcd.

#Entrada: a,b en R dominio euclídeo y como tercer parámetro recibe la operación de división del dominio
#Salida: mcd(a,b) (Notamos que el mcd no es único, ya que u⋅mcd es mcd para cualquier unidad u)
def Euclidean_Algorithm(a,b,f):
    r0=max(a,b)
    r1=min(a,b)
    while (r1 != 0):
        aux = r1
        _, r1 = f(r0, r1)
        r0 = aux
    return r0
︡ede51a08-9140-4923-8805-c69573ec5cbd︡{"done":true}
︠cce410e1-b4ef-4600-98f6-ed9ad2390a72︠
#OBSERVACIÓN: Comparamos el resultado obtenido con gcd de sage, hay que tener en cuenta la unicidad salvo unidades

#Prueba con enteros
print("Datos de entrada: {}, {}".format(10,5))
print("Resultado esperado: {}".format(gcd(10,5)))
print("Resultado obtenido: {}".format(ea_int(10,5)))
︡d236ec05-eb1a-4c85-93be-28a3eb96a97b︡{"stdout":"Datos de entrada: 10, 5\n"}︡{"stdout":"Resultado esperado: 5\n"}︡{"stdout":"Resultado obtenido: 5\n"}︡{"stdout":"5\n"}︡{"done":true}
︠b48d7742-bdb7-4376-8302-94e75d5f346as︠
#Prueba con polinomios sobre los racionales
R.<x> = QQ[]
print("Datos de entrada: {}, {}".format(x^2+x, x))
print("Resultado esperado: {}".format(gcd(x^2+x, x)))
print("Resultado obtenido: {}".format(ea_qx(x^2+x, x)))
︡1e86160a-440b-4a16-9eb7-ee7727b7f319︡{"stdout":"Datos de entrada: x^2 + x, x\n"}︡{"stdout":"Resultado esperado: x\n"}︡{"stdout":"Resultado obtenido: x\n"}︡{"done":true}
︠51736068-d10f-4126-857f-7276da6c65d4︠
#Prueba en anillos de polinomios sobre cuerpos finitos
R.<x> = PolynomialRing(FiniteField(7))
f = R.random_element(5)
g = R.random_element(3)
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(gcd(f,g)))
print("Resultado obtenido: {}".format(ea_fp(f, g)))
#Resultado multiplicado por la unidad 5
︡c6b5d88b-a509-49e5-9fbe-618773ac6b48︡{"stdout":"Datos de entrada: 4*x^5 + 2*x^4 + x^2 + 5*x + 6, 5*x^3 + x^2 + 3*x\n"}︡{"stdout":"Resultado esperado: x + 1\n"}︡{"stdout":"Resultado obtenido: 5*x + 5\n"}︡{"done":true}
︠4a70b56d-245a-490e-b810-dda1675da5bbs︠
#Prueba en el grupo cíclico Zp para p primo
p = 3 #Sea p cualquier primo
Zp = Integers(p)
f = ZZ.random_element(0, 5*p)
g = ZZ.random_element(0, 5*p)
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(gcd(f,g)%p))
print("Resultado obtenido: {}".format(ea_zp(f, g,p)))
︡4ef89671-c3fc-417f-8c26-304dad5372aa︡{"stdout":"Datos de entrada: 5, 0\n"}︡{"stdout":"Resultado esperado: 2\n"}︡{"stdout":"Resultado obtenido: 2\n"}︡{"done":true}
︠c27a4ec2-f3b6-492b-bb40-938633131b7ds︠
#Prueba para los enteros de Gauss
G = ZZ[I]
f = G.random_element()
g = G.random_element()
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(gcd(f,g)))
print("Resultado obtenido: {}".format(ea_g(f, g)))
︡0760d3c3-42af-45f5-bb39-c66f7568dd71︡{"stdout":"Datos de entrada: 0, 1\n"}︡{"stdout":"Resultado esperado: 1\n"}︡{"stdout":"Resultado obtenido: 1\n"}︡{"done":true}
︠b35e3a72-e138-4009-8c2d-c7e4922d26dc︠
#Otro ejemplo con enteros de Gauss
f = G(I-3)
g = G(I+5)
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(gcd(f,g)))
print("Resultado obtenido: {}".format(ea_g(f, g)))
#Resultado multiplicado por la unidad I
︡adafa5f9-03e0-478c-b818-f20ecebfb467︡{"stdout":"Datos de entrada: I - 3, I + 5\n"}︡{"stdout":"Resultado esperado: I + 1\n"}︡{"stdout":"Resultado obtenido: -I + 1\n"}︡{"done":true}
︠083e1388-9ab7-4585-89f6-5bc132d9f123︠









