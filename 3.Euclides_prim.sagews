︠4f8a9c78-2642-4272-833e-3e7c16ed6410︠
#Algoritmo 6.61 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard

#Método auxiliar que dado un polinomio devuelve su parte primitiva
def primitive_part(f, R):
    f = R(f)
    return R(f / f.content()) if f != 0 else 0

#ALGORITMO DE EUCLIDES PRIMITIVO
#Entrada: polinomios primitivos f,g en R[x] donde R es DFU de grados n >= m
#Salida: mcd(f,g) en R[x]

def euc_prim(f, g, R):
    n, m = f.degree(), g.degree()
    r0, r1 = f, g
    n0, n1 = n, m
    i = 1
    while r1!=0:
        #Pseudodivision
        #Para solucionar el problema de que la division no está necesariamente definida para todo DFU definimos la pseudodivisión, que consiste en multiplicar el dividendo por c^(1+m-n), donde c es el coeficiente principal del divisor, y m,n son los grados del dividendo y divisor respectivamente.
        #Así nos aseguramos poder dividir el coeficiente principal del dividendo entre el coeficiente principal del divisor
        a0 = r1.lc()^(1+n0-n1)*r0
        #Ahora dividimos de forma normal, quedándonos con la parte primitiva del resto
        q1, aux = a0.quo_rem(r1)
        r2 = primitive_part(aux, R)
        r2 = R(r2)
        n2 = r2.degree()
        #Actualizamos las variables
        i = i + 1
        r0, n0 = r1, n1
        r1 = r2
        #Devolvemos el resto normalizado, así hacemos que el mcd sea único
    return r0 // r0.lc()

#Para probar este método comparamos su salida con la del método gcd de sage

#Polinomios en una variable sobre el anillo de los enteros
R.<x> = ZZ[x]
f = R.random_element(4)
g = R.random_element(5)

print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(gcd(f,g)))
print("Resultado obtenido: {}".format(euc_prim(f, g, R)))
︡bf16e2bf-5ebe-485c-a8a1-72dc4069788e︡{"stdout":"Datos de entrada: -6*x^4 - 5*x^3 + x^2 - 2*x, -2*x^5 - x^4 - x^3 + x^2 - x\n"}︡{"stdout":"Resultado esperado: x\n"}︡{"stdout":"Resultado obtenido: x\n"}︡{"done":true}
︠0cf461f8-71a7-4fd6-9572-d26beab2a294s︠
f = (x^4+7*x+27)*4*(x^2+x+1)
g = (x^4+7*x+27)*(3*x+2)

print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(gcd(f,g)))
print("Resultado obtenido: {}".format(euc_prim(f, g, R)))
︡67e5ac28-10fa-4a44-a377-a0cb5c427e8f︡{"stdout":"Datos de entrada: 4*x^6 + 4*x^5 + 4*x^4 + 28*x^3 + 136*x^2 + 136*x + 108, 3*x^5 + 2*x^4 + 21*x^2 + 95*x + 54\n"}︡{"stdout":"Resultado esperado: x^4 + 7*x + 27\n"}︡{"stdout":"Resultado obtenido: x^4 + 7*x + 27\n"}︡{"done":true}
︠151e67cd-443e-4b15-ba8d-9eac33f49282︠










