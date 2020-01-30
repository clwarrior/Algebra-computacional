︠5abd3cef-da00-4ae8-bf8e-c3aaa4db75f4︠
#Algoritmo 2.5 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard

#Para polinomios con coeficientes en un DFU arbitrario, la división no está necesariamente definida, ya que necesitamos que en cada paso podamos dividir el coeficiente principal del dividendo entre el coeficiente principal del divisor.

#Para solucionar este problema definimos el siguiente algoritmo que emula la división 'manual' de polinomios.

#DIVISIÓN POLINÓMICA CON RESTO
#Entrada: polinomios a, b con grados n, m en R[x] siendo R un anillo conmutativo con unidad. El coeficiente director de b es una unidad y n>=m>=0
#Salida: q,r pertenecientes a R[x] que cumplen a = q·b + r y deg(r)<m

def div_pol(a,b,R):
    x = R['x'].0
    n, m = a.degree(), b.degree()
    #Inicialmente el 'resto' es todo el dividendo y el 'cociente' es cero
    r = a
    u = b.lc()^(-1)
    q = 0
    for i in range(n-m,-1,-1):
        #Cuando i es el grado correspondiente
        if r.degree() == m+i:
            #Actualizamos el cociente y el resto
            q_i = r.lc()*u
            q = q + q_i*x^i
            r = r - q_i*x^i*b
        else:
            q_i = 0
    return q,r



#EJEMPLOS
###########################################################
#Probamos el algoritmo de división comparando su resultado con el quo_rem dado por sage
R.<x> = ZZ['x']

f = x^2 +x + 1
g = x
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(f.quo_rem(g)))
print("Resultado obtenido: {}".format(div_pol(f, g, ZZ)))

︡cdbef3f4-6fd3-4dde-a9cb-cf8f06aa7b36︡{"stdout":"Datos de entrada: x^2 + x + 1, x\n"}︡{"stdout":"Resultado esperado: (x + 1, 1)\n"}︡{"stdout":"Resultado obtenido: (x + 1, 1)\n"}︡{"done":true}
︠74e51b56-ef10-461b-9cbb-43a8f533d42fs︠
f = R.random_element(5)
g = R.random_element(3)
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(f.quo_rem(g)))
print("Resultado obtenido: {}".format(div_pol(f, g, ZZ)))
︡fb4145b5-4a46-40f5-a67e-4d33e2d90a9c︡{"stdout":"Datos de entrada: -6*x^5 - x^3 + x^2 + 2*x, x^3 - x^2 + x + 1\n"}︡{"stdout":"Resultado esperado: (-6*x^2 - 6*x - 1, 12*x^2 + 9*x + 1)\n"}︡{"stdout":"Resultado obtenido: (-6*x^2 - 6*x - 1, 12*x^2 + 9*x + 1)\n"}︡{"done":true}
︠654a4a4f-650e-42af-a593-3fc4f72a25cfs︠
f = R.random_element(1)
g = R.random_element(2)
print("Datos de entrada: {}, {}".format(f,g))
print("Resultado esperado: {}".format(f.quo_rem(g)))
print("Resultado obtenido: {}".format(div_pol(f, g, ZZ)))
︡ce2563db-5d3f-4f61-a649-abd49c883039︡{"stdout":"Datos de entrada: -x, -4*x^2 - 3*x - 1\n"}︡{"stdout":"Resultado esperado: (0, -x)\n"}︡{"stdout":"Resultado obtenido: (0, -x)\n"}︡{"done":true}
︠8ba750c7-5814-43f2-ba8a-90577fccbcb0︠









