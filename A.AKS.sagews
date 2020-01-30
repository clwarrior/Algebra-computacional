︠c7f82e7a-a003-48b0-81b4-89dddfc9b117︠
#Método auxiliar que devuelve la longitud (el número de dígitos) de la representación binaria de un número dado
def len(a):
    if(a == 0): return 1
    else:
        return floor(log(abs(a), 2)) + 1

#Método auxiliar que calcula el orden de a módulo n. Es el k más pequeño para el cual a^k = 1 mod n
def orden(a, n):
    i = 1
    ini = a
    while(a != 1):
        a = a*ini
        a = a % n
        i = i + 1
    return i

#ALGORITMO DE PRIMALIDAD AKS
#Es un método determinista que dado un número natural n nos dice si es primo o compuesto
def AKS(n):
    #Sea a cada uno de los posibles factores de n
    for a in range(2, sqrt(n)+1):
        #Si n es una potencia de a, entonces es compuesto
        acum = a
        while(acum <= n):
            if(acum == n): return "Compuesto"
            acum = acum*acum
    #Buscamos el valor r más pequeño tal que se cumpla una de las siguientes:
    #     (i) gcd(n,r) > 1
    #     (ii) n tiene orden módulo r > 4 * log_2(n)^2
    r = 2
    while(gcd(n, r) <= 1 or (gcd(n, r) == 1 and orden(n % r, r) > 4*len(n)^2)):
        r = r + 1
    if(r == n): return "Primo"
    #Si r comparte un factor no trivial con n, n es compuesto
    if(gcd(n, r) > 1): return "Compuesto"
    #Definimo el anillo en que queremos comprobar la primalidad
    R.<x> = PolynomialRing(FiniteField(n))
    for j in range(1, 2*len(n)*floor(sqrt(r))+1):
        #Usando la generalización del pequeño teorema de Fermat:
        #    Si n y j son coprimos y n es primo, entonces se cumple
        #    (x + j)^n = x^n + j mod n
        f = R((x + j)^n)
        g = R(x^n + j)
        h = R(x^r - 1)
        _, modulo1 = f.quo_rem(h)
        _, modulo2 = g.quo_rem(h)
        #Si no se cumple la congruencia, podemos concluir que n no es primo
        if(modulo1 != modulo2): return "Compuesto"
    #Si tras hacer este proceso no hemos encontrado ningún factor, ni imcumplido ninguna congruencia, podemos afirmar que n es primo
    return "Primo"



#EJEMPLOS
####################################################################
#Para comprobar el resultado comparamos con el método is_prime() de sage

def es_primo(n):
    if is_prime(n): return "Primo"
    else: return "Compuesto"

n = 105929
print("Número: {}".format(n))
print("Resultado esperado: {}".format(es_primo(n)))
print("Resultado obtenido: {} \n".format(AKS(n)))

n = 524287
print("Número: {}".format(n))
print("Resultado esperado: {}".format(es_primo(n)))
print("Resultado obtenido: {} \n".format(AKS(n)))

n = 524288
print("Número: {}".format(n))
print("Resultado esperado: {}".format(es_primo(n)))
print("Resultado obtenido: {} \n".format(AKS(n)))

n = 524289
print("Número: {}".format(n))
print("Resultado esperado: {}".format(es_primo(n)))
print("Resultado obtenido: {} ".format(AKS(n)))
︡ec29879f-73cb-45ec-8a38-1461cac266c4︡{"stdout":"Número: 105929\n"}︡{"stdout":"Resultado esperado: Primo\n"}︡{"stdout":"Resultado obtenido: Primo \n\n"}︡{"stdout":"Número: 524287\n"}︡{"stdout":"Resultado esperado: Primo\n"}︡{"stdout":"Resultado obtenido: Primo \n"}︡{"stdout":"\n"}︡{"stdout":"Número: 524288\n"}︡{"stdout":"Resultado esperado: Compuesto\n"}︡{"stdout":"Resultado obtenido: Compuesto \n\n"}︡{"stdout":"Número: 524289\n"}︡{"stdout":"Resultado esperado: Compuesto\n"}︡{"stdout":"Resultado obtenido: Compuesto \n"}︡{"done":true}
︠3206d1b0-b3db-4fb9-acdf-91786038ac2f︠









