︠3e2576c9-aa93-4b74-8bad-2c264a052ea4s︠
#Algoritmo de primalidad de Miller-Rabin.

#Algoritmo para calcular rápidamente b^e mod m.
def fast_exp_mod(b,e,m):
    r = 1
    b = b%m
    while e>0:
        if e%2 == 1:
            r = (r*b)%m
        e = e//2
        b  = (b*b)%m
    return r


#Función que hace pasar al número, n, r rondas del test de primalidad
#Entrada: n>1 número impar, r número de vueltas del test

def probablemente_primo(n,rondas):
    if n == 2: return True
    if n<2 or n%2==0: return False

    #Calculamos k natural y m impar tales que: n-1 = (2^k)*m
    m = n-1
    k=0
    while m%2==0:
        k=k+1
        m=m//2

    for _ in range(rondas):
        #Tomamos un elemento aleatorio en [2, n-2)
        a = ZZ.random_element(2,n-2)
        #Calculamos a^m % n
        x = fast_exp_mod(a,m,n)
        if x==1 or x==n-1: continue #El test culmina con la conclusión de que n es probable primo y se sigue con la siguiente iteración para volver a realizar el test

        composite=True
        for _ in range(k-1):
            x = (x*x)%n
            if x==1: return False
            #Probablemente primo
            if x==n-1:
                composite=False
                break
        if composite: return False
    return True


print("¿{} es primo?: {}".format(15487457, probablemente_primo(15487457,3))) #Es primo
print("¿{} es primo?: {}".format(170141183460469231731687303715884105727, probablemente_primo(170141183460469231731687303715884105727,3))) #12º primo de mersenne
print("¿{} es primo?: {}".format(2635054262773128426267542653745695304430366239, probablemente_primo(2635054262773128426267542653745695304430366239,3))) #No es primo (producto de los dos anteriores)
print("¿{} es primo?: {}".format(3405514559213930613, probablemente_primo(3405514559213930613,3))) #3×29×41×311×397×4517×1711901
#Ahora vamos a ver cuántas veces falla el algoritmo aplicado al número de Carmichael, de un total de 10000 ejecuciones
a,b,c = 0,0,0

for i in range(10000):
    if probablemente_primo(3405514559213930613,1)==False: a=a+1
    if probablemente_primo(3405514559213930613,2)==False: b=b+1
    if probablemente_primo(3405514559213930613,3)==False: c=c+1
print("Fallos con 1 ronda: {}\nFallos con 2 rondas: {}\nFallos con 3 rondas: {}".format(10000 - a, 10000 - b, 10000 - c))
︡a0d5eefd-c2ec-4696-8860-0bd16ec769e4︡{"stdout":"¿15487457 es primo?: True\n"}︡{"stdout":"¿170141183460469231731687303715884105727 es primo?: True\n"}︡{"stdout":"¿2635054262773128426267542653745695304430366239 es primo?: False\n"}︡{"stdout":"¿3405514559213930613 es primo?: False\n"}︡{"stdout":"Fallos con 1 ronda: 0\nFallos con 2 rondas: 0\nFallos con 3 rondas: 0\n"}︡{"done":true}
︠473982c2-2876-4c14-8343-af8d98a50129︠









