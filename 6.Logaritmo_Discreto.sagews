︠e77da323-9be9-4968-aa4c-07f8e356d9e1︠
# Algoritmo para calcular el logaritmo discreto sobre Zq
def pollard_step(q, alpha, beta, x,a,b):
    opt = x%3
    if opt == 0:
        return (x*x)%q, (2*a)%(q-1), (2*b)%(q-1)
    elif opt == 1:
        return (alpha*x)%q, (a+1)%(q-1), b
    else:
        return (beta*x)%q, a, (b+1)%(q-1)

#Calcula el logaritmo discreto utilizando el algoritmo rho de pollard.
#Trabajamos en el grupo multiplicativo Zq*
#Devuelve un valor (sol) tal que base^sol = val mod q
def discrete_log_pollard(q, base, val):

    x1 = 1; a1=0; b1 = 0;
    x2 = 1; a2=0; b2 = 0;

    while(True):
        x1, a1, b1 = pollard_step(q, base, val, x1,a1,b1)
        x2, a2, b2 = pollard_step(q, base, val, x2,a2,b2)
        x2, a2, b2 = pollard_step(q, base, val, x2,a2,b2)
        if (x1 == x2): break
    if b1 == b2: return None
    #Ahora vamos a resolver la ecuación a+b*c == 0 módulo q-1
    var('c', domain='integer')
    a = a1 - a2; b = b1 - b2
    sol = solve_mod(a + b*c == 0, q-1)[0][0]
    return sol

#EJEMPLOS
#########################################################################
q = 1019
base = 2
valor = 5
assert GF(q)(base).order() == q
sol = discrete_log_pollard(q, base, valor)
print("log_{}({})={} mod {}".format(base, valor, sol, q))
print("({}^{} - {}) mod {} = {}".format(base, sol, valor, q, (base^sol - valor)%q))
︡56190717-8a0f-493a-bd56-377bad7dd757︡{"stdout":"log_2(5)=10 mod 1019\n"}︡{"stdout":"(2^10 - 5) mod 1019 = 0\n"}︡{"done":true}
︠9ac213f0-b079-438d-a6e0-beb110fb2309s︠

#Segundo ejemplo
q = 569
base = 5
valor = 262
assert GF(q)(base).order() == q
sol = discrete_log_pollard(q, base, valor)
print("log_{}({})={} mod {}".format(base, valor, sol, q))
print("({}^{} - {}) mod {} = {}".format(base, sol, valor, q, (base^sol - valor)%q))
︡be31e1b6-e6b9-4d0c-8722-29cdf9080070︡{"stdout":"log_5(262)=432 mod 569\n"}︡{"stdout":"(5^432 - 262) mod 569 = 0\n"}︡{"done":true}
︠4f3d7d86-e4a2-406b-975a-51a15515738a︠









