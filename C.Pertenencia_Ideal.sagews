︠891475ea-24a3-4cc2-a030-757ee4c441acs︠
#UTILIZAMOS COMO ALGORITMO AUXILIAR EL ALGORITMO DE DIVISIÓN EN VARIAS VARIABLES
#######################################################################
def divide(f, pols):
    r = 0
    p = f
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
#######################################################################

#Un elemento f pertenece a un ideal I, si al calcular la base de groebner del ideal
#f % base_groebner(I) == 0. (Pues eso quiere decir que se puede expresar en términos de los elementos de la base)
def pertenece(f, I):
    B = I.groebner_basis()
    _, r = divide(f, B)
    if r == 0: return True
    else: return False


#Ejemplo
#######################################################################
R.<x,y,z> = PolynomialRing(QQ,3)
I = Ideal([x^2+y+z-1,x+y^2+z-1,x+y+z^2-1])
pertenece(x^2+y+z-1, I)
︡86a8c143-46c1-4b85-a42f-1d9b53743047︡{"stdout":"True\n"}︡{"done":true}
︠82f26246-a5ea-4c84-b74f-232b6faf9058︠











