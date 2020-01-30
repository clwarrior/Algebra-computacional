︠c91c4d4e-b105-4f9b-acf2-b4fdfaff72bd︠
#Algoritmo 3.6 de Modern Computer Algebra por Joachim von zur Gathen y Jürgen Gerhard

#División en el anillo de los enteros de Gauss
def div_gauss(a, b):
    a, b = ZZ[i](a), ZZ[i](b)
    q, _ = a.quo_rem(b)
    Q = q.imag().round()*I + q.real().round()
    r = a - Q*b
    return Q, r

def div(a, b): return a.quo_rem(b)

#Métodos que nos dan como salida el mcd y los coeficientes de la Identidad de Bezout como una terna para distintos DE
def ea_ext_int(a, b): return Euclidean_Algorithm_Ext(a, b, div, Integer)
def ea_ext_zp(p, a, b):
    Zp = Integers(p)
    return Euclidean_Algorithm_Ext(Integer(Zp(a)), Integer(Zp(b)), div, Integer)
def ea_ext_fp(a, b, R): return Euclidean_Algorithm_Ext(a, b, div, R)
def ea_ext_qx(a, b): return Euclidean_Algorithm_Ext(a, b, div, QQ['x'])
def ea_ext_g(a, b): return Euclidean_Algorithm_Ext(a, b, div_gauss, ZZ[I])

#Métodos que nos imprimen un string con la Identidad de Bezout escrita de forma explícita para distintos DE
def main_int(a, b): return main(a, b, div, Integer)
def main_zp(p, a, b):
    Zp = Integers(p)
    aZp = Zp(a)
    bZp = Zp(b)
    r0, s0, t0 = Euclidean_Algorithm_Ext(Integer(aZp),Integer(bZp),div,Integer)
    print str(r0) + " = (" + str(s0) + ") * (" + str(a) + ") + (" + str(t0) + ") * (" + str(b) + ")"
def main_fp(a, b, R): return main(a, b, div, R)
def main_qx(a, b): return main(a, b, div, QQ['x'])
def main_g(a, b): return main(a, b, div_gauss, ZZ[I])

#El algoritmo de Euclides extendido se basa, al igual que el de Euclídes, en la igualdad mcd(a,b) = mcd(a, b%a). Ahora además vamos a ir obtener los coeficientes de la Identidad de Bezout. 
#Durante todo el algoritmo se cumple ri = a⋅si + b⋅ti para i = 0,1

#Entrada: a,b en R dominio euclídeo, como tercer parámetro recibe la operación de división del DE y como cuarto el DE
#Salida: r,s,t en R que cumple que r = mcd(a,b) y r = a*s + b*t
def Euclidean_Algorithm_Ext(a,b,f,K):
    #Inicialización, buscamos mcd(r0,r1)
    r0=a
    r1=b
    s0=1
    s1=0
    t0=0
    t1=1
    while (r1 != 0):
        #Dividimos r0//r1
        q,_ = f(r0, r1)
        #Guardamos los antiguos valores
        auxr1 = r1
        auxs1 = s1
        auxt1 = t1
        #Actualizamos los valores de r1,s1,t1
        r1=K(r0-q*r1)
        s1=s0-q*s1
        t1=t0-q*t1
        #Actualizamos los valores de r0,s0,t0
        r0 = auxr1
        s0 = auxs1
        t0 = auxt1
    return r0, s0, t0

#Método auxiliar que explicita la Identidad de Bezout
def main(a,b,f,K):
    r0, s0, t0 = Euclidean_Algorithm_Ext(a,b,f,K)
    print str(r0) + " = (" + str(s0) + ") * (" + str(a) + ") + (" + str(t0) + ") * (" + str(b) + ")"


#Para cada ejemplo usaremos el método tipo main que nos dará de forma explícita la Identidad de Bezout y el método que nos devuelve la terna de los tres números (lo hemos implementado porque puede ser de utilidad a la hora de implementar otros algoritmos)


#Enteros
print("Datos de entrada: {}, {}".format(15, 35))
main_int(15, 35)
ea_ext_int(15, 35)
︡71bbdc37-41e8-416c-9300-2b6b8ff98cf7︡{"stdout":"Datos de entrada: 15, 35\n"}︡{"stdout":"5 = (-2) * (15) + (1) * (35)\n"}︡{"stdout":"(5, -2, 1)\n"}︡{"done":true}
︠44b28448-35d1-4766-b51d-123129bd6571s︠

#Polinomios con coeficientes racionales
R.<x> = QQ[]
print("Datos de entrada: {}, {}".format(x^2+x, x))
main_qx(x^2+x, x)
ea_ext_qx(x^2+x, x)
︡b8ce76f6-1bdc-4014-9062-ff1e83aad12c︡{"stdout":"Datos de entrada: x^2 + x, x\n"}︡{"stdout":"x = (0) * (x^2 + x) + (1) * (x)\n"}︡{"stdout":"(x, 0, 1)\n"}︡{"done":true}
︠ca77ada7-9491-4213-aa04-b8ccd1bb8586s︠

#Elementos en anillos de polinomios sobre cuerpos finitos
R.<x> = PolynomialRing(FiniteField(7))
f = R.random_element(5)
g = R.random_element(3)
print("Datos de entrada: {}, {}".format(f,g))
main_fp(f, g, R)
ea_ext_fp(f, g, R)
︡dc32b17c-9d07-41e9-a426-6215bac6b688︡{"stdout":"Datos de entrada: 3*x^5 + x^4 + 5*x^3 + 6*x^2 + 4*x + 2, x^3 + 4*x^2 + 2*x + 6\n"}︡{"stdout":"2 = (x^2 + 6*x + 2) * (3*x^5 + x^4 + 5*x^3 + 6*x^2 + 4*x + 2) + (4*x^4 + 3*x^2 + 3*x + 2) * (x^3 + 4*x^2 + 2*x + 6)\n"}︡{"stdout":"(2, x^2 + 6*x + 2, 4*x^4 + 3*x^2 + 3*x + 2)\n"}︡{"done":true}
︠9b0de4ae-9dcf-4b07-84ac-65c5b47d606as︠

#Grupo cíclico Zp para p primo
p = 5 #Sea p cualquier primo
Zp = Integers(p)
aux1 = ZZ.random_element(0, 5*p)
aux2 = ZZ.random_element(0, 5*p)
print("Datos de entrada: {}, {}".format(aux1,aux2))
a = Zp(aux1)
b = Zp(aux2)
main_zp(p, aux1, aux2)
ea_ext_zp(p, a, b)
︡6cff01bb-1ddb-40fb-bbdc-b0926a1edcb8︡{"stdout":"Datos de entrada: 13, 8\n"}︡{"stdout":"3 = (0) * (13) + (1) * (8)\n"}︡{"stdout":"(3, 0, 1)\n"}︡{"done":true}
︠f2f0c4ef-1359-4508-856f-368a6e71f72bs︠
#Enteros de Gauss
G = ZZ[I]
a = G.random_element()
b = G.random_element()
print("Datos de entrada: {}, {}".format(a,b))
main_g(a, b)
ea_ext_g(a, b)
︡f68db6ea-52b8-45b5-a363-0d7aebff82a4︡{"stdout":"Datos de entrada: I, -4*I + 2\n"}︡{"stdout":"I = (1) * (I) + (0) * (-4*I + 2)"}︡{"stdout":"\n"}︡{"stdout":"(I, 1, 0)\n"}︡{"done":true}
︠b26bafa9-cba5-4b25-af1a-b8cf21d0d8b4s︠

#Otro ejemplo en los enteros de Gauss
aa = G(I-3)
bb = G(I+5)
print("Datos de entrada: {}, {}".format(aa,bb))
main_g(aa,bb)
ea_ext_g(aa,bb)
︡23ca0017-9f6e-48fc-8642-c96b886c900a︡{"stdout":"Datos de entrada: I - 3, I + 5\n"}︡{"stdout":"-I - 1 = (I - 2) * (I - 3) + (I - 1) * (I + 5)\n"}︡{"stdout":"(-I - 1, I - 2, I - 1)\n"}︡{"done":true}
︠bd9d6265-aaa2-4040-938a-101f340de4b4︠









