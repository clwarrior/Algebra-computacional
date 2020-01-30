︠767dc188-c818-4697-970a-ba6d202eb98bs︠
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

#Recibe un polinomio mónico libre de cuadrados, f, perteneciente a F_q[x], de grado n > 0. q es una potencia de un primo impar.
#Devuelve o bien un factor g de f, o "fallo"
#Basado en el algoritmo 14.31 de Modern Computer Algebra
def berlekamp(f,q):
    R=f.base_ring()
    RX=f.parent()
    n=f.degree()
    pol=1
    #Calculamos x^q mod f
    aux=fast_exp_mod(x,q,f)
    Q=matrix(R,n,n)

    #Rellenamos Q de forma que Q[i,j]:= coeficiente de (x^j) en (x^(qi) mod f)

    for i in range(n):
        for j in range(n): Q[i,j]=0

    #En cada iteración calculamos (x^(q*i) mod f). Los valores de la fila corresponden
    #A los coeficientes del polinomio
    #Para i = 0, sólo el primer elemento de la fila es != 0
    Q[0,0]=1
    for i in range(1,n):
        pol=(pol*aux)%f
        for j in range(pol.degree()+1): Q[i,j]=pol[j]

    #Calculamos el Kernel de Q-I , y una base para el kernel, B.
    for i in range(n): Q[i,i]=Q[i,i]-1
    K=Q.kernel()
    B=K.basis_matrix()
    dim=K.dimension()

    if dim==1: return f

    #Tomamos un vector de dim elementos aleatorios
    C=random_vector(R,dim)

    #a = c1*b1 + ... + c_dim*b_dim
    a=RX([0])
    for i in range(dim): a=a+C[i]*RX(B[i].list())

    g1=gcd(a,f)
    if g1!=1 and g1!=f: return g1

    #b = a^((q-1)/2) mod f
    b=fast_exp_mod(a, ((q-1)//2), f)
    g2=gcd(b-1,f)

    if g2!=1 and g2!=f: return g2
    else: return "failure"

#EJEMPLO
#################################################

p = 7
R.<x>=PolynomialRing(GF(p))
f=R.random_element((1, 10))
assert f.is_squarefree()
print "f = ", f
g=berlekamp(f,p)
if g!="failure":
    print "g = ",g
    h=f/g
    print "h = ", h
    print "g*h = ", g*h
else: print "failure"
︡57441c95-1122-4411-abe3-47e49c6c7941︡{"stdout":"f =  6*x^8 + 5*x^7 + 5*x^6 + x^5 + 6*x^4 + 4*x^3 + 3*x^2 + 3\n"}︡{"stdout":"g =  x^2 + 2*x + 3\nh =  6*x^6 + x^4 + 6*x^3 + 5*x^2 + 4*x + 1\ng*h =  6*x^8 + 5*x^7 + 5*x^6 + x^5 + 6*x^4 + 4*x^3 + 3*x^2 + 3\n"}︡{"done":true}
︠f88a6d24-96f4-4e37-af08-9fe45e0fe379s︠
q = 9
R.<x>=PolynomialRing(GF(q))
f=R.random_element((1, 10))
assert f.is_squarefree()
print "f = ", f
g=berlekamp(f,q)
if g!="failure":
    print "g = ",g
    h=f/g
    print "h = ", h
    print "g*h = ", g*h
else: print "failure"
︡1c4daf4d-0710-4ded-8f8b-0f8f7e183616︡{"stdout":"f =  2*x^10 + 2*x^9 + z2*x^8 + 2*x^7 + x^6 + (z2 + 1)*x^5 + 2*x^4 + 2*x^3 + (2*z2 + 1)*x^2 + 2*x + 2*z2 + 2\n"}︡{"stdout":"g =  x + z2 + 2\nh =  2*x^9 + (z2 + 1)*x^8 + 2*x^6 + z2*x^5 + z2*x^4 + x^3 + 2*z2*x^2 + (2*z2 + 2)*x + z2 + 2\ng*h =  2*x^10 + 2*x^9 + z2*x^8 + 2*x^7 + x^6 + (z2 + 1)*x^5 + 2*x^4 + 2*x^3 + (2*z2 + 1)*x^2 + 2*x + 2*z2 + 2\n"}︡{"done":true}
︠6bd6f9a5-45a2-4ca9-b3ff-8afa0272571c︠









