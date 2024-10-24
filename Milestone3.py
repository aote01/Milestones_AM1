from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
from scipy.optimize import newton

#Esquema numérico Euler U_n+1 = U_n + Δt * F_n

def Euler(U, dt, t, F):

    return U + dt * F(U, t)

#Esquema numérico RK4 U_n+1 = U_n + Δt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

def RK4(U, dt, t, F):

    k1 = F(U, t)
    k2 = F(U + dt * k1 / 2, t + dt / 2)
    k3 = F(U + dt * k2 / 2, t + dt / 2)
    k4 = F(U + dt * k3, t + dt)

    return U + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

#Esquema numérico de Euler inverso U_n+1 = U_n + Δt * F_n+1

def Euler_inverso(U, dt, t, F):

    def G(X):
       return X - U - dt * F(X, t) 
    
    return newton(G, U)

#Esquema numérico de Crank Nicolson U_n+1 = U_n + Δt / 2 * (F_n+1 + F_n)

def Crank_Nicolson(U, dt, t, F):

    def G(X):
       return X - U - dt / 2 * (F(U, t) + F(X, t))  
    
    return newton(G, U)

#Problema físico de Keppler

def Keppler(U, t):

    return array([U[2], U[3], - U[0] / (U[0] ** 2 + U[1] ** 2) ** 1.5, - U[1] / (U[0] ** 2 + U[1] ** 2) ** 1.5])

#Problema físico de un oscilador

def Oscilator(U,t):

    return array([U[1], - U[0]])


#Problema de Cauchy

def Cauchy(U0, t, Scheme, F):

    N = len(t)
    U = zeros((N, len(U0)))
    U[0, :] = U0
    for i in range(0, N-1):
        dt = t[i + 1] - t[i]
        U[i + 1, :] = Scheme(U[i, :], dt, t, F)

    return (U)

#Refinar malla
#Dada la partición t1 con N+1 'puntos' obtiene la partición de t2 que tiene con 2N+1 'puntos'
#Los nodos pares en t2 seguirán siendo los mismos de t1, los impares serán los puntos medios de los pares

def refinar_malla(t1):

    N = len(t1) - 1
    t2 = zeros(2 * N + 1)
    for i in range(0, N + 1):
        t2[2 * i] = t1[i] #pares
        t2[2 * i + 1] = (t1[i + 1] + t1[i] / 2)#impares
    t2[2 * N] = t1[N]

    return t2

#Hace una partición equiespaciada en N trozos de un segmento de la recta real entre a y b

def particion(a, b, N):

    t = zeros(N + 1)

    for i in range (0, N+1):
        t[i] = a + (b - a) / N * i

    return t

#Método de Richardson 

def Cauchy_error(F, Scheme, U0, t):

    N = len(t) - 1
    a = t[0]
    b = t[N]
    Error = zeros((N + 1, len(U0)))
    t1 = t
    t2 = particion(a, b, 2 * N)

    U1 = Cauchy(U0, t1, Scheme, F) 
    U2 = Cauchy(U0, t2, Scheme, F) 

    for i in range (0, N + 1):
        Error[i, :] = U2[2 * i, :] - U1[i, :]

    return U1, Error

t1 = particion(a = 0., b = 10, N = 1000)
U1, Error = Cauchy_error(F = Oscilator, Scheme = Euler, U0 = array([1, 0]), t = t1)

plt.plot(t1, U1[:, 0])


plt.plot(t1, Error[:, 0])

plt.show()

#a, b = 0., 1.
#N = 5

#t = linspace(a, b, N + 1)
#print(t)

#t1 = particion(a, b, N)
#print(t1)

#t2 = refinar_malla(t1)
#print(t2)

#t2 = particion(a, b, N)
#print(t2)

#