from numpy import array, zeros, log10, polyfit
from numpy.linalg import norm
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

def Cauchy_error_inicial(F, Scheme, U0, t):

    N = len(t) - 1
    a = t[0]
    b = t[N]
    F_Error = zeros((N + 1, len(U0)))
    t1 = t
    t2 = particion(a, b, 2 * N)

    U1 = Cauchy(U0, t1, Scheme, F) 
    U2 = Cauchy(U0, t2, Scheme, F) 
    for i in range(0, N + 1):
        F_Error[i, :] = U2[2 * i, :] - U1[i, :]

    return U1, F_Error

def Cauchy_error (F, Scheme, U0, t):
    N1 = len(t) - 1
    U1, F_Error1 = Cauchy_error_inicial(F, Scheme, U0, t)
    a = t[0]

    b = t[len(t) - 1]
    N2 = (len(t) -1) * 2
    t2 = particion(a, b, N2)
    U2, F_Error2 = Cauchy_error_inicial(F, Scheme, U0, t2)

    q = - (log10(norm(F_Error2[N2, :])) - log10(norm(F_Error1[N1, :])))/(log10(N2) - log10(N1))

    print(q)

    Error = F_Error1/(1 - 1/(2 ** q))

    return U1, Error

def Convergencia (U0, Scheme, F):
    a = 0.
    b = 10

    N2 = 8000
    N1 = 200

    N_array = zeros(((N2 - N1)//300))
    U = zeros(((N2 - N1)//300))

    for N in range (N1, N2, 300):
        t1 = particion(a, b, N)
        t2 = particion(a, b, 2 * N)
        
        U1 = Cauchy(U0, t1, Scheme, F)
        U2 = Cauchy(U0, t2, Scheme, F)
        
        U[(N - N1)//300] = log10(norm(U2[2 * N, :] - U1[N, :]))

        N_array[(N - N1)//300] = log10(N)

    y = U[U > -12]

    orden, b = - polyfit(N_array[0:len(y)], y, 1)
    
    return U, N_array, orden


t1 = particion(a = 0., b = 10, N = 1000)

U1, Error1 = Cauchy_error(F = Oscilator, Scheme = Euler, U0 = array([1, 0]), t = t1)

U2, Error2 = Cauchy_error(F = Oscilator, Scheme = RK4, U0 = array([1, 0]), t = t1)

U3, Error3 = Cauchy_error(F = Oscilator, Scheme = Euler_inverso, U0 = array([1, 0]), t = t1)

U4, Error4 = Cauchy_error(F = Oscilator, Scheme = Crank_Nicolson, U0 = array([1, 0]), t = t1)


U_convergencia, N_convergencia, orden1 = Convergencia (U0 = array([1, 0]), Scheme = Euler, F = Oscilator)

U_convergencia2, N_convergencia2, orden2 = Convergencia (U0 = array([1, 0]), Scheme = RK4, F = Oscilator)

U_convergencia3, N_convergencia3, orden3 = Convergencia (U0 = array([1, 0]), Scheme = Euler_inverso, F = Oscilator)

U_convergencia4, N_convergencia4, orden4 = Convergencia (U0 = array([1, 0]), Scheme = Crank_Nicolson, F = Oscilator)

# plt.plot(t1, U1[:, 0])

# plt.show()

plt.plot(t1, Error1[:, 0], '-b'  , lw = 1.0, label = r"Euler")
plt.plot(t1, Error2[:, 0] , '-r'  , lw = 1.0, label = r"RK4")
plt.plot(t1, Error3[:, 0] , '--m' , lw = 1.0, label = r'Euler inverso')
plt.plot(t1, Error4[:, 0] , '--c' , lw = 1.0, label = r'Crank Nicolson')
plt.xlabel(r't')
plt.ylabel(r'Error')
plt.legend(loc='lower left')
plt.title( r'Error en la integración del oscilador armónico para los diferentes esquemas temporales' )
plt.grid()
plt.show()
plt.show()


plt.plot(N_convergencia, U_convergencia, '-b'  , lw = 1.0, label = r"{} (q={})".format(Euler.__name__, orden1))
plt.plot(N_convergencia2, U_convergencia2 , '-r'  , lw = 1.0, label = r"{} (q={})".format(RK4.__name__, orden2) )
plt.plot(N_convergencia3, U_convergencia3 , '--m' , lw = 1.0, label = r'{} (q={})'.format(Euler_inverso.__name__, orden3) )
plt.plot(N_convergencia4, U_convergencia4 , '--c' , lw = 1.0, label = r'{} (q={})'.format(Crank_Nicolson.__name__, orden4) )
plt.xlabel(r'$log_{10}(N)$')
plt.ylabel(r'$log_{10}(U_2 - U_1)$')
plt.legend(loc='lower left')
plt.title( r'Ratio de convergencia para los diferentes esquemas temporales' )
plt.grid()
plt.show()
