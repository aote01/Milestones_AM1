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

#Esquema numérico de Euler inverso U_n+1 = U_n+1 + Δt * F_n

def Euler_inverso(U, dt, t, F):
    def G(X):
       return X - U - dt * F(U, t)   
    return newton(G, U)

#Esquema numérico de Euler inverso U_n+1 = U_n+1 + Δt / 2 * (F_n+1 + F_n)

def Crank_Nicolson(U, dt, t, F):
    def G(X):
       return X - U - dt / 2 * (F(U, t) + F(X, t))  
    return newton(G, U)

#Problema físico de Keppler

def Keppler(U, t):
    return array([U[2], U[3], - U[0] / (U[0] ** 2 + U[1] ** 2) ** 1.5, - U[1] / (U[0] ** 2 + U[1] ** 2) ** 1.5])

#Problema físico de un oscilador

def Oscilator(U,t):
    x = U[0]
    xdot = U[1]
    return array([xdot, - x])


#Problema de Cauchy

def Cauchy(U0, t, Scheme, F):
    N = len(t)
    U = zeros((N, len(U0)))
    U[0, :] = U0
    for i in range(0, N-1):
        dt = t[i + 1] - t[i]
        U[i + 1, :] = Scheme(U[i, :], dt, t, F)
    return (U)

#Se inicializa las diferentes variables y condiciones iniciales

N = 5000
dt = 20/N

t = linspace(0, N * dt, N+1)

Problema_fisico = Keppler

U0 = [1, 0, 0, 1]

U_Euler = Cauchy(U0, t, Euler, Problema_fisico)
U_RK4 = Cauchy(U0, t, RK4, Problema_fisico)
U_Euler_inverso = Cauchy(U0, t, Euler_inverso, Problema_fisico)
U_Crank_nicolson = Cauchy(U0, t, Crank_Nicolson, Problema_fisico)

plt.plot(U_Euler[:, 0], U_Euler[:, 1])
plt.plot(U_RK4[:, 0], U_RK4[:, 1])
plt.plot(U_Euler_inverso[:, 0], U_Euler_inverso[:, 1])
plt.plot(U_Crank_nicolson[:, 0], U_Crank_nicolson[:, 1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Keppler orbits trajectory')
plt.legend(["Euler", "RK4", "Euler inverso", "Crank-Nicolson"])
plt.show()
