from numpy import array, zeros, linspace
import matplotlib.pyplot as plt
#from scipy.optimize import newton

def Keppler(U, t):
    x = U[0]
    y = U[1]
    xdot = U[2]
    ydot = U[3]
    r = (x ** 2 + y ** 2) ** 1.5
    return array([xdot, ydot, - x / r, - y / r])

def Oscilator(U,t):
    x = U[0]
    xdot = U[1]
    return array([xdot, - x])

def Euler(U, dt, t, F):
    return U + dt * F(U, t)

# def Euler_inverso(U, dt, t, F):
#    def G(X):
#       return X - U - dt * F(X, t)   
#    return newton(G, U)

def Cauchy(U0, t, Scheme, F):
    N = len(t)
    U = zeros((N + 1, len(U0)))
    U[0, :] = U0
    for i in range(0, N-1):
        dt = t[i + 1] - t[i]
        U[i + 1, :] = Scheme(U[i, :], dt, t, F)
    return (U)

N = 800
dt = 0.1

U0 = array([0, 1])

t = linspace(0, N * dt, N)

U = Cauchy(U0, t, Euler, Oscilator)


y = t = linspace(0, N * dt, N)
print (len(U))
print (len(y))

plt.plot(U[:, 0], y)
plt.show()





