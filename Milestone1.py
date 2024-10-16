from numpy import array, zeros, linspace
import matplotlib.pyplot as plt

#N número de saltos temporales 

N = 2000

#Se declara el vector de estado U y se inicializa la componente inicial U_0

x = zeros(N + 1)
y = zeros(N + 1)
xdot = zeros(N + 1)
ydot = zeros(N + 1)
deltat = 0.1

x[0] = 1
y[0] = 0
xdot[0] = 0
ydot[0] = 1

#Esquema numérico Euler U_n+1 = U_n + Δt * F_n con Keppler

for i in range(1, N + 1):
    r = (x[i - 1] ** 2 + y[i - 1] ** 2) ** (3 / 2) 

    x[i] = x[i - 1] + deltat * xdot[i - 1] 
    y[i] = y[i - 1] + deltat * ydot[i - 1]

    xdot[i] = xdot[i - 1] + deltat * (- x[i - 1] / r)
    ydot[i] = ydot[i - 1] + deltat * (- y[i - 1] / r)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Keppler orbits trajectory')


#Esquema numérico RK U_n+1 = U_n + Δt / 6 * (k1 + 2 * k2 + 2 * k3 + k4) con Keppler

for i in range(1, N + 1):
    r = (x[i - 1] ** 2 + y[i - 1] ** 2) ** (3/2) 

    k1 = xdot[i-1]
    k2 = xdot[i-1] + k1 * deltat / 2
    k3 = xdot[i-1] + k2 * deltat / 2
    k4 = xdot[i-1] + k3 * deltat

    x[i] = x[i - 1] + deltat / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    k1 = ydot[i-1]
    k2 = ydot[i-1] + k1 * deltat / 2
    k3 = ydot[i-1] + k2 * deltat / 2
    k4 = ydot[i-1] + k3 * deltat

    y[i] = y[i - 1] + deltat / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    k1 = - x[i - 1] / r
    k2 = - (x[i - 1] + deltat / 2 * k1) / ((x[i-1] + deltat / 2 * k1) ** 2 + (y[i-1] + deltat / 2 * k1) ** 2) ** (3 / 2)
    k3 = - (x[i - 1] + deltat / 2 * k2) / ((x[i-1] + deltat / 2 * k2) ** 2 + (y[i-1] + deltat / 2 * k2) ** 2) ** (3 / 2)
    k4 = - (x[i - 1] + deltat * k3) / ((x[i-1] + deltat * k3) ** 2 + (y[i-1] + deltat * k3) ** 2) ** (3 / 2)

    xdot[i] = xdot[i - 1] + deltat / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    k1 = - y[i - 1] / r
    k2 = - (y[i - 1] + deltat / 2 * k1) / ((x[i-1] + deltat / 2 * k1) ** 2 + (y[i-1] + deltat / 2 * k1) ** 2) ** (3 / 2)
    k3 = - (y[i - 1] + deltat / 2 * k2) / ((x[i-1] + deltat / 2 * k2) ** 2 + (y[i-1] + deltat / 2 * k2) ** 2) ** (3 / 2)
    k4 = - (y[i - 1] + deltat * k3) / ((x[i-1] + deltat * k3) ** 2 + (y[i-1] + deltat * k3) ** 2) ** (3 / 2)

    ydot[i] = ydot[i - 1] + deltat / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Keppler orbits trajectory')
plt.legend(["Euler", "RK4"])
plt.show()