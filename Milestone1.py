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
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Keppler orbits trajectory')


#Esquema numérico RK4 U_n+1 = U_n + Δt / 6 * (k1 + 2 * k2 + 2 * k3 + k4) con Keppler

for i in range(0, N):
    r = (x[i] ** 2 + y[i] ** 2) ** (3/2) 

    k11 = xdot[i]
    k12 = ydot[i]
    k13 = - deltat * x[i] / r
    k14 = - deltat * y[i] / r

    k21 = xdot[i] + k13 * deltat / 2
    k22 = ydot[i] + k14 * deltat / 2
    k23 = - (x[i] + k11 * deltat / 2) / (((x[i] + k11 * deltat / 2) ** 2 + (y[i] + k12 * deltat / 2) ** 2) ** (3 / 2))
    k24 = - (y[i] + k12 * deltat / 2) / (((x[i] + k11 * deltat / 2) ** 2 + (y[i] + k12 * deltat / 2) ** 2) ** (3 / 2))
    
    k31 = xdot[i] + k23 * deltat / 2
    k32 = ydot[i] + k24 * deltat / 2
    k33 = - (x[i] + k21 * deltat / 2) / (((x[i] + k21 * deltat / 2) ** 2 + (y[i] + k22 * deltat / 2) ** 2) ** (3 / 2))
    k34 = - (y[i] + k22 * deltat / 2) / (((x[i] + k21 * deltat / 2) ** 2 + (y[i] + k22 * deltat / 2) ** 2) ** (3 / 2))
    
    k41 = xdot[i] + k33 * deltat
    k42 = ydot[i] + k34 * deltat
    k43 = - (x[i] + k31 * deltat) / (((x[i] + k31 * deltat) ** 2 + (y[i] + k32 * deltat) ** 2) ** (3 / 2))
    k44 = - (y[i] + k32 * deltat) / (((x[i] + k31 * deltat) ** 2 + (y[i] + k32 * deltat) ** 2) ** (3 / 2))

    x[i+1] = x[i] + (k11 + 2 * k21 + 2 * k31 + k41) * deltat / 6

    y[i+1] = y[i] + (k12 + 2 * k22 + 2 * k32 + k42) * deltat / 6

    xdot[i+1] = xdot[i] + (k13 + 2 * k23 + 2 * k33 + k43) * deltat / 6

    ydot[i+1] = ydot[i] + (k14 + 2 * k24 + 2 * k34 + k44) * deltat / 6

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Keppler orbits trajectory')
plt.legend(["Euler", "RK4"])
plt.show()