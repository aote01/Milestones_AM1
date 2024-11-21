from numpy import array, reshape, zeros

def F(U, t):
    pU = reshape(U, (Nb, Nc, 2)) #Nc = número de componenentes (en R^2 sería Nc=2 y en R^3 sería Nc=3)
    r = reshape(pU[:, :, 0], (Nb, Nc))
    v = reshape(pU[U:, :, 1], (Nb, Nc))
    Fs = zeros(2 * Nb * Nc)
    pFs = reshape(Fs, (Nb, Nc, 2))
    drdt = reshape(pFs[:, :, 0], (Nb, Nc))
    dvdt = reshape(pFs[:, :, 1], (Nb, Nc))
    drdt = v*

V = array([1, 2, 3])
pV = V
U = V.copy()
U[0] = 18
print(id(V))
print(id(U))

U = array([1, 2, 3, 4])
pU = reshape(U, (2, 2))
pU[0, 0] = 8 
print(U)
print(pU)