import numpy as np
import math
import matplotlib.pyplot as plt
N = 20
eps = 0.001
T = 14
dt = 0.1
L = 1
nu = 0.01
t = 0
h = L / N
#Создание матрицы для давлений:
def createA(N):
    A = np.zeros((N**2, N**2))
    for i in range(N):
        for j in range(N):
            ij = i*N+j
            count = 0

            if i!= 0:  # верхняя граница

                A[ij][ij - N] = -1
                A[ij][ij] += 1
                count += 1

            if j  != 0:  # левая

                A[ij][ij - 1] = -1
                count += 1

            if j!= N-1:  # правая

                A[ij][ij + 1] = -1
                count += 1

            if i!= N-1:  # нижняя

                A[ij][ij + N] = -1
                count += 1

            A[ij][ij] = count
    return(A)     
  
  def p_solve(A,b):
    return np.linalg.solve(A,b)
def diverg(u,v,h):
    b = np.zeros(N*N)
    for i in range(N):
        for j in range(N):
            b[i*N + j] = - h*(u[i][j+1]-u[i][j]+v[i+1][j]-v[i][j])
    return b
  
  
 def uv_solve(u_prev, v_prev, p_prev, N, dt, h, nu):
    uu = np.zeros((N, N + 1)) 
    vv = np.zeros((N + 1, N)) 
    for i in range (0, N):  
        for j in range (1, N): 
            uc = u_prev[i][j]
            vwn = v_prev[i][j - 1]
            vne = v_prev[i][j]
            vws = v_prev[i + 1][j - 1]
            ves = v_prev[i + 1][j]

            if i == 0:
                un = 2 - uc
            else:
                un = u_prev[i - 1][j]

            if i == (N-1):
                us = - uc
            else:
                us = u_prev[i + 1][j]

            if j == (0):
                uw = - uc

            else:
                uw = u_prev[i][j - 1]

            if j == (N):
                ue = - uc
            else:
                ue = u_prev[i][j + 1]
            
            pe = p_prev[i * N + j]
            pw = p_prev[i * N + j - 1]
                
            gradU = 0.25 / h * ((uc + ue) * (uc + ue) - (uw + uc) * (uw + uc) - (vwn + vne) * (un + uc) + (vws + ves) * (us + uc) )

            uu[i, j] = uc - dt * (gradU + nu / (h ** 2) * (4 * uc - uw - ue - us - un) + (pe - pw) * h)
#             uu[i][j] = uc - dt * (0.25 * (-(uw + uc)**2 + (uc + ue) ** 2 + (vws + ves) * (us + uc) - (vwn + vne) * (un + uc)) / h + nu * (4 * uc - ue - un - us - uw) / h / h + (p_prev[i*N + j] - p_prev[i*N + j - 1]) / h)

    for i in range (1, N):   
        for j in range (0, N):
            vc = v_prev[i][j]
            vn = v_prev[i-1][j]
            vs = v_prev[i+1][j]
            uwn = u_prev[i-1][j]
            uws = u_prev[i][j]
            une = u_prev[i-1][j+1]
            ues = u_prev[i][j+1]
            
            if (j == 0):
                vw = - vc
            else:
                vw = v_prev[i][j-1]
            if (j == N - 1):
                ve = - vc
            else:        
                ve = v_prev[i][j+1]
            pn = p_prev[(i - 1) * N + j]
            ps = p_prev[i * N + j]

#             vv[i][j] = vc - (h * 0.25 *((vn + vc)**2 - (vc + vs) ** 2 + (ve + vc) * (une + ues) - (vc + vw) * (uwn + uws))\
#                         + nu * (4 * vc - ve - vn - vs - vw) / h / h + (p_prev[i*N + j] - p_prev[(i-1)*N + j]) * h) * dt
            gradV = 0.25 / h * ((une + ues) * (ve + vc) - (uwn + uws) * (vc + vw) - ((vn + vc)) ** 2 + ((vs + vc)) ** 2 )
            vv[i, j] = vc - dt * (gradV + nu / h ** 2 * (4 * vc - vw - ve - vs - vn) + (ps - pn) * h)
#             vv[i][j] = vc - dt * (0.25 *((vn + vc)**2 - (vc + vs) ** 2 + (ve + vc) * (une + ues) - (vc + vw) * (uwn + uws)) / h  + nu * (4 * vc - ve - vn - vs - vw) / h / h + (p_prev[i*N + j] - p_prev[(i-1)*N + j]) / h)

    return(uu, vv)

def plot_solution(u, v, p, dt,streamplot=True):
        u = (u[:, :-1] + u[:, 1:]) / 2
        v = (v[1:, :] + v[:-1, :]) / 2
        print(u.shape)
        u = u[::-1, ::]
        v = -v[::-1, ::]
        p=p.reshape((N,N))[::-1,::]
        x = np.arange(h / 2, 1, h)
        y = np.arange(h / 2, 1, h)
        # y = np.arange(-1, -self.h/2, self.h)
        grid_x, grid_y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(10, 10))
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.streamplot(grid_x, grid_y, u, v, color='black')
        plt.contourf(grid_x, grid_y, p.reshape((N, N)))
        plt.title(f"Задача о каверне\nN = {N},E ={eps},Nu ={nu},T={T},dt ={dt},a = {alpha} ", fontsize=20)
#         plt.show()
        plt.savefig(f'images/five/nN = {N},E ={eps},Nu ={nu},T={T},dt ={dt},a = {alpha}.png')
  
  A = createA(N)
u_prev = np.zeros((N, N + 1))
v_prev = np.zeros((N + 1, N))
p_prev = np.zeros((N * N))
alpha = 0.8
# alpha_u = 0.85
while t < T:
    pp = p_prev
    iter = 0
    norm_b = 100
    while norm_b > eps:
        uu, vv = uv_solve(u_prev, v_prev, pp, N, dt, h, nu)
#         uu = alpha_u*uu+(1-alpha_u)*u_prev
#         vv = alpha_u*vv+(1-alpha_u)*v_prev
        b = diverg(uu, vv, h)
        norm_b = np.linalg.norm(b)
        print(norm_b)
        if norm_b > eps:
            p_delta = p_solve(A, b/dt)
            pp = pp + alpha*p_delta
            iter += 1

    u_prev, v_prev = uu, vv
    p_prev = pp

    print("iter ", iter)
    print('t: ', t)
    t += dt

plot_solution(u_prev, v_prev, p_prev, dt)  
