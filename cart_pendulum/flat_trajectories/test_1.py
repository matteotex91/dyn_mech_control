import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

if __name__=="__main__":
    L=1
    m=1
    M=10
    g=10
    dt=0.00001


    omegas=np.linspace(1.9734,1.9735,50)
    x_1_fin=[]
    F_trace=[]
    for o in tqdm(omegas):
        theta=-np.pi/2
        x=0
        theta_1=0
        x_1=0
        F_vec=[]
        # theta_vec=np.array([])
        # x_vec=np.array([])
        # theta_1_vec=np.array([])
        # x_1_vec=np.array([])
        t_fin=np.pi/o
        n_steps=int(t_fin/dt)
        for i in range(n_steps):
            angle=o*i*dt
            pi_cos_angle=np.cos(angle)
            pi_sin_angle=np.sin(angle)
            x2=((2*L*o*o*pi_cos_angle)/3+g*np.cos(pi_cos_angle/2))/np.sin(pi_cos_angle/2)
            x_1+=x2*dt
            x+=(0.5*x2+x_1)*dt
            F_vec.append((m+M)*x2-m*L*o*o*pi_cos_angle*np.sin(pi_cos_angle/2)/2-m*L*np.square(o*pi_sin_angle/2)*np.cos(pi_cos_angle/2))
        x_1_fin.append(x_1)
        F_trace.append(F_vec)
    plt.plot(omegas,x_1_fin)
    plt.xlabel("$\omega$")
    plt.ylabel("$x_{\omega}(T)$")
    plt.show()
    plt.pcolormesh(F_trace)
    plt.show()
            



