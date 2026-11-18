import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

if __name__=="__main__":
    L=1
    m=1
    M=10
    g=10
    dt=0.001


    omegas=np.linspace(1.9705,1.9735,100)
    x_1_fin=[]
    F_trace=[]
    for o in tqdm(omegas):
        theta=-np.pi/2
        x=0
        theta_1=0
        x_1=0
        F_vec=[]
        t_fin=np.pi/o
        n_steps=2000
        dt=t_fin/n_steps
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
    F_trace=np.array(F_trace)
    plt.plot(omegas,x_1_fin)
    plt.xlabel("$\omega$")
    plt.ylabel("$x_{\omega}(T)$")
    plt.show()
    plt.pcolormesh(F_trace)
    plt.show()
            



