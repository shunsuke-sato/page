import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

num_expand = 101

num_grid = 512
length = 1.0


xj = np.linspace(0.0, 2*length, num_grid+1)


wf_gs = np.sqrt(2/length)*np.sin(np.pi*xj/length)*np.heaviside(length -xj,0.0)


# Ground state wavefunction
wf = np.zeros(num_grid+1)

for m in range(1, num_expand):
    if(m == 2):
        cm = np.sqrt(2.0)/2.0
    elif(m%2 == 1):
        cm = np.sqrt(2.0)*(-1)**((m+3)//2)*(4/(4-m**2))/np.pi
    else:
        cm = 0

#    else:
#        cm = 4*np.sin(m*np.pi/2)/((4-m**2)*np.pi)
#        
#    cm = np.sqrt(2)*cm
    wf = wf + cm*np.sqrt(1.0/length)*np.sin(m*np.pi*xj/(2*length))


plt.title("Wavefunction")
plt.plot(xj, wf_gs, label ='Ground state (exact)')
plt.plot(xj, wf, label ='Ground state (expand)', linestyle='dashed')

plt.xlabel("x")
plt.ylabel("wavefunction")
plt.legend()

plt.savefig('wavefunction.pdf')

# Time propagation
omega = np.pi**2/(8*length**2) # Period of the oscillation
Tprop = 2.0*np.pi/omega
nt = 400
dt = Tprop/nt




# For loop for the time propagation
density_list = []

for it in range(nt+1):
    print(it)
    tt = it*dt

    wf = np.zeros(num_grid+1, dtype=complex)
    
    for m in range(1, num_expand):
        if(m == 2):
            cm = 1/2
        else:
            cm = 4*np.sin(m*np.pi/2)/((4-m**2)*np.pi)

        Em = np.pi**2*m**2/(8*length**2)
        zcm = np.sqrt(2)*cm*np.exp(-1j*tt*Em)
        wf = wf + zcm*np.sqrt(1.0/length)*np.sin(m*np.pi*xj/(2*length))

    rho = np.abs(wf)**2
    density_list.append(rho.copy())

        
# Define function to update plot for each frame of the animation
def update_plot(frame):
    plt.cla()
    plt.xlim([0.0, 2.0])
    plt.ylim([0.0, 3])

    plt.plot(xj, density_list[frame], label="$|\psi(x)|^2$ (calc.)")
    
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.legend(loc = 'upper right')


# Create the animation
fig = plt.figure()
ani = animation.FuncAnimation(fig, update_plot, frames=len(density_list), interval=100)
#ani.save('wavefunction_animation.gif', writer='imagemagick')
ani.save('density_animation.gif', writer='pillow')


