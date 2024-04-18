import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import sys
import time

def t(tmin, tmax, dt):
    ind = int((tmax-tmin)/dt) + int(1)
    time_ = np.zeros(ind)
    time_[0] = tmin
    for j1 in range(1, ind):
        time_[j1] = time_[j1-1] + dt
    #print'time: ', time, len(time)
    return time_

def Ga(n, eps, Delta, J, a, b, hbar):
    if (n == 0):
        v = (-1j/hbar)*(eps*a[n] +Delta*b[n])*(-0.500)
    if (n>0):
        v = (-1j/hbar)*( (eps*a[n] +Delta*b[n])*(-0.500) +J*a[n]*(np.conj(b[n-1])*b[n-1]) )
    return v

def Gb(n, eps, Delta, J, a, b, hbar):
    Nw = len(a)
    if (n<(Nw-1)):
        v = (-1j/hbar)*( (-eps*b[n] +Delta*a[n])*(-0.500) +J*b[n]*(np.conj(a[n+1])*a[n+1]) )
    if (n==(Nw-1)):
        v = (-1j/hbar)*( (-eps*b[n] +Delta*a[n])*(-0.500) )
    return v
    
sigma_x = np.matrix([[0.0, 1.0], [1.0, 0]])
sigma_y = np.matrix([[0.0, -1.0*1j], [1.0*1j, 0]])
sigma_z = np.matrix([[1.0, 0.0], [0.0, -1.0]])

def BS_coor(state):
    #sigma_x = np.matrix([[0.0, 1.0], [1.0, 0]])
    #sigma_y = np.matrix([[0.0, 1.0*1j], [1.0*1j, 0]])
    #sigma_z = np.matrix([[1.0, 0.0], [0.0, -1.0]])
    
    st = np.matrix(state)
    st_con = st.getH() 
    rho = st_con*st
    xi = np.trace(rho*sigma_x)
    yi = np.trace(rho*sigma_y)
    zi = np.trace(rho*sigma_z)
    return xi.real, yi.real, zi.real


#start main
start_time = time.time()
tmin = 0.0
tmax = 2000.0
dt = 0.0001*100.0
t = t(tmin, tmax, dt)
num_t = len(t)

hbar = 1.0 
N = 12 #number of doubel well
a = np.zeros(N, dtype = complex)
b = np.zeros(N, dtype = complex)

eps = 0.0 #0.1
Delta = 0.05 #0.05
J = 0.01 #0.05

#initial condition
#1
#######################
a[0] = 1.0 #np.cos(np.pi/8)#1.0
b[0] = 0.0 #np.sin(np.pi/8)#0.0
for n in range(1, N):
    a[n] = np.sqrt(1.0/2.0)
    b[n] = np.sqrt(1.0/2.0)
#######################

a_t = np.zeros([N, num_t], dtype = complex)
b_t = np.zeros([N, num_t], dtype = complex)

a_t[:,0] = a
b_t[:,0] = b

for l in range(1, num_t):
    tmp_a = a_t[:,l-1]
    tmp_b = b_t[:,l-1]

    for n in range(0, N):
        ### k update is incorrect (22/11/2021) ---> not edit now
        ### note on 13/03/2023 it can be used,
        ### bc it is evaluated at n^th component,
        ### input tmp_a and tmp_b as vector are trivial.  
        ka1 = Ga(n, eps, Delta, J, tmp_a, tmp_b, hbar)
        ka2 = Ga(n, eps, Delta, J, (tmp_a + 0.50*dt*ka1), (tmp_b + 0.50*dt*ka1), hbar)
        ka3 = Ga(n, eps, Delta, J, (tmp_a + 0.50*dt*ka2), (tmp_b + 0.50*dt*ka2), hbar)
        ka4 = Ga(n, eps, Delta, J, (tmp_a + dt*ka3), (tmp_b + dt*ka3), hbar)
        a_t[n, l] = tmp_a[n] +dt*(ka1 + 2.0*ka2 +2.0*ka3 +ka4)/6.0
        
        kb1 = Gb(n, eps, Delta, J, tmp_a, tmp_b, hbar) 
        kb2 = Gb(n, eps, Delta, J, (tmp_a + 0.50*dt*kb1), (tmp_b + 0.50*dt*kb1), hbar)
        kb3 = Gb(n, eps, Delta, J, (tmp_a + 0.50*dt*kb2), (tmp_b + 0.50*dt*kb2), hbar)
        kb4 = Gb(n, eps, Delta, J, (tmp_a + dt*kb3), (tmp_b + dt*kb3), hbar)
        b_t[n, l] = tmp_b[n] +dt*(kb1 + 2.0*kb2 +2.0*kb3 +kb4)/6.0

        
#Pa_t = np.zeros([N, num_t])
#Pb_t = np.zeros([N, num_t])

#for l1 in range(0, num_t):
#    for l2 in range(0, N):
#        Pa_t[l2, l1] =np.sqrt(1.0/2.0) (np.conj(a_t[l2,l1])*a_t[l2,l1]).real
#        Pb_t[l2, l1] = (np.conj(b_t[l2,l1])*b_t[l2,l1]).real

#xval_a = np.array(range(0,N))*2 
#xval_b = np.array(range(0,N))*2 +1
#for l in range(0, num_t, 50):
#    plt.figure(1)
#    #plt.bar(xval_a, Pa_t[:,l])
#    #plt.bar(xval_b, Pb_t[:,l])
#    plt.plot(Pa_t[:,l], '*-')
#    plt.plot(Pb_t[:,l], 'o-')
#    plt.ylim([0,1.2])
#    plt.show(block=False)
#    plt.pause(0.01)
#    plt.close(1)

#for l in range(0, N):
#    plt.plot(t,Pa_t[l,:]+l)
#    plt.plot(t,Pb_t[l,:]+l)
#plt.show()


sk = 100*2
num_sk = math.floor(num_t/sk)
a_t_sk = np.zeros([N, num_sk], dtype = complex)
b_t_sk = np.zeros([N, num_sk], dtype = complex)
t_sk = np.zeros(num_sk)

indx_sk =0
for l in range(0, num_sk):
    a_t_sk[:, l] = a_t[:, indx_sk]
    b_t_sk[:, l] = b_t[:, indx_sk]
    t_sk[l] = t[indx_sk]
    indx_sk = indx_sk + sk

Pa_t_sk = np.zeros([N, num_sk])
Pb_t_sk = np.zeros([N, num_sk])
x = np.zeros([N, num_sk])
y = np.zeros([N, num_sk])
z = np.zeros([N, num_sk])



for l1 in range(0, num_sk):
    for l2 in range(0, N):
        Pa_t_sk[l2, l1] = (np.conj(a_t_sk[l2,l1])*a_t_sk[l2,l1]).real
        Pb_t_sk[l2, l1] = (np.conj(b_t_sk[l2,l1])*b_t_sk[l2,l1]).real
        xi, yi, zi =BS_coor([a_t_sk[l2,l1], b_t_sk[l2,l1]])
        x[l2,l1] = xi
        y[l2,l1] = yi
        z[l2,l1] = zi

    
for l in range(0, N):
    plt.plot(t_sk,Pa_t_sk[l,:]+l, 'x')
    plt.plot(t_sk,Pb_t_sk[l,:]+l, 'o')
plt.show()

plt.plot(t_sk, x[0,:], '.r')
plt.plot(t_sk, y[0,:], '.g')
plt.plot(t_sk, z[0,:], '.b')
plt.show()


print('before skip: ', num_t)
print('frame: ', num_sk)
##############################

# Set the desired video resolution
video_width = 1280
video_height = 720

# Calculate the dpi based on the desired video resolution
dpi = 96  # Assuming a standard DPI for on-screen viewing

# Calculate the figsize based on the desired video resolution and dpi
fig_width = video_width / dpi
fig_height = video_height / dpi

# Create a 3D figure and axes                                                                
fig = plt.figure(figsize=(fig_width, fig_height))
ax = fig.add_subplot(111, projection='3d')
#ax = plt.axes(projection='3d')

# Set the limits and labels for the axes                                                     
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

# draw sphere
u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
xsp = np.cos(u)*np.sin(v)
ysp = np.sin(u)*np.sin(v)
zsp = np.cos(v)
xyzl_1 = np.array([-1.1, 1.1])
xyzl_0 = np.array([0, 0])

def animate(frm):
    ax.cla()  # Clear the previous frame
    ax.set_box_aspect((1,6,1))

    # Plot the sphere background
    ax.set_axis_off()
    #ax.plot_surface(xsp, ysp, zsp, color='lightblue', alpha=0.6, linewidth=0)
    
    sfd = 3
    for n in range(0, N):
        now_sf = sfd*n
        ax.plot_wireframe(xsp, ysp+now_sf, zsp, color="gray", alpha=0.1)
        ax.plot3D(xyzl_1, xyzl_0+now_sf, xyzl_0, 'r')
        ax.plot3D(xyzl_0, xyzl_1+now_sf, xyzl_0, 'g')
        ax.plot3D(xyzl_0, xyzl_0+now_sf, xyzl_1, 'b')
        ax.plot3D(x[n,:frm], y[n,:frm]+now_sf, z[n,:frm], '-', color='k')
        ax.plot3D(x[n,frm], y[n,frm]+now_sf, z[n,frm], 'o')
        ##ax.plot3D([0, x[n,frm]], [now_sf, y[n,frm]+now_sf], [0, z[n,frm]], '-m')

    # Set the limits and labels for the axes
    #ax.set_xlim([-1, 1])
    #ax.set_ylim([-1, 1])
    #ax.set_zlim([-1, 1])
    #ax.set_xlabel('X')
    #ax.set_ylabel('Y')
    #ax.set_zlabel('Z')

    ax.view_init(elev=10.0, azim=0.0)
    ##ax.view_init(elev=30, azim=30)

# Set the number of frames and the interval between frames                                   
num_frames = num_sk
interval = 50  # milliseconds 

# Create the animation                                                                       
animation = FuncAnimation(fig, animate, frames=num_frames, interval=interval)

# Set the filename for the output
filename = '00_test_v2.mp4'

#filename = '03_edit02_output_mp4/05_initial_state4/05_1_J_lessthan_Delta.mp4'
#filename = '03_edit02_output_mp4/05_initial_state4/05_2_J_equal_Delta.mp4'
#filename = '03_edit02_output_mp4/05_initial_state4/05_3_J_morethan_Delta.mp4'

#filename = '03_edit02_output_mp4/04_1qubit/04_6_sigmaxz_init_c.mp4'

# Save the animation as an MP4 file                                                          
animation.save(filename, writer='ffmpeg', dpi=dpi)

print('total time', time.time()-start_time)
plt.show()

