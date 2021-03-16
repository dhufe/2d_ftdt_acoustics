#
#
#
#
#

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def circle(indices, xm, ym, r):

    x = r - 1
    y = 0
    dx = 1
    dy = 1
    err = dx - (r * 2)


    while ( x >= y):
        indices [ xm + x , ym + y] = True
        indices [ xm + y , ym + x] = True
        indices [ xm - y , ym + x] = True
        indices [ xm - x , ym + y] = True
        indices [ xm - x , ym - y] = True
        indices [ xm - y , ym - x] = True
        indices [ xm + y , ym - x] = True
        indices [ xm + x , ym - y] = True

        if err <= 0:
            y+=1
            err += dy
            dy += 2

        if err > 0:
            x-=1
            dy += 2
            err += dx - ( 2*r )


## Computing grid and resolutions
# Grid size x-direction
NX = 500
# Grid size y-direction
NY = 500

## Properties of the acoustic excitation
# Frequency
freq = 60.0e3           # Hz
# wave speed
cmax = 343              # m/s
# wave length
lamda = cmax / freq
# Spatial resolution
# spatial stability critria: dx must be smaller or equal than lambda_min / 20
# where lambda_min is the shortest wavelength in the model!
dx = lamda/40 #.01
# Time domain resolution
# time stability criteria: dt must be smaller or equal than dx / ( sqrt(2) * c_max )
# where c_max is the highest wave speed in the model!
dt = dx/(cmax*np.sqrt(2) ) # 20.0e-6
# mesh grid
x = np.arange ( 0, (NX)*dx, dx )
y = np.arange ( 0, (NY)*dx, dx )
xx, yy = np.meshgrid( y, x )

## Properties of the fluid like density and viscosity
# Density
rho = 1.241
# Bulk viscosity: c^2 x rho
kappa = np.power(cmax, 2.0) * rho #142.0e3

## Computing magnitutes for two dimensional particle velocity and pressure
Vx = np.zeros ( ( NX + 1, NY     ) )
Vy = np.zeros ( ( NX    , NY + 1 ) )
P  = np.zeros ( ( NX    , NY     ) )

NFrames = 3000

SourceWidth  = 20
SourceHeight = 10

# setup indices
ind = np.full(( NX, NY), False, dtype=bool)

# setup indices
ind = np.full(( NX, NY), False, dtype=bool)

Pmx = int(np.floor(NX/4 - 10 + 1 ))
Pmy = int(np.floor(NY/2 + 1 ))

circle ( ind, Pmx, Pmy, 20 )

Pmx = int(np.floor(NX/2 - 10 + 1 ))
Pmy = int(np.floor(NY/2 + 1 ))

circle ( ind, Pmx, Pmy, 20 )

Pmx = int(np.floor(3*NX/4 - 10  + 1 ))
Pmy = int(np.floor(  NY/2 + 1 ))

circle ( ind, Pmx, Pmy, 20 )

Pmx = int(np.floor(NX - 30  + 1 ))
Pmy = int(np.floor(NY/2 + 1 ))

circle ( ind, Pmx, Pmy, 20 )

# Pxs = int(np.floor(NX/2 - SourceWidth/2  + 1 ))
# Pxe = int(np.floor(NX/2 + SourceWidth/2  + 1 ))

# Pys = int(np.floor(NY/2 - SourceHeight/2 + 1 ))
# Pye = int(np.floor(NY/2 + SourceHeight/2 + 1 ))

# ind[ Pxs:Pxe, Pys:Pys ] = True

## Visual stuff
# Colormap
colormap = 'RdBu'
# Plot creation


fig = plt.figure(figsize=(12.8, 7.2))
ax  = fig.add_subplot(1,1,1)
cax  = ax.pcolormesh( xx, yy, P, vmin=-1, vmax=1, cmap=colormap)
ax.set_xlabel ( r'Position $x$ / $m$' )
ax.set_ylabel ( r'Position $y$ / $m$' )
ax.set_xlim   ( y[0], y[-1] )
ax.set_ylim   ( x[0], x[-1] )

fig.colorbar(cax)

fig.tight_layout()

image_step = 400

## help variables
dt_over_rho_x_dx = dt / ( rho * dx )
kappa_x_dt_over_dx = kappa * dt / dx

print ( 'Spatial stepsize  %00.3f mm.' % ( dx*1e3 ) )
print ( 'Time stepsize     %00.3f us.' % ( dt*1e6 ) )
print ( 'Volume elasticity %03.3f .' % ( kappa*1e-3 ) )
print ( 'Pulse width         %d    ' % ( int((1.0/freq)/dt) ) )

n = 0

def updatefig (n):
    # Updating particle velocities
    for i in range (2,NX):
        for j in range ( 1, NY ):
            Vx[i,j] -=  dt_over_rho_x_dx * ( P[i,j] - P[i-1,j] )

    for i in range (1,NX):
        for j in range ( 2, NY):
            Vy[i,j] -= dt_over_rho_x_dx * ( P[i,j] - P[i,j-1] )


    # Update sound pressure
    for i in range (1, NX):
        for j in range (1,NY):
            P[i,j] -=  ( ( Vx[i+1,j] - Vx[i,j] ) + ( Vy[i,j+1] - Vy[i,j] ) )

    # Acoustic source ( during one period)
    if (n < (1.0/freq)/dt):
        P[ ind ] += .33*(1.0-np.cos(2.0*np.pi*freq*n*dt))/2.0 * np.sin(2.0*np.pi*freq*n*dt)

    if (n < (1.0/(.5*freq))/dt):
        P[ ind ] += .33*(1.0-np.cos(2.0*np.pi*.5*freq*n*dt))/2.0 * np.sin(2.0*np.pi*.5*freq*n*dt)

    if (n < (1.0/(.25*freq))/dt):
        P[ ind ] += .33*(1.0-np.cos(2.0*np.pi*.25*freq*n*dt))/2.0 * np.sin(2.0*np.pi*.25*freq*n*dt)


    if (( n + 1 ) % 100 == 1) or ( n == 0):
        print ( '--- processing step %03d ---' % ( n ))

    # Updating the current calculation step
    n += 1

    # Updating data
    cax.set_array ( P.flatten() )
    #cax2.set_array ( np.abs(P).flatten() )
    #ax.set_title("circular membrane: l={}, m={}-Mode".format(l+1,m))
    ax.set_title("Time step {} ms".format( int((n*dt*1e3*10))/10.0 ) )
    #ax[1].set_title("Time step {} ms".format( int((n*dt*1e3*10))/10.0 ) )
    return cax,

anim = animation.FuncAnimation(fig, updatefig, frames=NFrames-1,interval=image_step, blit=True)

anim
Writer = animation.writers['ffmpeg']
writer = Writer ( fps=100, metadata=dict(artist='dkotscha' ), bitrate=6000)
anim.save('ftdt_acoustic_one_source_2d.mp4', writer=writer )
plt.close()

