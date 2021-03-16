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

def circle(indices, ym, xm, r):

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

def SourceRect( indices, xs , ys, width, height ):
    indices[ ys : ys + height, xs : xs + width ] = True

## Computing grid and resolutions
# Grid size x-direction
NX = 500
# Grid size y-direction
NY = 1000
NSources = 8
NFrames = 15000

## Properties of the acoustic excitation
# Frequency
freq = 250.0e3           # Hz
# wave speed
cmax = 343              # m/s
# wave length
lamda = cmax / freq
# Spatial resolution
# spatial stability critria: dx must be smaller or equal than lambda_min / 20
# where lambda_min is the shortest wavelength in the model!
dx = lamda/10 #.01
# Time domain resolution
# time stability criteria: dt must be smaller or equal than dx / ( sqrt(2) * c_max )
# where c_max is the highest wave speed in the model!
dt = dx/(cmax*np.sqrt(2) ) # 20.0e-6
# mesh grid
x = np.arange ( 0, (NX)*dx, dx )
y = np.arange ( 0, (NY)*dx, dx )
xx, yy = np.meshgrid( x, y )

## Properties of the fluid like density and viscosity
# Density
rho = 1.241
# Bulk viscosity: c^2 x rho
kappa = np.power(cmax, 2.0) * rho #142.0e3

## Computing magnitutes for two dimensional particle velocity and pressure
Vx = np.zeros ( ( NY + 1, NX     ) )
Vy = np.zeros ( ( NY    , NX + 1 ) )
P  = np.zeros ( ( NY    , NX     ) )

sigma_x = np.ones ( ( NY, NX ) )
sigma_y = np.ones ( ( NY, NX ) )
Excitation = np.full(( NSources, NY, NX), False, dtype=bool)

SourceWidth  = 40
SourceHeight = 10

PMLWidth = 15

sigma_max = 5

for i in range ( 0, PMLWidth):
    sigma = np.exp ( - 2 * sigma_max * (i/PMLWidth) )

    sigma_x [: , PMLWidth - i    -1          ] = sigma 
    sigma_x [: , NX - PMLWidth + i           ] = sigma

    sigma_y [PMLWidth - i - 1, :             ] = sigma
    sigma_y [NY - PMLWidth + i ,:            ] = sigma

# sigma_y [0:PMLWidth, :  ] = 

# setup indices
ind = np.full(( NX, NY), False, dtype=bool)


dxStep = NX // (NSources + 1)
dyStep = NY // (NSources + 1)

for iSource in range(0, NSources ):
    Pmy = 100

    if iSource == 0:
        Pmx = int( dxStep // 2)
    else: 
        Pmx = int( ( iSource * dxStep + dxStep // 2 ) ) 

    print ( 'Source (%d): Px %f, Py %f, width = %f, height = %f.' % ( iSource, Pmx*dx , Pmy*dx, SourceWidth*dx, SourceHeight*dx )  )
    SourceRect (  Excitation[iSource][:][:], Pmx, Pmy, SourceWidth, SourceHeight )    


## Visual stuff
# Colormap
colormap = 'RdBu'
# Plot creation

fig = plt.figure(figsize=(12.8, 7.2))
ax  = fig.add_subplot(1,1,1)
cax  = ax.pcolormesh( xx, yy, P, vmin=-1, vmax=1, cmap=colormap, shading='auto')
ax.set_xlabel ( r'Position $x$ / $m$' )
ax.set_ylabel ( r'Position $y$ / $m$' )
ax.set_xlim   ( y[0], y[-1] )
ax.set_ylim   ( x[0], x[-1] )

fig.colorbar(cax)

fig.tight_layout()

image_step = 200

## help variables
dt_over_rho_x_dx = dt / ( rho * dx )
kappa_x_dt_over_dx = kappa * dt / dx
phaseshift = np.pi*.125
n = 0

print ( 'Spatial stepsize  %00.3f mm.' % ( dx*1e3 ) )
print ( 'Time stepsize     %00.3f us.' % ( dt*1e6 ) )
print ( 'Volume elasticity %03.3f   .' % ( kappa*1e-3 ) )
print ( 'Frames            %0003d   .' % ( NFrames ))
print ( 'Pulse width         %03d us.' % ( int((1.0/freq)/dt) ) )
print ( 'Phaseshift         %3.3f us.' % ( phaseshift *1e6/ (2*np.pi*freq )  ) )

def updatefig ( n ):
    # Updating particle velocities
    for i in range (2,NY):
        for j in range ( 1, NX ):
            Vx[i,j] -=  dt_over_rho_x_dx * sigma_x[i,j] * ( P[i,j] - P[i-1,j] )

    for i in range (1,NY):
        for j in range ( 2, NX):
            Vy[i,j] -= dt_over_rho_x_dx * sigma_y[i,j] * ( P[i,j] - P[i,j-1] )


    # Update sound pressure
    for i in range (1, NY):
        for j in range (1,NX):
            P[i,j] -=  ( ( Vx[i+1,j] - Vx[i,j] ) + ( Vy[i,j+1] - Vy[i,j] ) )

    # Acoustic source ( during one period)
    for iSource in range(0, NSources):
        ind = Excitation[iSource][:][:]
        if ( n * dt >= ( 5/freq + iSource * phaseshift)):
            P[ ind ] += (1.0-np.cos(2.0*np.pi*freq*n*dt + iSource * phaseshift))/2.0 * np.sin(2.0*np.pi*freq*n*dt + iSource * phaseshift )

    if (( n + 1 ) % 100 == 1) and ( n != 0):
        print ( '--- processing step %03d / %3.1f ms ---' % ( n , int((n*dt*1e3*10))/10.0  ))

    # Updating the current calculation step
    # n += 1

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
writer = Writer ( fps=100, metadata=dict(artist='dhufschl' ), bitrate=6000)
anim.save('ftdt_acoustic_4_sources_phase_2d.mp4', writer=writer )
plt.close()

