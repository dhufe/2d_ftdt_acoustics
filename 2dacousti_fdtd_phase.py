import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import seaborn as sns

def set_style():
    # This sets reasonable defaults for font size for
    # a figure that will go in a paper
    sns.set_context("paper")

    # Set the font to be serif, rather than sans
    sns.set(font='serif',style="ticks")

    # Make the background white, and specify the
    # specific font family
    sns.set_style("white", {
        "font.family": "serif",
        "font.serif": ["Times", "Palatino", "serif"]
    })
    sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

def set_size(fig):
    fig.set_size_inches(8, 5)
    fig.tight_layout()

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
NSources = 4
NFrames = 1000
# Density
rho = 1.241
# Frequency
freq = 250.0e3           # Hz
# wave speed
cmax = 343              # m/s
# single source width
SourceWidth  = 40
# single source height
SourceHeight = 10
# PML width
PMLWidth = 25
# maximum damping of the PML
sigma_max = 20
# colormap
colormap = 'RdYlBu'
# time between two frames
image_step = 200
# angle of incident 
alpha = 15

# wave length
lamda = cmax / freq
# Spatial resolution
# spatial stability critria: dx must be smaller or equal than lambda_min / 20
# where lambda_min is the shortest wavelength in the model!
dx = lamda/20
# Time domain resolution
# time stability criteria: dt must be smaller or equal than dx / ( sqrt(2) * c_max )
# where c_max is the highest wave speed in the model!
dt = dx/(cmax*np.sqrt(2) )

# independent phase shift of each source
phaseshift = (2 * np.pi / lamda ) * SourceWidth * dx * np.sin ( np.pi * alpha / 180)

# mesh grid
x = np.arange ( 0, (NX)*dx, dx )
y = np.arange ( 0, (NY)*dx, dx )
xx, yy = np.meshgrid( x, y )

## Properties of the fluid like density and viscosity
# Bulk viscosity: c^2 x rho
kappa = np.power(cmax, 2.0) * rho #142.0e3

## Computing magnitutes for two dimensional particle velocity and pressure
# particle velocity in x direction
Vx = np.zeros ( ( NY + 1, NX     ) )
# particle velocity in y direction
Vy = np.zeros ( ( NY    , NX + 1 ) )
# pressure
P  = np.zeros ( ( NY    , NX     ) )
# numeric absorption (PML) in x direction
sigma_x = np.ones ( ( NY, NX ) )
# numeric absorption (PML) in x direction
sigma_y = np.ones ( ( NY, NX ) )
# source configuration
Excitation = np.full(( NSources, NY, NX), False, dtype=bool)

### create the PML
for i in range ( 0, PMLWidth):
    sigma = 1 + sigma_max * (i/PMLWidth)
    sigma_x [: , PMLWidth - i    -1          ] = sigma 
    sigma_x [: , NX - PMLWidth + i           ] = sigma
    sigma_y [PMLWidth - i - 1, :             ] = sigma
    sigma_y [NY - PMLWidth + i ,:            ] = sigma


dxStep = 10
dyStep = NY // (NSources + 1)

offset = int ( NX//2 - ( dxStep * ( .5 + ( NSources//2 - 1 ) ) + SourceWidth * NSources//2 )  )

### create the acoustic sources
for iSource in range(0, NSources ):
    Pmy = 100

    if iSource == 0:
        Pmx = offset
    else: 
        Pmx = int( offset + iSource * (dxStep + SourceWidth ) )
    
    print ( 'Source (%d): Px %f, Py %f, width = %f, height = %f.' % ( iSource, Pmx*dx , Pmy*dx, SourceWidth*dx, SourceHeight*dx )  )
    SourceRect (  Excitation[iSource, :, :], Pmx, Pmy, SourceWidth, SourceHeight )    


sigma_x[ 0:Pmy + SourceHeight - 1, : ] = sigma_max
sigma_y[ 0:Pmy + SourceHeight - 1, : ] = sigma_max

# Plot creation
set_style()
fig = plt.figure(figsize=(12.8, 7.2))
ax  = fig.add_subplot(1,1,1)
cax  = ax.pcolormesh( xx, yy, P, vmin=-1, vmax=1, cmap=colormap, shading='auto')
ax.set_xlabel ( r'Position $x$ / $m$' )
ax.set_ylabel ( r'Position $y$ / $m$' )
ax.set_ylim   ( y[0], y[-1] )
ax.set_xlim   ( x[0], x[-1] )
fig.colorbar(cax)
fig.tight_layout()

## help variables for computation
dt_over_rho_x_dx = dt / ( rho * dx )
kappa_x_dt_over_dx = kappa * dt / dx

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
            Vx[i,j] -=  dt_over_rho_x_dx * (1/sigma_x[i,j]) * ( P[i,j] - P[i-1,j] )

    for i in range (1,NY):
        for j in range ( 2, NX):
            Vy[i,j] -= dt_over_rho_x_dx * (1/sigma_y[i,j]) * ( P[i,j] - P[i,j-1] )


    # Update sound pressure
    for i in range (1, NY):
        for j in range (1,NX):
            P[i,j] -=  ( ( Vx[i+1,j] - Vx[i,j] ) + ( Vy[i,j+1] - Vy[i,j] ) )

    # Acoustic source ( during n period)
    for iSource in range(0, NSources):
        ind = Excitation[iSource][:][:]
        if ( n * dt <= ( 2/freq + iSource * phaseshift)):
            P[ ind ] += (1.0-np.cos(2.0*np.pi*freq*n*dt + iSource * phaseshift))/2.0 * np.sin(2.0*np.pi*freq*n*dt + iSource * phaseshift )

    if (( n + 1 ) % 100 == 1) and ( n != 0):
        print ( '--- processing step %03d / %3.1f ms ---' % ( n , int((n*dt*1e3*10))/10.0  ))

    # Updating the current calculation step
    # n += 1

    # Updating data
    cax.set_array ( P.flatten() )
    #cax2.set_array ( np.abs(P).flatten() )
    ax.set_title("Time step {} ms".format( int((n*dt*1e3*10))/10.0 ) )
    return cax,

anim = animation.FuncAnimation(fig, updatefig, frames=NFrames-1, blit=True)

anim
Writer = animation.writers['ffmpeg']
writer = Writer ( fps=100, metadata=dict(artist='dhufschl' ), bitrate=6000)
anim.save('ftdt_acoustic_8_pml_2d.mp4', writer=writer )
plt.close()

