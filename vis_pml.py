import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

## Computing grid and resolutions
# Grid size x-direction
NX = 500
# Grid size y-direction
NY = 500

## Properties of the acoustic excitation
# Frequency
freq = 50.0e3           # Hz
# wave speed
cmax = 343              # m/s
# wave length
lamda = cmax / freq
# Spatial resolution
# spatial stability critria: dx must be smaller or equal than lambda_min / 20
# where lambda_min is the shortest wavelength in the model!
dx = lamda/10 #.01

def SourceRect( indices, xs , ys, width, height ):
    indices[ ys : ys + height, xs : xs + width ] = True

x = np.arange ( 0, (NX)*dx, dx )
y = np.arange ( 0, (NY)*dx, dx )
ind = np.full(( NX, NY), False, dtype=bool)

sigma_x = np.ones ( ( NX , NY   ) )
sigma_y = np.ones ( ( NX   , NY ) )

SourceWidth  = 40
SourceHeight = 10

PMLWidth = 20

sigma_max = 5

for i in range ( 0, PMLWidth):
    sigma = np.exp ( - 2 * sigma_max * (i/PMLWidth) )

    sigma_x [: , PMLWidth - i    -1          ] = sigma 
    sigma_x [: , NX - PMLWidth + i           ] = sigma

    sigma_y [PMLWidth - i - 1, :             ] = sigma
    sigma_y [NY - PMLWidth + i ,:            ] = sigma

xx, yy = np.meshgrid( y, x )

fig, ax = plt.subplots(1)
#ax.plot ( sigma_x[ NY//2, : ] )
ax.pcolormesh ( sigma_x + sigma_y , shading = 'auto' )


plt.savefig ( 'pml.png' , dpi = 300)
plt.close()