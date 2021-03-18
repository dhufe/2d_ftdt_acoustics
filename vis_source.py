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
freq = 250.0e3           # Hz
# wave speed
cmax = 343              # m/s
# wave length
lamda = cmax / freq
# Spatial resolution
# spatial stability critria: dx must be smaller or equal than lambda_min / 20
# where lambda_min is the shortest wavelength in the model!
dx = lamda/20 #.01

def SourceRect( indices, xs , ys, width, height ):
    indices[ ys : ys + height, xs : xs + width ] = True

x = np.arange ( 0, (NX)*dx, dx )
y = np.arange ( 0, (NY)*dx, dx )
ind = np.full(( NX, NY), False, dtype=bool)

NSources = 4
SourceWidth  = 25
SourceHeight = 40
dxStep = 10
dyStep = NY // (NSources + 1)

Wc = NSources * ( SourceWidth + dxStep ) - dxStep 

print ( Wc ) 
print ( dxStep ) 

offset = int ( NX//2 - ( dxStep * ( .5 + ( NSources // 2 - 1 ) ) + SourceWidth * NSources//2 )  )

print ( offset ) 
for iSource in range(0, NSources ):
    Pmy = 100

    if iSource == 0:
        Pmx = offset
    else: 
        Pmx = int( offset + iSource * (dxStep + SourceWidth ) )
    
    print ( 'Source (%d): Px %f, Py %f, width = %f, height = %f.' % ( iSource, Pmx*dx , Pmy*dx, SourceWidth*dx, SourceHeight*dx )  )
    SourceRect (  ind, Pmx, Pmy, SourceWidth, SourceHeight )    

xx, yy = np.meshgrid( y, x )

fig, ax = plt.subplots(1)
ax.pcolormesh ( xx, yy, ind, shading = 'auto' )

plt.savefig ( 'source.png' , dpi = 300)
plt.close()
