__author__ = "Navid Rezazadeh"
__copyright__ = "Copyright (C) 2022 Navid Rezazadeh"
__license__ = "Public Domain"
__version__ = "1.0"

# A function to plot a radiation pattern in a 3d polar plot

############################################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib import colors
from matplotlib import cm

def plot_3d_polar(F:np.ndarray, theta:np.ndarray=np.arange(0, 180 , 1), 
                    phi:np.ndarray=np.arange(0, 360 , 1),
                    levels:np.ndarray=np.arange(-30,30,5),  
                    view_ang:list=[40, 30]):
    """ 3d polar radiation pattern where
    The color map is proportional to R rather than Z.


    Inputs:
    F(theta, phi) is the far field matrix over the entire angular space of the spherical coordinate system,
    i.e. 0<theta<180, 0<phi<360. F matrix rows correpsond to theta and columns correspond to phi,
    it has shape (m, n), where m is the length of theta and n is the length of phi 
    
    levels sets the range to display, its min is set as coordinate origin and the max is the maximum radius shown

    view_ang is [el, az] in deg
    """
    
    cmap = cm.afmhot # colormap
    
    # Setting a minimum value for the pattern that corresponds to radius = 0
    F_min = levels[0]
    threshold_indices = F[:, :] < F_min 
    R = F
    R[threshold_indices] = F_min # Force this as the minimum of the pattern
    R = R - F_min
    # Setting a max value for radius of plot
    R_max = levels[-1]-levels[0] 

    fig = plt.figure(figsize=(5.5, 5), dpi=200, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1, projection='3d', frame_on = False)
    ax.grid(False)
    ax.set_xlim(-(abs(R_max)+1), abs(R_max)+1)
    ax.set_ylim(-(abs(R_max)+1), abs(R_max)+1)
    ax.set_zlim(-5, abs(R_max)+1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_axis_off()
    
    # Add coordinate axes
    marg = 5
    x = np.linspace(0, 1, 40) 
    y = np.zeros(len(x))
    X, Y  = np.meshgrid(x, y)
    X = (R_max + marg) * X
    Y = (R_max + marg) * Y
    Z = 0 * X
    ax.plot_wireframe(X, Y, Z, linewidth=0.5, rstride=3, cstride=3,
                        color = 'k', alpha=1)
    
    y = np.linspace(0, 1, 40)
    x = np.zeros(len(y))
    X, Y  = np.meshgrid(x, y)
    X = (R_max + marg) * X
    Y = (R_max + marg) * Y
    Z = 0 * X
    ax.plot_wireframe(X, Y, Z, linewidth=0.5, rstride=3, cstride=3,
                        color = 'k', alpha=1)
    
    x = 0*np.linspace(0, 1, 40)
    y = 0*np.linspace(0, 1, 40)
    X, Y  = np.meshgrid(x, y)
    z = (R_max + marg*2) * np.linspace(0, 1, 40)
    z = np.expand_dims(z, axis=0)
    Z = np.repeat(z, len(x), axis=0)
    ax.plot_wireframe(X, Y, Z, linewidth=0.5, rstride=3, cstride=3,
                        color = 'k', alpha=1)

    ax.text(R_max+8, 0, 0, 'x', 'x', color = 'k')
    ax.text(0, R_max+5, 0, 'y', 'y', color = 'k')
    ax.text(-10, -10, R_max+5, 'z', 'x', color = 'k')

    # Spherical to Cartesian transformation
    theta_mesh, phi_mesh = np.meshgrid(theta, phi, indexing='ij')
    THETA = np.deg2rad(theta_mesh)
    PHI = np.deg2rad(phi_mesh)
    X = R * np.sin(THETA) * np.cos(PHI)
    Y = R * np.sin(THETA) * np.sin(PHI)
    Z = R * np.cos(THETA)

    rgb = np.ones((Z.shape[0], Z.shape[1], 4))
    color_amp = R/R_max # Color is the normalized amplitude of the radial component
    rgb = cm.afmhot(color_amp)
    light = LightSource(azdeg=45, altdeg=45)
    rgb = light.shade_rgb(rgb, Z, blend_mode='overlay', fraction=1)
    
    ax.plot_surface(X, Y, Z,
                    rstride = 1, cstride = 1, facecolors=rgb,
                    antialiased=False, linewidth=0)
    
    # Add a color bar which maps values to colors
    norm = colors.Normalize(vmin=F_min, vmax=R_max+F_min)
    cbax = fig.add_axes([0.7, 0.05, 0.15, 0.8], frame_on=False)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=cbax, shrink=1, aspect=15)
    cbax.set_xticks([])
    cbax.set_yticks([])
    cb.set_label('Directivity (dBi)')

    # Adjusting the view
    ax.view_init(view_ang[0], view_ang[1])
    ax.set_xlim(-(abs(R_max)+1), abs(R_max)+1)
    ax.set_ylim(-(abs(R_max)+1), abs(R_max)+1)
        
    return fig

if __name__ == '__main__':
    
    theta = np.arange(0, 180, 5)
    phi = np.arange(0, 361, 1)
    theta_mesh, phi_mesh = np.meshgrid(theta, phi, indexing='ij')
    F = abs(20*np.cos(np.pi*theta_mesh/180))
    plot_3d_polar(F, theta, phi, levels=np.arange(5, 25), view_ang=[40, 30])
    plt.show()