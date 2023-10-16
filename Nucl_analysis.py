import matplotlib.pyplot as plt
import LabelLib as ll 
import MDAnalysis as mda
import nglview as nv
import numpy as np
import scipy as sp
from scipy import stats
from scipy.stats import binom
from tqdm.auto import tqdm, trange

def annot_max(x,y, ax=None):
    xmax = x[np.argmax(y)]
    ymax = y.max()
    text= "x={:.3f}".format(xmax)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",)
    kw = dict(xycoords='data',
              arrowprops=arrowprops, 
              bbox=bbox_props)
    ax.annotate(text, xy=(xmax, ymax), xytext=(xmax, ymax-0.5), **kw)

def gaus( mu, sigma,flac = True, start = 0.01, stop = 0.99, step = 0.01): 
    '''
    The function constructing the Gaussian distribution
    flac -- True to plot the function, False to return the x and y coordinate arrays. 
    '''
    if flac == True:
        x_data = np.arange(start, stop, step)
        y_data = stats.norm.pdf(x_data, mu, sigma)
        
        plt.plot(x_data, y_data)
        annot_max(x_data, y_data)
    else:
        x_data = np.arange(start, stop, step)[::-1]
        x_data = R0*((1-x_data)/x_data)**(1/6)
        mu = R0*((1-mu)/mu)**(1/6)
        sigma = R0*((1-sigma)/sigma)**(1/6)
        # print(x_data)## y-axis as the gaussian
        y_data = stats.norm.pdf(x_data, mu, sigma)
        ## plot data
        # plt.plot(x_data, y_data)

    return x_data, y_data
def get_vdv_radiuses(atoms_mda): #прнимает выборку обьект атомов из mda
    atomtypes = atoms_mda.types # 
    radii_dict=mda.topology.tables.vdwradii #
    radii=[radii_dict[atomtype] if atomtype in radii_dict else 0.0 for atomtype in atomtypes ]
    return radii

def savePqr_2(fileName, grid, weights):
    '''
    Function for saving pqr files
    '''
    points = grid.points()
    template = 'ATOM{0: 7}   AV  AV{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n'
    r = grid.discStep * 0.5
    weights/=np.max(np.unique(weights))
    with open(fileName, "w") as out:
        for i, (p, w_p) in enumerate(zip(points.T[:,:3], weights)):
            x, y, z = p
            w =  w_p
            sz = template.format(i + 1, 1, x, y, z, w, r)
            out.write(sz)

def setings_label(u, chain_Dye='L',attach_n_Dye=97):
    u.atoms.positions.T
    atoms=np.vstack([u.select_atoms(f'not ((segid {chain_Dye} and resnum {attach_n_Dye}))').atoms.positions.T, # or (segid {chain_Cy5} and resnum {attach_n_Cy5})
                     get_vdv_radiuses(u.select_atoms(f'not ((segid {chain_Dye} and resnum {attach_n_Dye}))').atoms)]) 

    anchor_atom_sel='((name C5 and resname DC DT) or (name C8 and resname DA DG) or (name CA))'# выбираем конкретный атом от которого считаем радиусы
    attachDy=u.select_atoms(f'segid {chain_Dye} and resnum {attach_n_Dye} and {anchor_atom_sel}')

    Dye_attachment_coordinat=attachDy[0].position 

    return atoms, Dye_attachment_coordinat 

def create_cloud_AV1(
    atoms,
    Dye_attachment_coordinat,
    linker_length_Dye=20,
    linker_width=1,
    dye_radius=7,
    simulation_grid_spacing=1,
    Emean_list=[],
    meanDist_list=[],
):
    av1 = ll.dyeDensityAV1(
        atoms,
        Dye_attachment_coordinat,
        linker_length_Dye,
        linker_width,
        dye_radius,
        simulation_grid_spacing,
    )
    if av1.points().size == 0:
        Emean = None
        meanDist = None
        Emean_list.append(Emean)
        meanDist_list.append(meanDist)
    else:
        pass
    return av1


def labellib_cloud_av1(
    u,
    chain_Dye="L",
    attach_n_Dye=15,
    linker_length_Dye=20,
    dye_radius=7,
    linker_width=1,
    simulation_grid_spacing=1,
    Emean_list=[],
    meanDist_list=[],
):
    atoms, Dye_attachment_coordinat = setings_label(
        u, chain_Dye=chain_Dye, attach_n_Dye=attach_n_Dye
    )

    av1 = create_cloud_AV1(
        atoms=atoms,
        Dye_attachment_coordinat=Dye_attachment_coordinat,
        linker_length_Dye=linker_length_Dye,
        linker_width=linker_width,
        dye_radius=dye_radius,
        simulation_grid_spacing=simulation_grid_spacing,
        Emean_list=Emean_list,
        meanDist_list=meanDist_list,
    )
    #     if av1.points().size == 0:
    #         av1 = 0

    return av1, Dye_attachment_coordinat


def labellib_length_calculating_av1(
    universe,
    chain_Dy="L",
    attach_n_Dy=24,
    linker_length_Dy=16,
    linker_width=1,
    dye_radius=7,
    simulation_grid_spacing=1,
    pqr1_name=None,
):
    """
    universe -- pass command with pdb file + dynamics file
    chain_Dy -- chain names 
    attach_n_Dy -- the number of nucleotide to which we attach the tag on the chain_Dy chain.
    linker_length_Dy -- linker length
    linker_width -- linker width
    dye_radius -- radii of our tags
    simulation_grid_spacing -- grid step
    """
    u = universe

    av1, Dy_attachment_point = labellib_cloud_av1(
        u,
        chain_Dye=chain_Dy,
        attach_n_Dye=attach_n_Dy,
        linker_length_Dye=linker_length_Dy,
        dye_radius=dye_radius,
        linker_width=linker_width,
        simulation_grid_spacing=simulation_grid_spacing,
    )
    if not pqr1_name is None:
        savePqr_2(
            f"{pqr1_name}" + ".pqr",
            grid=av1,
            weights=weights,
        )
    return av1