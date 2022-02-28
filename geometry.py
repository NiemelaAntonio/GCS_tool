import numpy as np
from astropy.coordinates import SkyCoord
from numpy import pi, sin, cos, tan, arcsin, sqrt
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
from sunpy.coordinates import frames, sun


def skeleton(alpha, distjunc, straight_vertices, front_vertices, k):
    """
    Compute the axis of the GCS CME model. Based on IDL version shellskeleton.pro
    https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/shellskeleton.pro

    Parameters
    ----------
    alpha: width of CME (half angle)
    distjunc: CME height (in length units, e.g. solar radii)
    straight_vertices: number of vertices along straight edges
    front_vertices: number of vertices along front
    k: GCS ratio

    Returns
    -------

    """
    # compute entire loop lenght
    loop_length = distjunc * (1. + (alpha +pi / 2) * tan(alpha))
    # copute circular part half length
    circle_length = distjunc * np.tan(alpha) * (2. * alpha +pi) / 2

    gamma = arcsin(k)

    # calculate the points of the straight line part
    pRstart = np.array([0, sin(alpha), cos(alpha)])  # start on the limb
    pLstart = np.array([0, -sin(alpha), cos(alpha)])
    pslR = np.outer(np.linspace(0, distjunc, straight_vertices), np.array([0, sin(alpha), cos(alpha)]))
    pslL = np.outer(np.linspace(0, distjunc, straight_vertices), np.array([0, -sin(alpha), cos(alpha)]))
    rsl = tan(gamma) * norm(pslR, axis=1)
    casl = np.full(straight_vertices, -alpha)

    # calculate the points of the circular part
    beta = np.linspace(-alpha, pi / 2, front_vertices)
    hf = distjunc
    h = hf / cos(alpha)
    rho = hf * tan(alpha)

    X0 = (rho + h * k ** 2 * sin(beta)) / (1 - k ** 2)
    rc = sqrt((h ** 2 * k ** 2 - rho ** 2) / (1 - k ** 2) + X0 ** 2)
    cac = beta

    pcR = np.array([np.zeros(beta.shape), X0 * cos(beta), h + X0 * sin(beta)]).T
    pcL = np.array([np.zeros(beta.shape), -X0 * cos(beta), h + X0 * sin(beta)]).T

    r = np.concatenate((rsl, rc[1:], np.flipud(rc)[1:], np.flipud(rsl)[1:]))
    ca = np.concatenate((casl, cac[1:], pi-np.flipud(cac)[1:], pi-np.flipud(casl)[1:]))
    p = np.concatenate((pslR, pcR[1:], np.flipud(pcL)[1:], np.flipud(pslL)[1:]))

    return p, r, ca


def gcs_mesh(alpha, height, straight_vertices, front_vertices, circle_vertices, k):
    """
    Calculate GCS model mesh. Based on IDL version cmecloud.pro.
    https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/idl/scraytrace/cmecloud.pro

    Parameters
    ----------
    alpha: width of CME (half angle, in radians)
    height: CME height (in length units, e.g. solar radii)
    straight_vertices: number of vertices along straight edges
    front_vertices: number of vertices along front
    circle_vertices: number of vertices along each circle
    k: GCS ratio

    Returns
    -------orking on python 3.7.3. Make sure it's available on 2.7
    GCS mesh (in length units, e.g. solar radii) in the following coordinate system (Heliographic Stonyhurst):
    - Z axis: Sun-Earth line projected onto Sun's equatorial plane
    - Y axis: Sun's north pole
    - Origin: center of the Sun

    """
    # calculate position

    # calculate height of skeleton depending on height of leading edge
    distjunc = height *(1-k)*cos(alpha)/(1.+sin(alpha))
    p, r, ca = skeleton(alpha, distjunc, straight_vertices, front_vertices, k)

    theta = np.linspace(0, 2 * pi, circle_vertices)
    pspace = np.arange(0, (front_vertices + straight_vertices) * 2 - 3)
    u, v = np.meshgrid(theta, pspace)
    u, v = u.flatten(), v.flatten()

    mesh = r[v, np.newaxis] * np.array([cos(u), sin(u) * cos(ca[v]), sin(u) * sin(ca[v])]).T + p[v]

    # mesh = np.concatenate([
    #    r * np.array([costheta, sintheta * cos(ca), sintheta * sin(ca)]).T + np.array([x, y, 0])
    #    for (x, y), r, ca in zip(p, r, ca)
    # ])

    return mesh, u, v


def rotate_mesh(mesh, neang):
    """
    Rotates the GCS mesh into the correct orientation

    Parameters
    ----------
    mesh GCS mesh, returned by gcs_mesh
    neang angles: longitude, latitude, and rotation angles (in radians)

    Returns
    -------
    rotated GCS mesh

    """
    return Rotation.from_euler('zyx', [neang[2], neang[1], neang[0]]).apply(mesh)


def gcs_mesh_rotated(alpha, height, straight_vertices, front_vertices, circle_vertices, k, lat, lon, tilt):
    """
    Calculate GCS model mesh and rotate it into the correct orientation.
    Convenience function that combines gcs_mesh and rotate_mesh.

    Parameters
    ----------
    alpha: width of CME (half angle, in radians)
    height: CME height (in length units, e.g. solar radii)
    straight_vertices: number of vertices along straight edges
    front_vertices: number of vertices along front
    circle_vertices: number of vertices along each circle
    k: GCS ratio
    lat: CME latitude in radians
    lon: CME longitude in radians (0 = Earth-directed)
    tilt: CME tilt angle in radians

    Returns
    -------
    rotated GCS mesh
    """
    mesh, u, v = gcs_mesh(alpha, height, straight_vertices, front_vertices, circle_vertices, k)
    mesh = rotate_mesh(mesh, [lon, lat, tilt])
    return mesh, u, v


def gcs_mesh_sunpy(date, alpha, height, straight_vertices, front_vertices, circle_vertices, k, lat, lon, tilt):
    """
    Provides the GCS model mesh in SunPy SkyCoord format. This can be directly plotted using SunPy.

    Parameters
    ----------
    date: date and time of observation (Python datetime instance)
    alpha: width of CME (half angle, in radians)
    height: CME height (in length units, e.g. solar radii)
    straight_vertices: number of vertices along straight edges
    front_vertices: number of vertices along front
    circle_vertices: number of vertices along each circle
    k: GCS ratio
    lat: CME latitude in radians
    lon: CME longitude in radians (0 = Earth-directed)
    tilt: CME tilt angle in radians

    Returns
    -------
    GCS mesh as SunPy SkyCoord
    """
    mesh, u, v = gcs_mesh_rotated(alpha, height, straight_vertices, front_vertices, circle_vertices, k, lat, lon, tilt)
    m = mesh.T[[2, 1, 0], :] * sun.constants.radius
    m[1, :] *= -1
    mesh_coord = SkyCoord(*(m), frame=frames.HeliographicStonyhurst,
                          obstime=date, representation_type='cartesian')
    return mesh_coord


def apex_radius(alpha, height, k):
    """
    Calculates the cross-section radius of the flux rope at the apex, based on GCS parameters alpha, height and kappa.

    Parameters
    ----------
    alpha: width of CME (half angle, in radians)
    height: CME height (in length units, e.g. solar radii)
    k: GCS ratio

    Returns
    -------
    apex radius (in length units, e.g. solar radii)
    """
    h = height * (1 - k) * cos(alpha) / (1. + sin(alpha))
    b = h / cos(alpha)
    rho = h * tan(alpha)
    return k * (b + rho) / (1 - k ** 2)




"""
This module provides dictionaries for generating
`~matplotlib.colors.LinearSegmentedColormap`, and a dictionary of these
dictionaries.
"""
import pathlib

import matplotlib.colors as colors
import numpy as np

import astropy.units as u

__all__ = [
    'aia_color_table', 'sswidl_lasco_color_table', 'eit_color_table',
    'sxt_color_table', 'xrt_color_table', 'trace_color_table',
    'sot_color_table', 'hmi_mag_color_table', 'suvi_color_table',
    'rhessi_color_table', 'std_gamma_2', 'euvi_color_table', 'solohri_lya1216_color_table',
]


CMAP_DATA_DIR = pathlib.Path(__file__).parent.absolute() / 'data'


def create_cdict(r, g, b):
    """
    Create the color tuples in the correct format.
    """
    i = np.linspace(0, 1, r.size)
    cdict = {name: list(zip(i, el / 255.0, el / 255.0))
             for el, name in [(r, 'red'), (g, 'green'), (b, 'blue')]}
    return cdict


def _cmap_from_rgb(r, g, b, name):
    cdict = create_cdict(r, g, b)
    return colors.LinearSegmentedColormap(name, cdict)


def cmap_from_rgb_file(name, fname):
    """
    Create a colormap from a RGB .csv file.

    The .csv file must have 3  equal-length columns of integer data, with values
    between 0 and 255, which are the red, green, and blue values for the colormap.

    Parameters
    ----------
    name : str
        Name of the colormap.
    fname : str
        Filename of data file. Relative to the sunpy colormap data directory.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    data = np.loadtxt(CMAP_DATA_DIR / fname, delimiter=',')
    if data.shape[1] != 3:
        raise RuntimeError(f'RGB data files must have 3 columns (got {data.shape[1]})')
    return _cmap_from_rgb(data[:, 0], data[:, 1], data[:, 2], name)


def get_idl3():
    # The following values describe color table 3 for IDL (Red Temperature)
    return np.loadtxt(CMAP_DATA_DIR / 'idl_3.csv', delimiter=',')


def solohri_lya1216_color_table():
    solohri_lya1216 = get_idl3()
    solohri_lya1216[:, 2] = solohri_lya1216[:, 0] * np.linspace(0, 1, 256)
    return _cmap_from_rgb(*solohri_lya1216.T, 'SolO EUI HRI Lyman Alpha')



def create_aia_wave_dict():
    idl_3 = get_idl3()
    r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]

    c0 = np.arange(256, dtype='f')
    c1 = (np.sqrt(c0) * np.sqrt(255.0)).astype('f')
    c2 = (np.arange(256)**2 / 255.0).astype('f')
    c3 = ((c1 + c2 / 2.0) * 255.0 / (c1.max() + c2.max() / 2.0)).astype('f')

    aia_wave_dict = {
        1600*u.angstrom: (c3, c3, c2),
        1700*u.angstrom: (c1, c0, c0),
        4500*u.angstrom: (c0, c0, b0 / 2.0),
        94*u.angstrom: (c2, c3, c0),
        131*u.angstrom: (g0, r0, r0),
        171*u.angstrom: (r0, c0, b0),
        193*u.angstrom: (c1, c0, c2),
        211*u.angstrom: (c1, c0, c3),
        304*u.angstrom: (r0, g0, b0),
        335*u.angstrom: (c2, c0, c1)
    }
    return aia_wave_dict


@u.quantity_input
def aia_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SDO AIA images.

    Based on aia_lct.pro part of SDO/AIA on SSWIDL written by Karel
    Schrijver (2010/04/12).

    Parameters
    ----------
    wavelength : `~astropy.units.quantity`
        Wavelength for the desired AIA color table.
    """
    aia_wave_dict = create_aia_wave_dict()
    try:
        r, g, b = aia_wave_dict[wavelength]
    except KeyError:
        raise ValueError("Invalid AIA wavelength. Valid values are "
                         "1600,1700,4500,94,131,171,193,211,304,335.")

    return _cmap_from_rgb(r, g, b, 'SDO AIA {:s}'.format(str(wavelength)))



@u.quantity_input
def eit_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SOHO EIT images.
    """
    # SOHO EIT Color tables
    # EIT 171 IDL Name EIT Dark Bot Blue
    # EIT 195 IDL Name EIT Dark Bot Green
    # EIT 284 IDL Name EIT Dark Bot Yellow
    # EIT 304 IDL Name EIT Dark Bot Red
    try:
        color = {171*u.angstrom: 'dark_blue', 195*u.angstrom: 'dark_green',
                 284*u.angstrom: 'yellow', 304*u.angstrom: 'dark_red',
                 }[wavelength]
    except KeyError:
        raise ValueError("Invalid EIT wavelength. Valid values are "
                         "171, 195, 284, 304.")

    return cmap_from_rgb_file('SOHO EIT {:s}'.format(str(wavelength)), f'eit_{color}.csv')



def sswidl_lasco_color_table(number):
    """
    Returns one of the SSWIDL-defined color tables for SOHO LASCO images.

    This function is included to allow users to access the SSWIDL-
    defined LASCO color tables provided by SunPy. It is recommended to
    use the function 'lasco_color_table' to obtain color tables for use
    with LASCO data and Helioviewer JP2 images.
    """
    try:
        return cmap_from_rgb_file(f'SOHO LASCO C{number}', f'lasco_c{number}.csv')
    except OSError:
        raise ValueError("Invalid LASCO number. Valid values are 2, 3.")



# Translated from the JP2Gen IDL SXT code lct_yla_gold.pro.  Might be better
# to explicitly copy the numbers from the IDL calculation.  This is a little
# more compact.
sxt_gold_r = np.concatenate((np.linspace(0, 255, num=185,
                                         endpoint=False), 255 * np.ones(71)))
sxt_gold_g = 255 * (np.arange(256)**1.25) / (255.0**1.25)
sxt_gold_b = np.concatenate((np.zeros(185), 255.0 * np.arange(71) / 71.0))

grayscale = np.arange(256)
grayscale.setflags(write=False)


def sxt_color_table(sxt_filter):
    """
    Returns one of the fundamental color tables for Yokhoh SXT images.
    """
    try:
        r, g, b = {
            'al': (sxt_gold_r, sxt_gold_g, sxt_gold_b),
            'wh': (grayscale, grayscale, grayscale)
        }[sxt_filter]
    except KeyError:
        raise ValueError("Invalid SXT filter type number. Valid values are "
                         "'al', 'wh'.")
    return _cmap_from_rgb(r, g, b, 'Yohkoh SXT {:s}'.format(sxt_filter.title()))



def xrt_color_table():
    """
    Returns the color table used for all Hinode XRT images.
    """
    idl_3 = get_idl3()
    r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]
    return _cmap_from_rgb(r0, g0, b0, 'Hinode XRT')



def cor_color_table(number):
    """
    Returns one of the fundamental color tables for STEREO coronagraph images.
    """
    # STEREO COR Color tables
    if number not in [1, 2]:
        raise ValueError("Invalid COR number. Valid values are " "1, 2.")

    return cmap_from_rgb_file(f'STEREO COR{number}', f'stereo_cor{number}.csv')


def trace_color_table(measurement):
    """
    Returns one of the standard color tables for TRACE JP2 files.
    """
    if measurement == 'WL':
        return cmap_from_rgb_file(f'TRACE {measurement}', 'grayscale.csv')

    try:
        return cmap_from_rgb_file(f'TRACE {measurement}', f'trace_{measurement}.csv')
    except OSError:
        raise ValueError(
            "Invalid TRACE filter waveband passed.  Valid values are "
            "171, 195, 284, 1216, 1550, 1600, 1700, WL")



def sot_color_table(measurement):
    """
    Returns one of the standard color tables for SOT files (following osdc
    convention).

    The relations between observation and color have been defined in
    hinode.py
    """
    idl_3 = get_idl3()
    r0, g0, b0 = idl_3[:, 0], idl_3[:, 1], idl_3[:, 2]
    try:
        r, g, b = {
            'intensity': (r0, g0, b0),
        }[measurement]
    except KeyError:
        raise ValueError(
            "Invalid (or not supported) SOT type. Valid values are: "
            "intensity")

    return _cmap_from_rgb(r, g, b, f'Hinode SOT {measurement:s}')



def iris_sji_color_table(measurement, aialike=False):
    """
    Return the standard color table for IRIS SJI files.
    """
    # base vectors for IRIS SJI color tables
    c0 = np.arange(0, 256)
    c1 = (np.sqrt(c0) * np.sqrt(255)).astype(np.uint8)
    c2 = (c0**2 / 255.).astype(np.uint8)
    c3 = ((c1 + c2 / 2.) * 255. / (np.max(c1) + np.max(c2) / 2.)).astype(
        np.uint8)
    c4 = np.zeros(256).astype(np.uint8)
    c4[50:256] = (1 / 165. * np.arange(0, 206)**2).astype(np.uint8)
    c5 = ((1 + c1 + c3.astype(np.uint)) / 2.).astype(np.uint8)

    rr = np.ones(256, dtype=np.uint8) * 255
    rr[0:176] = np.arange(0, 176) / 175. * 255.
    gg = np.zeros(256, dtype=np.uint8)
    gg[100:256] = np.arange(0, 156) / 155. * 255.
    bb = np.zeros(256, dtype=np.uint8)
    bb[150:256] = np.arange(0, 106) / 105. * 255.
    agg = np.zeros(256, dtype=np.uint8)
    agg[120:256] = np.arange(0, 136) / 135. * 255.
    abb = np.zeros(256, dtype=np.uint8)
    abb[190:256] = np.arange(0, 66) / 65. * 255.

    if aialike:
        color_table = {
            '1330': (c1, c0, c2),
            '1400': (rr, agg, abb),
            '2796': (rr, c0, abb),
            '2832': (c3, c3, c2),
        }
    else:
        color_table = {
            '1330': (rr, gg, bb),
            '1400': (c5, c2, c4),
            '2796': (c1, c3, c2),
            '2832': (c0, c0, c2),
        }

    color_table.update({
        '1600': (c1, c0, c0),
        '5000': (c1, c1, c0),
        'FUV': (rr, gg, bb),
        'NUV': (c1, c3, c2),
        'SJI_NUV': (c0, c0, c0)
    })

    try:
        r, g, b = color_table[measurement]
    except KeyError:
        raise ValueError("Invalid IRIS SJI waveband.  Valid values are \n" +
                         str(list(color_table.keys())))

    return _cmap_from_rgb(r, g, b, f'IRIS SJI {measurement:s}')


def hmi_mag_color_table():
    """
    Returns an alternate HMI Magnetogram color table; from Stanford
    University/JSOC.

    Examples
    --------
    >>> # Example usage for NRT data:
    >>> import sunpy.map
    >>> import matplotlib.pyplot as plt
    >>> hmi = sunpy.map.Map('fblos.fits')  # doctest: +SKIP
    >>> hmi.plot_settings['cmap'] = plt.get_cmap('hmimag')  # doctest: +SKIP
    >>> hmi.peek(vmin=-1500.0, vmax=1500.0)  # doctest: +SKIP

    >>> # OR (for a basic plot with pixel values on the axes)
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import sunpy.map
    >>> hmi = sunpy.map.Map('fblos.fits')  # doctest: +SKIP
    >>> plt.imshow(np.clip(hmi.data, -1500.0, 1500.0), cmap=plt.get_cmap('hmimag'), origin='lower')  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

    >>> # Example usage for science (Level 1.0) data:
    >>> import numpy as np
    >>> import sunpy.map
    >>> hmi = sunpy.map.Map('hmi.m_45s.2014.05.11_12_00_45_TAI.magnetogram.fits')  # doctest: +SKIP
    >>> hmir = hmi.rotate()  # doctest: +SKIP
    >>> hmir.plot_settings['cmap'] = plt.get_cmap('hmimag')  # doctest: +SKIP
    >>> hmir.peek(vmin=-1500.0, vmax=1500.0)  # doctest: +SKIP

    References
    ----------
    * `Stanford Colortable (pdf) <http://jsoc.stanford.edu/data/hmi/HMI_M.ColorTable.pdf>`_
    """
    return cmap_from_rgb_file('SDO HMI magnetogram', 'hmi_mag.csv')



def stereo_hi_color_table(camera):
    if camera not in [1, 2]:
        raise ValueError("Valid HI cameras are 1 and 2")
    return cmap_from_rgb_file(f'STEREO HI{camera}', f'hi{camera}.csv')


@u.quantity_input
def suvi_color_table(wavelength: u.angstrom):
    """
    Returns one of the fundamental color tables for SUVI images.
    SUVI uses AIA color tables.
    """
    aia_wave_dict = create_aia_wave_dict()
    try:
        if wavelength == 195*u.angstrom:
            r, g, b = aia_wave_dict[193*u.angstrom]
        elif wavelength == 284*u.angstrom:
            r, g, b = aia_wave_dict[335*u.angstrom]
        else:
            r, g, b = aia_wave_dict[wavelength]
    except KeyError:
        raise ValueError(
            "Invalid SUVI wavelength. Valid values are "
            "94, 131, 171, 195, 284, 304."
        )
    return _cmap_from_rgb(r, g, b, 'GOES-R SUVI {:s}'.format(str(wavelength)))



def rhessi_color_table():
    return cmap_from_rgb_file("rhessi", "rhessi.csv")



def std_gamma_2():
    return cmap_from_rgb_file("std_gamma_2", "std_gamma_2.csv")



def euvi_color_table(wavelength: u.angstrom):
    try:
        return cmap_from_rgb_file(f'EUVI {str(wavelength)}', f'euvi_{int(wavelength.value)}.csv')
    except OSError:
        raise ValueError(
            "Invalid EUVI wavelength. Valid values are "
            "171, 195, 284, 304."
        )


linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as (0, (1, 1)) or ':'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')]  # Same as '-.'

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]


    
