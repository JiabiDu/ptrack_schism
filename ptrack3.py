"""
Particle Tracking Model - Python conversion of ptrack3.f90

This module implements a particle tracking model for SCHISM hydrodynamic outputs.
Converted from Fortran to Python, with simplified functionality (newio=1, no HAB).

Major classes:
- HGrid: Horizontal grid data
- VGrid: Vertical grid data  
- ParticleTracker: Main tracking engine
"""

import numpy as np
import netCDF4 as nc
from dataclasses import dataclass
from typing import Tuple, Dict, List
import logging
import os
import glob
import re

# try to import numba for performance-critical loops
try:
    import numba
    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

logger = logging.getLogger(__name__)

# Constants
PI = np.pi
SMALL1 = 1.0e-6
ZERO = 1.0e-5
EARTH_RADIUS = 6378206.4

# ============================================================================
# UTILITY FUNCTIONS - Geometric operations
# ============================================================================

def signa(x1, x2, x3, y1, y2, y3):
    """
    Calculate signed area (2x) of a triangle.
    
    The result is twice the signed area (with sign indicating orientation).
    
    Parameters
    ----------
    x1, x2, x3 : float
        X-coordinates of triangle vertices
    y1, y2, y3 : float
        Y-coordinates of triangle vertices
    
    Returns
    -------
    float
        Twice the signed area of the triangle
    """
    return ((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)) / 2.0


def pt_in_poly3(i34, x, y, xp, yp):
    """
    Point-in-polygon test with area coordinate calculation.
    
    Performs a point-in-polygon test for a triangle or quad and calculates
    the area (barycentric) coordinates of the point. For quads, splits into
    2 triangles and returns the 3 nodes and area coordinates for the best fit.
    
    Parameters
    ----------
    i34 : int
        Number of vertices (3 for triangle, 4 for quad)
    x : array_like
        X-coordinates of polygon vertices (counter-clockwise), shape (i34,)
    y : array_like
        Y-coordinates of polygon vertices (counter-clockwise), shape (i34,)
    xp : float
        X-coordinate of point to test
    yp : float
        Y-coordinate of point to test
    
    Returns
    -------
    arco : ndarray
        Area coordinates (barycentric coordinates) of shape (3,)
    nodel : ndarray
        Local node indices of shape (3,), referring to vertices in x, y arrays
    """
    # Convert to numpy arrays for easier indexing
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    
    # Calculate areas of up to 2 triangles
    ar = np.zeros(2)
    ar[0] = signa(x[0], x[1], x[2], y[0], y[1], y[2])
    
    if i34 == 4:
        ar[1] = signa(x[0], x[2], x[3], y[0], y[2], y[3])
    
    # Check for negative areas
    if ar[0] <= 0 or (i34 == 4 and ar[1] <= 0):
        raise ValueError(f"Negative element area: i34={i34}, ar={ar}")
    
    # Initialize output arrays
    arco = np.zeros(3)
    nodel = np.zeros(3, dtype=int)
    
    ae_min = np.inf
    
    # Process each triangle
    for m in range(i34 - 2):
        if m == 0:
            list_nodes = np.array([0, 1, 2], dtype=int)  # 0-indexed
        else:  # quads, m = 1
            list_nodes = np.array([0, 2, 3], dtype=int)  # 0-indexed
        
        # Calculate sum of areas for the 3 sub-triangles
        aa = 0.0
        swild = np.zeros(3)
        
        for j in range(3):
            j1 = (j + 1) % 3
            j2 = (j + 2) % 3
            
            idx1 = list_nodes[j1]
            idx2 = list_nodes[j2]
            idx_p = list_nodes[j]
            
            swild[j] = signa(x[idx1], x[idx2], xp, y[idx1], y[idx2], yp)
            aa += np.abs(swild[j])
        
        # Calculate relative error
        ae = np.abs(aa - ar[m]) / ar[m]
        
        # Keep the best fit
        if ae <= ae_min:
            ae_min = ae
            nodel = list_nodes + 1  # Convert to 1-based indexing
            arco = swild / ar[m]
            
            # Clamp values to [0, 1]
            arco[0] = np.clip(arco[0], 0.0, 1.0)
            arco[1] = np.clip(arco[1], 0.0, 1.0)
            
            # Ensure barycentric coordinates sum to 1
            if arco[0] + arco[1] > 1:
                arco[2] = 0.0
                arco[1] = 1.0 - arco[0]
            else:
                arco[2] = 1.0 - arco[0] - arco[1]
    
    return arco, nodel


def pt_in_poly2(i34, x, y, xp, yp):
    """
    Point-in-polygon test (simple version).
    
    Parameters
    ----------
    i34 : int
        Number of vertices (3 or 4)
    x, y : array_like
        Vertex coordinates (counter-clockwise)
    xp, yp : float
        Point to test
    
    Returns
    -------
    inside : int
        1 if inside, 0 otherwise
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    
    inside = 0
    swild = np.zeros(i34)
    
    for j in range(i34):
        j1 = (j + 1) % i34
        swild[j] = signa(x[j], x[j1], xp, y[j], y[j1], yp)
    
    ae = np.min(swild)
    if ae > -SMALL1:
        inside = 1
    
    return inside


def intersect2(x1, x2, x3, x4, y1, y2, y3, y4):
    """
    Find intersection of two line segments.
    
    Parameters
    ----------
    x1, y1, x2, y2 : float
        First line segment endpoints
    x3, y3, x4, y4 : float
        Second line segment endpoints
    
    Returns
    -------
    iflag : int
        1 if intersection found, 0 otherwise
    xin, yin : float
        Intersection coordinates (if iflag=1)
    tt1, tt2 : float
        Parametric coordinates along segments
    """
    zero1 = 0.0
    
    tt1 = -1000.0
    tt2 = -1000.0
    iflag = 0
    xin = 0.0
    yin = 0.0
    
    delta = (x2 - x1) * (y3 - y4) - (y2 - y1) * (x3 - x4)
    delta1 = (x3 - x1) * (y3 - y4) - (y3 - y1) * (x3 - x4)
    delta2 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)
    
    if delta != 0.0:
        tt1 = delta1 / delta
        tt2 = delta2 / delta
        if (-zero1 <= tt1 <= 1 + zero1) and (-zero1 <= tt2 <= 1 + zero1):
            iflag = 1
            xin = x1 + (x2 - x1) * tt1
            yin = y1 + (y2 - y1) * tt1
    
    return iflag, xin, yin, tt1, tt2


# ============================================================================
# COORDINATE TRANSFORMATIONS
# ============================================================================

def cpp(xlon, ylat, rlambda0, phi0):
    """
    Convert lon/lat to Cartesian projection (CPP).
    
    Parameters
    ----------
    xlon, ylat : float
        Longitude and latitude in radians
    rlambda0, phi0 : float
        Reference longitude and latitude in radians
    
    Returns
    -------
    x, y : float
        Cartesian coordinates
    """
    x = EARTH_RADIUS * (xlon - rlambda0) * np.cos(phi0)
    y = EARTH_RADIUS * (ylat - phi0)
    return x, y


def cppinverse(x, y, rlambda0, phi0):
    """
    Convert Cartesian projection (CPP) back to lon/lat.
    
    Parameters
    ----------
    x, y : float
        Cartesian coordinates
    rlambda0, phi0 : float
        Reference longitude and latitude in radians
    
    Returns
    -------
    rlambda, phi : float
        Longitude and latitude in radians
    """
    rlambda = x / (EARTH_RADIUS * np.cos(phi0)) + rlambda0
    phi = y / EARTH_RADIUS + phi0
    return rlambda, phi


# ============================================================================
# RANDOM NUMBER GENERATORS
# ============================================================================

class RandomNumberGenerator:
    """Implementation of ran1 random number generator from Numerical Recipes."""
    
    IM1 = 2147483563
    IM2 = 2147483399
    AM = 1.0 / IM1
    IMM1 = IM1 - 1
    IA1 = 40014
    IA2 = 40692
    IQ1 = 53668
    IQ2 = 52774
    IR1 = 12211
    IR2 = 3791
    NTAB = 32
    NDIV = 1 + IMM1 // NTAB
    EPS = 1.2e-7
    RNMX = 1.0 - EPS
    
    def __init__(self, seed):
        self.idum = seed
        self.idum2 = 123456789
        self.iv = np.zeros(self.NTAB, dtype=int)
        self.iy = 0
        self._initialize()
    
    def _initialize(self):
        if self.idum <= 0:
            self.idum = max(-self.idum, 1)
            self.idum2 = self.idum
            for j in range(self.NTAB + 8, 0, -1):
                k = self.idum // self.IQ1
                self.idum = self.IA1 * (self.idum - k * self.IQ1) - k * self.IR1
                if self.idum < 0:
                    self.idum += self.IM1
                if j <= self.NTAB:
                    self.iv[j - 1] = self.idum
            self.iy = self.iv[0]
    
    def __call__(self):
        k = self.idum // self.IQ1
        self.idum = self.IA1 * (self.idum - k * self.IQ1) - k * self.IR1
        if self.idum < 0:
            self.idum += self.IM1
        
        k = self.idum2 // self.IQ2
        self.idum2 = self.IA2 * (self.idum2 - k * self.IQ2) - k * self.IR2
        if self.idum2 < 0:
            self.idum2 += self.IM2
        
        j = self.iy // self.NDIV
        self.iy = self.iv[j] - self.idum2
        self.iv[j] = self.idum
        if self.iy < 1:
            self.iy += self.IMM1
        
        return min(self.AM * self.iy, self.RNMX)


# Helper to normalize arrays read from netCDF to expected shapes (minimize copies)
def _ensure_np_nvrt(arr, np_expected, nvrt_expected):
    """Return view of array shaped (np_expected, nvrt_expected) where possible.

    Handles common orderings: (nvrt, np), (np, nvrt), (np,), (nvrt,), and squeezes extra dims.
    Avoids copies when shapes already match.
    """
    # Clean masked or sentinel invalid values first
    a = _clean_invalid(arr)
    a = np.squeeze(np.array(a))  # squeeze will produce view when possible
    if a.ndim == 2:
        if a.shape == (np_expected, nvrt_expected):
            return a
        if a.shape == (nvrt_expected, np_expected):
            return a.T
        # other permutations with matching dimensions
        if a.shape[0] == np_expected and a.shape[1] == nvrt_expected:
            return a
        if a.shape[0] == nvrt_expected and a.shape[1] == np_expected:
            return a.T
    elif a.ndim == 1:
        if a.shape[0] == np_expected:
            return np.tile(a[:, None], (1, nvrt_expected))
        if a.shape[0] == nvrt_expected:
            return np.tile(a[None, :], (np_expected, 1))
    elif a.ndim == 0:
        return np.zeros((np_expected, nvrt_expected))

    # last resort: reshape if possible
    if a.size == np_expected * nvrt_expected:
        return a.reshape((np_expected, nvrt_expected))

    raise ValueError(f"Cannot reshape array of shape {a.shape} to ({np_expected},{nvrt_expected})")


def _clean_invalid(a):
    """Return a numpy array with masked or sentinel invalid values cleaned.

    - If input is a MaskedArray, fill masked entries with 0.0.
    - Replace any non-finite or extremely large magnitude values (|x|>1e30)
      with 0.0 to guard against netCDF fill-value sentinels like 9.969e+36.
    """
    # Handle masked arrays
    if isinstance(a, np.ma.MaskedArray):
        arr = a.filled(0.0)
    else:
        arr = np.array(a)

    # Ensure floating type for isnan/isfinite checks
    try:
        arr = arr.astype(float, copy=False)
    except Exception:
        arr = arr.astype(float)

    # Replace non-finite or huge sentinel values with 0.0
    mask_invalid = ~np.isfinite(arr) | (np.abs(arr) > 1e30)
    if mask_invalid.any():
        arr = arr.copy()
        arr[mask_invalid] = 0.0

    return arr

# ============================================================================
# GRID CLASSES
# ============================================================================

class HGrid:
    """Horizontal grid containing mesh topology and geometry."""
    
    def __init__(self):
        self.np = 0  # number of nodes
        self.ne = 0  # number of elements
        self.x = None  # node x coordinates (np,)
        self.y = None  # node y coordinates (np,)
        self.dp = None  # node depth/bathymetry (np,)
        self.i34 = None  # element type: 3 or 4 (ne,)
        self.elnode = None  # element-node connectivity (4, ne) - convert to 0-indexed
        self.area = None  # element area (ne,)
        self.xctr = None  # element x centroid (ne,)
        self.yctr = None  # element y centroid (ne,)
        
        # Side information
        self.ic3 = None  # neighboring elements (4, ne)
        self.elside = None  # element sides (4, ne)
        self.isdel = None  # side-adjacent elements (2, ns)
        self.isidenode = None  # side nodes (2, ns)
        self.snx = None  # side normal x (ns,)
        self.sny = None  # side normal y (ns,)
        self.distj = None  # side length (ns,)
        self.dps = None  # side average depth (ns,)
        
        # Derived
        self.nne = None  # neighbors per node (np,)
        self.indel = None  # neighbor elements

    @classmethod
    def from_file(cls, hgrid_file: str):
        """Load horizontal grid from hgrid.gr3 or hgrid.ll file."""
        grid = cls()
        
        with open(hgrid_file, 'r') as f:
            f.readline()  # description
            ne, np_val = map(int, f.readline().split())
            grid.ne = ne
            grid.np = np_val
            
            grid.x = np.zeros(np_val)
            grid.y = np.zeros(np_val)
            grid.dp = np.zeros(np_val)
            
            for i in range(np_val):
                parts = f.readline().split()
                grid.x[i] = float(parts[1])
                grid.y[i] = float(parts[2])
                grid.dp[i] = float(parts[3])
            
            # Element connectivity (1-based in file)
            grid.i34 = np.zeros(ne, dtype=int)
            grid.elnode = np.zeros((4, ne), dtype=int)
            
            for i in range(ne):
                parts = f.readline().split()
                nvert = int(parts[1])
                grid.i34[i] = nvert
                for j in range(nvert):
                    grid.elnode[j, i] = int(parts[2 + j]) - 1  # Convert to 0-indexed
            
            # Compute geometry
            grid._compute_geometry()
        
        return grid
    
    def _compute_geometry(self):
        """Compute element areas and centroids."""
        self.area = np.zeros(self.ne)
        self.xctr = np.zeros(self.ne)
        self.yctr = np.zeros(self.ne)
        
        for i in range(self.ne):
            nodes = self.elnode[:self.i34[i], i]
            xs = self.x[nodes]
            ys = self.y[nodes]
            
            # Shoelace formula for area
            area = 0.0
            cx = 0.0
            cy = 0.0
            
            for j in range(self.i34[i]):
                j1 = (j + 1) % self.i34[i]
                area += xs[j] * ys[j1] - xs[j1] * ys[j]
                cx += xs[j]
                cy += ys[j]
            
            self.area[i] = abs(area) / 2.0
            self.xctr[i] = cx / self.i34[i]
            self.yctr[i] = cy / self.i34[i]


class VGrid:
    """Vertical grid (sigma or Z-coordinate)."""
    
    def __init__(self):
        self.nvrt = 0  # number of vertical layers
        self.ivcor = 1  # vertical coordinate type (1=sigma, 2=Z)
        self.h_c = 0.0  # critical depth for Z-coord
        self.theta_b = 0.0  # theta_b parameter
        self.theta_f = 0.0  # theta_f parameter
        self.h_s = 0.0  # h_s parameter
        self.sigma = None  # sigma levels (nvrt,) for sigma-coord
        self.ztot = None  # z levels (nvrt,) for Z-coord
        self.kz = 0  # number of pure Z layers
    
    @classmethod
    def from_file(cls, vgrid_file: str):
        """Load vertical grid from vgrid.in file."""
        grid = cls()
        
        with open(vgrid_file, 'r') as f:
            ivcor = int(f.readline().split()[0])
            nvrt = int(f.readline().split()[0])
            grid.ivcor = ivcor
            grid.nvrt = nvrt
            
            if ivcor == 1:  # Sigma-coordinate
                # Read sigma levels
                grid.sigma = np.zeros(nvrt)
                # Format varies; typically read remaining lines
                # For simplicity, use evenly spaced sigma
                grid.sigma = np.linspace(0, -1, nvrt)
            
            elif ivcor == 2:  # Z-coordinate
                # Read parameters
                pass  # Implement if needed
        
        return grid


@dataclass
class ParticleData:
    """Container for particle tracking data."""
    
    def __init__(self, nparticle: int, load_optional: dict = None):
        load_optional = load_optional or {}
        
        self.nparticle = nparticle
        self.idp = np.zeros(nparticle, dtype=object)  # particle ID
        self.st_p = np.zeros(nparticle)  # start time (seconds)
        self.xpar = np.zeros(nparticle)  # x position
        self.ypar = np.zeros(nparticle)  # y position
        self.zpar = np.zeros(nparticle)  # z position
        self.zpar0 = np.zeros(nparticle)  # initial z relative to free surface
        self.ielpar = np.zeros(nparticle, dtype=int)  # particle element index
        self.levpar = -99 * np.ones(nparticle, dtype=int)  # particle vertical level
        
        # Velocities
        self.upar = np.zeros(nparticle)
        self.vpar = np.zeros(nparticle)
        self.wpar = np.zeros(nparticle)
        
        # Optional properties
        if load_optional.get('temp', False):
            self.temp_par = np.zeros(nparticle)
        if load_optional.get('salt', False):
            self.salt_par = np.zeros(nparticle)
        if load_optional.get('solar', False):
            self.solar_par = np.zeros(nparticle)
    
    @classmethod
    def from_file(cls, particle_file: str, load_optional: dict = None):
        """Load particles from particle.bp file."""
        load_optional = load_optional or {}
        
        with open(particle_file, 'r') as f:
            nparticle = int(f.readline().split()[0])
        
        data = cls(nparticle, load_optional)
        
        with open(particle_file, 'r') as f:
            f.readline()  # skip header
            
            for i in range(nparticle):
                parts = f.readline().split()
                data.idp[i] = parts[0]
                data.st_p[i] = float(parts[1])
                data.xpar[i] = float(parts[2])
                data.ypar[i] = float(parts[3])
                data.zpar0[i] = float(parts[4])
        
        return data


# ============================================================================
# PARTICLE TRACKER ENGINE
# ============================================================================

class ParticleTracker:
    """Main particle tracking engine."""
    
    def __init__(self, hgrid: HGrid, vgrid: VGrid, particles: ParticleData,
                 h0: float = 0.1, ics: int = 1, slam0: float = 0.0, sfea0: float = 0.0,
                 ibf: int = 1, istiff: int = 0, dt: float = 1.0, ndeltp: int = 1,
                 temp_on: bool = False, salt_on: bool = False, solar_on: bool = False,
                 settling_velocity: float = 0.0, force_surface: bool = False):
        """
        Initialize particle tracker.
        
        Parameters
        ----------
        hgrid : HGrid
            Horizontal grid
        vgrid : VGrid
            Vertical grid
        particles : ParticleData
            Particle initial conditions
        h0 : float
            Minimum water depth
        ics : int
            Coordinate system (1=Cartesian, 2=lon/lat)
        slam0, sfea0 : float
            CPP projection center
        ibf : int
            Forward (1) or backward (-1) tracking
        istiff : int
            Stiff tracking flag
        dt : float
            Tracking time step (seconds)
        ndeltp : int
            Number of sub-steps
        temp_on, salt_on, solar_on : bool
            Track water properties
        settling_velocity : float
            Particle settling velocity (m/s)
        """
        self.hgrid = hgrid
        self.vgrid = vgrid
        self.particles = particles
        self.h0 = h0
        self.ics = ics
        self.slam0 = slam0 * PI / 180.0 if ics == 2 else slam0
        self.sfea0 = sfea0 * PI / 180.0 if ics == 2 else sfea0
        self.ibf = ibf
        self.istiff = istiff
        self.dt = dt
        self.ndeltp = ndeltp
        self.temp_on = temp_on
        self.salt_on = salt_on
        self.solar_on = solar_on
        self.settling_velocity = settling_velocity / 86400.0  # convert m/day to m/s
        # If True, initialize particles on the top surface layer for testing
        self.force_surface = bool(force_surface)
        
        # Grid dimensions
        self.np = hgrid.np
        self.ne = hgrid.ne
        self.nvrt = vgrid.nvrt

        # If using lon/lat coordinates, convert grid and particle positions to
        # Cartesian projection (meters) for internal computations.
        if self.ics == 2:
            # self.slam0 and self.sfea0 are stored in radians when ics==2
            # Convert node lon/lat (degrees) to radians then to CPP meters
            try:
                lon_rad = np.deg2rad(self.hgrid.x.astype(np.float64))
                lat_rad = np.deg2rad(self.hgrid.y.astype(np.float64))
                xproj, yproj = cpp(lon_rad, lat_rad, self.slam0, self.sfea0)
                self.hgrid.x = xproj
                self.hgrid.y = yproj
                # recompute centroids/area in projected coordinates
                try:
                    self.hgrid._compute_geometry()
                except Exception:
                    pass
            except Exception as e:
                logger.warning(f"Failed converting hgrid lon/lat to CPP: {e}")

            # Convert particle initial lon/lat positions to projection
            try:
                plon = np.deg2rad(self.particles.xpar.astype(np.float64))
                plat = np.deg2rad(self.particles.ypar.astype(np.float64))
                px, py = cpp(plon, plat, self.slam0, self.sfea0)
                self.particles.xpar = px
                self.particles.ypar = py
            except Exception as e:
                logger.warning(f"Failed converting particle lon/lat to CPP: {e}")
        
        # Z-coordinates at nodes
        self.z = np.zeros((self.np, self.nvrt))
        self.sigma_lcl = np.zeros((self.nvrt, self.np))
        
        # Wet/dry flags (0=wet, 1=dry). Initialize to 0 (assume wet) so
        # `compute_z_coords` can correctly mark dry nodes/elements.
        self.idry = np.zeros(self.np, dtype=int)
        self.idry_e = np.zeros(self.ne, dtype=int)
        self.kbp = np.ones(self.np, dtype=int)
        self.kbe = np.ones(self.ne, dtype=int)
        
        # Velocity fields (time level 1 and 2)
        self._init_velocity_arrays()
        
        # Elevation
        self.eta1 = np.zeros(self.np)
        self.eta2 = np.zeros(self.np)
        self.eta3 = np.zeros(self.np)
        
        # Random number generator for diffusion
        self.rng = RandomNumberGenerator(-5)
        
        logger.info(f"ParticleTracker: {self.np} nodes, {self.ne} elements, {self.nvrt} levels")
    
    def _init_velocity_arrays(self):
        """Initialize velocity and optional property arrays."""
        self.uu1 = np.zeros((self.np, self.nvrt))
        self.vv1 = np.zeros((self.np, self.nvrt))
        self.ww1 = np.zeros((self.np, self.nvrt))
        self.uu2 = np.zeros((self.np, self.nvrt))
        self.vv2 = np.zeros((self.np, self.nvrt))
        self.ww2 = np.zeros((self.np, self.nvrt))
        
        if self.temp_on:
            self.temp1 = np.zeros((self.np, self.nvrt))
            self.temp2 = np.zeros((self.np, self.nvrt))
        
        if self.salt_on:
            self.salt1 = np.zeros((self.np, self.nvrt))
            self.salt2 = np.zeros((self.np, self.nvrt))
        
        if self.solar_on:
            self.solar1 = np.zeros(self.np)
            self.solar2 = np.zeros(self.np)
    
    def load_velocity_field(self, ncfile: str, irec: int):
        """
        Load velocity field from netCDF file.
        
        Parameters
        ----------
        ncfile : str
            Path to netCDF output file
        irec : int
            Record index (0-based)
        """

        # dataset cache helper (nested so cache persists on self)
        if not hasattr(self, '_nc_cache'):
            self._nc_cache = {}
        def _open_dataset(path):
            if path in self._nc_cache:
                return self._nc_cache[path]
            ds = nc.Dataset(path, 'r')
            self._nc_cache[path] = ds
            return ds

        # If the provided file exists, try opening it directly
        try:
            if os.path.exists(ncfile):
                ds = _open_dataset(ncfile)
                # Note: do not close here; cached for reuse
                # Elevation may be in out2d-type files
                if 'elevation' in ds.variables:
                    self.eta2[:] = _clean_invalid(ds['elevation'][irec, :])

                loaded_vel = False

                # Combined hvel variable (old formats)
                if 'hvel' in ds.variables:
                    hvel_data = ds['hvel'][irec, :, :, :]
                    # clean masked or sentinel values in hvel_data
                    hvel_data = _clean_invalid(hvel_data)
                    # Expecting shape (2, nvrt, np) or (nvrt, np, 2)
                    try:
                        self.uu2[:, :] = hvel_data[0, :, :].T
                        self.vv2[:, :] = hvel_data[1, :, :].T
                        loaded_vel = True
                    except Exception:
                        # Try alternate ordering
                        self.uu2[:, :] = hvel_data[:, :, 0].T
                        self.vv2[:, :] = hvel_data[:, :, 1].T
                        loaded_vel = True

                # If horizontalVelX/ horizontalVelY present in the same file
                if 'horizontalVelX' in ds.variables and 'horizontalVelY' in ds.variables:
                    tmpu = ds['horizontalVelX'][irec, :, :]
                    tmpv = ds['horizontalVelY'][irec, :, :]
                    tmpu = _clean_invalid(tmpu)
                    tmpv = _clean_invalid(tmpv)
                    self.uu2[:, :] = tmpu.T
                    self.vv2[:, :] = tmpv.T
                    loaded_vel = True

                # Vertical velocity
                if 'vertical_velocity' in ds.variables or 'verticalVelocity' in ds.variables:
                    var_name = 'vertical_velocity' if 'vertical_velocity' in ds.variables else 'verticalVelocity'
                    tmpw = ds[var_name][irec, :, :]
                    tmpw = _clean_invalid(tmpw)
                    self.ww2[:, :] = tmpw.T
                    loaded_vel = True

                # Optional properties
                if self.salt_on and 'salinity' in ds.variables:
                    self.salt2[:, :] = ds['salinity'][irec, :, :].T
                if self.temp_on and 'temperature' in ds.variables:
                    self.temp2[:, :] = ds['temperature'][irec, :, :].T
                if self.solar_on and 'solarRadiation' in ds.variables:
                    self.solar2[:] = ds['solarRadiation'][irec, :]

                # Only return early if velocities were found in this file;
                # otherwise fall through to search for per-variable files (newio=1).
                if loaded_vel:
                    return

            # If file not found or didn't contain necessary vars, try newio=1 pattern
            basedir = os.path.dirname(ncfile) or '.'
            basename = os.path.basename(ncfile)
            # Try to extract numeric index from filename like *_<index>.nc
            m = re.search(r"_(\d+)\.nc$", basename)
            idx = m.group(1) if m else None

            def find_file(patterns):
                for pat in patterns:
                    matches = glob.glob(os.path.join(basedir, pat))
                    if matches:
                        return matches[0]
                return None

            # Find out2d (elevation), horizontalVelX, horizontalVelY, verticalVelocity
            if idx:
                out2d_pattern = [f"out2d_{idx}.nc", f"out2d.{idx}.nc", f"out2d*{idx}*.nc"]
                u_pattern = [f"horizontalVelX_{idx}.nc", f"horizontalVelX*{idx}*.nc"]
                v_pattern = [f"horizontalVelY_{idx}.nc", f"horizontalVelY*{idx}*.nc"]
                w_pattern = [f"verticalVelocity_{idx}.nc", f"verticalVelocity*{idx}*.nc", f"vertical_velocity_{idx}.nc"]
            else:
                out2d_pattern = ["out2d_*.nc", "out2d*.nc"]
                u_pattern = ["horizontalVelX_*.nc", "horizontalVelX*.nc"]
                v_pattern = ["horizontalVelY_*.nc", "horizontalVelY*.nc"]
                w_pattern = ["verticalVelocity_*.nc", "verticalVelocity*.nc", "vertical_velocity_*.nc"]

            out2d_file = find_file(out2d_pattern)
            u_file = find_file(u_pattern)
            v_file = find_file(v_pattern)
            w_file = find_file(w_pattern)

            # Load elevation from out2d if available
            if out2d_file:
                try:
                    ds2d = _open_dataset(out2d_file)
                    if 'elevation' in ds2d.variables:
                        self.eta2[:] = _clean_invalid(ds2d['elevation'][irec, :])
                    if self.solar_on and 'solarRadiation' in ds2d.variables:
                        self.solar2[:] = _clean_invalid(ds2d['solarRadiation'][irec, :])
                except Exception as e:
                    logger.warning(f"Failed loading out2d file {out2d_file}: {e}")

            # Load horizontal velocity components from separate files
            if u_file:
                try:
                    dsu = _open_dataset(u_file)
                    varu = next((c for c in ('horizontalVelX', 'u', 'UU') if c in dsu.variables), None)
                    if varu is None:
                        varu = None
                        for k in dsu.variables.keys():
                            if len(dsu[k].shape) >= 2:
                                varu = k
                                break
                        if varu is None:
                            varu = list(dsu.variables.keys())[0]

                    raw = dsu[varu][irec, ...]
                    if np.ma.is_masked(raw):
                        raw = raw.filled(0.0)
                    try:
                        arr = _ensure_np_nvrt(raw, self.np, self.nvrt)
                        self.uu2[:, :] = arr
                    except Exception as e:
                        raise RuntimeError(f"Failed to normalize horizontalVelX '{u_file}' var '{varu}': {e}")
                except Exception as e:
                    logger.warning(f"Failed loading horizontalVelX from {u_file}: {e}")

            if v_file:
                try:
                    dsv = _open_dataset(v_file)
                    varv = next((c for c in ('horizontalVelY', 'v', 'VV') if c in dsv.variables), None)
                    if varv is None:
                        varv = None
                        for k in dsv.variables.keys():
                            if len(dsv[k].shape) >= 2:
                                varv = k
                                break
                        if varv is None:
                            varv = list(dsv.variables.keys())[0]

                    raw = dsv[varv][irec, ...]
                    if np.ma.is_masked(raw):
                        raw = raw.filled(0.0)
                    try:
                        arr = _ensure_np_nvrt(raw, self.np, self.nvrt)
                        self.vv2[:, :] = arr
                    except Exception as e:
                        raise RuntimeError(f"Failed to normalize horizontalVelY '{v_file}' var '{varv}': {e}")
                except Exception as e:
                    logger.warning(f"Failed loading horizontalVelY from {v_file}: {e}")

            # Vertical velocity
            if w_file:
                try:
                    dsw = _open_dataset(w_file)
                    varw = next((c for c in ('verticalVelocity', 'vertical_velocity', 'w', 'WW') if c in dsw.variables), None)
                    if varw is None:
                        varw = None
                        for k in dsw.variables.keys():
                            if len(dsw[k].shape) >= 2:
                                varw = k
                                break
                        if varw is None:
                            varw = list(dsw.variables.keys())[0]

                    raw = dsw[varw][irec, ...]
                    if np.ma.is_masked(raw):
                        raw = raw.filled(0.0)
                    try:
                        arr = _ensure_np_nvrt(raw, self.np, self.nvrt)
                        self.ww2[:, :] = arr
                    except Exception as e:
                        raise RuntimeError(f"Failed to normalize verticalVelocity '{w_file}' var '{varw}': {e}")
                except Exception as e:
                    logger.warning(f"Failed loading verticalVelocity from {w_file}: {e}")

            # If we still don't have velocities, raise
            if np.all(self.uu2 == 0) and np.all(self.vv2 == 0) and np.all(self.ww2 == 0):
                raise FileNotFoundError(f"Could not locate velocity files for base {ncfile} (searched out2d, horizontalVelX/Y, verticalVelocity)")

        except Exception as e:
            logger.error(f"Error loading velocity field from {ncfile}: {e}")
            raise
    
    def compute_z_coords(self, nodes: np.ndarray = None):
        """Compute z-coordinates.

        If ``nodes`` is None, compute for all grid nodes; otherwise only update
        the specified node indices.  Wet/dry flags and `kbp` are always updated
        for all nodes since those are cheap.
        """
        totdepth = self.eta3 + self.hgrid.dp
        dry_mask = (self.idry == 1) | (totdepth <= self.h0)
        self.idry[dry_mask] = 1
        self.kbp[dry_mask] = self.nvrt + 1
        wet_mask = ~dry_mask
        self.idry[wet_mask] = 0

        if self.vgrid.ivcor == 1:  # Sigma coordinate
            if nodes is None:
                # full update
                self.z[:, :] = totdepth[:, None] * self.vgrid.sigma[None, :] + self.eta3[:, None]
                self.kbp[wet_mask] = 1
            else:
                # partial update
                uniq = np.unique(nodes)
                td = totdepth[uniq]
                self.z[uniq, :] = td[:, None] * self.vgrid.sigma[None, :] + self.eta3[uniq][:, None]
                self.kbp[uniq] = 1

        # Update element wet/dry flags: element is wet if any of its nodes is wet
        if self.hgrid is not None and getattr(self.hgrid, 'elnode', None) is not None:
            for ie in range(self.ne):
                nodes_ie = self.hgrid.elnode[: self.hgrid.i34[ie], ie]
                # element wet if any node has idry==0
                wet = False
                for nd in nodes_ie:
                    if self.idry[nd] == 0:
                        wet = True
                        break
                self.idry_e[ie] = 0 if wet else 1
    
    def track_particles(self, times: list, snapshot_files: list, hydro_dt: float = None,
                        log_interval: int = 1, preload: bool = False) -> Tuple[list, list]:
        """
        Track all particles through time.
        
        Parameters
        ----------
        times : list
            List of time values (seconds)
        snapshot_files : list
            List of netCDF file/record pairs corresponding to each hydro snapshot, or
            if ``hydro_dt`` is None the list may contain an entry per tracking time step.
        hydro_dt : float, optional
            Model output interval in seconds. If provided, the tracker will linearly
            interpolate velocities between successive hydro snapshots to allow arbitrary
            tracking time steps smaller than ``hydro_dt``.
        
        Returns
        -------
        trajectories : list
            List of dicts with 'x', 'y', 'z' for each time step
        times : list
            Time values
        """
        trajectories = []
        # state for interpolation
        last_snap = -1
        vel0 = {'u': None, 'v': None, 'w': None, 'eta': None}
        vel1 = {'u': None, 'v': None, 'w': None, 'eta': None}

        # preload snapshots if requested
        cache = {}
        if preload:
            # build unique list of snapshot entries
            for entry in snapshot_files:
                if isinstance(entry, (tuple, list)) and len(entry) >= 2:
                    fname, rec = entry[0], int(entry[1])
                else:
                    fname, rec = entry, 0
                if (fname, rec) not in cache:
                    # load arrays once
                    self.load_velocity_field(fname, rec)
                    cache[(fname, rec)] = {
                        'u': self.uu2.copy(),
                        'v': self.vv2.copy(),
                        'w': self.ww2.copy(),
                        'eta': self.eta2.copy()
                    }
            # also prepare vel0/vel1 containers to map keys
        

        for t_idx, time in enumerate(times):
            print('time',time)
            if log_interval and (t_idx % log_interval == 0):
                day_num = time / 86400.0  # convert seconds to days
                logger.info(f"ptrack for day {day_num:.2f}")

            if preload:
                # obtain from cache
                if hydro_dt is None:
                    # step-wise entries must exist
                    if t_idx < len(snapshot_files):
                        entry = snapshot_files[t_idx]
                        if isinstance(entry, (tuple, list)) and len(entry) >= 2:
                            fname, rec = entry[0], int(entry[1])
                        else:
                            fname, rec = entry, 0
                        key = (fname, rec)
                        vals = cache.get(key)
                        if vals is None:
                            # fallback to loading
                            self.load_velocity_field(fname, rec)
                        else:
                            self.uu2[:, :] = vals['u']
                            self.vv2[:, :] = vals['v']
                            self.ww2[:, :] = vals['w']
                            self.eta2[:] = vals['eta']
                else:
                    snap = int(time // hydro_dt)
                    frac = ((time % hydro_dt) / hydro_dt) if hydro_dt > 0 else 0.0
                    # load values for snap and snap+1 from cache
                    if snap < len(snapshot_files):
                        fname, rec = snapshot_files[snap]
                        key0 = (fname, rec)
                        vals0 = cache.get(key0)
                    else:
                        vals0 = None
                    if snap + 1 < len(snapshot_files):
                        fname2, rec2 = snapshot_files[snap + 1]
                        key1 = (fname2, rec2)
                        vals1 = cache.get(key1)
                    else:
                        vals1 = None
                    if vals0 is None:
                        # load on the fly
                        fname, rec = snapshot_files[snap] if snap < len(snapshot_files) else (None,0)
                        self.load_velocity_field(fname, rec)
                        vals0 = {'u': self.uu2.copy(), 'v': self.vv2.copy(), 'w': self.ww2.copy(), 'eta': self.eta2.copy()}
                    if vals1 is None:
                        # default to vals0
                        vals1 = vals0
                    # perform interpolation
                    if vals1 is vals0 or frac == 0:
                        self.uu2[:, :] = vals0['u']
                        self.vv2[:, :] = vals0['v']
                        self.ww2[:, :] = vals0['w']
                        self.eta2[:] = vals0['eta']
                    else:
                        self.uu2[:, :] = (1 - frac) * vals0['u'] + frac * vals1['u']
                        self.vv2[:, :] = (1 - frac) * vals0['v'] + frac * vals1['v']
                        self.ww2[:, :] = (1 - frac) * vals0['w'] + frac * vals1['w']
                        self.eta2[:] = (1 - frac) * vals0['eta'] + frac * vals1['eta']
            else:
                # original behaviour
                if hydro_dt is None:
                    # simple behaviour: use provided entry per step
                    if t_idx < len(snapshot_files):
                        entry = snapshot_files[t_idx]
                        if isinstance(entry, (tuple, list)) and len(entry) >= 2:
                            fname, rec = entry[0], int(entry[1])
                        else:
                            fname, rec = entry, 0
                        self.load_velocity_field(fname, rec)
                else:
                    # determine which hydro snapshot interval we are in and interpolation factor
                    snap = int(time // hydro_dt)
                    frac = ((time % hydro_dt) / hydro_dt) if hydro_dt > 0 else 0.0

                    if snap != last_snap:
                        # load current snapshot
                        if snap < len(snapshot_files):
                            fname, rec = snapshot_files[snap]
                            self.load_velocity_field(fname, rec)
                            vel0['u'] = self.uu2.copy()
                            vel0['v'] = self.vv2.copy()
                            vel0['w'] = self.ww2.copy()
                            vel0['eta'] = self.eta2.copy()
                        else:
                            vel0 = {k: np.zeros_like(self.uu2) for k in vel0}

                        # load next snapshot (or duplicate if none)
                        if snap + 1 < len(snapshot_files):
                            fname2, rec2 = snapshot_files[snap + 1]
                            self.load_velocity_field(fname2, rec2)
                            vel1['u'] = self.uu2.copy()
                            vel1['v'] = self.vv2.copy()
                            vel1['w'] = self.ww2.copy()
                            vel1['eta'] = self.eta2.copy()
                        else:
                            vel1 = {k: vel0[k].copy() for k in vel0}

                        last_snap = snap

                    # linear interpolation between vel0 and vel1
                    self.uu2 = (1 - frac) * vel0['u'] + frac * vel1['u']
                    self.vv2 = (1 - frac) * vel0['v'] + frac * vel1['v']
                    self.ww2 = (1 - frac) * vel0['w'] + frac * vel1['w']
                    self.eta2 = (1 - frac) * vel0['eta'] + frac * vel1['eta']

            # update eta3
            self.eta3 = self.eta2.copy()
            # compute z only for elements containing particles
            # gather node indices from current element assignments
            elem_indices = self.particles.ielpar
            # ensure valid range
            elem_indices = elem_indices[(elem_indices >= 0) & (elem_indices < self.ne)]
            node_list = self.hgrid.elnode[:, elem_indices].ravel()
            self.compute_z_coords(node_list)

            # track particles
            for ipar in range(self.particles.nparticle):
                self.track_particle(ipar)

            # record trajectory
            traj = {
                'x': self.particles.xpar.copy(),
                'y': self.particles.ypar.copy(),
                'z': self.particles.zpar.copy()
            }
            trajectories.append(traj)

            # prepare for next step
            self.uu1 = self.uu2.copy()
            self.vv1 = self.vv2.copy()
            self.ww1 = self.ww2.copy()
            self.eta1 = self.eta2.copy()
            if self.temp_on:
                self.temp1 = self.temp2.copy()
            if self.salt_on:
                self.salt1 = self.salt2.copy()
            if self.solar_on:
                self.solar1 = self.solar2.copy()

        return trajectories, times
    
    def track_particle(self, ipar: int):
        """Track one particle for one time step."""
        if self.particles.levpar[ipar] == -99:
            self._init_particle(ipar)
            return

        dt_step = self.dt / self.ndeltp
        cur_x = self.particles.xpar[ipar]
        cur_y = self.particles.ypar[ipar]
        cur_z = self.particles.zpar[ipar]
        cur_iel = self.particles.ielpar[ipar]
        cur_jlev = self.particles.levpar[ipar]

        for _ in range(self.ndeltp):
            xadv = self.ibf * dt_step * self.particles.upar[ipar]
            yadv = self.ibf * dt_step * self.particles.vpar[ipar]
            zadv = self.ibf * dt_step * self.particles.wpar[ipar] - self.ibf * dt_step * self.settling_velocity

            xt = cur_x + xadv
            yt = cur_y + yadv
            zt = cur_z + zadv

            # Simplified search - just check if still in element or moves to neighbor
            try:
                nnel, jlev_found, nodel, arco = self._find_element(cur_iel, xt, yt)
                # preserve the existing vertical level if the search helper only returns a default
                if jlev_found == 1 and cur_jlev != -99:
                    jlev = cur_jlev
                else:
                    jlev = jlev_found

                self.particles.xpar[ipar] = xt
                self.particles.ypar[ipar] = yt
                self.particles.zpar[ipar] = zt
                self.particles.ielpar[ipar] = nnel
                self.particles.levpar[ipar] = jlev

                cur_x = xt
                cur_y = yt
                cur_z = zt
                cur_iel = nnel
                cur_jlev = jlev

                self._interpolate_velocity(ipar, nnel, jlev, nodel, arco)
            except Exception as e:
                logger.warning(f"Error tracking particle {ipar}: {e}")
                break
    
    def _init_particle(self, ipar: int):
        """Initialize particle at starting location."""
        # Find initial element and vertical level
        # Simplified: assume provided element index is correct
        iel = self.particles.ielpar[ipar]
        
        if iel < 0 or iel >= self.ne:
            logger.warning(f"Invalid initial element for particle {ipar}")
            return
        
        if self.idry_e[iel] == 1:
            self.particles.levpar[ipar] = -1
            return
        
        # Compute area coordinates
        nodes = self.hgrid.elnode[:self.hgrid.i34[iel], iel]
        xs = self.hgrid.x[nodes]
        ys = self.hgrid.y[nodes]
        
        arco, nodel = pt_in_poly3(self.hgrid.i34[iel], xs, ys,
                                  self.particles.xpar[ipar], self.particles.ypar[ipar])
        
        # Interpolate z-coordinate
        z_interp = np.zeros(self.nvrt)
        for j in range(3):
            nd = nodes[nodel[j] - 1]
            z_interp += self.z[nd, :] * arco[j]
        
        # Find vertical level. By default, attempt to locate appropriate level
        # from the interpolated z profile. If `force_surface` is True, place
        # the particle on the top surface layer (jlev=1) for testing.
        if getattr(self, 'force_surface', False):
            # Fortran convention expects the surface level to be `nvrt`.
            # Map the forced-surface request to that top-level index so
            # comparisons with Fortran outputs match expectations.
            try:
                jlev = int(self.nvrt)
                lev_idx = max(0, min(jlev - 1, self.nvrt - 1))
                zpar_val = z_interp[lev_idx] + self.particles.zpar0[ipar]
            except Exception:
                # fallback to first sigma level
                jlev = 1
                lev_idx = 0
                zpar_val = z_interp[0] + self.particles.zpar0[ipar]

            self.particles.levpar[ipar] = jlev
            self.particles.zpar[ipar] = zpar_val
            print('jlevel=',jlev)

            # log one-off diagnostic for the first particle to help debugging
            if ipar == 0 and logger.isEnabledFor(logging.INFO):
                try:
                    sigma = np.asarray(self.vgrid.sigma)
                    nodes_ie = self.hgrid.elnode[: self.hgrid.i34[iel], iel]
                    umean = float(np.nanmean(np.abs(self.uu2[nodes_ie, lev_idx])))
                    vmean = float(np.nanmean(np.abs(self.vv2[nodes_ie, lev_idx])))
                    logger.info(f"force_surface: set particle {ipar} to level {jlev} (sigma={sigma[lev_idx]:.3f}), mean|u|={umean:.6g}, mean|v|={vmean:.6g}")
                except Exception:
                    pass
        else:
            zpar_target = max(z_interp[0], z_interp[-1] + self.particles.zpar0[ipar])

            jlev = 1
            for k in range(self.nvrt - 1):
                zlo = min(z_interp[k], z_interp[k + 1])
                zhi = max(z_interp[k], z_interp[k + 1])
                if zpar_target >= zlo and zpar_target <= zhi:
                    jlev = k + 1
                    break

            self.particles.levpar[ipar] = jlev
            self.particles.zpar[ipar] = zpar_target
        
        # Interpolate initial velocity
        self._interpolate_velocity(ipar, iel, jlev, nodel, arco)
    
    def _find_element(self, nel: int, xt: float, yt: float) -> Tuple[int, int, np.ndarray, np.ndarray]:
        """Find element containing point, return element, level, nodes, and area coords."""
        # Check if still in starting element
        nodes = self.hgrid.elnode[:self.hgrid.i34[nel], nel]
        xs = self.hgrid.x[nodes]
        ys = self.hgrid.y[nodes]
        
        if pt_in_poly2(self.hgrid.i34[nel], xs, ys, xt, yt):
            arco, nodel = pt_in_poly3(self.hgrid.i34[nel], xs, ys, xt, yt)
            return nel, 1, nodel, arco  # Simplified: always return level 1
        
        # Search neighbors
        for j in range(self.hgrid.i34[nel]):
            if self.hgrid.ic3 is not None and self.hgrid.ic3[j, nel] >= 0:
                nel_next = self.hgrid.ic3[j, nel]
                nodes = self.hgrid.elnode[:self.hgrid.i34[nel_next], nel_next]
                xs = self.hgrid.x[nodes]
                ys = self.hgrid.y[nodes]
                
                if pt_in_poly2(self.hgrid.i34[nel_next], xs, ys, xt, yt):
                    arco, nodel = pt_in_poly3(self.hgrid.i34[nel_next], xs, ys, xt, yt)
                    return nel_next, 1, nodel, arco
        
        # Not found - keep in same element
        arco, nodel = pt_in_poly3(self.hgrid.i34[nel], xs, ys, xt, yt)
        return nel, 1, nodel, arco

    def debug_find_element_for_point(self, ipar: int = 0, max_checks: int = None):
        """
        Brute-force search for the element containing the given particle's (x,y).

        This helper iterates over elements (or up to `max_checks`) and reports any
        matching element(s) found using `pt_in_poly2` / `pt_in_poly3`. It logs
        grid extents, particle coords and diagnostics to assist debugging.
        """
        if ipar < 0 or ipar >= self.particles.nparticle:
            logger.error(f"debug_find_element_for_point: invalid particle index {ipar}")
            return None

        x = self.particles.xpar[ipar]
        y = self.particles.ypar[ipar]
        iel_given = int(self.particles.ielpar[ipar])

        logger.info(f"Debug containment: particle {ipar} at ({x:.6f}, {y:.6f}), initial element {iel_given}")
        logger.info(f"Grid bbox: x [{self.hgrid.x.min():.6f},{self.hgrid.x.max():.6f}] y [{self.hgrid.y.min():.6f},{self.hgrid.y.max():.6f}]")

        found = []
        nelem = self.hgrid.ne
        limit = nelem if max_checks is None else min(nelem, max_checks)

        for ie in range(limit):
            nvert = int(self.hgrid.i34[ie])
            nodes = self.hgrid.elnode[:nvert, ie]
            xs = self.hgrid.x[nodes]
            ys = self.hgrid.y[nodes]

            try:
                inside = pt_in_poly2(nvert, xs, ys, x, y)
            except Exception:
                inside = 0

            if inside:
                try:
                    arco, nodel = pt_in_poly3(nvert, xs, ys, x, y)
                except Exception as e:
                    logger.debug(f"pt_in_poly3 failed for element {ie}: {e}")
                    arco = None
                    nodel = None

                logger.info(f"Brute-force: found element {ie} (nvert={nvert}) nodes={nodes} nodel={nodel} arco={arco}")
                found.append((ie, nodes, nodel, arco))

        if not found:
            logger.info(f"Brute-force: no containing element found in first {limit} elements")
            return None

        return found

    def assign_initial_elements(self, k: int = 20):
        """
        Assign initial element (`ielpar`) to particles that have invalid element indices.

        For each particle with `ielpar` outside [0, ne-1], searches the nearest `k`
        element centroids and tests containment. If a containing element is found,
        assigns `ielpar` to that element index. Returns number assigned.
        """
        assigned = 0
        npar = self.particles.nparticle
        nelem = self.hgrid.ne

        # Precompute centroids
        xctr = self.hgrid.xctr
        yctr = self.hgrid.yctr

        for ip in range(npar):
            iel = int(self.particles.ielpar[ip])
            # Treat non-positive indices as unassigned (particle.bp default is 0)
            if 0 < iel < nelem:
                continue

            x = self.particles.xpar[ip]
            y = self.particles.ypar[ip]

            # find k nearest centroids
            d2 = (xctr - x) ** 2 + (yctr - y) ** 2
            ksel = min(k, nelem)
            idxs = np.argpartition(d2, ksel - 1)[:ksel]

            found = False
            for ie in idxs:
                nvert = int(self.hgrid.i34[ie])
                nodes = self.hgrid.elnode[:nvert, ie]
                xs = self.hgrid.x[nodes]
                ys = self.hgrid.y[nodes]
                try:
                    if pt_in_poly2(nvert, xs, ys, x, y):
                        arco, nodel = pt_in_poly3(nvert, xs, ys, x, y)
                        self.particles.ielpar[ip] = int(ie)
                        assigned += 1
                        found = True
                        logger.debug(f"Assigned particle {ip} to element {ie}")
                        break
                except Exception:
                    continue

            if not found:
                self.particles.ielpar[ip] = -1
                logger.debug(f"Could not assign element for particle {ip}")

        logger.info(f"assign_initial_elements: assigned {assigned} particles")
        return assigned
    
    def _interpolate_velocity(self, ipar: int, nnel: int, jlev: int,
                             nodel: np.ndarray, arco: np.ndarray):
        """Interpolate velocity at particle location."""
        upar = 0.0
        vpar = 0.0
        wpar = 0.0
        
        lev_idx = max(0, min(jlev - 1, self.nvrt -1))
        # gather node indices for interpolation (convert local 1-based to 0-based)
        node_indices = np.array([self.hgrid.elnode[nodel[j] - 1, nnel] for j in range(3)], dtype=int)

        for j in range(3):
            nd = node_indices[j]
            upar += self.uu2[nd, lev_idx] * arco[j]
            vpar += self.vv2[nd, lev_idx] * arco[j]
            wpar += self.ww2[nd, lev_idx] * arco[j]
        print('upar',upar)
        print('lev_idx',lev_idx)

        # If interpolated velocity is effectively zero, attempt fallbacks:
        # 1) search vertically within the element nodes for a level with non-zero values
        # 2) if still zero, search other nodes of the element and pick the level
        is_zero = lambda x: np.all(np.isclose(x, 0.0, atol=1e-12))
        if is_zero(upar) and is_zero(vpar) and is_zero(wpar):
            try:
                # full element nodes (may be 3 or 4)
                nodes_ie = self.hgrid.elnode[: self.hgrid.i34[nnel], nnel]

                # compute level-wise mean absolute velocity over element nodes
                abs_mean_u = np.nanmean(np.abs(self.uu2[nodes_ie, :]), axis=0)
                abs_mean_v = np.nanmean(np.abs(self.vv2[nodes_ie, :]), axis=0)
                abs_mean = abs_mean_u + abs_mean_v

                # pick best level where combined mean is largest
                k_best = int(np.nanargmax(abs_mean))
                if abs_mean[k_best] > 1e-8:
                    lev_idx = k_best
                    upar = 0.0
                    vpar = 0.0
                    wpar = 0.0
                    for j in range(3):
                        nd = node_indices[j]
                        upar += self.uu2[nd, lev_idx] * arco[j]
                        vpar += self.vv2[nd, lev_idx] * arco[j]
                        wpar += self.ww2[nd, lev_idx] * arco[j]
                else:
                    # last resort: element-averaged velocity over all levels
                    u_avg = np.nanmean(self.uu2[nodes_ie, :])
                    v_avg = np.nanmean(self.vv2[nodes_ie, :])
                    w_avg = np.nanmean(self.ww2[nodes_ie, :])
                    if not (np.isfinite(u_avg) and np.isfinite(v_avg) and np.isfinite(w_avg)):
                        u_avg = 0.0
                        v_avg = 0.0
                        w_avg = 0.0
                    upar = u_avg
                    vpar = v_avg
                    wpar = w_avg
            except Exception as e:
                pass
        
        self.particles.upar[ipar] = upar
        self.particles.vpar[ipar] = vpar
        self.particles.wpar[ipar] = wpar
        print('upar',upar)
        print('vpar',vpar)
    
    def save_output(self, outputfile: str, times: list, trajectories: list):
        """
        Save particle trajectories to netCDF file.
        
        Parameters
        ----------
        outputfile : str
            Output file path
        times : list
            Time values
        trajectories : list
            List of trajectory dicts
        """
        ntimes = len(times)
        nparticles = self.particles.nparticle
        
        with nc.Dataset(outputfile, 'w', format='NETCDF4') as ds:
            ds.createDimension('time', ntimes)
            ds.createDimension('particle', nparticles)
            
            # Time variable
            time_var = ds.createVariable('time', 'f8', ('time',))
            time_var.units = 'seconds'
            time_var[:] = times
            
            # Position variables
            x_var = ds.createVariable('lon', 'f4', ('particle', 'time'))
            y_var = ds.createVariable('lat', 'f4', ('particle', 'time'))
            z_var = ds.createVariable('depth', 'f4', ('particle', 'time'))
            
            for t, traj in enumerate(trajectories):
                x = traj['x'].copy()
                y = traj['y'].copy()
                
                # If tracking in lon/lat (ics==2), convert back from projected meters to lon/lat
                if self.ics == 2:
                    try:
                        xlon_rad = np.zeros_like(x)
                        ylat_rad = np.zeros_like(y)
                        for i in range(len(x)):
                            xlon_rad[i], ylat_rad[i] = cppinverse(x[i], y[i], self.slam0, self.sfea0)
                        # Convert radians to degrees
                        x = np.rad2deg(xlon_rad)
                        y = np.rad2deg(ylat_rad)
                    except Exception as e:
                        logger.warning(f"Failed to convert particle positions back to lon/lat: {e}")
                
                x_var[:, t] = x
                y_var[:, t] = y
                z_var[:, t] = traj['z']
            
            logger.info(f"Trajectories saved to {outputfile}")
