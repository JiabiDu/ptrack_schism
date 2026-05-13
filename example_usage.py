#!/usr/bin/env python
"""
Example usage of the ptrack3 Python particle tracking module.

This demonstrates how to:
1. Load grid and vertical coordinate data
2. Load particle initial conditions
3. Track particles through time using netCDF hydrodynamic outputs
4. Save results to netCDF
"""

import numpy as np
from ptrack3 import HGrid, VGrid, ParticleData, ParticleTracker
import logging
import os

# Configure logging; allow override with PYTHONLOGGING environment variable
log_level = os.getenv('PYTHONLOGGING', '').upper()
if log_level:
    try:
        level = getattr(logging, log_level)
    except AttributeError:
        try:
            level = int(log_level)
        except ValueError:
            level = logging.INFO
else:
    level = logging.INFO

logging.basicConfig(level=level)
logger = logging.getLogger(__name__)


def main():
    """Main particle tracking workflow."""
    
    # Configuration parameters
    hydroDir = '/scratch/user/jdu/GoM/run19q3/'  # Directory containing hgrid, vgrid, and particle files
    h0 = 0.1  # Minimum water depth (m)
    ibf = 1  # Forward tracking (1) or backward (-1)
    istiff = 0  # Stiff coordinate flag (0=moving free surface, 1=fixed below free surface)
    dt =  7200.0  # Output time step (seconds) from hydrodynamic model
    ndeltp = 30  # Sub-steps for particle tracking
    rnday = 2.0  # Total simulation days
    ics = 2  # Coordinate system (1=Cartesian, 2=lon/lat)
    slam0 = -95  # CPP projection center longitude
    sfea0 = 29.5  # CPP projection center latitude
    
    # Optional tracers
    temp_on = False  # Track temperature
    salt_on = False  # Track salinity
    solar_on = False  # Track solar radiation
    settling_velocity = 0.0  # m/day
    
    # File paths
    hgrid_file = f'{hydroDir}/hgrid.gr3'
    vgrid_file = f'{hydroDir}/vgrid.in'
    particle_file = 'particle.bp'
    output_file = 'particle_tracks.nc'
    
    # =========================================================================
    # STEP 1: Load Grid Data
    # =========================================================================
    logger.info("Loading grid data...")
    
    hgrid = HGrid.from_file(hgrid_file)
    logger.info(f"  Horizontal grid: {hgrid.np} nodes, {hgrid.ne} elements")
    
    vgrid = VGrid.from_file(vgrid_file)
    logger.info(f"  Vertical grid: {vgrid.nvrt} layers, ivcor={vgrid.ivcor}")
    
    # =========================================================================
    # STEP 2: Load Particle Initial Conditions
    # =========================================================================
    logger.info("Loading particles...")
    
    optional_props = {
        'temp': temp_on,
        'salt': salt_on,
        'solar': solar_on
    }
    
    particles = ParticleData.from_file(particle_file, optional_props)
    logger.info(f"  {particles.nparticle} particles loaded")
    
    # =========================================================================
    # STEP 3: Initialize Particle Tracker
    # =========================================================================
    logger.info("Initializing particle tracker...")
    
    tracker = ParticleTracker(
        hgrid=hgrid,
        vgrid=vgrid,
        particles=particles,
        h0=h0,
        ics=ics,
        slam0=slam0,
        sfea0=sfea0,
        ibf=ibf,
        istiff=istiff,
        dt=dt,
        ndeltp=ndeltp,
        temp_on=temp_on,
        salt_on=salt_on,
        solar_on=solar_on,
        settling_velocity=settling_velocity,
        force_surface=True
    )
    
    # Quick containment debug for the first particle to help locate element/search issues
    try:
        logger.info("Running containment debug for first particle...")
        tracker.debug_find_element_for_point(0, max_checks=5000)
    except Exception as e:
        logger.warning(f"Containment debug failed: {e}")

    # Assign initial elements for particles if not provided in particle.bp
    try:
        n_assigned = tracker.assign_initial_elements(k=50)
        logger.info(f"Assigned initial elements to {n_assigned} particles")
    except Exception as e:
        logger.warning(f"assign_initial_elements failed: {e}")

    # Update z-coordinates and wet/dry flags after assigning elements
    try:
        tracker.compute_z_coords()
    except Exception as e:
        logger.warning(f"compute_z_coords after assignment failed: {e}")

    # =========================================================================
    # STEP 4: Run Particle Tracking
    # =========================================================================
    logger.info("Tracking particles through time...")
    
    # Prepare time steps and corresponding netCDF files
    ntime = int(rnday * 86400 / dt)
    times = np.arange(ntime) * dt

    # Hydrodynamic output cadence and records-per-file
    hydro_dt = 7200          # hydrodynamic output interval in seconds (2 hours)
    recs_per_file = 12       # number of records stored per file (e.g., 1 day per file)
    seconds_per_file = hydro_dt * recs_per_file

    # logging every few steps so we can see progress; set to 0 to disable
    # If you are just testing the conversion you can limit the number of steps.
    # Set test_mode=False to run the full simulation.
    test_mode = True
    max_steps = 864
    if test_mode and len(times) > max_steps:
        times = times[:max_steps]

    # Build list of hydro snapshots (file, record) for the (possibly truncated) times
    snaps_needed = sorted(set(int(t // hydro_dt) for t in times))
    snapshot_files = []
    for s in snaps_needed:
        t_snap = s * hydro_dt
        file_idx = int(t_snap // seconds_per_file) + 1
        record_in_file = int((t_snap % seconds_per_file) // hydro_dt)
        fname = f"{hydroDir}/outputs/out2d_{file_idx}.nc"
        snapshot_files.append((fname, record_in_file))

    # Run tracking with interpolation enabled
    trajectories, times = tracker.track_particles(times, snapshot_files,
                                                 hydro_dt=hydro_dt,
                                                 log_interval=10,
                                                 preload=True)
    
    # =========================================================================
    # STEP 5: Save Results
    # =========================================================================
    logger.info(f"Saving results to {output_file}...")
    
    tracker.save_output(output_file, times, trajectories)
    
    logger.info("Completed successfully!")


if __name__ == '__main__':
    main()
