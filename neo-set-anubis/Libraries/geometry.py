import numpy as np
from constants import *
from helpers import *

## Functions to determine ANUBIS acceptance
        
# Check if particle decays outside ATLAS, within cavern
def in_cavern_edit(r, eta, phi):
    decay_pos = list(to_cartesian(r, eta, phi))

    # Shift by ATLAS offset
#    decay_pos[0] += origin[0]

    # Check if within ATLAS
    if np.sqrt(np.square(decay_pos[0] - origin[0]) + np.square(decay_pos[1])) < atlas_t_max and \
       atlas_z_min < decay_pos[2] < atlas_z_max:
        return False
    
    # Check if within cavern rectangle bounds
    if np.abs(decay_pos[0]) < cavern_x_max and \
       np.abs(decay_pos[1]) < cavern_y_max and \
       np.abs(decay_pos[2]) < cavern_z_max:
        return True
    
    # Check if within curved part of cavern
    if np.abs(decay_pos[0]) < cavern_x_max and \
       np.abs(decay_pos[2]) < cavern_z_max and \
       np.sqrt(np.square(decay_pos[0]) +
               np.square(decay_pos[1] -
                         (cavern_sagitta + cavern_y_max -
                          cavern_r_curvature))) < cavern_r_curvature:
        # Veto if in-shaft and passed straight tracking stations
        if np.sqrt(np.square(decay_pos[0]) +
                   np.square(decay_pos[2] - (anubis_z_min + anubis_radius))) < anubis_radius and \
           decay_pos[1] > anubis_ts_height[0]:
            return False
        if np.sqrt(np.square(decay_pos[0]) +
                   np.square(decay_pos[2] - (px16_z_min - px16_radius))) < px16_radius and \
           decay_pos[1] > anubis_ts_height[0]:
            return False

        return True

    # Must be out of cavern
    return False

# Check if particle decays outside ATLAS, within cavern
def in_cavern(r, eta, phi):
    decay_pos = list(to_cartesian(r, eta, phi))

    # Shift by ATLAS offset
    decay_pos[0] += origin[0]

    # Check if within ATLAS
    if np.sqrt(np.square(decay_pos[0] - origin[0]) + np.square(decay_pos[1])) < atlas_t_max and \
       atlas_z_min < decay_pos[2] < atlas_z_max:
        return False

    # Check if within cavern rectangle bounds
    if np.abs(decay_pos[0]) < cavern_x_max and \
       np.abs(decay_pos[1]) < cavern_y_max and \
       np.abs(decay_pos[2]) < cavern_z_max:
        return True

    # Check if within curved part of cavern
    if np.abs(decay_pos[0]) < cavern_x_max and \
       np.abs(decay_pos[2]) < cavern_z_max and \
       np.sqrt(np.square(decay_pos[0]) +
               np.square(decay_pos[1] -
                         (cavern_sagitta + cavern_y_max -
                          cavern_r_curvature))) < cavern_r_curvature:
        # Veto if in-shaft and passed straight tracking stations
        if np.sqrt(np.square(decay_pos[0]) +
                   np.square(decay_pos[2] - (anubis_z_min + anubis_radius))) < anubis_radius and \
           decay_pos[1] > anubis_ts_height[0]:
            return False
        if np.sqrt(np.square(decay_pos[0]) +
                   np.square(decay_pos[2] - (px16_z_min - px16_radius))) < px16_radius and \
           decay_pos[1] > anubis_ts_height[0]:
            return False

        return True

    # Must be out of cavern
    return False

# Check if particle decays within PX14 shaft (ANUBIS)
def in_shaft(r, eta, phi):
    decay_pos = list(to_cartesian(r, eta, phi))

    # Shift by ATLAS offset
    decay_pos[0] += origin[0]

    # Check height
    if decay_pos[1] < anubis_ts_height[0] or \
       decay_pos[1] > anubis_ts_height[-1]:
        return False

    # Check lateral
    if np.sqrt(np.square(decay_pos[0]) +
               np.square(decay_pos[2] - (anubis_z_min + anubis_radius))) > anubis_radius:
        return False

    # Check pipe cut-off
    if decay_pos[0] < anubis_x_min:
        return False

    # Must be in-shaft
    return True

# Compute point of intersection of jet particle in local coordinates of ceiling station
def ceil_intersection(eta, phi, decay_pos):
    # Compute y coordinate of curvature center
    y_curve_cent = cavern_y_max + cavern_sagitta - cavern_r_curvature

    # Compute slope of jet particle
    m = np.tan(phi)

    # Pre-compute often used quantity
    cydydxm = y_curve_cent - decay_pos[1] + (decay_pos[0] * m)

    # Solve for intersection between curve of cavern and line of jet particle
    # Take correct intersection
    if phi < 0 or phi > np. pi:
        x = ((np.square(m) * cydydxm) - \
             np.sqrt(np.square(m) * (np.square(cavern_r_curvature) * (1 + np.square(m)) -
                                     np.square(cydydxm)))) / \
            (m + np.power(m, 3))

        y = ((decay_pos[1] - (decay_pos[0] * m) + (y_curve_cent * np.square(m))) - \
             np.sqrt(np.square(m) * (np.square(cavern_r_curvature) * (1 + np.square(m)) -
                                     np.square(cydydxm)))) / \
            (1 + np.square(m))
    else:
        x = ((np.square(m) * cydydxm) + \
             np.sqrt(np.square(m) * (np.square(cavern_r_curvature) * (1 + np.square(m)) -
                                     np.square(cydydxm)))) / \
             (m + np.power(m, 3))

        y = ((decay_pos[1] - (decay_pos[0] * m) + (y_curve_cent * np.square(m))) + \
             np.sqrt(np.square(m) * (np.square(cavern_r_curvature) * (1 + np.square(m)) -
                                     np.square(cydydxm)))) / \
            (1 + np.square(m))

    # Find distance between intersection point and vertical center
    disp = np.sqrt(np.square(x) + np.square(y - (cavern_y_max + cavern_sagitta)))

    # Find angle between cavern ceiling and intersection point, with correct sign
    alpha = np.arccos((2 * np.square(cavern_r_curvature) - np.square(disp)) /
                      (2 * np.square(cavern_r_curvature))) * (-1 if x < 0 else 1)

    # Get arc length given angle, reverse coordinates
    arc_length = cavern_r_curvature * alpha

    # Get ceiling height at intersection point
    ceil_height = np.sqrt(np.square(cavern_r_curvature) - np.square(x)) + y_curve_cent

    # Find z distance covered to travel to cavern ceiling
    z = (ceil_height - decay_pos[1]) / np.tan(to_theta(eta))
    z += decay_pos[2]

    return arc_length, y, z

# Check if jet particle intersects PX14 tracking station
def ts_stats(eta, phi, decay_pos):
    num_hits = np.full(len(anubis_ts_height), 0)

    # Veto if particle passes through concrete
    d_to_px14 = anubis_ts_height[0] - decay_pos[1]
    if d_to_px14 > 0:
        # Get jet position at first TS
        r = d_to_px14 / (np.sin(to_theta(eta)) *
                         np.sin(phi))
        x, _, z = to_cartesian(r, eta, phi)
        x += decay_pos[0]
        z += decay_pos[2]

        # Check if outside first TS
        if np.sqrt(np.square(x) +
                   np.square(z -
                             (anubis_z_min + anubis_radius))) > anubis_radius or \
           x < anubis_x_min:
            # Return no hits
            return num_hits

    # Iterate over tracking stations
    for i, ts_height in enumerate(anubis_ts_height):
        # Find y difference
        dy = ts_height - decay_pos[1]

        # Ensure particle has not already passed TS
        if dy < 0:
            continue

        # Find jet position at TS
        r = dy / (np.sin(to_theta(eta)) *
                  np.sin(phi))
        x, _, z = to_cartesian(r, eta, phi)
        x += decay_pos[0]
        z += decay_pos[2]

        # Veto backwards-going particles
        if phi < 0 and phi > -np.pi:
            continue

        # Check if decay is within TS bounds
        if np.sqrt(np.square(x) +
                   np.square(z -
                             (anubis_z_min + anubis_radius))) < anubis_radius and \
           x > anubis_x_min:
            num_hits[i] += 1

    return num_hits

# Check if jet particle intersects ceiling tracking station
def ceiling_stat(eta, phi, decay_pos):
    # Ensure decay occurs within cavern
    llp_r, llp_eta, llp_phi = to_spherical(decay_pos[0] - origin[0],
                                           decay_pos[1],
                                           decay_pos[2])
    if not in_cavern(llp_r, llp_eta, llp_phi):
        return False

    # Find local coordinates of intersection with ceiling
    jet_x, jet_y, jet_z = ceil_intersection(eta, phi, decay_pos)

    # Check if intersection is below cavern
    if jet_y < 0:
        return False

    # Check if intersection is within z-bounds
    if jet_z > cavern_z_max or jet_z < -cavern_z_max:
        return False

    # Check if intersection is wihtin x-bounds

    # Get arc length of ceiling
    ceil_angle = np.arctan2(cavern_x_max, cavern_r_curvature - cavern_sagitta)
    ceil_length = cavern_r_curvature * ceil_angle

    if jet_x > ceil_length or jet_x < -ceil_length:
        return False

    # Must be wihtin the ceiling bounds
    return True
