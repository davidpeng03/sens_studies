import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from geometry import *
from constants import *
from helpers import *

# Plot histogram
def plot_histogram(ax, var, bins, min_range, max_range, col, e_col, xlabel, title, flow):
    # Bin histogram
    y, binBounds = np.histogram(var, bins=bins, range=(min_range, max_range))
    err = np.sqrt(y)
    binCents = (binBounds[1:] + binBounds[:-1]) / 2
    binWidth = binBounds[1] - binBounds[0]
    ax.step(binCents, y, where='mid', linewidth=2, color=col)
    ax.errorbar(binCents, y, yerr=err, ls='none', linewidth=2, color=e_col, zorder=100)

    # Under/overflow
    if flow:
        # Calculate underflow and overflow
        under = np.sum(np.array(var) < min_range)
        over  = np.sum(np.array(var) > max_range)
        ax.text(0.98, 0.95, 'Underflow = ' + str(under),
                ha='right', va='top', transform = ax.transAxes, fontsize='large')
        ax.text(0.98, 0.9, 'Overflow = ' + str(over),
                ha='right', va='top', transform = ax.transAxes, fontsize='large')

    # Style axes
    ax.set_xlabel(xlabel, fontsize='large')
    ax.tick_params(axis='x', labelsize='large')

    ax.set_ylabel('Count', fontsize='large')
    ax.tick_params(axis='y', labelsize='large')

    ax.set_title(title, fontsize='x-large')

# Plot cutflow
def plot_cutflow(ax, cut_titles, cut_values, c, label, title):
    # Plot and format
    ax.step(cut_titles, cut_values, where='mid', linewidth=2, color=c, label=label)

    # Style plot
    ax.tick_params('x', labelrotation=-80, labelsize='x-large')
    ax.tick_params('y', labelsize='x-large')
    ax.set_yscale('log')
    ax.set_ylabel('Fraction of Events', fontsize='x-large')
    ax.set_title(title, fontsize='x-large')
    ax.legend()

# Plot ATLAS and ANUBIS
def plot_detectors(ax1, ax2):
    ## YZ projection

    # Set plot limits and labels
    ax1[0].set_xlim((-80, 80))
    ax1[0].set_xlabel('z (m)', fontsize='large')

    ax1[0].set_ylim((-20, 100))
    ax1[0].set_ylabel('y (m)', fontsize='large')

    # Plot ATLAS bounds

    # Plot vertexing limit
    ax1[0].plot([atlas_z_min, atlas_z_max], [atlas_t_max, atlas_t_max], c=colors[1], ls='dotted')
    ax1[0].plot([atlas_z_min, atlas_z_max], [-atlas_t_max, -atlas_t_max], c=colors[1], ls='dotted')
    ax1[0].plot([atlas_z_min, atlas_z_min], [-atlas_t_max, atlas_t_max], c=colors[1], ls='dotted')
    ax1[0].plot([atlas_z_max, atlas_z_max], [-atlas_t_max, atlas_t_max], c=colors[1], ls='dotted')

    # Plot physical limit
    ax1[0].plot([atlas_z_min, atlas_z_max], [atlas_full_max, atlas_full_max], c=colors[0])
    ax1[0].plot([atlas_z_min, atlas_z_max], [-atlas_full_max, -atlas_full_max], c=colors[0])
    ax1[0].plot([atlas_z_min, atlas_z_min], [-atlas_full_max, atlas_full_max], c=colors[0])
    ax1[0].plot([atlas_z_max, atlas_z_max], [-atlas_full_max, atlas_full_max], c=colors[0])
    ax1[0].text(0, atlas_full_max + 1, 'ATLAS', ha='center', va='bottom')

    # Plot cavern bounds
    ax1[0].plot([-cavern_z_max, cavern_z_max],
                [cavern_y_max + cavern_sagitta, cavern_y_max + cavern_sagitta],
                c=colors[1])
    ax1[0].text(px16_z_min - (2 * px16_radius) - 1, cavern_y_max + cavern_sagitta + 1,
                'CAVERN CEILING TS', ha='right', va='bottom')
    ax1[0].plot([-cavern_z_max, cavern_z_max], [-cavern_y_max, -cavern_y_max], c=colors[2], ls='dotted')
    ax1[0].plot([-cavern_z_max, -cavern_z_max], [-cavern_y_max, cavern_y_max + cavern_sagitta],
                c=colors[2], ls='dotted')
    ax1[0].plot([cavern_z_max, cavern_z_max], [-cavern_y_max, cavern_y_max + cavern_sagitta],
                c=colors[2], ls='dotted')
    ax1[0].text(-cavern_z_max + 1, cavern_y_max + cavern_sagitta - 1, 'CAVERN', ha='left', va='top')

    # Plot IP
    ax1[0].scatter([0], [0], c=colors[1])
    ax1[0].text(1, -1, 'IP', ha='left', va='top')

    # Plot PX14
    y_max = 100

    ax1[0].plot([anubis_z_min, anubis_z_min], [cavern_y_max + cavern_sagitta, y_max], c=colors[2], ls='dotted')
    ax1[0].plot([anubis_z_min + (2 * anubis_radius), anubis_z_min + (2 * anubis_radius)],
                [cavern_y_max + cavern_sagitta, y_max], c=colors[2], ls='dotted')
    ax1[0].text(anubis_z_min + 1, y_max - 1, 'PX14', ha='left', va='top')

    # Plot PX16
    ax1[0].plot([px16_z_min, px16_z_min], [cavern_y_max + cavern_sagitta, y_max], c=colors[2], ls='dotted')
    ax1[0].plot([px16_z_min - (2 * px16_radius), px16_z_min - (2 * px16_radius)],
                [cavern_y_max + cavern_sagitta, y_max], c=colors[2], ls='dotted')
    ax1[0].text(px16_z_min - (2 * px16_radius) + 1, y_max - 1, 'PX16', ha='left', va='top')

    # Plot ANUBIS
    for i in range(1, 5):
        ax1[0].plot([anubis_z_min, anubis_z_min + (2 * anubis_radius)],
                    [anubis_ts_height[i - 1], anubis_ts_height[i - 1]], c=colors[1])
        ax1[0].text(anubis_z_min + 1, anubis_ts_height[i - 1] + 1,
                    'TS ' + str(i), ha='left', va='bottom')

    ## Plot YZ

    # Set plot limits and labels
    ax1[1].set_xlim((-80, 80))
    ax1[1].set_xlabel('x (m)', fontsize='large')

    ax1[1].set_ylim((-20, 100))
    ax1[1].set_ylabel('y (m)', fontsize='large')

    # Plot ATLAS bounds

    # Plot vertexing limit
    xs = np.linspace(-atlas_t_max, atlas_t_max, num=100)
    ax1[1].plot(xs + origin[0], np.sqrt(np.square(atlas_t_max) - np.square(xs)), c=colors[1], ls='dotted')
    ax1[1].plot(xs + origin[0], -np.sqrt(np.square(atlas_t_max) - np.square(xs)), c=colors[1], ls='dotted')

    # Plot physical limit
    xs = np.linspace(-atlas_full_max, atlas_full_max, num=100)
    ax1[1].plot(xs + origin[0], np.sqrt(np.square(atlas_full_max) - np.square(xs)), c=colors[0])
    ax1[1].plot(xs + origin[0], -np.sqrt(np.square(atlas_full_max) - np.square(xs)), c=colors[0])
    ax1[1].text(origin[0], atlas_full_max + 1, 'ATLAS', ha='center', va='bottom')

    # Plot cavern bounds

    # Plot rectangular
    ax1[1].plot([-cavern_x_max, cavern_x_max], [-cavern_y_max, -cavern_y_max], c=colors[2], ls='dotted')
    ax1[1].plot([-cavern_x_max, -cavern_x_max], [-cavern_y_max, cavern_y_max], c=colors[2], ls='dotted')
    ax1[1].plot([cavern_x_max, cavern_x_max], [-cavern_y_max, cavern_y_max], c=colors[2], ls='dotted')
    ax1[1].text(-cavern_x_max - 1, cavern_y_max - 1, 'CAVERN', ha='right', va='top')

    # Plot curved
    xs_left = np.linspace(-cavern_x_max, -anubis_radius, num=100)
    ax1[1].plot(xs_left, np.sqrt(np.square(cavern_r_curvature) - np.square(xs_left)) +
                         cavern_y_max -
                         np.sqrt(np.square(cavern_r_curvature) - np.square(cavern_x_max)),
                c=colors[1])
    xs_right = np.linspace(anubis_radius, cavern_x_max, num=100)
    ax1[1].plot(xs_right, np.sqrt(np.square(cavern_r_curvature) - np.square(xs_right)) +
                          cavern_y_max -
                          np.sqrt(np.square(cavern_r_curvature) - np.square(cavern_x_max)),
                c=colors[1])
    ax1[1].text(anubis_radius + 5, anubis_ts_height[0] - 2, 'CAVERN CEILING TS', ha='left', va='top')

    # Plot IP
    ax1[1].scatter([origin[0]], [0], c=colors[1])
    ax1[1].text(origin[0] + 1, -1, 'IP', ha='left', va='top')

    ax1[1].plot([-anubis_radius, -anubis_radius], [px14_height, y_max], c=colors[2], ls='dotted')
    ax1[1].plot([anubis_radius, anubis_radius], [px14_height, y_max], c=colors[2], ls='dotted')
    ax1[1].text(-anubis_radius + 1, y_max - 1, 'PX14', ha='left', va='top')

    # Plot ANUBIS
    for i in range(1, 5):
        ax1[1].plot([-anubis_radius, anubis_radius],
                    [anubis_ts_height[i - 1], anubis_ts_height[i - 1]], c=colors[1])
        ax1[1].text(-anubis_radius + 1, anubis_ts_height[i - 1] + 1,
                    'TS ' + str(i), ha='left', va='bottom')

    ## Plot XZ four times

    ax2[0].set_ylabel('z (m)', fontsize='large')
    for i in range(1, len(anubis_ts_height) + 1):
        # Set plot limits and labels
        ax2[i].set_xlim((-30, 30))
        ax2[i].set_xlabel('x (m)', fontsize='large')

        ax2[i].set_ylim((-45, 45))

        ax2[i].set_title('Tracking Station %i' % (i), fontsize='large')

        # Plot ATLAS bounds

        # Plot physical limit only
        ax2[i].plot([-atlas_full_max + origin[0], -atlas_full_max + origin[0]],
                    [atlas_z_min, atlas_z_max], c=colors[0])
        ax2[i].plot([-atlas_full_max + origin[0], atlas_full_max + origin[0]],
                    [atlas_z_max, atlas_z_max], c=colors[0])
        ax2[i].plot([atlas_full_max + origin[0], atlas_full_max + origin[0]],
                    [atlas_z_max, atlas_z_min], c=colors[0])
        ax2[i].plot([atlas_full_max + origin[0], -atlas_full_max + origin[0]],
                    [atlas_z_min, atlas_z_min], c=colors[0])
        ax2[i].text(-atlas_full_max + origin[0] + 1, -atlas_z_max + 1, 'ATLAS', ha='left', va='bottom')

        # Plot cavern bounds
        ax2[i].plot([-cavern_x_max, -cavern_x_max], [-cavern_z_max, cavern_z_max], c=colors[2], ls='dotted')
        ax2[i].plot([cavern_x_max, cavern_x_max], [-cavern_z_max, cavern_z_max], c=colors[2], ls='dotted')
        ax2[i].plot([-cavern_x_max, cavern_x_max], [-cavern_z_max, -cavern_z_max], c=colors[2], ls='dotted')
        ax2[i].plot([-cavern_x_max, cavern_x_max], [cavern_z_max, cavern_z_max], c=colors[2], ls='dotted')
        ax2[i].text(-cavern_x_max + 1, -cavern_z_max - 1, 'CAVERN', ha='left', va='top')

        # Plot PX16
        xs = np.linspace(-px16_radius, px16_radius, num=100)
        ax2[i].plot(xs, np.sqrt(np.square(px16_radius) - np.square(xs)) + px16_z_min - px16_radius,
                    c=colors[2], ls='dotted')
        ax2[i].plot(xs, -np.sqrt(np.square(px16_radius) - np.square(xs)) + px16_z_min - px16_radius,
                    c=colors[2], ls='dotted')
        ax2[i].text(0, px16_z_min - px16_radius, 'PX16', ha='center', va='center')

        # Plot IP
        ax2[i].scatter([origin[0]], [0], c=colors[1])
        ax2[i].text(1 + origin[0], -1, 'IP', ha='left', va='top')

        xs = np.linspace(anubis_x_min, anubis_radius, num=100)
        zs1 = np.sqrt(np.square(anubis_radius) - np.square(xs)) + anubis_z_min + anubis_radius
        zs2 = -np.sqrt(np.square(anubis_radius) - np.square(xs)) + anubis_z_min + anubis_radius
        ax2[i].plot(xs, zs1, c=colors[0])
        ax2[i].plot(xs, zs2, c=colors[0])
        ax2[i].plot([anubis_x_min, anubis_x_min], [zs1[0], zs2[0]], c=colors[0])
        ax2[i].text(0, anubis_z_min + (2 * anubis_radius) + 1, 'PX14', ha='center', va='bottom')

    ## Initialize ceiling unroll

    # Set plot limits and labels
    ax2[0].set_xlim((-30, 30))
    ax2[0].set_xlabel('Distance along ceiling curve (m)', fontsize='large')
    ax2[0].set_ylim((-45, 45))
    ax2[0].set_title('Ceiling Station Unrolled', fontsize='large')

    # Plot IP
    ax2[0].scatter([ip_x_unrolled], [0], c=colors[1])
    ax2[0].text(ip_x_unrolled + 1, -1, 'IP', ha='left', va='top')

    # Get arc length of ceiling
    ceil_angle_tot = np.pi - 2 * (np.arctan2(cavern_r_curvature - cavern_sagitta, cavern_x_max))
    ceil_length = cavern_r_curvature * ceil_angle_tot / 2

    # Plot ceiling bounds
    ax2[0].plot([-ceil_length, -ceil_length], [-cavern_z_max, cavern_z_max], c=colors[2], ls='dotted')
    ax2[0].plot([ceil_length, ceil_length], [-cavern_z_max, cavern_z_max], c=colors[2], ls='dotted')
    ax2[0].plot([-ceil_length, ceil_length], [-cavern_z_max, -cavern_z_max], c=colors[2], ls='dotted')
    ax2[0].plot([-ceil_length, ceil_length], [cavern_z_max, cavern_z_max], c=colors[2], ls='dotted')
    ax2[0].text(-ceil_length + 1, -cavern_z_max - 1, 'CAVERN CEILING', ha='left', va='top')

    # Plot ovals for shafts

    # PX16
    xs = np.linspace(-px16_radius, px16_radius, num=100)
    zs1 = np.sqrt(np.square(px16_radius) - np.square(xs)) + px16_z_min - px16_radius
    zs2 = -np.sqrt(np.square(px16_radius) - np.square(xs)) + px16_z_min - px16_radius

    # Skew xs
    height = anubis_ts_height[0] - (cavern_y_max + cavern_sagitta - cavern_r_curvature)
    xs = np.arctan2(xs, height) * cavern_r_curvature

    ax2[0].plot(xs, zs1,
                c=colors[2], ls='dotted')
    ax2[0].plot(xs, zs2,
                c=colors[2], ls='dotted')
    ax2[0].text(0, px16_z_min - px16_radius, 'PX16', ha='center', va='center')

    # PX14
    xs = np.linspace(anubis_x_min, anubis_radius, num=100)
    zs1 = np.sqrt(np.square(anubis_radius) - np.square(xs)) + anubis_z_min + anubis_radius
    zs2 = -np.sqrt(np.square(anubis_radius) - np.square(xs)) + anubis_z_min + anubis_radius

    # Skew xs
    xs = np.arctan2(xs, height) * cavern_r_curvature
    x_min_skewed = np.arctan2(anubis_x_min, height) * cavern_r_curvature

    ax2[0].plot(xs, zs1, c=colors[2], ls='dotted')
    ax2[0].plot(xs, zs2, c=colors[2], ls='dotted')
    ax2[0].plot([x_min_skewed, x_min_skewed], [zs1[0], zs2[0]], c=colors[2], ls='dotted')
    ax2[0].text(0, anubis_z_min + anubis_radius, 'PX14', ha='center', va='center')

    # Add custom legend
    leg_els = [Line2D([0], [0], marker='x', color=colors[3], lw=0, markersize=7,
                      markeredgewidth=2, label='Observed'),
               Line2D([0], [0], marker='x', color=colors[4], lw=0, markersize=7,
                      markeredgewidth=2, label='Not observed')]
    ax2[0].legend(handles=leg_els, loc='upper left', fontsize='medium')

# Plot LLP event
def plot_llp(ax1, ax2, decay_pos, ceil_stat, px14_stats, jets, title=""):
    # Set title
    ax1[0].set_title(title)

    ax1[0].plot([0, decay_pos[2]], [0, decay_pos[1]], c=colors[2], linestyle='dashed')
    ax1[1].plot([origin[0], decay_pos[0]], [0, decay_pos[1]], c=colors[2], linestyle='dashed')

    # Add legend
    leg_lines = [Line2D([0], [0], color=colors[2], lw=2, linestyle='dashed'),
                 Line2D([0], [0], color=colors[3], lw=2)]
    ax1[0].legend(leg_lines, ['LLP', 'Jet Particle'], loc="upper left")

    # Display decay statistics on plot
    ax1[1].text(-75, 95, 'Number of Jet Particles Observed', ha='left', va='top')
    for i, stat in enumerate(px14_stats):
        ax1[1].text(-75, 90 - (5 * (i + 1)), 'TS %i: %i' % (i + 1, stat), ha='left', va='top')
    ax1[1].text(-75, 90, 'Ceiling: %i' % ceil_stat, ha='left', va='top')

    # Plot jet particles in cross-section
    for eta, phi, _ in jets:

        # Get position sufficiently far away
        jet_x, jet_y, jet_z = to_cartesian(100000, eta, phi)

        # Add offset
        jet_x += decay_pos[0]
        jet_y += decay_pos[1]
        jet_z += decay_pos[2]

        # Plot particle
        ax1[0].plot([decay_pos[2], jet_z], [decay_pos[1], jet_y], c=colors[3])
        ax1[1].plot([decay_pos[0], jet_x], [decay_pos[1], jet_y], c=colors[3])

    # Plot jet particles in bird's eye
    for i, ts_height in enumerate(anubis_ts_height):

        # Get displacement to tracking station
        dy = ts_height - decay_pos[1]

        if dy < 0 or decay_pos[1] < 0:
            continue

        # Iterate over jet particles
        for eta, phi, stats in jets:

            # Remove backwards-going particles
            if phi < 0 and phi > -np.pi:
                continue

            # Get dr to tracking station
            dr = dy / (np.sin(to_theta(eta)) * np.sin(phi))
            jet_x, jet_y, jet_z = to_cartesian(dr, eta, phi)

            # Add offset
            jet_x += decay_pos[0]
            jet_y += decay_pos[1]
            jet_z += decay_pos[2]

            # Get color depending on intersection
            c = colors[3] if stats[i + 1] else colors[4]

            # Plot intersection with TS
            ax2[i + 1].scatter([jet_x], [jet_z], marker='x', c=c)

    # Get decay position of particle
    llp_r, llp_eta, llp_phi = to_spherical(decay_pos[0] - origin[0], decay_pos[1], decay_pos[2])
    decays_in_cavern = in_cavern(llp_r, llp_eta, llp_phi)

    # Unroll ceiling
    for eta, phi, stats in jets:

        # Check if particle decays beyond ceiling
        if not decays_in_cavern:
            continue

        jet_x, jet_y, jet_z = ceil_intersection(eta, phi, decay_pos)

        # Don't include if intersects below cavern
        if jet_y < 0:
            continue

        # Get color
        c = colors[3] if stats[0] else colors[4]

        # Plot
        ax2[0].scatter([jet_x], [jet_z], marker='x', c=c)
