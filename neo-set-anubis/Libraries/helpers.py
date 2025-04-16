import numpy as np

## Helper functions

# Calculate boost (gamma)
def calculate_boost(pt, eta, m):
    p = pt * np.cosh(eta)
    e = np.sqrt(p * p + m * m)
    return e / m

# Calculate beta
def calculate_beta(gamma):
    return np.sqrt(1 - (1 / np.square(gamma)))

# Convert eta to theta
def to_theta(eta):
    return 2 * np.arctan(np.exp(-eta))

# Convert theta to eta
def to_eta(theta):
    return -1 * np.log(np.tan(theta / 2))

# Scale phi to [-pi, pi]
def scale_phi(phi):
    return ((phi + np.pi) % (2 * np.pi)) - np.pi

# Exponential decay of particle with certain boost and ctau
def expo_decay(boost, beta, ctau):
    return np.random.exponential(scale=(boost * beta * ctau))

# Convert r, eta, phi to Cartesian coordinates
def to_cartesian(r, eta, phi):
    x = r * np.sin(to_theta(eta)) * np.cos(phi)
    y = r * np.sin(to_theta(eta)) * np.sin(phi)
    z = r * np.cos(to_theta(eta))

    return x, y, z

# Convert Cartesian coordinates to r, eta, phi
def to_spherical(x, y, z):
    # Perform conversion
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arctan2(np.sqrt(np.square(x) + np.square(y)), z)
    if x > 0:
        phi = np.arctan(y / x)
    elif x < 0 and y >= 0:
        phi = np.arctan(y / x) + np.pi
    elif x < 0 and y < 0:
        phi = np.arctan(y / x) - np.pi
    elif x == 0 and y > 0:
        phi = np.pi / 2
    elif x == 0 and y < 0:
        phi = -np.pi / 2
    else:
        phi = np.nan

    # Scale phi
    phi = scale_phi(phi)

    return r, to_eta(theta), phi

# Rotate about x-axis by alpha
def rotate_x(r, theta, phi, alpha):
    # Convert to Cartesian
    x, y, z = to_cartesian(r, to_eta(theta), phi)

    # Apply rotation
    new_x = x
    new_y = y * np.cos(alpha) - z * np.sin(alpha)
    new_z = y * np.sin(alpha) + z * np.cos(alpha)

    # Convert back to spherical
    new_r, new_eta, new_phi = to_spherical(new_x, new_y, new_z)

    # Return result
    return new_r, to_theta(new_eta), new_phi

# Rotate about z-axis by alpha
def rotate_z(r, theta, phi, alpha):
    # Convert to Cartesian
    x, y, z = to_cartesian(r, to_eta(theta), phi)

    # Apply rotation
    new_x = x * np.cos(alpha) - y * np.sin(alpha)
    new_y = x * np.sin(alpha) + y * np.cos(alpha)
    new_z = z

    # Convert back to spherical
    new_r, new_eta, new_phi = to_spherical(new_x, new_y, new_z)

    # Return result
    return new_r, to_theta(new_eta), new_phi
