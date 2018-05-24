import numpy as np

def r_hat(theta,phi):
    rhx = np.cos(phi) * np.cos(theta)
    rhy = np.sin(phi) * np.cos(theta)
    rhz = np.sin(theta)
    return np.stack((rhx, rhy, rhz))

def theta_hat(theta, phi):
    thx = -np.cos(phi) * np.sin(theta)
    thy = -np.sin(phi) * np.sin(theta)
    thz = np.cos(theta)
    return np.stack((thx, thy, thz))

def phi_hat(theta, phi):
    phx = -np.sin(phi)
    phy = np.cos(phi)
    phz = np.zeros_like(phi)
    return np.stack((phx, phy, phz))

def rotation_matrix(axis, angle):
    """
    Rodrigues' rotation matrix formula
    """
    K = np.array([
        [0.,-axis[2], axis[1]],
        [axis[2],0.,-axis[0]],
        [-axis[1],axis[0],0.]
    ])

    I = np.identity(3)

    R = I + np.sin(angle) * K + (1. - np.cos(angle)) * np.dot(K,K)

    return R

def spherical_coordinates_map(R, theta, phi):
    """
    Returns the spherical coordinates of the point specified by p = R . q,
    where q is the 3D position vector of the point specified by (theta,phi) and
    R is the 3D rotation matrix that relates two coordinate charts.
    """
    q_hat_1 = np.cos(phi) * np.cos(theta)
    q_hat_2 = np.sin(phi) * np.cos(theta)
    q_hat_3 = np.sin(theta)
    q_hat = np.stack((q_hat_1,q_hat_2,q_hat_3))
    p_hat = np.einsum('ab...,b...->a...',R ,q_hat)
    beta = np.arcsin(p_hat[-1,:])
    alpha = np.arctan2(p_hat[1,:],p_hat[0,:])
    return (beta,alpha)

def spherical_basis_transformation_components(theta, phi, R):

    if np.allclose(R, np.eye(3)) == True:
        cosX = np.ones_like(theta)
        sinX = np.zeros_like(theta)

    else:
        beta, alpha = spherical_coordinates_map(R, theta, phi)

        th = theta_hat(theta, phi)
        ph = phi_hat(theta, phi)

        bh = np.einsum('ab...,b...->a...', R.T, theta_hat(beta, alpha))
        ah = np.einsum('ab...,b...->a...', R.T, phi_hat(beta, alpha))

        cosX = np.einsum('a...,a...', bh, th)
        sinX = np.einsum('a...,a...', bh, ph)

    return cosX, sinX
