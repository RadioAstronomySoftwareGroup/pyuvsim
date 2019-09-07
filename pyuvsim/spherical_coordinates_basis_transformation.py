# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np


def r_hat(theta, phi):
    rhx = np.cos(phi) * np.sin(theta)
    rhy = np.sin(phi) * np.sin(theta)
    rhz = np.cos(theta)
    return np.stack((rhx, rhy, rhz))


def theta_hat(theta, phi):
    thx = np.cos(phi) * np.cos(theta)
    thy = np.sin(phi) * np.cos(theta)
    thz = -np.sin(theta)
    return np.stack((thx, thy, thz))


def phi_hat(theta, phi):
    phx = -np.sin(phi)
    phy = np.cos(phi)
    phz = np.zeros_like(phi)
    return np.stack((phx, phy, phz))


def spherical_coordinates_map(R, theta, phi):
    """
    Returns the spherical coordinates of the point specified by p = R . q,
    where q is the 3D position vector of the point specified by (theta,phi) and
    R is the 3D rotation matrix that relates two coordinate charts.
    """
    # This is NOT written to be vectorized for multiple (theta, phi)

    # Replace with function call?
    q_hat_1 = np.cos(phi) * np.sin(theta)
    q_hat_2 = np.sin(phi) * np.sin(theta)
    q_hat_3 = np.cos(theta)
    q_hat = np.stack((q_hat_1, q_hat_2, q_hat_3))

    p_hat = np.einsum('ab...,b...->a...', R, q_hat)
    # Should test for shape of p_hat

    # Should write a function to do this as well, i.e., pull back angles from
    # a vector (or use healpix?)
    beta = np.arccos(p_hat[2])
    alpha = np.arctan2(p_hat[1], p_hat[0])
    if alpha < 0:
        alpha += 2. * np.pi
    # alpha[alpha < 0] += 2. * np.pi

    return (beta, alpha)


def spherical_basis_transformation_components(theta, phi, R):
    beta, alpha = spherical_coordinates_map(R, theta, phi)

    th = theta_hat(theta, phi)
    ph = phi_hat(theta, phi)

    bh = np.einsum('ab...,b...->a...', R.T, theta_hat(beta, alpha))
    ah = np.einsum('ab...,b...->a...', R.T, phi_hat(beta, alpha))

    cosX = np.einsum('a...,a...', bh, th)
    sinX = np.einsum('a...,a...', bh, ph)

    return cosX, sinX, th, ph, bh, ah, beta, alpha

def spherical_basis_transformation_components2(beta, alpha, theta, phi, R):
    #beta, alpha = spherical_coordinates_map(R, theta, phi)

    th = theta_hat(theta, phi)
    ph = phi_hat(theta, phi)

    bh = np.einsum('ab...,b...->a...', R.T, theta_hat(beta, alpha))
    ah = np.einsum('ab...,b...->a...', R.T, phi_hat(beta, alpha))

    cosX = np.einsum('a...,a...', bh, th)
    sinX = np.einsum('a...,a...', bh, ph)

    return cosX, sinX, th, ph, bh, ah, beta, alpha

def spherical_basis_transformation_components_two_points(beta, alpha, theta, phi, R):
    #beta, alpha = spherical_coordinates_map(R, theta, phi)

    th = theta_hat(theta, phi)
    ph = phi_hat(theta, phi)

    bh = np.einsum('ab...,b...->a...', R.T, theta_hat(beta, alpha))
    ah = np.einsum('ab...,b...->a...', R.T, phi_hat(beta, alpha))

    cosX = np.einsum('a...,a...', bh, th)
    sinX = np.einsum('a...,a...', bh, ph)

    return cosX, sinX #, th, ph, bh, ah, beta, alpha

# From spherical_coordinates_basis_transformation in cst2ijones
def axis_angle_rotation_matrix(axis, angle):
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

def verify_orthogonal(R,tol=1e-15):
    return np.allclose(np.matmul(R,R.T), np.eye(3), atol=tol)

def verify_unit(v,tol=1e-15):
    return np.allclose(np.dot(v,v), 1, atol=tol)

def pts2rot(theta1, phi1, theta2, phi2):
    """ Given two points on the sphere, find the rotation matrix that connects them.  
    Probably done somewhere else better """
    r1 = scbt.r_hat(theta1, phi1)
    r2 = scbt.r_hat(theta2, phi2)
    n = np.cross(r1,r2)
    # Note that Psi is between 0 and pi
    sinPsi = np.sqrt(np.dot(n,n))
    n_hat = n/sinPsi # Trouble lurks if Psi = 0.
    cosPsi = np.dot(r1,r2)
    Psi = np.arctan2(sinPsi, cosPsi)
    rotation = axis_angle_rotation_matrix(n_hat, Psi)
    try:
        assert(verify_unit(r1))
        assert(verify_unit(r2))
        assert(verify_unit(n_hat))
        assert(verify_orthogonal(rotation))
    except:
        print('r1', r1)
        print('r2', r2) 
        print('n', n) 
        print('len n', np.dot(n,n)) 
        print('n_hat', n_hat) 
        print('len n_hat', np.dot(n_hat,n_hat)) 
        print('sinPsi', sinPsi) 
        print('R', rotation) 
        raise
    return rotation

def vecs2rot(r1, r2):
    """ Given two vectors on the sphere, find the rotation matrix that connects them.  
    Probably done somewhere else better """
    #r1 = scbt.r_hat(theta1, phi1)
    #r2 = scbt.r_hat(theta2, phi2)
    n = np.cross(r1,r2)
    # Note that Psi is between 0 and pi
    sinPsi = np.sqrt(np.dot(n,n))
    n_hat = n/sinPsi # Trouble lurks if Psi = 0.
    cosPsi = np.dot(r1,r2)
    Psi = np.arctan2(sinPsi, cosPsi)
    rotation = axis_angle_rotation_matrix(n_hat, Psi)
    try:
        assert(verify_unit(r1))
        assert(verify_unit(r2))
        assert(verify_unit(n_hat))
        assert(verify_orthogonal(rotation))
    except:
        print('r1', r1)
        print('r2', r2) 
        print('n', n) 
        print('len n', np.dot(n,n)) 
        print('n_hat', n_hat) 
        print('len n_hat', np.dot(n_hat,n_hat)) 
        print('sinPsi', sinPsi) 
        print('R', rotation) 
        raise
    return rotation
