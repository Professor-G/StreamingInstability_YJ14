#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 12:01:01 2024

@author: daniel
"""
import numpy as np
import pkg_resources
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def dsharp_model(q, wavelength, grain_sizes, bin_approx=False):
    """
    Compute dust absorption and scattering opacities using the DSHARP model.

    This function supports both monodisperse and polydisperse grain size distributions,
    allowing for either a full grain size integration or a binned approximation (`bin_approx`).
    
    Opacities are computed using the DSHARP dust opacity tables (Birnstiel et al. 2018) and
    are valid for wavelengths in the range 1e-5 to 10 cm and grain sizes from 1e-5 to 100 cm.
    It assumes a constant internal density across all grain sizes of 1.675 (g/cm³), as per the adopted dust model.

    Notes
    -----
    When `bin_approx=True`, the grain size distribution is discretized into bins such that:
    - The first bin always starts at 1e-5 cm (minimum DSHARP grain size).
    - Each bin's `a_max` is given by the input `grain_sizes`.
    - Each subsequent bin begins slightly above the previous bin edge (offset by 1e-10 cm)
      to avoid overlap and mass overestimation.

    If `bin_approx=False`, a single grain size is assumed (monodisperse case).

    Warning
    -------
    Opacities cannot be computed for exactly 1e-5 cm unless a small positive offset is applied
    (e.g., `1e-5 + 1e-10` cm), due to our binning strategy.

    Parameters
    ----------
    q : float
        Power-law index of the grain size distribution; n(a) ∝ a^(-q).
    wavelength : float
        Wavelength at which to evaluate the opacities [cm].
    grain_sizes : list or ndarray
        Array of maximum grain sizes for each bin [cm].
    bin_approx : bool, optional
        If True, enables the binned approximation for polydisperse models.
        Default is False (monodisperse case).

    Returns
    -------
    absorption_opacities : np.ndarray
        Absorption opacities (cm²/g) for each grain bin.
    scattering_opacities : np.ndarray
        Scattering opacities (cm²/g) for each grain bin.
    grain_size_bins : np.ndarray
        Array of (a_min, a_max) tuples for each grain size bin [cm].
    """

    # This is the data extracted from the DSHARP code released as part of Birnstiel+18
    resource_package = __name__

    resource_path = '/'.join(('data', 'all_grain_sizes.txt'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    a_orig = np.loadtxt(file)
    
    resource_path = '/'.join(('data', 'all_wavelengths.txt'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    lam = np.loadtxt(file)

    resource_path = '/'.join(('data', 'kappa_abs.npy'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    k_abs_orig = np.load(file)

    resource_path = '/'.join(('data', 'kappa_sca.npy'))
    file = pkg_resources.resource_filename(resource_package, resource_path)
    k_sca_orig = np.load(file)

    # Transpose to have shape (n_nu, n_a)
    kappa_a_nu_orig, kappa_a_nu_orig_sca = k_abs_orig.T, k_sca_orig.T

    # Input checks
    #
    if not isinstance(q, (int, float)): raise TypeError(f"q must be a single numeric value (int or float), but got {type(q)}.")
    if not isinstance(wavelength, (int, float)): raise TypeError(f"wavelength must be a single numeric value (int or float), but got {type(wavelength)}.")
    if not (lam.min() <= wavelength <= lam.max()): raise ValueError(f"Input wavelength {wavelength} cm is outside the valid range [{lam.min()} cm, {lam.max()} cm].")
    if not isinstance(bin_approx, bool): raise TypeError(f"bin_approx must be a boolean, but got {type(bin_approx)}.")
    if isinstance(grain_sizes, (list, np.ndarray)):
        if not all(a_orig.min() <= size <= a_orig.max() for size in grain_sizes): raise ValueError(f"One or more grain sizes are outside the valid range [{a_orig.min()} cm, {a_orig.max()} cm].")
    elif isinstance(grain_sizes, (int, float)):
        if grain_sizes <= 1e-5: raise ValueError('Grain size must be greater than 1e-5 cm!')
        if not (a_orig.min() <= grain_sizes <= a_orig.max()): raise ValueError(f"Grain size {grain_sizes} cm is outside the valid range [{a_orig.min()} cm, {a_orig.max()} cm].")
    else:
        raise TypeError("Grain sizes must be a single value (int or float) or a list/array of values.")
    
    #
    #

    # Find the index jnu corresponding to where lambda = wavelength
    jnu = np.argmin(np.abs(lam - wavelength))
    lambda_jnu = lam[jnu]

    #
    # Need to ensure the grain_size array is sorted! This routine will not work unless the grains are descending order (small to big)
    # As the shearing_box class allows for the grain sizes to be input in any order, the grains will be sorted
    # For use in this routine only after which they will be re-ordered using inverse permutation
    if isinstance(grain_sizes, (list, np.ndarray)):
        grain_sizes = np.array(grain_sizes)
        _original_grain_sizes_ = grain_sizes.copy()
        sort_index = np.argsort(grain_sizes)
        grain_sizes = grain_sizes[sort_index]
    else:
        grain_sizes = np.array([grain_sizes])
        sort_index = [0]

    #
    # Now run loop but get the opacity at the desired wavelength only (thus can only do one wavelength at a time)
    #

    # To store opacities and bins
    num_bins = len(grain_sizes) 
    opacity_abs = np.zeros(num_bins)
    opacity_sca = np.zeros(num_bins)
    bin_arrays = []

    for kappa, kappa_a in zip(['abs', 'sca'], [kappa_a_nu_orig[jnu, :], kappa_a_nu_orig_sca[jnu, :]]): # Kappa will track whether it is absorption or scattering
        #
        a_min = 1e-5 # Limits from the DSHARP dust model
        a_max = max(grain_sizes) 

        # Need to interpolate for full control over parameter space
        a_new = np.logspace(np.log10(a_min), np.log10(a_max), 5000)

        # Ensure that the bin edges are included in a_new, so combine and sort
        a_new = np.unique(np.concatenate((a_new, grain_sizes)))
        a_new.sort()

        # To have full control over param space we need to interp kappa_a to a_new, in log-log space thus first need to remove zeros
        mask = (kappa_a > 0) & (a_orig > 0)
        a_valid, kappa_a_valid= a_orig[mask], kappa_a[mask]

        # Interpolate log(kappa_a) vs log(a)
        interp_func = interp1d(np.log10(a_valid), np.log10(kappa_a_valid), kind='linear', fill_value='extrapolate')
        log_kappa_a_new = interp_func(np.log10(a_new))
        kappa_a_new = 10 ** log_kappa_a_new

        # Compute na, ma, da for a_new
        na = a_new ** (-q)
        ma = a_new ** 3  # assuming constant density therefore m(a) da scales with radius only
        #da = np.gradient(a_new)

        # Define bin edges
        bin_edges_lower = []
        bin_edges_upper = grain_sizes.copy()

        for i in range(num_bins):
            #
            if i == 0:
                bin_edges_lower.append(a_min)  # Start at a_min
            else:
                # Left edge is slightly greater than previous amax to avoid overlap on the left edge
                bin_edges_lower.append(grain_sizes[i - 1] + 1e-10)

        #print('Bin edges (cm):')
        #for i in range(num_bins): print(f'Bin {i+1}: lower edge = {bin_edges_lower[i]}, upper edge = {bin_edges_upper[i]}')
        #
        #
        for i in range(num_bins):
            #
            idx_min_full = 0
            idx_max_full = np.searchsorted(a_new, bin_edges_upper[i], side='right') - 1
            #
            idx_min_bin = np.searchsorted(a_new, bin_edges_lower[i], side='left')
            idx_max_bin = idx_max_full
            #
            if bin_approx is False:
                # Assuming full grain size distribution going from a_min=1e-5 to a_max
                if idx_max_full < idx_min_full or idx_max_full >= len(a_new):
                    if kappa == 'abs':
                        opacity_abs[i] = np.nan
                    else:
                        opacity_sca[i] = np.nan
                else:
                    numerator_full = np.trapz(
                        kappa_a_new[idx_min_full:idx_max_full + 1] *
                        na[idx_min_full:idx_max_full + 1] *
                        ma[idx_min_full:idx_max_full + 1],
                        x=a_new[idx_min_full:idx_max_full + 1]
                    )
                    denominator_full = np.trapz(
                        na[idx_min_full:idx_max_full + 1] *
                        ma[idx_min_full:idx_max_full + 1],
                        x=a_new[idx_min_full:idx_max_full + 1]
                    )
                    if kappa == 'abs':
                        opacity_abs[i] = numerator_full / denominator_full if denominator_full != 0 else np.nan
                    else:
                        opacity_sca[i] = numerator_full / denominator_full if denominator_full != 0 else np.nan
                    #
                    if kappa == 'abs': bin_arrays.append(a_new[idx_min_full:idx_max_full + 1])  #Only append bins the first time
            #
            elif bin_approx:
                # Binned approximation
                if idx_max_bin < idx_min_bin or idx_max_bin >= len(a_new):
                    if kappa == 'abs':
                        opacity_abs[i] = np.nan  # No grains in bin
                    else:
                        opacity_sca[i] = np.nan
                else:
                    numerator_bin = np.trapz(
                        kappa_a_new[idx_min_bin:idx_max_bin + 1] *
                        na[idx_min_bin:idx_max_bin + 1] *
                        ma[idx_min_bin:idx_max_bin + 1],
                        x=a_new[idx_min_bin:idx_max_bin + 1]
                    )
                    denominator_bin = np.trapz(
                        na[idx_min_bin:idx_max_bin + 1] *
                        ma[idx_min_bin:idx_max_bin + 1],
                        x=a_new[idx_min_bin:idx_max_bin + 1]
                    )
                    #print('bin', a_new[idx_min_bin:idx_max_bin+1])
                    if kappa == 'abs':
                        opacity_abs[i] = numerator_bin / denominator_bin if denominator_bin != 0 else np.nan
                    else:
                        opacity_sca[i] = numerator_bin / denominator_bin if denominator_bin != 0 else np.nan
                    #
                    if kappa == 'abs': bin_arrays.append(a_new[idx_min_bin:idx_max_bin+1])  #Only append bins the first time
                    
    # Reorder outputs to match the original order of grain_sizes
    inverse_sort_index = np.argsort(sort_index)
    opacity_abs = opacity_abs[inverse_sort_index]
    opacity_sca = opacity_sca[inverse_sort_index]
    bin_arrays = [bin_arrays[i] for i in inverse_sort_index]

    return opacity_abs, opacity_sca, bin_arrays
    