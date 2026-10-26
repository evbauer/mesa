#!/usr/bin/env python3
"""Write a minimal conforming WDM_AV processed HDF5 table for testing the Fortran
loader, following eos/eosWDM_AV_builder/FORMAT.md.

The nodal data are the analytic "ideal" model used by test_interp.f90:
fully ionized ideal ions + ideal nondegenerate electrons, no radiation.  All
constants match MESA const_def so the Fortran comparison matches to roundoff.

Usage:  python3 -I make_test_h5.py <out_dir>

Writes <out_dir>/eosWDM_AV_he4.h5          (valid, format_version=1)
       <out_dir>/eosWDM_AV_badver.h5       (format_version=99, otherwise valid)

Datasets are written with h5py shape (ntau, nxi) so that forum/Fortran receive
them dimensioned (nxi, ntau).  String attributes are written as fixed-length
ASCII (forum reads character attributes with a fixed-length string type).
"""

import sys
import os
import math
import numpy as np
import h5py

# constants, matching MESA const_def exactly
PI = 3.1415926535897932384626433832795028841971693993751
AVO = 6.02214076e23
AMU = 1.0 / AVO
KERG = 1.380649e-16
PLANCK_H = 6.62607015e-27
ME = 9.1093837015e-28


def model_constants(A, Z):
    K = KERG / AMU
    B = (1.0 + Z) / A
    Ci = (2.0 * PI * (A * AMU) * KERG / PLANCK_H**2) ** 1.5
    Ce = (2.0 * PI * ME * KERG / PLANCK_H**2) ** 1.5
    a_i = math.log(Ci) + math.log(A * AMU) + 1.0
    a_e = math.log(2.0 * Ce) - math.log(Z) + math.log(A * AMU) + 1.0
    A0 = (a_i + Z * a_e) / A
    return K, B, A0


def nodal(A, Z, xi, tau):
    """return (f, f_xi, f_tau, f_xixi, f_tautau, f_xitau,
               f_xixitau, f_xitautau, f_xixitautau)"""
    K, B, A0 = model_constants(A, Z)
    E = math.exp(tau)
    P1 = A0 + 1.5 * B * tau - B * xi
    f = -K * E * P1
    f_xi = K * B * E
    f_tau = -K * E * (P1 + 1.5 * B)
    f_xixi = 0.0
    f_tautau = -K * E * (P1 + 3.0 * B)
    f_xitau = K * B * E
    f_xixitau = 0.0
    f_xitautau = K * B * E
    f_xixitautau = 0.0
    return (f, f_xi, f_tau, f_xixi, f_tautau, f_xitau,
            f_xixitau, f_xitautau, f_xixitautau)


def write_table(path, fmt_version, A, Z, isotope, element,
                nxi, ntau, xi_min, dxi, tau_min, dtau, order=5):
    names = ['f', 'f_xi', 'f_tau', 'f_xixi', 'f_tautau', 'f_xitau',
             'f_xixitau', 'f_xitautau', 'f_xixitautau']
    # h5py shape (ntau, nxi)
    data = {nm: np.zeros((ntau, nxi), dtype='f8') for nm in names}
    for j in range(ntau):
        tau = tau_min + j * dtau
        for i in range(nxi):
            xi = xi_min + i * dxi
            vals = nodal(A, Z, xi, tau)
            for k, nm in enumerate(names):
                data[nm][j, i] = vals[k]

    with h5py.File(path, 'w') as f:
        f.attrs['format_version'] = np.int32(fmt_version)
        f.attrs.create('isotope', np.bytes_(isotope))
        f.attrs.create('element', np.bytes_(element))
        f.attrs['Z'] = np.int32(Z)
        f.attrs['A'] = np.float64(A)
        f.attrs['interp_order'] = np.int32(order)
        f.attrs['nxi'] = np.int32(nxi)
        f.attrs['ntau'] = np.int32(ntau)
        f.attrs['xi_min'] = np.float64(xi_min)
        f.attrs['dxi'] = np.float64(dxi)
        f.attrs['tau_min'] = np.float64(tau_min)
        f.attrs['dtau'] = np.float64(dtau)
        f.attrs['includes_radiation'] = np.int32(0)
        f.attrs['e_offset'] = np.float64(0.0)
        f.attrs['f_offset'] = np.float64(0.0)
        f.attrs['has_Zbar'] = np.int32(0)
        f.attrs['has_eta'] = np.int32(0)
        for nm in names:
            f.create_dataset(nm, data=data[nm])
    print(f'wrote {path}  ({nxi} x {ntau}, order {order}, format_version {fmt_version})')


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out_dir, exist_ok=True)

    ln10 = math.log(10.0)
    # logRho in [-10, 6], 10 pts/decade; logT in [3, 9], 10 pts/decade
    nxi = 161
    ntau = 61
    xi_min = ln10 * (-10.0)
    dxi = ln10 * 0.1
    tau_min = ln10 * 3.0
    dtau = ln10 * 0.1

    A_he4 = 4.0026
    Z_he4 = 2

    write_table(os.path.join(out_dir, 'eosWDM_AV_he4.h5'), 1,
                A_he4, Z_he4, 'he4', 'He', nxi, ntau, xi_min, dxi, tau_min, dtau)
    write_table(os.path.join(out_dir, 'eosWDM_AV_badver.h5'), 99,
                A_he4, Z_he4, 'he4', 'He', nxi, ntau, xi_min, dxi, tau_min, dtau)


if __name__ == '__main__':
    main()
