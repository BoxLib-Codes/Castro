#!/bin/env python3

import glob
import os
import operator

import numpy as np
import matplotlib.pyplot as plt

import detonation as dt

import yt
yt.funcs.mylog.setLevel(50)


def basic_plot(dset, label):
    """label should have format specifiers like {0.cfl},
    and we'll fill them in from the Detonation object"""

    fig, axs = plt.subplots(2, 1, figsize=(7, 10), constrained_layout=True)

    ax1 = axs.flatten()[0]
    ax2 = axs.flatten()[1]


    for p in dset:
        ax1.plot(p.data.x, p.data.T, label=label.format(p))
        ax2.plot(p.data.x, np.abs(p.data.enuc), label=label.format(p))

        ax1.legend(frameon=False, fontsize="small")
        ax1.set_ylabel("T [K]")
        ax1.set_yscale("log")

        ax2.legend(frameon=False, fontsize="small")
        ax2.set_xlabel("x [cm]")
        ax2.set_ylabel("enuc [erg/g/s]")
        ax2.set_yscale("log")
        ax2.set_ylim(1.e14)

    return fig


if __name__ == "__main__":

    # get all the data
    run_dirs = glob.glob("det_[rst]*")
    runs = []
    for run in run_dirs:
        try:
            if not os.path.isfile(os.path.join(run, "Backtrace.0")):
                print(run)
                runs.append(dt.Detonation(run))
        except IndexError:
            # the run didn't produce output -- it might still be running?
            print(f"run {run} didn't produce output")

    runs.sort()

    nzones = {q.nzones for q in runs}
    dtnuces = {q.dtnuce for q in runs}
    CFLs = {q.cfl for q in runs}

    # for a single resolution show SDC profile vs # of iters

    nz = sorted(list(nzones))[1]
    dset = [q for q in runs if q.integrator == "simplified-SDC" and q.nzones == nz]

    print("making master SDC plot")
    print(dset)

    fig = basic_plot(dset, "CFL = {0.cfl}; dtnuc_e = {0.dtnuce}; iters = {0.niters}")

    fig.suptitle(f"simplified-SDC, nzones = {nz}")
    fig.savefig(f"profile_sdc_iter_compare_{nz}.png")

    # make a plot of T, enuc vs. x for different SDC CFL
    # group by different dtnuc_e

    nz = sorted(list(nzones))[1]
    for dtnuce in dtnuces:
        dset = [q for q in runs if q.integrator == "simplified-SDC" and
                q.nzones == nz and q.dtnuce == dtnuce and q.niters == 2]

        fig = basic_plot(dset, "CFL = {0.cfl}")

        fig.suptitle(f"simplified-SDC, 2 iterations, nzones = {nz}, dtnuc_e = {dtnuce}")
        fig.savefig(f"profile_sdc_2iter_{nz}_dtnuce_{dtnuce}.png")


    # make a plot of T, enuc vs. x for Strang and SDC
    nz = sorted(list(nzones))[1]
    CFL = sorted(list(CFLs))[-1]

    for dtnuce in dtnuces:
        dset = [q for q in runs if q.nzones == nz and q.dtnuce == dtnuce and q.cfl == CFL]
        
        fig = basic_plot(dset, "{0.integrator} {0.niter_str}")

        fig.suptitle(f"nzones = {nz}, CFL = {CFL}, dtnuc_e = {dtnuce}")
        fig.savefig(f"profile_intergrator_comp_{nz}_cfl_{CFL}_dtnuce_{dtnuce}.png")
