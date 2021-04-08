import glob
import os
import operator
import tqdm

import numpy as np
import matplotlib.pyplot as plt

import yt

yt.funcs.mylog.setLevel(50)

class Profile:
    """read a plotfile using yt and store the 1d profile for T and enuc"""

    def __init__(self, plotfile):
        ds = yt.load(plotfile)

        time = float(ds.current_time)
        ad = ds.all_data()

        # Sort the ray values by 'x' so there are no discontinuities
        # in the line plot
        srt = np.argsort(ad['x'])
        x_coord = np.array(ad['x'][srt])
        temp = np.array(ad['Temp'][srt])
        enuc = np.array(ad['enuc'][srt])

        self.time = time
        self.x = x_coord
        self.T = temp
        self.enuc = enuc

    def find_x_for_T(self, T_0=1.e9):
        """ given a profile x(T), find the x_0 that corresponds to T_0 """

        # our strategy here assumes that the hot ash is in the early
        # part of the profile.  We then find the index of the first
        # point where T drops below T_0
        try:
            idx = np.where(self.T < T_0)[0][0]
        except IndexError:
            idx = len(self.T)-1

        T1 = self.T[idx-1]
        x1 = self.x[idx-1]

        T2 = self.T[idx]
        x2 = self.x[idx]

        slope = (x2 - x1)/(T2 - T1)

        return x1 + slope*(T_0 - T1)


class Detonation:
    def __init__(self, name, precompute=True):
        self.name = name

        self.cfl = None
        self.nzones = None
        self.integrator = None
        self.niters = None
        self.order = None
        self.quadrature = None
        self.dtnuce = None
        self.has_started = False
        self.crashed = False

        # read the meta data
        with open(os.path.join(name, "run.meta")) as mf:
            for line in mf.readlines():
                k, v = [q.strip() for q in line.split("=")]
                if k == "cfl":
                    self.cfl = float(v)
                elif k == "nzones":
                    self.nzones = int(v)
                elif k == "integrator":
                    self.integrator = v.strip()
                elif k == "niters":
                    self.niters = int(v)
                elif k == "dtnuc_e":
                    self.dtnuce = float(v)
                elif k == "order":
                    self.order = int(v)
                elif k == "quadrature":
                    self.quadrature = int(v)

        if self.niters is not None:
            self.niter_str = f"(iters = {self.niters})"
        else:
            self.niter_str = ""

        cwd = os.getcwd()
        os.chdir(name)

        # did we crash?
        if os.path.isfile("Backtrace.0"):
            self.crashed = True

        # find all the output (plot) files
        self.files = glob.glob("*plt?????")
        self.files.sort()
        os.chdir(cwd)

        if precompute:

            # precompute the velocity and the data profiles
            if len(self.files) >= 3:
                self.v, self.v_sigma = self.get_velocity()
            else:
                self.v, self.v_sigma = 0.0, 0.0

            if len(self.files) >= 1:
                self.data = self.get_data()

                self.end_time = self.data.time
                self.nsteps = int(self.files[-1].split("plt")[1])
                self.has_started = True
            else:
                self.end_time = -1.0
                self.nsteps = -1.0
                self.has_started = False

    def __repr__(self):
        return self.name

    def __lt__(self, other):
        """sort by CFL number and resolution and then # of SDC
        iterations"""

        # first sort on integrator
        if self.integrator == other.integrator:

            # next sort on iterators
            if self.niters == other.niters:

                # next sort on order
                if self.order == other.order:

                    # next sort on CFL
                    if self.cfl == other.cfl:

                        # next sort on dtnuce
                        if self.dtnuce == other.dtnuce:

                            # finally sort on the number of zones
                            return self.nzones < other.nzones

                        else:
                            return self.dtnuce < other.dtnuce

                    else:
                        return self.cfl < other.cfl

                else:
                    return self.order < other.order

            else:
                return self.niters < other.niters

        else:
            # strang comes first
            if self.integrator == "Strang":
                return True
            return False

    def get_velocity(self):
        """look at the 3 pairings of the last 3 plotfiles and estimate the
        velocity by finite-differencing

        """

        vs = []
        pairs = [(-2, -1), (-3, -1), (-3, -2)]

        for i1, i2 in pairs:
            f1 = self.files[i1]
            p1 = Profile(os.path.join(self.name, f1))

            f2 = self.files[i2]
            p2 = Profile(os.path.join(self.name, f2))

            # we'll do this by looking at 3 different temperature
            # thresholds and averaging
            T_ref = [2.e9, 3.e9, 4.e9]

            for T0 in T_ref:
                x1 = p1.find_x_for_T(T0)
                x2 = p2.find_x_for_T(T0)
                vs.append((x1 - x2)/(p1.time - p2.time))

        vs = np.array(vs)
        v = np.mean(vs)
        v_sigma = np.std(vs)
        return v, v_sigma

    def get_data(self):
        """get the temperature and energy generation rate from the last
        plotfile"""

        return Profile(os.path.join(self.name, self.files[-1]))

