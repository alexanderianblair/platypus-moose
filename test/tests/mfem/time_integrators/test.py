import mms
import unittest
from mooseutils import fuzzyAbsoluteEqual

# Expected temporal convergence order of each time integration scheme, keyed by the command line
# argument selecting it, with the initial timestep of a refinement study whose last refinement is
# in the asymptotic regime
SCHEMES = {
    "Executioner/scheme=implicit-euler": (1, 0.1),
    "Executioner/scheme=explicit-euler": (1, 0.05),
    "Executioner/scheme=crank-nicolson": (2, 0.1),
    "Executioner/scheme=explicit-midpoint": (2, 0.025),
    "Executioner/scheme=dirk": (2, 0.1),
    "Executioner/scheme=explicit-tvd-rk-2": (2, 0.025),
    "Executioner/mfem_scheme=implicit-midpoint": (2, 0.1),
    "Executioner/mfem_scheme=sdirk23": (3, 0.00625),
    "Executioner/mfem_scheme=sdirk33": (3, 0.00625),
    "Executioner/mfem_scheme=sdirk34": (4, 0.003125),
    "Executioner/mfem_scheme=esdirk32": (2, 0.1),
    "Executioner/mfem_scheme=esdirk33": (3, 0.0125),
    "Executioner/mfem_scheme=generalized-alpha": (2, 0.1),
    "Executioner/mfem_scheme=rk3-ssp": (3, 0.0125),
    "Executioner/mfem_scheme=rk4": (4, 0.0125),
}

# The reaction kernel and convective boundary condition are treated explicitly, and the diffusion
# kernel implicitly
IMEX_SPLIT = ["Kernels/reaction/implicit=false", "BCs/right/implicit=false"]
IMEX_SCHEMES = {
    "Executioner/mfem_scheme=imex-euler": (1, 0.1),
    "Executioner/mfem_scheme=imex-ars222": (2, 0.025),
    "Executioner/mfem_scheme=imex-ars232": (2, 0.05),
    "Executioner/mfem_scheme=imex-ars343": (3, 0.025),
}


def run(scheme, dt, *args):
    name = scheme.split("=")[-1]
    df = mms.run_temporal(
        "mms.i", 3, scheme, *args, file_base=name + "_{}", dt=dt, y_pp=["l2_error"]
    )
    fig = mms.ConvergencePlot(xlabel=r"$\Delta$t", ylabel="$L_2$ Error")
    fig.plot(df, label=name, marker="o", markersize=8, slope_precision=2, num_fitted_points=2)
    return fig.label_to_slope[name]


class TestSchemes(unittest.TestCase):
    def test(self):
        for scheme, (order, dt) in SCHEMES.items():
            with self.subTest(scheme=scheme):
                slope = run(scheme, dt)
                print("%s, %f" % (scheme, slope))
                self.assertTrue(fuzzyAbsoluteEqual(slope, order, 0.1))


class TestIMEXSchemes(unittest.TestCase):
    def test(self):
        for scheme, (order, dt) in IMEX_SCHEMES.items():
            with self.subTest(scheme=scheme):
                slope = run(scheme, dt, *IMEX_SPLIT)
                print("%s, %f" % (scheme, slope))
                self.assertTrue(fuzzyAbsoluteEqual(slope, order, 0.1))
