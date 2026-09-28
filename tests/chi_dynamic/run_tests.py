"""Build and run the frequency-dependent Casida numerics against LAPACK references."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fc", default="gfortran", help="Fortran compiler")
    parser.add_argument("--lapack", default="-llapack -lblas", help="LAPACK/BLAS link flags")
    parser.add_argument("--openmp", action="store_true", help="compile with -fopenmp (threaded frequency loop)")
    args = parser.parse_args()
    tests = Path(__file__).resolve().parent
    repo = tests.parents[1]
    with tempfile.TemporaryDirectory(prefix="chi-dynamic-") as temporary:
        build = Path(temporary)
        program = build / "test_chi_dynamic.exe"
        subprocess.run(
            [
                args.fc,
                "-cpp",
                "-ffree-form",
                "-ffree-line-length-none",
                "-fcheck=all",
                *(["-fopenmp"] if args.openmp else []),
                str(tests / "pars.f90"),
                str(repo / "src/tddft/TDDFT_Chi_dynamic.F"),
                str(tests / "test_dynamic.f90"),
                "-o",
                str(program),
                *args.lapack.split(),
            ],
            cwd=build,
            check=True,
        )
        subprocess.run([str(program)], cwd=build, check=True)


if __name__ == "__main__":
    main()
