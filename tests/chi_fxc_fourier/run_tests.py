"""Build and run the actual TDDFT_Chi_Fourier module against a small grid."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fc", default="gfortran", help="Fortran compiler")
    args = parser.parse_args()
    tests = Path(__file__).resolve().parent
    repo = tests.parents[1]
    with tempfile.TemporaryDirectory(prefix="chi-fxc-fourier-") as temporary:
        build = Path(temporary)
        program = build / "test_chi_fourier.exe"
        subprocess.run(
            [
                args.fc,
                "-cpp",
                "-ffree-form",
                "-ffree-line-length-none",
                "-fcheck=all",
                str(tests / "pars.f90"),
                str(repo / "src/tddft/TDDFT_Chi_Fourier.F"),
                str(tests / "test_fourier.f90"),
                "-o",
                str(program),
            ],
            cwd=build,
            check=True,
        )
        subprocess.run([str(program)], cwd=build, check=True)


if __name__ == "__main__":
    main()
