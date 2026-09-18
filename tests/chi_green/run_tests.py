"""Compile the production numerical module without a Yambo/NetCDF installation.

Usage: python tests/chi_green/run_tests.py --fc gfortran
Full synthetic pipeline: add --lapack-source /path/to/Reference-LAPACK.
Host-array precision: add --single-precision (default: double).
Real MPI bubble: add --mpifc mpifort --mpiexec mpiexec.
No downloads or writes outside a temporary build directory.
"""
import argparse
import os
import re
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
parser.add_argument('--single-precision', action='store_true', help='Use single-precision Yambo host arrays')
parser.add_argument('--lapack-source', type=Path,
                    help='Optional Reference-LAPACK source tree for the native pipeline test')
parser.add_argument('--mpifc', help='Optional GNU-compatible MPI Fortran wrapper for real MPI bubble tests')
parser.add_argument('--mpiexec', default='mpiexec', help='MPI launcher for --mpifc')
parser.add_argument('--mpi-ranks', nargs='+', type=int, default=[2, 3],
                    help='MPI rank counts; default includes an empty synthetic rank')
args = parser.parse_args()
if any(ranks < 1 for ranks in args.mpi_ranks):
    parser.error('--mpi-ranks must be positive')
root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='yambo-chi-green-') as directory:
    build = Path(directory)
    (build / 'pars.f90').write_text(
        'module pars\n use iso_fortran_env, only:real64\n'
        ' integer,parameter::DP=real64\nend module\n', encoding='utf-8')
    executable = build / ('numerics.exe' if os.name == 'nt' else 'numerics')
    subprocess.run([args.fc, '-cpp', '-ffree-form', '-std=f2008', '-O0',
                    '-Wall', '-Wextra', '-fcheck=all',
                    str(build / 'pars.f90'),
                    str(root / 'src/modules/mod_Chi_Green.F'),
                    str(root / 'tests/chi_green/test_numerics.F90'),
                    '-o', str(executable)], cwd=build, check=True)
    subprocess.run([str(executable)], cwd=build, check=True)
    if args.lapack_source or args.mpifc:
        # Production free-form files use Yambo memory macros. This fixture
        # removes only accounting; allocation/deallocation semantics remain.
        (build / 'y_memory.h').write_text(
            ' implicit none\n#define YAMBO_ALLOC(x,s) allocate(x s)\n'
            '#define YAMBO_FREE(x) if(allocated(x)) deallocate(x)\n', encoding='utf-8')
    if args.lapack_source:
        sources = [root / 'tests/chi_green/pipeline_mocks.F90',
                   root / 'src/modules/mod_QP.F',
                   root / 'src/modules/mod_Chi.F',
                   root / 'src/modules/mod_Chi_Green.F',
                   root / 'src/qp/QP_Green_Function.F',
                   root / 'src/pol_function/Chi_G_bubble.F',
                   root / 'src/pol_function/Chi_fxc.F',
                   root / 'tests/chi_green/pipeline_io.F90',
                   root / 'tests/chi_green/test_pipeline.F90']
        objects = []
        for index, source in enumerate(sources):
            obj = build / f'host_{index}.o'
            precision = ['-D_TEST_SINGLE'] if args.single_precision else []
            subprocess.run([args.fc, '-cpp', '-D_TEST_REAL_QP_MODULE', *precision,
                            '-ffree-form', '-ffree-line-length-none',
                            '-I', str(build), '-fcheck=all', '-O0', '-c', str(source),
                            '-o', str(obj)], cwd=build, check=True)
            objects.append(obj)
        lapack = args.lapack_source.resolve()
        available = {}
        for folder in ['BLAS/SRC', 'SRC', 'INSTALL']:
            for source in (lapack / folder).glob('*'):
                if source.suffix.lower() in ['.f', '.f90']:
                    available.setdefault(source.stem.lower(), source)
        executable = build / ('pipeline.exe' if os.name == 'nt' else 'pipeline')
        compiled = set()
        for iteration in range(30):
            link = subprocess.run([args.fc, *map(str, objects), '-o', str(executable)],
                                  cwd=build, capture_output=True, text=True)
            if link.returncode == 0:
                break
            names = set(re.findall(r"undefined reference to [`']([a-zA-Z0-9_]+)_['’]", link.stderr))
            needed = names - compiled
            if not needed or any(name not in available for name in needed):
                raise RuntimeError(link.stderr)
            for name in sorted(needed):
                source = available[name]
                obj = build / f'lapack_{name}.o'
                form = '-ffree-form' if source.suffix.lower() == '.f90' else '-ffixed-form'
                subprocess.run([args.fc, form, '-w', '-O0', '-c', str(source), '-o', str(obj)],
                               cwd=build, check=True)
                objects.append(obj)
                compiled.add(name)
        else:
            raise RuntimeError('Could not resolve LAPACK dependency closure')
        subprocess.run([str(executable)], cwd=build, check=True)
        for scenario, message in [('conditioning', 'ill-conditioned'),
                                  ('missing_state', 'cover every response band'),
                                  ('legacy', 'legacy spectra are unsupported'),
                                  ('reference', 'starting energies differ'),
                                  ('static_wrong_kind', 'actual COHSEX ndb.QP'),
                                  ('static_residue', 'requires unit residues'),
                                  ('static_linewidth', 'cannot have a QP linewidth'),
                                  ('static_crossing', 'cross the fixed KS chemical potential'),
                                  ('static_missing', 'cover every response band'),
                                  ('static_reference', 'starting energies differ'),
                                  ('static_offdiag', 'off-diagonal COHSEX'),
                                  ('static_duplicate', 'duplicate COHSEX state'),
                                  ('mpi_band', 'use k-only MPI'),
                                  ('mpi_q', 'use k-only MPI'),
                                  ('mpi_g', 'use k-only MPI')]:
            rejected = subprocess.run([str(executable), scenario], cwd=build,
                                      capture_output=True, text=True)
            if rejected.returncode == 0 or message not in rejected.stdout:
                raise RuntimeError(f'Pipeline did not reject {scenario}: {rejected.stdout}')
            print(f'PASS: {scenario} rejected')
    if args.mpifc:
        # Separate module directory prevents serial/MPI fixture contamination.
        mpi_build = build / 'mpi'
        mpi_build.mkdir()
        sources = [root / 'tests/chi_green/pipeline_mocks.F90',
                   root / 'src/modules/mod_QP.F',
                   root / 'src/modules/mod_Chi.F',
                   root / 'src/modules/mod_Chi_Green.F',
                   root / 'src/pol_function/Chi_G_bubble.F',
                   root / 'tests/chi_green/pipeline_io.F90',
                   root / 'tests/chi_green/test_mpi.F90']
        objects = []
        for index, source in enumerate(sources):
            obj = mpi_build / f'mpi_{index}.o'
            precision = ['-D_TEST_SINGLE'] if args.single_precision else []
            subprocess.run([args.mpifc, '-cpp', '-D_TEST_REAL_QP_MODULE', '-D_TEST_MPI', *precision,
                            '-ffree-form', '-ffree-line-length-none',
                            '-I', str(build), '-fcheck=all', '-O0', '-c', str(source),
                            '-o', str(obj)], cwd=mpi_build, check=True)
            objects.append(obj)
        executable = mpi_build / ('mpi.exe' if os.name == 'nt' else 'mpi')
        subprocess.run([args.mpifc, *map(str, objects), '-o', str(executable)], cwd=mpi_build, check=True)
        for ranks in args.mpi_ranks:
            subprocess.run([args.mpiexec, '-n', str(ranks), str(executable)], cwd=mpi_build, check=True)
