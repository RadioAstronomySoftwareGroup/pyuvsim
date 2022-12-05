# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2019 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Testing environment setup and teardown for pytest."""
import os
import pickle as pkl
import re
import warnings
from subprocess import (DEVNULL, CalledProcessError, TimeoutExpired,
                        check_output)

import pytest
from _pytest._code.code import ExceptionChainRepr
from _pytest.reports import TestReport
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.utils import iers
from pyuvdata import UVBeam
from pyuvdata.data import DATA_PATH

from pyuvsim import mpi
from pyuvsim.astropy_interface import MoonLocation, hasmoon

issubproc = os.environ.get('TEST_IN_PARALLEL', 0)
try:
    issubproc = bool(int(issubproc))
except ValueError:
    pass


def pytest_collection_modifyitems(session, config, items):
    # Enforce that the profiler test is run last.

    if len(items) <= 1:
        return
    for ii, it in enumerate(items):
        if 'profiler' in it.name:
            profiler_index = ii
            break
    items.append(items.pop(profiler_index))     # Move profiler tests to the end.


def pytest_configure(config):
    """Register an additional marker."""
    config.addinivalue_line(
        "markers",
        "parallel(n, timeout=70): mark test to run in n parallel mpi processes."
        " Optionally, set a timeout in seconds."
    )


def pytest_addoption(parser):
    parser.addoption("--nompi", action="store_true", help="skip mpi-parallelized tests.")


def pytest_runtest_setup(item):
    if 'parallel' in item.keywords:
        if mpi is None:
            pytest.skip("Need mpi4py to run parallelized tests.")
        elif item.config.getoption('nompi', False):
            pytest.skip("Skipping parallelized tests with --nompi option.")


@pytest.hookimpl(hookwrapper=True)
def pytest_exception_interact(node, call, report):
    if issubproc:
        from pyuvsim import mpi  # noqa
        if report.failed:
            pth = f"/tmp/mpitest_{report.head_line}"
            try:
                os.makedirs(pth)
            except OSError:
                pass
            with open(os.path.join(pth, f"report_rank{mpi.rank}.pkl"), 'wb') as ofile:
                pkl.dump(report, ofile)
            raise call.excinfo.value
    yield


def pytest_runtest_call(item):
    # If a test is to be run in parallel, spawn a subprocess that runs it in parallel.
    parmark = item.get_closest_marker("parallel")
    if parmark is None:
        return
    nproc = 1
    if len(parmark.args) >= 1:
        try:
            nproc = int(parmark.args[0])
        except ValueError:
            raise ValueError(f"Invalid number of processes: {parmark.args[0]}")

    timeout = 70
    if 'timeout' in parmark.kwargs:
        timeout = float(parmark.kwargs['timeout'])

    call = ['mpirun', '--host', 'localhost:10', '-n', str(nproc),
            "python", "-m", "pytest", "{:s}::{:s}".format(str(item.fspath), str(item.name))]
    if not issubproc:
        try:
            envcopy = os.environ.copy()
            envcopy['TEST_IN_PARALLEL'] = '1'
            check_output(call, env=envcopy, stderr=DEVNULL, timeout=timeout)
        except (TimeoutExpired, CalledProcessError) as err:
            message = f"Parallelized test {item.name} failed"
            if isinstance(err, TimeoutExpired):
                message += (f" after {timeout} seconds due to timeout. \nThe timeout may be set"
                            " via the ``timeout`` keyword in the ``parallel`` decorator. \n"
                            "A stalled test may be caused by an inconsistent MPI state. Check"
                            " that blocking operations are reached by all processes")
            message += "."
            raise AssertionError(message) from err

        # If passing, do not run after this function.
        def blank_func(*args, **kwargs):
            return
        item.obj = blank_func


class _MPIExceptionChainRepr(ExceptionChainRepr):
    """
    Extension of _pytest._code.code.ExceptionChainRepr.

    Changes the toterminal method to better represent exceptions
    coming from multiple MPI processes.
    """

    mpi_ranks = []

    def __init__(self, exception_chain_repr, mpi_ranks):
        if isinstance(mpi_ranks, int):
            mpi_ranks = [mpi_ranks]
        self.mpi_ranks.extend(mpi_ranks)
        super().__init__(exception_chain_repr.chain)

    def append(self, exception_chain_repr, mpi_rank):
        self.mpi_ranks.append(mpi_rank)
        self.chain.append(exception_chain_repr.chain[0])

    def toterminal(self, tw):
        for ii, element in enumerate(self.chain):
            rank = self.mpi_ranks[ii]
            tw.line("")
            tw.sep('— ', f"MPI Rank = {rank}", cyan=True)
            element[0].toterminal(tw)
            if element[2] is not None:
                tw.line("")
                tw.line(element[2], yellow=True)
        super(ExceptionChainRepr, self).toterminal(tw)


def pytest_runtest_makereport(item, call):
    report = TestReport.from_item_and_call(item, call)
    parmark = item.get_closest_marker("parallel")

    # Is a parallel test but not currently within the subprocess.
    if report.failed and (not issubproc) and (parmark is not None) and (report.when == 'call'):
        pth = f"/tmp/mpitest_{report.head_line}"
        if os.path.exists(pth):
            for ii, repf in enumerate(os.listdir(pth)):
                if not repf.endswith(".pkl"):
                    continue
                rank = int(re.findall('(?<=rank)[0-9]+', repf)[0])
                rank_report = pkl.load(open(os.path.join(pth, repf), 'rb'))
                os.remove(os.path.join(pth, repf))
                if ii == 0:
                    report = rank_report
                    report.longrepr = _MPIExceptionChainRepr(report.longrepr, rank)
                else:
                    report.longrepr.append(rank_report.longrepr, rank)
            os.rmdir(pth)
        nproc = int(parmark.args[0])
        outcome = call.excinfo.exconly()
        mpitest_info = f"Run with {nproc} processes."
        mpitest_info += f"\nFailed with message: {outcome}"
        report.longrepr.addsection("MPI Info", mpitest_info)
    return report


@pytest.fixture(autouse=True, scope="session")
def _setup_and_teardown_package():
    # Do a calculation that requires a current IERS table. This will trigger
    # automatic downloading of the IERS table if needed, including trying the
    # mirror site in python 3 (but won't redownload if a current one exists).
    # If there's not a current IERS table and it can't be downloaded, turn off
    # auto downloading for the tests and turn it back on once all tests are
    # completed (done by extending auto_max_age).
    try:
        t1 = Time.now()
        t1.ut1
    except (Exception):
        iers.conf.auto_max_age = None

    # yield to allow tests to run
    yield

    iers.conf.auto_max_age = 30


@pytest.fixture(autouse=True)
def _ignore_deprecation():
    warnings.filterwarnings('ignore', message='Achromatic gaussian beams will not be supported',
                            category=PendingDeprecationWarning)
    warnings.filterwarnings('ignore', message='"initialize_catalog_from_params will not return'
                                              ' recarray by default in the future.',
                            category=PendingDeprecationWarning)


@pytest.fixture(scope='session')
def cst_beam():
    beam = UVBeam()

    freqs = [150e6, 123e6]

    cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
    beam_files = [os.path.join(DATA_PATH, 'NicCSTbeams', f) for f in cst_files]
    beam.read_cst_beam(
        beam_files, beam_type='efield', frequency=freqs,
        telescope_name='HERA', feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
        model_name='E-field pattern - Rigging height 4.9m', model_version='1.0'
    )
    beam.x_orientation = 'east'
    beam.freq_interp_kind = 'linear'
    beam.interpolation_function = 'az_za_simple'
    beam.peak_normalize()
    return beam


@pytest.fixture(scope='session')
def hera_loc():
    return EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)


@pytest.fixture(scope='session')
def apollo_loc():
    if hasmoon:
        return MoonLocation(lat=0.6875, lon=24.433, height=0)
    else:
        return None
