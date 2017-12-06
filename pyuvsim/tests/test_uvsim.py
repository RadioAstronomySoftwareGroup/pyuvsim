import nose.tools as nt
import uvsim

def test_source():
       sim = uvsim.Source()
       attributes = ['ra', 'dec', 'polarization_angle', 'flux_calc_done', 'update_done', 'epoch']
       for a in attributes:
              nt.ok_(hasattr(sim,a))
