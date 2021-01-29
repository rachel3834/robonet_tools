import numpy as np
import compare_xmatches
from pyDANDIA import metadata
from pyDANDIA import logs
from astropy.table import Table, Column

def test_crossmatch_pixel_positions():

    log = logs.start_stage_log( '.', 'test_compare_crossmatch' )

    meta1 = metadata.MetaData()
    nstars1 = 20000    # Must be greater than nstars2
    data = [Column(name='x', data=np.arange(0,nstars1,1)),
            Column(name='y', data=np.arange(0,nstars1,1))]
    setattr(meta1,'star_catalog',([],Table(data)))


    meta2 = metadata.MetaData()
    nstars2 = 1000
    data = [Column(name='x', data=np.arange(0.01,(nstars2+0.01),1)),
            Column(name='y', data=np.arange(0.01,(nstars2+0.01),1))]
    setattr(meta2,'star_catalog',([],Table(data)))

    threshold = 0.02

    # TEST 1: Find closest-matching stars
    (idx1, idx2) = compare_xmatches.crossmatch_pixel_positions(meta1, meta2, log)

    # As the number of stars in the two catalogs isn't necessarily the same,
    # but the pixel positions are a monotonically increasing sequence,
    # the nearest match to the later stars in the longer catalog will be the
    # last star in the shorter catalog.
    assert( idx1[0:nstars2] == idx2[0:nstars2] ).all()
    assert( idx2[nstars2:] == nstars2-1 ).all()

    # TEST 2: Require less than allowed separation for a match:
    (idx1, idx2) = compare_xmatches.crossmatch_pixel_positions(meta1, meta2, log,
                                                                threshold)
    assert( idx1[0:nstars2] == idx2[0:nstars2] ).all()
    assert(len(idx1) == nstars2)

    logs.close_log(log)

def test_calculate_separations_on_sky():

    log = logs.start_stage_log( '.', 'test_compare_crossmatch' )

    nstars = 10
    id_start = 405644951584364160

    data = [ Column(name='x1', data=np.arange(1.0,float(nstars+1),1.0)),
             Column(name='y1', data=np.arange(1.0,float(nstars+1),1.0)),
             Column(name='ra1', data=np.linspace(250.0, 265.0, nstars)),
             Column(name='dec1', data=np.linspace(-25.0, -28.0, nstars)),
             Column(name='gaia_id1', data=np.arange(id_start, id_start+nstars, 1)),

             Column(name='x2', data=np.arange(2.0,float(nstars+2),1.0)),
             Column(name='y2', data=np.arange(2.0,float(nstars+2),1.0)),
             Column(name='ra2', data=np.linspace(250.001, 265.001, nstars)),
             Column(name='dec2', data=np.linspace(-25.001, -28.001, nstars)),
             Column(name='gaia_id2', data=np.arange(id_start, id_start+nstars, 1)),

             Column(name='separation', data=np.zeros(nstars)) ]
    matched_data = Table(data)

    matched_data = compare_xmatches.calculate_separations_on_sky(matched_data, log)

    assert( (matched_data['separations'] != 0.0).all() )

    logs.close_log(log)

if __name__ == '__main__':
    #test_crossmatch_pixel_positions()
    test_calculate_separations_on_sky()
