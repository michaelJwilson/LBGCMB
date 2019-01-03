'''
desisim.specsim
===============

DESI wrapper functions for external specsim classes.
'''

from    __future__  import absolute_import, division, print_function
import  desiutil.log

## Cached simulators, keyed by config string
_simulators  = dict()

## Cached defaults after loading a new simulator, to be used to reset a
## simulator back to a reference state before returning it as a cached copy
_simdefaults = dict()

log          = desiutil.log.get_logger()

def get_simulator(config='desi', num_fibers=1):
    '''
    Returns new or cached specsim.simulator.Simulator object;
    Also adds placeholder for BGS fiberloss if that isn't already in the config.
    '''
    
    key = (config, num_fibers)
    
    if key in _simulators:
        log.debug('Returning cached {} simulator'.format(config))
        
        qsim                                  = _simulators[key]
        defaults                              = _simdefaults[key]

        qsim.source.focal_xy                  = defaults['focal_xy']
                     
        qsim.atmosphere.airmass               = defaults['airmass']
        qsim.atmosphere.moon.moon_phase       = defaults['moon_phase']
        qsim.atmosphere.moon.separation_angle = defaults['moon_angle']
        qsim.atmosphere.moon.moon_zenith      = defaults['moon_zenith']

        qsim.observation.exposure_time        = defaults['exposure_time']

    else:
        log.debug('Creating new {} simulator'.format(config))

        ##  New config; create Simulator object
        import  specsim.simulator

        qsim = specsim.simulator.Simulator(config, num_fibers)

        ## TODO FIXME HACK: desimodel/specsim doesn't have BGS fiberloss yet,
        ## so scale from LRG
        if 'bgs' not in qsim.instrument.fiber_acceptance_dict.keys():
            log.warning('Treating BGS fiberloss = 0.5 * LRG fiberloss')

            qsim.instrument.fiber_acceptance_dict['bgs'] = 0.5 * qsim.instrument.fiber_acceptance_dict['lrg']

        ## Cache defaults to reset back to original state later.
        defaults                  = dict()
        
        defaults['focal_xy']      = qsim.source.focal_xy

        defaults['airmass']       = qsim.atmosphere.airmass
        defaults['exposure_time'] = qsim.observation.exposure_time

        defaults['moon_phase']    = qsim.atmosphere.moon.moon_phase
        defaults['moon_angle']    = qsim.atmosphere.moon.separation_angle
        defaults['moon_zenith']   = qsim.atmosphere.moon.moon_zenith

        _simulators[key]          = qsim
        _simdefaults[key]         = defaults

    return qsim

