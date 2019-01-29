import  json
import  sfdmap
import  sqlite3
import  numpy               as  np
import  pylab               as  pl
import  matplotlib.pyplot   as  plt

from    astropy             import units as u
from    astropy.coordinates import SkyCoord

def radec2project(ra, dec):
    return (np.radians(ra) - np.pi, np.radians(dec))

##  https://arxiv.org/pdf/1809.01669.pdf
##  Appendix C1. 

##  https://www.lsst.org/scientists/simulations/opsim/opsim-v335-benchmark-surveys
##  Knutt Olsen LSST notebooks.


verbose       =   False

##  http://ops2.lsst.org/docs/current/architecture.html
conn          =   sqlite3.connect('/global/cscratch1/sd/mjwilson/LSST/minion_1016_sqlite.db')
c             =   conn.cursor()

if verbose:
  for Table in ['ObsHistory', 'Field', 'Proposal', 'ObsHistory_Proposal', 'Proposal_Field']:
    c.execute('PRAGMA TABLE_INFO({})'.format(Table))

    print([tup[1] for tup in c.fetchall()])

'''                                                                                                                                                         
##  Pg. 27 of https://github.com/LSSTScienceCollaborations/ObservingStrategy/blob/pdf/whitepaper/LSST_Observing_Strategy_White_Paper.pdf                   

(52, u'conf/survey/GalacticPlaneProp.conf', u'WL', 4367689744, u'minion', 1016)                 --                                                          
(53, u'conf/survey/SouthCelestialPole-18.conf', u'WL', 4367691152, u'minion', 1016)             --                                                          
(54, u'conf/survey/Universal-18-0824B.conf', u'WLTSS', 4367690832, u'minion', 1016)             --  WFD.                                                    
(55, u'conf/survey/NorthEclipticSpur-18c.conf', u'WLTSS', 4367690768, u'minion', 1016)          --                                                          
(56, u'conf/survey/DDcosmology1.conf', u'WLTSS', 4367691088, u'minion', 1016)                   --  Deep Drilling fields                                     
'''

##  A many-to-many relationship table that stores the fields fieldID from the Field table that maps to the field centers or specified by proposal propID.
c.execute('SELECT {coi1},{coi2} FROM {tn} WHERE {cn}=54'.format(coi1='Field_fieldID', coi2='proposal_field_id', tn='Proposal_Field', cn='Proposal_propID')) 

##  WFD:  Field_fieldID, proposal_field_id.
WFD_FieldID, WFD_PropFieldID = np.split(np.array(c.fetchall()), 2, 1)

##  This table maps visits to a field to the proposal or proposals that requested it.
c.execute('SELECT {coi1},{coi2} FROM {tn} WHERE {cn} = 54'.\
          format(coi1='ObsHistory_obsHistID', coi2='obsHistory_propID', tn='ObsHistory_Proposal', cn='Proposal_propID'))

##  (29, 21064471)
##  (30, 21064477)
_ObsHistID, _PropIDs = np.split(np.array(c.fetchall()), 2, 1)

##  Field properties.
c.execute('SELECT {coi1},{coi2},{coi3},{coi4} FROM {tn}'.\
          format(coi1='fieldID', coi2='fieldFov', coi3='fieldRA', coi4='fieldDec', tn='Field'))

##  Max field is 5292.
_FieldID, _FieldFOV, _FieldRa, _FieldDec = np.split(np.array(c.fetchall()), 4, 1)

##  This table keeps a record of each visit made by the telescope during a simulated survey.
##  Note:  all visits (2 x 15s exposures) where filter is u.
c.execute('SELECT {coi1},{coi2},{coi3},{coi4},{coi5} FROM {tn} WHERE {cn}="u"'.\
          format(coi1='Field_fieldID', coi2='expMJD', coi3='visitExpTime', coi4='fiveSigmaDepth', coi5='obsHistID', tn='ObsHistory', cn='filter'))

FieldIDs, MJDs, ExpTimes, FiveSigs, ObsHistID = np.split(np.array(c.fetchall()), 5, 1)
ObsHistID   = ObsHistID.astype(np.int)

print(MJDs)

assert np.all(np.diff(MJDs)       >= 0.)
assert np.all(np.diff(ObsHistID)  >= 0.)
assert np.all(np.diff(_ObsHistID) >= 0.)

##  Closing the connection to the database file                                                                                                                                     
conn.close()
  
WFDFIELDIDS = []

##  By default, a scaling of 0.86 is applied to the map values to reflect the recalibration by Schlafly & Finkbeiner (2011)                                   
##  https://github.com/kbarbary/sfdmap.                                                                                                                      
sfd         = sfdmap.SFDMap('/global/homes/m/mjwilson/SFD/')

YEAR        = 5
_ObsHistID  = _ObsHistID[:np.round(YEAR * 0.1 * len(_ObsHistID)).astype(np.int)]

for i, x in enumerate(ObsHistID):
  if x in _ObsHistID:
    WFDFIELDIDS.append(FieldIDs[i])

    print(100. * i / len(ObsHistID))
                                                                                     
for i, id in enumerate(_FieldID):                                                                                                                         
  if id in WFDFIELDIDS:
    ##  RA, DEC [degrees].
    ebv       = sfd.ebv(_FieldRa[i], _FieldDec[i])

    c         = SkyCoord(ra=_FieldRa[i] * u.degree, dec= _FieldDec[i] * u.degree, frame='icrs')
    cg        = c.galactic
    l, b      = cg.l, cg.b

    if (ebv > 0.3) | (np.abs(b) < 10.0 * u.degree):
      ##  Remove regions of high SFD extinction.  Remove the first matching value. 
      WFDFIELDIDS.remove(id)


WFDFIELDIDS    = np.array(WFDFIELDIDS).astype(np.int)
uWFDIDs, uVis  = np.unique(WFDFIELDIDS, return_counts=True)

##  Visits per WFD field.
hitmap         = dict(zip(uWFDIDs, uVis))

##  Current expected performance
single_m5      = {'u': 23.98, 'g': 24.91, 'r': 24.42, 'i': 23.97, 'z': 23.38, 'y': 22.47}

for i, x in enumerate(hitmap.keys()):  
  ##  Coadded 5 sig. depth =  One. Visit. 5. Sigma. Depth. + 2.5 * np.log10(np.sqrt(Vists. Per. Filter.))
  hitmap[x] = {'CoAdd5sig': single_m5['u'] + 2.5 * np.log10(np.sqrt(hitmap[x])), 'Visits': hitmap[x]}    

for i, y in enumerate(_FieldID.astype(np.int)):
  y   = y[0]

  if y in list(hitmap.keys()):
      hitmap[y] = {'CoAdd5sig': hitmap[y]['CoAdd5sig'], 'Visits': hitmap[y]['Visits'], 'RA': _FieldRa[i][0],\
                   'DEC': _FieldDec[i][0], 'FOV': _FieldFOV[i][0]}
                
with open('lsst_u_fiveyr.json', 'w') as outfile:
  json.dump(hitmap, outfile)

print('\n\nDone.\n\n')
