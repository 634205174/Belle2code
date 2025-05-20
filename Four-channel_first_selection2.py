#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################
#  This script performs the study of: 
#   B0 -> J/psi K_S0
#   B+ -> J/psi K+
#  Based on Belle's data, using b2bii module
#  and Belle II data 
#  Including both J/psi -> mu+mu- and J/psi -> e+e- decay modes
# 
##########################################

import sys
import os
import ROOT
from ROOT import Belle2

import basf2 as b2
import modularAnalysis as ma
import variables as v
import variables.collections as vc
import variables.utils as vu
import vertex as vt
import flavorTagger
import stdV0s

import b2biiConversion as b2c
import b2biiMonitors as b2m

from basf2 import B2ERROR, B2WARNING
from b2bii import isB2BII
from modularAnalysis import cutAndCopyList, reconstructDecay, applyCuts, copyLists
from vertex import treeFit, kFit

from stdCharged import stdPi, stdK, stdE, stdMu
from stdV0s import stdLambdas
from stdPhotons import stdPhotons
from stdPi0s import stdPi0s
from stdHyperons import stdXi, stdXi0, stdOmega, goodXi, goodXi0, goodOmega

from variables.MCGenTopo import mc_gen_topo
from pdg import add_particle
from basf2 import statistics
from basf2 import set_nprocesses
#===============Begin analysis script===============================#
main = b2.create_path()
set_nprocesses(4)
#===============initial setting and read data===begin==========#
if len(sys.argv) != 4:
   sys.exit('Must provide following input parameters: [mc/data] [input_Belle_MDST_file] [output_BelleII_ROOT_file_dir] [isbsub].\n'
            'Check for more Belle MC or Data at http://bweb3.cc.kek.jp/ \n'
            'Example: http://bweb3/montecarlo.php?ex=33&rs=1&re=9999&ty=evtgen-mixed&dt=on_resonance&bl=caseB&st=Any')

mc_or_data = sys.argv[1].lower()
isMC = {'mc':True, 'data':False}.get(mc_or_data, None)
if isMC is None:
   sys.exit('First parameter must be "mc" or "data" to indicate whether MC or real data')

outputBelle2ROOTFile_dir = sys.argv[2]
isbsub = sys.argv[3]

if isbsub == "signal":
        outputFile = outputBelle2ROOTFile_dir
elif isbsub == "bsub":
        outputFile = outputBelle2ROOTFile_dir
else:
        outputFile = outputBelle2ROOTFile_dir + '/test.root'

ma.inputMdstList(['/gpfs/group/belle2/users2022/zepengxu/run/physics_analysis/Beta_angle/MC_gen/Jpsi_Ks0_ee/mdst0001.root'],environmentType='default',path=main)

#===============initial setting and read data===end==========#
b2.conditions.append_globaltag(ma.getAnalysisGlobaltag())
b2.conditions.prepend_globaltag("pid_nn_release08_v1")
b2.conditions.prepend_globaltag('user_liuyux_KsFinder_dev')
#===============some explanation=============================#
"""
Study of B0 -> J/psi K_S0 and B+ -> J/psi K+
J/psi-> mu+mu- and J/psi -> e+e-
For isospin asymmetry measurement
"""

#==============Begin: set globaltag=============================#
#==============End: set globaltag=============================#

#===============Begin: event selections and output==========#
if isMC:
    mc_vars = vc.mc_truth+ vc.mc_variables + vc.pid
else:
    mc_vars =  []
# Defined standard list for MC truth variables including isSignal
mc_truth_vars = ['isSignal', 'mcErrors', 'mcPDG', 'genMotherPDG'] if isMC else []
def add_ks_selection_aliases():
    """创建KS选择所需的变量别名"""
    v.variables.addAlias('M_lambda_p', 'useAlternativeDaughterHypothesis(M, 0:p+)')
    v.variables.addAlias('M_lambda_antip', 'useAlternativeDaughterHypothesis(M, 1:anti-p-)')
    v.variables.addAlias('daughtersDeltaZ', 'daughterDiffOf(0, 1, dz)')
    v.variables.addAlias('cosVertexMomentum', 'cosAngleBetweenMomentumAndVertexVector')
    v.variables.addAlias('pip_nPXDHits', 'daughter(0,nPXDHits)')
    v.variables.addAlias('pin_nPXDHits', 'daughter(1,nPXDHits)')
    v.variables.addAlias('pip_nSVDHits', 'daughter(0,nSVDHits)')
    v.variables.addAlias('pin_nSVDHits', 'daughter(1,nSVDHits)')
    v.variables.addAlias('daughterAngleDiffInMother', 'useRestFrame(daughterAngle(0, 1))')
    v.variables.addAlias('pip_p', 'daughter(0,p)')
    v.variables.addAlias('pin_p', 'daughter(1,p)')
    v.variables.addAlias('pip_dr', 'daughter(0,dr)')
    v.variables.addAlias('pin_dr', 'daughter(1,dr)')
    v.variables.addAlias('pip_cosTheta', 'daughter(0,cosTheta)')
    v.variables.addAlias('pin_cosTheta', 'daughter(1,cosTheta)')
    v.variables.addAlias('pip_protonID', 'daughter(0,protonID)')
    v.variables.addAlias('pin_protonID', 'daughter(1,protonID)')
    v.variables.addAlias('pipPDG', 'daughter(0,mcPDG)')
    v.variables.addAlias('pinPDG', 'daughter(1,mcPDG)')

# 定义KS选择的变量集合
def add_ks_variable_collection():
    """添加KS选择器变量集合"""
    add_ks_selection_aliases()
    inputVariablesList = [
        'cosVertexMomentum',
        'flightDistance',
        'significanceOfDistance',
        'cosHelicityAngleMomentum',
        'ImpactXY',
        'decayAngle(0)',
        'decayAngle(1)',
        'daughterAngleDiffInMother',
        'daughtersDeltaZ',
        'pip_nSVDHits', 'pip_nPXDHits',
        'pin_nSVDHits', 'pin_nPXDHits',
        'pip_dr', 'pin_dr',
        'pip_protonID', 'pin_protonID',
        'M_lambda_p', 'M_lambda_antip',
        'pip_p', 'pin_p',
        'pip_cosTheta', 'pin_cosTheta',
        'ArmenterosLongitudinalMomentumAsymmetry',
        'pinPDG', 'pipPDG',
    ]
    vu.add_collection(inputVariablesList, 'ks_selector_info')
#===============Begin: generator level match================#
v.variables.addAlias("decaymode", "extraInfo(decayModeID)")
if isMC:
    # Fill basic MC particle lists
    ma.fillParticleListsFromMC([('gamma:gen', ''),('e+:gen', ''),('mu+:gen', ''),('pi+:gen', ''),('K+:gen', '')],path=main)
    ma.fillParticleListsFromMC([('K_S0:gen','')],skipNonPrimary=True,addDaughters=True,skipNonPrimaryDaughters=False,path=main)
    
    # Find both J/psi decay modes in MC
    ma.findMCDecay('J/psi:genmm','J/psi -> mu+ mu-',path=main)
    ma.findMCDecay('J/psi:genee','J/psi -> e+ e-',path=main)
    ma.findMCDecay('K_S0:gen_2pi',  'K_S0 -> pi+ pi-',path=main)
    
    # Fill particle lists for all decay products
    ma.fillParticleListsFromMC([
        ('pi+:gen2pi', 'isDescendantOfList(K_S0:gen_2pi)'),
        ('mu+:gen2mu', 'isDescendantOfList(J/psi:genmm)'),
        ('e+:gen2e', 'isDescendantOfList(J/psi:genee)')
    ],path=main)
    
    # Reconstruct MC J/psi from both decay modes
    ma.reconstructMCDecay('J/psi:gen2mm -> mu+:gen2mu mu-:gen2mu','',dmID=1,path=main)
    ma.reconstructMCDecay('J/psi:gen2ee -> e+:gen2e e-:gen2e','',dmID=2,path=main)
    
    # Combine both J/psi decay modes
    copyLists('J/psi:gen',[
        'J/psi:gen2mm',
        'J/psi:gen2ee'
    ],path=main)
    
    # Reconstruct K_S0 and B0 at MC level
    #ma.fillParticleListFromMC('K_S0:gen2pi -> pi+:gen2pi pi-:gen2pi ','',dmID=1,path=main)
    ma.fillParticleListFromMC('K_S0:gen', cut='',skipNonPrimary=True, addDaughters=True, skipNonPrimaryDaughters=True, path=main)
    #copyLists('K_S0:gen',['K_S0:genpi'],path=main)
    
    # Reconstruct B0 MC decay
    ma.reconstructMCDecay('B0:genJ/psi -> J/psi:gen K_S0:gen','',dmID=1,path=main)
    copyLists('B0:gen',['B0:genJ/psi'],path=main)
    
    # Reconstruct B+ MC decay
    ma.reconstructMCDecay('B+:genJ/psi -> J/psi:gen K+:gen','',dmID=1,path=main)
    copyLists('B+:gen',['B+:genJ/psi'],path=main)
    
    # Set event extra info for B0
    v.variables.addAlias('JpsiDM',    'daughter(0,decaymode)')
    ma.variablesToEventExtraInfo('B0:gen',{
        'isInList(B0:gen)':'isGenB0',
        'charge':'genB0charge',
        'decaymode':'genB0DM',
        'JpsiDM':'genB0JpsiDM'
    },option=1, path=main)
    
    # Set event extra info for B+
    ma.variablesToEventExtraInfo('B+:gen',{
        'isInList(B+:gen)':'isGenBp',
        'charge':'genBpcharge',
        'decaymode':'genBpDM',
        'JpsiDM':'genBpJpsiDM'
    },option=1, path=main)
    
    # Define aliases
    v.variables.addAlias('isgenB0',   'eventExtraInfo(isGenB0)')
    v.variables.addAlias('genB0c',    'eventExtraInfo(genB0charge)')
    v.variables.addAlias('genB0',     'eventExtraInfo(genB0DM)')
    v.variables.addAlias('genB0Jpsi', 'eventExtraInfo(genB0JpsiDM)')
    v.variables.addAlias('ngenB0',    'nParticlesInList(B0:gen)')
    
    v.variables.addAlias('isgenBp',   'eventExtraInfo(isGenBp)')
    v.variables.addAlias('genBpc',    'eventExtraInfo(genBpcharge)')
    v.variables.addAlias('genBp',     'eventExtraInfo(genBpDM)')
    v.variables.addAlias('genBpJpsi', 'eventExtraInfo(genBpJpsiDM)')
    v.variables.addAlias('ngenBp',    'nParticlesInList(B+:gen)')
    
    gen_inform = ['isgenB0','genB0c','genB0','genB0Jpsi','ngenB0',
                  'isgenBp','genBpc','genBp','genBpJpsi','ngenBp']
    
    mc_detailed_vars = [
    # B kinematics
    'useCMSFrame(p)', 'useCMSFrame(pt)', 'useCMSFrame(pz)', 'useCMSFrame(E)',

    # J/psi kinematics
    'daughter(0,useCMSFrame(p))', 'daughter(0,useCMSFrame(pt))',
    'daughter(0,useCMSFrame(pz))', 'daughter(0,useCMSFrame(E))',

    # K_S0/K+ kinematics
    'daughter(1,useCMSFrame(p))', 'daughter(1,useCMSFrame(pt))',
    'daughter(1,useCMSFrame(pz))', 'daughter(1,useCMSFrame(E))',

    # Decay vertices and angles
    'mcFlightDistance', 'mcDecayVertexX', 'mcDecayVertexY', 'mcDecayVertexZ',
    'daughter(1,mcFlightDistance)', 'daughterAngle(0,1)'
    ]

    # B0 4-momentum
    v.variables.addAlias('b0_E', 'E')
    v.variables.addAlias('b0_p', 'p')
    v.variables.addAlias('b0_px', 'px')
    v.variables.addAlias('b0_py', 'py')
    v.variables.addAlias('b0_pz', 'pz')
    
    # B0 4-momentum in CMS frame
    v.variables.addAlias('b0_cms_E', 'useCMSFrame(E)')
    v.variables.addAlias('b0_cms_p', 'useCMSFrame(p)')
    v.variables.addAlias('b0_cms_px', 'useCMSFrame(px)')
    v.variables.addAlias('b0_cms_py', 'useCMSFrame(py)')
    v.variables.addAlias('b0_cms_pz', 'useCMSFrame(pz)')
    

    gen_vars = vc.deltae_mbc + mc_truth_vars + mc_detailed_vars + [
    # Basic properties
    'charge', 'decaymode', 'JpsiDM',
    
    # B 4-momentum
    'b0_E', 'b0_p', 'b0_px', 'b0_py', 'b0_pz',
    'b0_cms_E', 'b0_cms_p', 'b0_cms_px', 'b0_cms_py', 'b0_cms_pz',
    
    # Additional useful variables
    'mcFlightDistance', 'daughter(1, mcFlightDistance)',
    'daughterAngle(0,1)'  # Angle between J/psi and K_S0/K+
    ] #+  ['PDG', 'genmotherPDG', 'isSignal']
    
    # Add these to final gen_vars
    ma.variablesToNtuple("B0:gen", gen_vars, filename=outputFile, treename='genB0', path=main)
    ma.variablesToNtuple("B+:gen", gen_vars, filename=outputFile, treename='genBp', path=main)
else:
    gen_inform = []
#===============End:generator level match==================#

#===============Begin: event level variables======================#
event_shape_vars = ['foxWolframR1','foxWolframR2','foxWolframR3','foxWolframR4','thrust','thrustAxisCosTheta']
ma.applyEventCuts('[nTracks >= 4]',path=main)
ma.fillParticleList(decayString='pi+:goodtracks',cut='pt> 0.1',path=main)
ma.fillParticleList(decayString='mu+:goodtracks',cut='pt> 0.1',path=main)
ma.fillParticleList(decayString='e+:goodtracks',cut='pt> 0.1',path=main)
ma.buildEventShape(['pi+:goodtracks', 'mu+:goodtracks', 'e+:goodtracks'],path=main)

# Add variable for continuum suppression
v.variables.addAlias('r2', 'foxWolframR2')

evtLel = vc.reco_stats + ['Ecms','eventRandom','isMC'] + gen_inform + event_shape_vars + ['r2']
#===============End: event level variables======================#

#===============Begin: event selections=====================#

#==============Begin: tracks selection=====================#



b2.conditions.prepend_globaltag("pid_nn_release08_v1_expert")
#expert version for PID, it can compare performance between different PID version
payloadName = 'pid_nn_release08_v1_expert'
for p, pdg in zip(['pion', 'kaon', 'p', 'electron', 'muon', 'deuteron'],
                  ['211',  '321', '2212', '11',     '13',   '1000010020']):
    v.variables.addAlias(f'{p}IDNN_my', f'pidNeuralNetworkValueExpert({pdg},{payloadName})')
kaon_pid  = 'kaonIDNN > 0.01'
e_pid  = 'electronIDNN_my > 0.01'
mu_pid  = 'muonIDNN_my > 0.01'
#mu_pid = 'muonID_noSVD>0.99'
#e_pid = 'electronID_noTOP > 0.9'

mu_mom_cut = 'p > 0.1'
e_mom_cut = 'p > 0.1'

goodtrack = 'dr < 0.5 and abs(dz) < 2 and thetaInCDCAcceptance'

# Create standard particle lists
stdPi('all', path=main)
stdK('all', path=main)
stdMu('all', path=main)
stdE('all', path=main)


# Match the base lists BEFORE saving ntuples and applying cuts
if isMC:
    ma.matchMCTruth('pi+:all', path=main)
    ma.matchMCTruth('K+:all', path=main)
    ma.matchMCTruth('mu+:all', path=main)
    ma.matchMCTruth('e+:all', path=main)

# -- Save Track Ntuples (CHANGED to save :all lists with cut variables) --
# Save ntuples from the ':all' lists, including variables needed for cuts + isSignal
track_base_vars = ['p', 'pt', 'theta', 'phi', 'charge', 'dr', 'dz', 'nCDCHits', 'thetaInCDCAcceptance']
track_pid_vars = ['pionIDNN_my', 'kaonIDNN_my', 'muonIDNN_my', 'electronIDNN_my']
track_vars_to_save = track_base_vars + track_pid_vars + mc_truth_vars 
# Make sure mc_truth_vars (containing isSignal) is added
#ma.variablesToNtuple('pi+:all',track_vars_to_save, filename=outputFile, treename='cut_pions_all', path=main)
#ma.variablesToNtuple('K+:all', track_vars_to_save, filename=outputFile, treename='cut_kaons_all', path=main)
#ma.variablesToNtuple('mu+:all', track_vars_to_save, filename=outputFile, treename='cut_muons_all', path=main)
#ma.variablesToNtuple('e+:all', track_vars_to_save, filename=outputFile, treename='cut_electrons_all', path=main)

# -- Apply full track selections (Defined combined cut strings) --
pion_full_cut = f'{goodtrack}'
kaon_full_cut = f'{goodtrack} and {kaon_pid}'
muon_full_cut = f'{goodtrack} and {mu_pid} and {mu_mom_cut}'
electron_full_cut = f'{goodtrack} and {e_pid} and {e_mom_cut}'


# Apply cuts to create good track lists
cutAndCopyList('pi+:mypi', 'pi+:all', pion_full_cut, path=main)
cutAndCopyList('K+:myK', 'K+:all', kaon_full_cut, path=main)
cutAndCopyList('mu+:mymu', 'mu+:all', muon_full_cut, path=main)
cutAndCopyList('e+:myeuncorrected', 'e+:all', electron_full_cut, path=main)
trk_vars = mc_truth_vars + ['pt', 'p', 'E', 'charge', 'PDG']
#==============End: tracks selection=====================#
# apply Bremsstrahlung correction to electrons
v.variables.addAlias(
    "goodFWDGamma", "passesCut(clusterReg == 1 and clusterE > 0.075)"
)
v.variables.addAlias(
    "goodBRLGamma", "passesCut(clusterReg == 2 and clusterE > 0.05)"
)
v.variables.addAlias(
    "goodBWDGamma", "passesCut(clusterReg == 3 and clusterE > 0.1)"
)
v.variables.addAlias(
    "goodGamma", "passesCut(goodFWDGamma or goodBRLGamma or goodBWDGamma)"
)
ma.fillParticleList("gamma:brems", "goodGamma", path=main)
ma.correctBrems("e+:mye", "e+:myeuncorrected", "gamma:brems", path=main)
v.variables.addAlias("isBremsCorrected", "extraInfo(bremsCorrected)")

if isMC:
    # Match J/psi candidates
    ma.matchMCTruth("K+:myK", path=main)

#==============Begin: variables for composite particles=====================#
#variables for basic kinematics
v.variables.addAlias('cmsp', 'useCMSFrame(p)')
v.variables.addAlias('daughterAngle','daughterAngle(0,1)')
v.variables.addAlias('decayAngle','decayAngle(0)')
v.variables.addAlias('decayAngle1','decayAngle(1)')
kinematics_vars = ['InvM','M','decaymode','charge','p','cmsp','pt','daughterAngle','decayAngle','decayAngle1']
kinematics_vars += mc_truth_vars 

#PID variable
#PID_vars = ['kaonIDNN', 'muonIDNN_my', 'electronIDNN_my', 'p']

#v.variables.addAlias('jpsi_mu_plus_ID', 'daughter(0, daughter(0, muonIDNN_my))') 
#v.variables.addAlias('jpsi_mu_minus_ID', 'daughter(0, daughter(1, muonIDNN_my))') 
#v.variables.addAlias('jpsi_e_plus_ID', 'daughter(0, daughter(0, electronIDNN_my))')  
#v.variables.addAlias('jpsi_e_minus_ID', 'daughter(0, daughter(1, electronIDNN_my))')
#
#
#v.variables.addAlias('jpsi_mu_plus_p', 'daughter(0, daughter(0, p))')
#v.variables.addAlias('jpsi_mu_minus_p', 'daughter(0, daughter(1, p))')
#v.variables.addAlias('jpsi_e_plus_p', 'daughter(0, daughter(0, p))')
#v.variables.addAlias('jpsi_e_minus_p', 'daughter(0, daughter(1, p))')

v.variables.addAlias('jpsi_mu_plus_ID', 'daughter(0, muonIDNN_my)')
v.variables.addAlias('jpsi_mu_minus_ID', 'daughter(1, muonIDNN_my)')
v.variables.addAlias('jpsi_e_plus_ID', 'daughter(0, electronIDNN_my)')
v.variables.addAlias('jpsi_e_minus_ID', 'daughter(1, electronIDNN_my)')


v.variables.addAlias('jpsi_mu_plus_p', 'daughter(0, p)')
v.variables.addAlias('jpsi_mu_minus_p', 'daughter(1, p)')
v.variables.addAlias('jpsi_e_plus_p', 'daughter(0, p)')
v.variables.addAlias('jpsi_e_minus_p', 'daughter(1, p)')

v.variables.addAlias('ks_pi_plus_pionID', 'daughter(0, pionIDNN_my)')
v.variables.addAlias('ks_pi_minus_pionID', 'daughter(1, pionIDNN_my)')

v.variables.addAlias('ks_pi_plus_p', 'daughter(0, p)')
v.variables.addAlias('ks_pi_minus_p', 'daughter(1, p)')
ks_pion_vars = [
    'ks_pi_plus_pionID', 'ks_pi_minus_pionID',
    'ks_pi_plus_p', 'ks_pi_minus_p',
]


#K+
v.variables.addAlias('k_plus_ID', 'daughter(1, kaonIDNN)')
v.variables.addAlias('k_plus_p', 'daughter(1, p)')

#Jpsi daughter
Jpsi_daughter_alias_vars = ['jpsi_mu_plus_ID', 'jpsi_mu_minus_ID', 'jpsi_e_plus_ID', 'jpsi_e_minus_ID', 'k_plus_ID',
                  'jpsi_mu_plus_p', 'jpsi_mu_minus_p', 'jpsi_e_plus_p', 'jpsi_e_minus_p', 'k_plus_p']

#variables before mass fit
extraM={'M':'extraM','InvM':'extraInvM'}
v.variables.addAlias("M0", "extraInfo(extraM)")
v.variables.addAlias("InvM0", "extraInfo(extraInvM)")
extraM_vars = ['M0','InvM0']
kinematics_vars += extraM_vars

extraB={'Mbc':'extraMbc','deltaE':'extraDeltaE'}
v.variables.addAlias("Mbc0", "extraInfo(extraMbc)")
v.variables.addAlias("deltaE0", "extraInfo(extraDeltaE)")
extraB_vars = ['Mbc0','deltaE0']

#variables for fit
v.variables.addAlias('chisq','extraInfo(chiSquared)')
v.variables.addAlias('chisq2ndf','formula(extraInfo(chiSquared)/extraInfo(ndf))')
fit_vars = ['chisq','chiProb','chisq2ndf']
#variables for vertex
v.variables.addAlias('distSign','formula(flightDistance/flightDistanceErr)')
v.variables.addAlias('cosvp','cosAngleBetweenMomentumAndVertexVector')
v.variables.addAlias('cosvpxy','cosAngleBetweenMomentumAndVertexVectorInXYPlane')
vertex_vars = ['distance','distSign','flightDistance','flightDistanceErr','cosvp','cosHelicityAngleMomentum']

#variables for B meson
v.variables.addAlias('cosHeli','cosHelicityAngle(0,0)')
v.variables.addAlias('cosHeliBeam','cosHelicityAngleBeamMomentum(0)')

# Add J/psi decay mode variable (1=mu+mu-, 2=e+e-)
v.variables.addAlias('jpsi_decaymode', 'daughter(0, decaymode)')

bplus_vars = vc.deltae_mbc + extraB_vars +['cosHeli','cosHeliBeam','cosToThrustOfEvent','jpsi_decaymode']



v.variables.addAlias('jpsi_mm_mass', 'ifInfo(decaymode,1,daughter(0, M),0)')
v.variables.addAlias('jpsi_ee_mass', 'ifInfo(decaymode,2,daughter(0, M),0)')
#==============End: variables for composite particles=====================#

#FOM condition
# event
v.variables.addAlias('evt_nTracks', 'nTracks')

# track
v.variables.addAlias('trk_dr', 'dr')
v.variables.addAlias('trk_dz', 'abs(dz)')
v.variables.addAlias('trk_inCDC', 'thetaInCDCAcceptance')

# PID

# J/psi
v.variables.addAlias('jpsi_M', 'daughter(0, M)')
v.variables.addAlias('jpsi_decaymode', 'daughter(0, decaymode)')
v.variables.addAlias('jpsi_chiProb', 'daughter(0, chiProb)')

# K_S0
v.variables.addAlias('ks_M', 'daughter(1, M)')
v.variables.addAlias('ks_flightDist', 'daughter(1, flightDistance)')
v.variables.addAlias('ks_distSign', 'daughter(1, distSign)')
v.variables.addAlias('ks_chiProb', 'daughter(1, chiProb)')

# B meson
v.variables.addAlias('b_chiProb', 'chiProb')
v.variables.addAlias('b_chisq2ndf', 'chisq2ndf')

# ntuple
event_cuts = ['evt_nTracks', 'r2']
track_cuts = ['trk_dr', 'trk_dz', 'trk_inCDC']
pid_cuts = ['kaonIDNN', 'muonIDNN_my', 'electronIDNN_my'] + Jpsi_daughter_alias_vars
jpsi_cuts = ['jpsi_M', 'jpsi_decaymode', 'jpsi_chiProb']
ks_cuts = ['ks_M', 'ks_flightDist', 'ks_distSign', 'ks_chiProb']
b_cuts = ['Mbc', 'deltaE', 'b_chiProb', 'b_chisq2ndf']

# save cut
all_cuts = event_cuts + b_cuts + jpsi_cuts + ks_cuts + pid_cuts

#========================#

# 
cutAndCopyList('pi+:loosepi', 'pi+:mypi', '', path=main)
cutAndCopyList('K+:looseK', 'K+:myK', '', path=main)
cutAndCopyList('mu+:loosemu', 'mu+:mymu', '', path=main)
cutAndCopyList('e+:loosee', 'e+:mye', '', path=main)

# Defined variable list including fit vars and mc_truth_vars
jpsi_base_vars = ['InvM', 'M', 'p', 'pt', 'E', 'charge', 'decaymode']
jpsi_fit_vars = ['chiProb', 'chisq2ndf'] # Aliases defined below
jpsi_vars_to_save = jpsi_base_vars + jpsi_fit_vars + Jpsi_daughter_alias_vars + mc_truth_vars 


# Loose J/psi
reconstructDecay('J/psi:looseJpsimm -> mu+:loosemu mu-:loosemu', '2.8 < M < 3.4', path=main)
reconstructDecay('J/psi:looseJpsiee -> e+:loosee e-:loosee ?addbrems', '2.8 < M < 3.4', path=main)
copyLists('J/psi:looseJpsi', ['J/psi:looseJpsimm', 'J/psi:looseJpsiee'], path=main)

# -- Apply Vertex Fit to LOOSE Candidates 
vt.treeFit('J/psi:looseJpsimm', 0.0, path=main)
vt.treeFit('J/psi:looseJpsiee', 0.0, path=main)

# -- Define Aliases for Fit Variables --
v.variables.addAlias('chisq','extraInfo(chiSquared)')
v.variables.addAlias('chisq2ndf','formula(extraInfo(chiSquared)/extraInfo(ndf))')

# -- MC Match LOOSE J/psi (Ensured it's before ntuple) --
if isMC:
    ma.matchMCTruth('J/psi:looseJpsimm', path=main)
    ma.matchMCTruth('J/psi:looseJpsiee', path=main)

# -- Save LOOSE J/psi Ntuples (CHANGED to include fit vars, renamed tree) --
ma.variablesToNtuple('J/psi:looseJpsimm', jpsi_vars_to_save,
                     filename=outputFile, treename='cut_jpsi_mm_loose', path=main)
ma.variablesToNtuple('J/psi:looseJpsiee', jpsi_vars_to_save,
                     filename=outputFile, treename='cut_jpsi_ee_loose', path=main)


stdV0s.stdKshorts(path=main)
# Loose K_S0
cutAndCopyList('K_S0:looseKs', 'K_S0:merged', '0.45 < M < 0.55', path=main)

# Loose B meson
reconstructDecay('B0:loose -> J/psi:looseJpsi K_S0:looseKs', 'Mbc > 5.2 and abs(deltaE) < 0.3', path=main)
reconstructDecay('B+:loose -> J/psi:looseJpsi K+:looseK', 'Mbc > 5.2 and abs(deltaE) < 0.3', path=main)

#============ MC MATCHING ============#
if isMC:
    ma.matchMCTruth('K_S0:looseKs', path=main)
    ma.matchMCTruth('K+:looseK', path=main)
    ma.matchMCTruth('B0:loose', path=main)
    ma.matchMCTruth('B+:loose', path=main)

#============save cut value to ntuple ============#


# save cut
ma.variablesToNtuple('pi+:loosepi', track_cuts + ['pionIDNN', 'p'] + mc_truth_vars,
                    filename=outputFile, treename='cut_pions', path=main)
ma.variablesToNtuple('K+:looseK', track_cuts + ['kaonIDNN', 'p'] + mc_truth_vars,
                    filename=outputFile, treename='cut_kaons', path=main)
ma.variablesToNtuple('mu+:loosemu', track_cuts + ['muonIDNN_my', 'p'] + mc_truth_vars,
                    filename=outputFile, treename='cut_muons', path=main)
ma.variablesToNtuple('e+:loosee', track_cuts + ['electronIDNN_my', 'p'] + mc_truth_vars, 
                    filename=outputFile, treename='cut_electrons', path=main)

#save K_S0 cut
ks_vars = ['M', 'flightDistance', 'distSign', 'chiProb', 'cosAngleBetweenMomentumAndVertexVector',
           'daughter(0, p)', 'daughter(1, p)'] + mc_truth_vars 
#save B0 cut
b0_vars = ['Mbc', 'deltaE', 'chiProb', 'chisq2ndf', 'r2', 'jpsi_M', 'jpsi_chiProb',
           'ks_flightDist', 'ks_distSign'] + pid_cuts + mc_truth_vars 
# save B+ cut
bp_vars = ['Mbc', 'deltaE', 'chiProb', 'chisq2ndf', 'r2', 'jpsi_M', 'jpsi_chiProb',
           'k_plus_ID', 'k_plus_p'] + pid_cuts + mc_truth_vars


ma.variablesToNtuple('K_S0:looseKs', ks_vars,
                    filename=outputFile, treename='cut_ks', path=main)

ma.variablesToNtuple('B0:loose', b0_vars,
                    filename=outputFile, treename='cut_b0', path=main)

ma.variablesToNtuple('B+:loose', bp_vars,
                    filename=outputFile, treename='cut_bp', path=main)





#==============Begin: J/psi selection=====================#
# J/psi -> mu+mu- reconstruction
reconstructDecay('J/psi:myJpsimm -> mu+:mymu mu-:mymu', '2.8 < M < 3.4', dmID=1, path=main)
if isMC:
    # Match J/psi candidates
    ma.matchMCTruth("J/psi:myJpsimm", path=main)
#vt.treeFit('J/psi:myJpsimm', 0.0,massConstraint=['J/psi'], path=main)

vt.treeFit('J/psi:myJpsimm', 0.0, massConstraint=['J/psi'], path=main)
applyCuts('J/psi:myJpsimm', '', path=main)

# J/psi -> e+e- reconstruction
#reconstructDecay('J/psi:myJpsiee -> e+:mye e-:mye', '2.91 < M < 3.19', dmID=2, path=main)  # Slightly wider window for e+e-
ma.reconstructDecay(
    "J/psi:myJpsiee -> e+:mye e-:mye ?addbrems",
    cut="2.8 < M < 3.4", dmID=2,
    path=main
)
if isMC:
    # Match J/psi candidates
    ma.matchMCTruth("J/psi:myJpsiee", path=main)
#vt.treeFit('J/psi:myJpsiee', 0.0, massConstraint=['J/psi'], path=main)
vt.treeFit('J/psi:myJpsiee', 0.0,  massConstraint=['J/psi'], path=main)
applyCuts('J/psi:myJpsiee', '', path=main)

# Combine both J/psi lists
copyLists('J/psi:myJpsi', ['J/psi:myJpsimm', 'J/psi:myJpsiee'], path=main)
ma.variablesToNtuple('J/psi:myJpsimm', jpsi_vars_to_save,
                     filename=outputFile, treename='jpsi_mm', path=main)
ma.variablesToNtuple('J/psi:myJpsiee', jpsi_vars_to_save,
                     filename=outputFile, treename='jpsi_ee', path=main)
#==============End: J/psi selection=====================#

#==============Begin: K_S0 selection=====================#
# Get K_S0 candidates with standard reconstruction
#stdV0s.stdKshorts(path=main)
add_ks_selection_aliases()


# Apply additional quality cuts to K_S0 - using standard cuts
#cutAndCopyList('K_S0:goodKs', 'K_S0:merged',
#               'flightDistance > 0.1 and distSign > 3.0', path=main)

ma.cutAndCopyList('K_S0:preKs', 'K_S0:merged',
                    '0.45 < M < 0.55',
                    path=main)


main.add_module('MVAMultipleExperts',
              listNames=['K_S0:preKs'],
              extraInfoNames=['KsSelector_V0Selector', 'KsSelector_LambdaVeto'],
              identifiers=['Ks_LGBM_V0Selector', 'Ks_LGBM_LambdaVeto'])


# Define variables to save for Ks pre-MVA cut, including MVA scores
ks_standard_vars = ['M']
ks_preMVA_vars = ['M',
                  'extraInfo(KsSelector_V0Selector)', # Save MVA score 1
                  'extraInfo(KsSelector_LambdaVeto)']  # Save MVA score 2
#ks_MVA_vars =['extraInfo(KsSelector_V0Selector)', # Save MVA score 1
#                  'extraInfo(KsSelector_LambdaVeto)']  # Save MVA score 2
v.variables.addAlias('ks_v0', 'daughter(1, extraInfo(KsSelector_V0Selector))')
v.variables.addAlias('ks_lambda', 'daughter(1, extraInfo(KsSelector_LambdaVeto))')

# Update ks_MVA_vars to use the new aliases
ks_MVA_vars = ['ks_v0', 'ks_lambda']
if isMC:
    # Match kaon candidates
    ma.matchMCTruth("K_S0:merged", path=main)
    ma.matchMCTruth("K_S0:preKs", path=main)

# Save the ntuple BEFORE applying the MVA cut, but AFTER calculating scores
ma.variablesToNtuple('K_S0:merged', ks_standard_vars + mc_truth_vars,
                     filename=outputFile, treename='ks_standard', path=main)
#ma.variablesToNtuple('K_S0:preKs', ks_preMVA_vars + mc_truth_vars + ks_pion_vars,
#                     filename=outputFile, treename='cut_ks_preMVA', path=main)

#ma.cutAndCopyList('K_S0:goodKs', 'K_S0:preKs',
#                     'extraInfo(KsSelector_V0Selector) > 0.90 and '
#                     'extraInfo(KsSelector_LambdaVeto) > 0.11',
#                     path=main)
ma.cutAndCopyList('K_S0:goodKs', 'K_S0:preKs',
                     '',
                     path=main)
if isMC:
    ma.matchMCTruth("K_S0:goodKs", path=main)

ma.variablesToNtuple('K_S0:goodKs', ks_preMVA_vars + mc_truth_vars + ks_pion_vars, filename=outputFile, treename='ks_MVA', path=main)
#ma.variablesToNtuple('J/psi:myJpsiee', jpsi_vars_to_save,
#                   filename=outputFile, treename='jpsi_ee', path=main)
#==============End: K_S0 selection=====================#i

#==============Begin: B0 selection=====================#
# Reconstruct B0 candidates with combined J/psi list and improved K_S0 selection
ma.reconstructDecay(
    "B0 -> J/psi:myJpsi K_S0:goodKs",
    cut="Mbc > 5.2 and abs(deltaE) < 0.3",
    path=main,
)


v.variables.addAlias('nKsPreMVA', 'nParticlesInList(K_S0:preKs)')
# Count particles AFTER the MVA cut is applied (in the final selected list)
v.variables.addAlias('nKsGood', 'nParticlesInList(K_S0:goodKs)')

# Apply vertex fit to B0 candidates
#vt.treeFit('B0', 0.0, massConstraint=['J/psi'], path=main)  # 0.0 means no IP constraint
vt.treeFit('B0', 0.0, path=main)  # 0.0 means no IP constraint
applyCuts('B0', 'chiProb > 0.001 and chisq2ndf < 10', path=main)

# Apply continuum suppression to B0
#applyCuts('B0', 'r2 < 0.7', path=main)

# Apply additional quality cuts to B0 candidates
applyCuts('B0', 
         'Mbc > 5.2 and abs(deltaE) < 0.3 and ' + 
         'daughter(0, chiProb) > 0.001',
         path=main)
#==============End: B0 selection=====================#

#==============Begin: B+ selection=====================#
# Reconstruct B+ candidates with combined J/psi list and K+ selection
ma.reconstructDecay(
    "B+ -> J/psi:myJpsi K+:myK",
    cut="Mbc > 5.2 and abs(deltaE) < 0.3",
    path=main,
)

# Apply vertex fit to B+ candidates
#vt.treeFit('B+', 0.0, massConstraint=['J/psi'], path=main)  # 0.0 means no IP constraint
vt.treeFit('B+', 0.0, path=main)  # 0.0 means no IP constraint
applyCuts('B+', 'chiProb > 0.001 and chisq2ndf < 10', path=main)

# Apply continuum suppression to B+
#applyCuts('B+', 'r2 < 0.7', path=main)

# Apply additional quality cuts to B+ candidates
applyCuts('B+', 
         'Mbc > 5.2 and abs(deltaE) < 0.3 and ' + 
         'daughter(0, chiProb) > 0.001',
         path=main)
#==============End: B+ selection=====================#

#==============Begin: MC matching (if MC)=====================#
# If MC, add detailed MC matching
if isMC:
    # Match B candidates
    ma.matchMCTruth("B0", path=main)
    ma.matchMCTruth("B+", path=main)
#==============End: MC matching=====================#

#==============Begin: BCS =====================#
# Best Candidate Selection - keep only the best B candidate per event
ma.rankByLowest("B0", 'abs(deltaE)', numBest=1, path=main)
ma.rankByLowest("B+", 'abs(deltaE)', numBest=1, path=main)
#==============End: BCS =====================#

#==============Begin: outputs =====================#
# Update the variable list to include all the new variables
B_vars = evtLel + kinematics_vars + bplus_vars + vertex_vars + fit_vars 
basic_vars = ['nCDCHits', 'dr', 'dz']
# Variables for B0 specific components
B0_jpsi_vars = vu.create_aliases_for_selected(fit_vars + vertex_vars + kinematics_vars, 'B0 -> ^J/psi K_S0')
B0_ks_vars = vu.create_aliases_for_selected(fit_vars + vertex_vars + kinematics_vars, 'B0 -> J/psi ^K_S0')
B0_kaon_vars = vu.create_aliases_for_selected(trk_vars, 'B0 -> J/psi ^K_S0')
B0_ks_mva_vars = vu.create_aliases_for_selected(['extraInfo(KsSelector_V0Selector)', 
                                              'extraInfo(KsSelector_LambdaVeto)'], 
                                             'B0 -> J/psi ^K_S0')
# Variables for B+ specific components
Bp_jpsi_vars = vu.create_aliases_for_selected(fit_vars + vertex_vars + kinematics_vars, 'B+ -> ^J/psi K+')
Bp_k_vars = vu.create_aliases_for_selected(fit_vars + vertex_vars + kinematics_vars, 'B+ -> J/psi ^K+')
Bp_kaon_vars = vu.create_aliases_for_selected(trk_vars, 'B+ -> J/psi ^K+')

# Output combined analysis for B0 and B+
ma.variablesToNtuple("B0",basic_vars + B_vars + B0_jpsi_vars + B0_ks_vars + B0_kaon_vars + Jpsi_daughter_alias_vars + mc_truth_vars + ks_MVA_vars + B0_ks_mva_vars, 
                    filename=outputFile, treename='rec_B0', path=main)

ma.variablesToNtuple("B+",basic_vars + B_vars + Bp_jpsi_vars + Bp_k_vars + Bp_kaon_vars + mc_truth_vars, 
                    filename=outputFile, treename='rec_Bp', path=main)

# Separate B0 candidates by J/psi decay mode
ma.cutAndCopyList('B0:mumu', 'B0', 'daughter(0, decaymode) == 1', path=main)
ma.cutAndCopyList('B0:ee', 'B0', 'daughter(0, decaymode) == 2', path=main)

# Separate B+ candidates by J/psi decay mode
ma.cutAndCopyList('B+:mumu', 'B+', 'daughter(0, decaymode) == 1', path=main)
ma.cutAndCopyList('B+:ee', 'B+', 'daughter(0, decaymode) == 2', path=main)

# Output separated by J/psi decay mode for B0
ma.variablesToNtuple('B0:mumu', basic_vars + B_vars + B0_jpsi_vars + B0_ks_vars + B0_kaon_vars + mc_truth_vars + ks_MVA_vars, 
                    filename=outputFile, treename='rec_B0_mumu', path=main)
ma.variablesToNtuple('B0:ee', basic_vars + B_vars + B0_jpsi_vars + B0_ks_vars + B0_kaon_vars + mc_truth_vars + ks_MVA_vars + B0_ks_mva_vars, 
                    filename=outputFile, treename='rec_B0_ee', path=main)

# Output separated by J/psi decay mode for B+
ma.variablesToNtuple('B+:mumu', basic_vars + B_vars + Bp_jpsi_vars + Bp_k_vars + Bp_kaon_vars + mc_truth_vars, 
                    filename=outputFile, treename='rec_Bp_mumu', path=main)
ma.variablesToNtuple('B+:ee', basic_vars + B_vars + Bp_jpsi_vars + Bp_k_vars + Bp_kaon_vars + mc_truth_vars, 
                    filename=outputFile, treename='rec_Bp_ee', path=main)
#==============End: outputs =====================#

# put by jwpark
#main.add_module("ParticlePrinter",listName="B0:gen", variables= ["M", "p", "E"])
#ma.printMCParticles(onlyPrimaries=True, path=main)
ma.summaryOfLists(
    ['pi+:gen', 'pi+:all', 'pi+:mypi', 'pi+:loosepi'],
    path=main
)
# Process the path
b2.process(main)
print(b2.statistics)
