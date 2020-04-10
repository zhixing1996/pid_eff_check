#!/usr/bin/env python
"""
Parameter Base
"""

__author__ = "JING Maoqiang <jingmq@ihep.ac.cn>"
__copyright__ = "Copyright (c) JING Maoqiang"
__created__ = "[2020-03-29 Sun 17:08]" 

import sys 
import os
import ROOT 


# ---------------------------------------------
# Parameters
# ---------------------------------------------

def pid_eff(var):
    if var == 'costheta':
        xmin = -1.0
        xmax = 1.0
        xbins = 20
    if var == 'ptrk':
        xmin = 0.2
        xmax = 1.3
        xbins = 22
    if var == 'phi':
        xmin = -3.14
        xmax = 3.14
        xbins = 20
    return xmin, xmax, xbins
def data_path(particle):
    if particle == 'e':
        data_path = '/besfs/groups/cal/dedx/jingmq/example/qiujf/hadron_track/electron/electron.root'
    if particle == 'mu':
        data_path = '/besfs/groups/cal/dedx/jingmq/example/qiujf/hadron_track/muon/muon.root'
    if particle == 'pi':
        data_path = '/besfs/groups/cal/dedx/jingmq/example/qiujf/hadron_track/pion/pion.root'
    if particle == 'K':
        data_path = '/besfs/groups/cal/dedx/jingmq/example/qiujf/hadron_track/kaon/kaon.root'
    if particle == 'p':
        data_path = '/besfs/groups/cal/dedx/jingmq/example/qiujf/hadron_track/proton/proton.root'
    return data_path
def com_patches(var, particle):
    file_list = []
    file_list.append('/besfs/groups/cal/dedx/jingmq/bes/pid_eff_check/python/files/dedx_pideff_' + var + '_' + particle + '_705_data.root')
    # file_list.append('/besfs/groups/cal/dedx/jingmq/bes/pid_eff_check/python/files/test_dedx_pideff_' + var + '_' + particle + '_705_data.root')
    file_list.append('/besfs/groups/cal/dedx/jingmq/bes/pid_eff_check/python/files/test_dedx_pideff_' + var + '_K_705_data.root')
    file_list.append('/besfs/groups/cal/dedx/jingmq/bes/pid_eff_check/python/files/test_dedx_pideff_' + var + '_p_705_data.root')
    return file_list
