#!/usr/bin/env python
"""
Plot pid efficiency
"""

__author__ = "Maoqiang JING <jingmq@ihep.ac.cn>"
__copyright__ = "Copyright (c) Maoqiang JING"
__created__ = "[2020-03-27 Fri 10:53]"

import sys, os
import logging
from math import *
from ROOT import *
from params import pid_eff, data_path
from array import array
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')
gStyle.SetOptTitle(0)
gStyle.SetOptTitle(0)
gStyle.SetCanvasDefH(900)
gStyle.SetCanvasDefW(450)
gStyle.SetNdivisions(505,"xy")
gStyle.SetStatFontSize(0.05)

def usage():
    sys.stdout.write('''
NAME
    plot_pid_eff.py

SYNOPSIS
    ./plot_pid_eff.py [var] [particle] [ptrk_min] [ptrk_max]

AUTHOR
    Maoqiang JING <jingmq@ihep.ac.cn>

DATE
    March 2020
\n''')

def set_canvas_style(mbc):
    mbc.SetFillColor(0)
    mbc.SetLeftMargin(0.15)
    mbc.SetRightMargin(0.15)
    mbc.SetTopMargin(0.1)
    mbc.SetBottomMargin(0.15)

def plot(var, particle, ptrk_min, ptrk_max):
    try:
        chain = TChain('n103')
        chain.Add(data_path(particle))
    except:
        logging.error(data_path(particle) + ' is invalid!')
        sys.exit()

    if var == 'costheta' or var == 'phi':
        cut_ptrkp = '('+ptrk_min+'<ptrk<'+ptrk_max+')&&'
        cut_ptrkm = '('+ptrk_min+'<ptrk<'+ptrk_max+')&&'
    elif var == 'ptrk':
        cut_ptrkp = ''
        cut_ptrkm = ''
    mbc = TCanvas('mbc', 'mbc')
    set_canvas_style(mbc)
    mbc.Divide(2, 3)
    xmin, xmax, xbins = pid_eff(var)
    h_plus_nume = TH1F('h_plus_nume', ' ', xbins, xmin, xmax)
    h_plus_deno = TH1F('h_plus_deno', ' ', xbins, xmin, xmax)
    h_minus_nume = TH1F('h_minus_nume', ' ', xbins, xmin, xmax)
    h_minus_deno = TH1F('h_minus_deno', ' ', xbins, xmin, xmax)
    if particle == 'e':
        cut_prob = '(prob[0]>prob[1]&&prob[0]>prob[2]&&prob[0]>prob[3]&&prob[0]>prob[4])'
    if particle == 'mu':
        cut_prob = '(prob[1]>prob[0]&&prob[1]>prob[2]&&prob[1]>prob[3]&&prob[1]>prob[4])'
    if particle == 'pi':
        cut_prob = '(prob[2]>prob[3]&&prob[2]>prob[4])'
    if particle == 'K':
        cut_prob = '(prob[3]>prob[2]&&prob[3]>prob[4])'
    if particle == 'p':
        cut_prob = '(prob[4]>prob[2]&&prob[4]>prob[3])'
    mbc.cd(1)
    chain.Draw(var + '>>h_plus_nume', cut_ptrkp + '(charge==1)&&' + cut_prob)
    mbc.cd(2)
    chain.Draw(var + '>>h_plus_deno', cut_ptrkp + '(charge==1)')
    mbc.cd(3)
    chain.Draw(var + '>>h_minus_nume', cut_ptrkm + '(charge==-1)&&' + cut_prob)
    mbc.cd(4)
    chain.Draw(var + '>>h_minus_deno', cut_ptrkm + '(charge==-1)')

    mbc.cd(5)
    h_plus_nume.Sumw2()
    h_plus_deno.Sumw2()
    h_plus_pid_eff = TH1F('h_plus_eff', ' ', xbins, xmin, xmax)
    h_plus_pid_eff.Divide(h_plus_nume, h_plus_deno)
    if particle == 'e':
        frame1 = TH2F('frame1', 'Eff of e^{+} for data', xbins, xmin, xmax, 22, 0, 1.1)
    if particle == 'mu':
        frame1 = TH2F('frame1', 'Eff of #mu^{+} for data', xbins, xmin, xmax, 22, 0, 1.1)
    if particle == 'pi':
        frame1 = TH2F('frame1', 'Eff of #pi^{+} for data', xbins, xmin, xmax, 22, 0, 1.1)
    elif particle == 'K':
        frame1 = TH2F('frame1', 'Eff of K^{+} for data', xbins, xmin, xmax, 22, 0, 1.1)
    elif particle == 'p':
        frame1 = TH2F('frame1', 'Eff of p for data', xbins, xmin, xmax, 22, 0, 1.1)
    frame1.Draw()
    h_plus_pid_eff.SetLineColor(2)
    h_plus_pid_eff.Draw('same')

    mbc.cd(6)
    h_minus_nume.Sumw2()
    h_minus_deno.Sumw2()
    h_minus_pid_eff = TH1F('h_minus_eff', ' ', xbins, xmin, xmax)
    h_minus_pid_eff.Divide(h_minus_nume, h_minus_deno)
    if particle == 'e':
        frame2 = TH2F('frame2', 'Eff of e^{-} for data', xbins, xmin, xmax, 22, 0, 1.1)
    if particle == 'mu':
        frame2 = TH2F('frame2', 'Eff of #mu^{-} for data', xbins, xmin, xmax, 22, 0, 1.1)
    if particle == 'pi':
        frame2 = TH2F('frame2', 'Eff of #pi^{-} for data', xbins, xmin, xmax, 22, 0, 1.1)
    elif particle == 'K':
        frame2 = TH2F('frame2', 'Eff of K^{-} for data', xbins, xmin, xmax, 22, 0, 1.1)
    elif particle == 'p':
        frame2 = TH2F('frame2', 'Eff of #bar{p} for data', xbins, xmin, xmax, 22, 0, 1.1)
    frame2.Draw()
    h_minus_pid_eff.SetLineColor(2)
    h_minus_pid_eff.Draw('same')

    if not os.path.exists('./files/'):
        os.makedirs('./files/')
    out_root = TFile('./files/dedx_pideff_' + var + '_' + particle + '_705_data.root', 'RECREATE')
    h_plus_pid_eff.Write()
    h_minus_pid_eff.Write()
    out_root.Write()

    if not os.path.exists('./figs/'):
        os.makedirs('./figs/')
    mbc.SaveAs('./figs/pid_eff_' + var + '_' + particle + '.pdf')
    raw_input('Enter anything to end...')

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args)<4:
        usage()
        sys.exit()
    var = args[0]
    particle = args[1]
    ptrk_min = args[2]
    ptrk_max = args[3]

    plot(var, particle, ptrk_min, ptrk_max)
