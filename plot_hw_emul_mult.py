#!/bin/env python

# plot the validation histograms made with skimL1.py
#
# davide.gerbaudo@gmail.com
# Feb 2016

import os

import ROOT as R
R.gROOT.SetBatch(1)

def main():
    input_file = R.TFile.Open('/tmp/tmptrig.root')
    output_dir = mkdir_if_needed('/tmp/gerbaudo/plots')

    can = R.TCanvas('can')

    can.cd()
    hnames = ['h_l1jpsi_emall_mult_match',
              'h_l1btag_mu4ab_mult_match',
              'h_l1btag_cj15ab_mult_match']
    titles = ['L1_JPSI-1M5-EM7',
              'L1_BTAG_MU4J15',
              'L1_BTAG_MU4J15']

    for hn, title in zip(hnames, titles):
        can.Clear()
        h = input_file.Get(hn)
        h.SetTitle(title)
        h.Draw('colz text')
        can.Update()
        stats = h.FindObject('stats')
        st_width  = abs(stats.GetX1NDC()-stats.GetX2NDC())
        st_heigth = abs(stats.GetY1NDC()-stats.GetY2NDC())
        stats.SetX1NDC(0.5)
        stats.SetY1NDC(0.5)
        stats.SetX2NDC(0.5+st_width)
        stats.SetY2NDC(0.5+st_heigth)
        can.Update()
        can.SaveAs(os.path.join(output_dir, hn+'.png'))

def mkdir_if_needed(dirname) :
    dest_dir = None
    if os.path.exists(dirname) and os.path.isdir(dirname) :
        dest_dir = dirname
    elif not os.path.exists(dirname) :
        os.makedirs(dirname)
        dest_dir = dirname
    return dest_dir

if __name__=='__main__':
    main()
