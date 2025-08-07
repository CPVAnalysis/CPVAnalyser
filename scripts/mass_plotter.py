import os
import sys
import ROOT
from glob import glob
sys.path.append('../')
from nanoTools import NanoTools
from tools import Tools
from array import array
from mva_tools import MVATools
import matplotlib.pyplot as plt


do_plot_mass = False
do_create_tree = False
do_create_mva_tree = True


glob_path = '/eos/user/a/anlyon/CPVGen/V00_*/*/*/*/*/flat'
filename = 'flat_bparknano.root'

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)


# get data
chain = ROOT.TChain('signal_tree')
list_files = [f for f in glob('{}/{}'.format(glob_path, filename))]
for f in list_files:
   print f
   if NanoTools().checkFlatFile(f, True, branch_check=True, branchname='bs_mass'): 
      chain.Add(f)
   else:
      print 'file not valid'
print chain.GetEntries()

#for entry in chain:
#   print entry.bs_mass


# get signal
f_signal = ROOT.TFile.Open('/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/flat_bparknano.root')
tree_signal = f_signal.Get('signal_tree')
      
if do_plot_mass:
   canv = ROOT.TCanvas('canv', 'canv', 900, 800)
   pad_up = ROOT.TPad("pad_up","pad_up",0.02,0,1,1)
   pad_up.Draw()
   pad_up.cd()

   hist = ROOT.TH1D('hist', 'hist', 50, 5.1, 5.6)
   hist_sel = ROOT.TH1D('hist_sel', 'hist_sel', 50, 5.1, 5.6)
   hist_signal = ROOT.TH1D('hist_signal', 'hist_signal', 50, 5.1, 5.6)

   # data
   #chain.Draw("bs_mass_corr>>hist")
   #chain.Draw("bs_mass_corr>>hist_sel", "phi1_pt>4 && phi2_pt>3 && deltar_min<0.07 && deltar_max<1.3 && bs_lxysig>7 && bs_cos2d>0.95")
   chain.Draw("bs_mass_corr>>hist_sel", "phi1_pt> 6.5 && phi2_pt > 5 && bs_lxysig > 15 && bs_cos2d > 0.995") # cutbased_1
   #chain.Draw("bs_mass_corr>>hist_sel", "bs_lxysig > 12 && bs_sv_prob>0.1 && bs_cos2d > 0.995 && bs_pt>7") # cutbased_GK_1
   #chain.Draw("bs_mass_corr>>hist_sel", "bs_lxysig > 12 && bs_sv_prob>0.1 && bs_cos2d > 0.995 && bs_pt>7 && deltar_max<1") # cutbased_GK_2
   #chain.Draw("bs_mass_corr>>hist_sel", "phi1_pt> 6.5 && phi2_pt > 5 && deltar_min < 0.05 && deltar_max < 0.75 && bs_lxysig > 15 && bs_cos2d > 0.995") # cutbased_2
   #chain.Draw("bs_mass_corr>>hist_sel", "phi1_pt> 6.5 && phi2_pt > 5 && deltar_min < 0.05 && deltar_max < 0.75 && bs_lxysig > 15 && bs_cos2d > 0.995 && bs_sv_prob>0.4") # cutbased_3
   #chain.Draw("bs_mass_corr>>hist_sel", "phi1_pt>7 && phi2_pt>5 && bs_lxysig > 25") # cutbased_4

   hist_sel.SetLineWidth(3)
   hist_sel.SetLineColor(ROOT.kBlue+1)

   # signal
   tree_signal.Draw("bs_mass_corr>>hist_signal", "phi1_pt>4 && phi2_pt>3 && deltar_min<0.07 && deltar_max<1.3 && bs_lxysig>7 && bs_cos2d>0.95")

   # normalise signal events to number in data
   n_data = hist_sel.Integral()
   n_signal = hist_signal.Integral()
   hist_signal.Scale(0.75 * (5.45-5.3)/(5.7-5.05) * n_data/n_signal)
         
   hist_signal.SetLineWidth(2)
   hist_signal.SetLineStyle(ROOT.kDashed)
   hist_signal.SetLineColor(ROOT.kOrange+1)

   # create frame  
   frame = hist_sel.Clone('frame')

   frame.SetTitle('')
   frame.GetXaxis().SetTitle('#it{m}(K^{+}K^{-}K^{+}K^{-}) (GeV)')
   frame.GetXaxis().SetLabelSize(0.04)
   frame.GetXaxis().SetTitleSize(0.047)
   frame.GetXaxis().SetTitleOffset(1.1)

   frame.GetYaxis().SetTitle('Events / Bin')
   frame.GetYaxis().SetLabelSize(0.04)
   frame.GetYaxis().SetTitleSize(0.047)
   frame.GetYaxis().SetTitleOffset(1.5)
   frame.GetYaxis().SetRangeUser(0, 1.3*hist_sel.GetMaximum())
       

   frame.Draw()
   hist_sel.Draw('hist same')
   #hist_signal.Draw('hist same')

   Tools().printInnerCMSTag(pad_up, 'Preliminary', True, x_pos=0.17, y_pos=0.835, size=0.55)
   Tools().printLumiTag(pad_up, 33.6, size=0.43, offset=0.01)

   ROOT.gStyle.SetOptStat(0)

   canv.cd()
   canv.SaveAs('./myPlots/bs_mass_corr.png')
   canv.SaveAs('./myPlots/bs_mass_corr.pdf')


if do_create_tree:

   label = 'cutbased1'
   root_filename = './myPlots/selection/tree_{}.root'.format(label)
   out_file = ROOT.TFile(root_filename, 'RECREATE')

   treename = 'signal_tree'
   tree = ROOT.TTree(treename, treename)
       
   the_bs_mass_corr = array('d', [0])
   tree.Branch('bs_mass_corr', the_bs_mass_corr, 'bs_mass_corr/D')
    
   for ientry, entry in enumerate(chain):
      if ientry%1000 == 0:
        print '{}% events processed'.format(round(float(ientry)/float(chain.GetEntries())*100, 1))

      #if entry.bs_lxysig > 12 and entry.bs_sv_prob>0.1 and entry.bs_cos2d > 0.995 and entry.bs_pt>7 and entry.deltar_max<1: # GK2
      if entry.phi1_pt> 6.5 and entry.phi2_pt > 5 and entry.bs_lxysig > 15 and entry.bs_cos2d > 0.995:  # cutbased_1
         the_bs_mass_corr[0] = entry.bs_mass_corr
         tree.Fill()

   tree.Write()
   out_file.Close()

   print '--> {} created'.format(root_filename)


if do_create_mva_tree:

   #label = 'mva'
   #files = list_files
   #training_label = 'test_2025Mar11_12h42m51s' # does not look good 
   #training_label = 'test_2025Mar11_14h18m43s' # is unbalanced (D1_0)
   #training_label = 'test_2025Mar11_14h55m14s' # balanced (C1_0)
   training_label = 'test_2025Mar11_16h33m47s' # removed invmass --> ok
   #training_label = 'test_2025Mar11_16h54m22s/' # balanced
   #training_label = 'test_2025Mar11_17h59m34s' # sig weighted by 20

   # for data
   #for i in range(1, 7):
   #  print i
   #
   #  glob_path = '/eos/user/a/anlyon/CPVGen/V00_*{}/*/*/*/*/flat'.format(i)
   #  filename = 'flat_bparknano.root'

   #  files = [f for f in glob('{}/{}'.format(glob_path, filename))]
   #  label = i

   #  MVATools(path_mva='./outputs').getFileWithScore(files=files, training_label=training_label, label=label, force_overwrite=True)

   # for signal
   files = ['/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/flat_bparknano.root']

   label = 'signal'
   MVATools(path_mva='./outputs').getFileWithScore(files=files, training_label=training_label, label=label, force_overwrite=True)

