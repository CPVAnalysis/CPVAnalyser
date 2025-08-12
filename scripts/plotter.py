import os
import ROOT
import math
from os import path

import sys
sys.path.append('../objects')
from tools import Tools
from quantity import quantities

#TODO fetch qte from PhiToKK from Bs (e.g. deltaR)
#TODO add min and max deltaR between k1 and (k3, k4)
#TODO add min and max deltaR between k2 and (k3, k4)


class Plotter(Tools):
  def __init__(self, quantity='', data_filenames='', signal_filenames='', outdirlabel='', isnano=False):
    self.quantity = quantity
    self.data_filenames = data_filenames
    self.signal_filenames = signal_filenames
    self.outdirlabel = outdirlabel

    self.tools = Tools()
    if isnano:
        self.treename = 'Events'
        self.branchname = 'nano'
    else:
        self.treename = 'signal_tree'
        self.branchname = 'flat'

    if isnano:
        self.selection_sig = 'BsToPhiPhiTo4K_isMatched==1 && abs(BsToPhiPhiTo4K_Bs_ct_2D_cm/BsToPhiPhiTo4K_gen_Bs_ct -1)<1'
        self.selection_bkg = 'fabs(BsToPhiPhiTo4K_Bs_mass - 5.367) > 0.2'
    else:
        self.selection_sig = 'ismatched == 1'
        self.selection_bkg = 'fabs(bs_mass_corr - 5.367) > 0.2'

    self.do_shape = True
    self.add_overflow = True

    self.outputdir = './myPlots/{}'.format(self.outdirlabel)
    if not path.exists(self.outputdir):
      os.system('mkdir -p {}'.format(self.outputdir))


  def getMaxRangeY(self, hist1, hist2, do_log=False, use_sig=False):
    margin = 0.15 if do_log==False else 0.5
    if use_sig: 
      max_hist1 = [hist.GetMaximum() for hist in hist1]
      the_max_hist1 = max(max_hist1)
      max_range = the_max_hist1+margin*the_max_hist1 if the_max_hist1>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    else:
      max_range = hist1.GetMaximum()+margin*hist1.GetMaximum() if hist1.GetMaximum()>hist2.GetMaximum() else hist2.GetMaximum()+margin*hist2.GetMaximum()
    return max_range


  def plot(self):

    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    # create the canvas
    canv_name = 'canv_{}_{}_{}_{}'.format(self.quantity.label, outdirlabel.replace('/', '_'), self.quantity.do_log, self.do_shape)
    canv = self.tools.createTCanvas(name=canv_name, dimx=1200, dimy=1000)
    canv.cd()

    # define the pads
    pad_up = ROOT.TPad("pad_up","pad_up",0.02,0,1,1)
    if self.quantity.do_log: 
      pad_up.SetLogy()
    pad_up.Draw()
    canv.cd()

    # prepare the legend
    leg = self.tools.getRootTLegend(xmin=0.55, ymin=0.61, xmax=0.8, ymax=0.82, size=0.038)

    pad_up.cd()

    # background
    hist_data_tot = ROOT.TH1D('hist_data_tot', 'hist_data_tot', self.quantity.nbins, self.quantity.bin_min, self.quantity.bin_max)
    hist_data_tot.Sumw2()

    int_data_tot = 0.
    overflow_data_tot = 0
    error_overflow_data_tot = 0.

    hist_data_stack = ROOT.THStack('hist_data_stack', '')
    data_label = ''

    for idata, data_filename in enumerate(self.data_filenames):
      f_data = self.tools.getRootFile(data_filename)
      tree_data = self.tools.getTree(f_data, self.treename)

      selection_data = self.selection_bkg

      hist_data_name = 'hist_data_{}_{}_{}_{}'.format(self.quantity, outdirlabel.replace('/', '_'), self.quantity.do_log, self.do_shape)
      hist_data = self.tools.createHisto(tree_data, self.quantity, hist_name=hist_data_name, branchname=self.branchname, selection=selection_data)
      hist_data.Sumw2()

      if self.add_overflow:
        overflow_data_tot = (hist_data.GetBinContent(hist_data.GetNbinsX()) + hist_data.GetBinContent(hist_data.GetNbinsX()+1))
        error_overflow_data_tot += math.sqrt(math.pow(hist_data.GetBinError(hist_data.GetNbinsX()), 2) + math.pow(hist_data.GetBinError(hist_data.GetNbinsX()+1), 2)) 
        hist_data.SetBinContent(hist_data.GetNbinsX(), overflow_data_tot)
        hist_data.SetBinError(hist_data.GetNbinsX(), error_overflow_data_tot)
        hist_data.SetBinContent(hist_data.GetNbinsX()+1, 0)
        hist_data.SetBinError(hist_data.GetNbinsX()+1, 0)

      if self.do_shape: 
        int_data_tot += hist_data.Integral()

      hist_data_tot.Add(hist_data)

      if self.do_shape and hist_data.Integral() != 0: 
        hist_data.Scale(1./hist_data.Integral())
        hist_data_stack.Add(hist_data)

    if self.do_shape and int_data_tot != 0.: hist_data_tot.Scale(1./int_data_tot)

    leg.AddEntry(hist_data_tot, 'Background', 'elpf')

    # set the style
    hist_data_tot.SetFillColor(ROOT.kBlue-3)
    hist_data_tot.SetFillStyle(3005)

    # signal
    signal_hists = []
    for signal_filename in self.signal_filenames:
      f_sig = self.tools.getRootFile(signal_filename)
      tree_sig = self.tools.getTree(f_sig, self.treename)

      hist_signal_name = 'hist_signal_{}_{}_{}_{}'.format(self.quantity, self.outdirlabel.replace('/', '_'), self.quantity.do_log, self.do_shape)

      selection_signal = self.selection_sig

      weight_sig = '(1)'
      #if add_weight_hlt : weight_sig +=0 * ({})'.format(weight_hlt)
      #if add_weight_pu : weight_sig +=0 * ({})'.format(weight_pusig)
      #if add_weight_muid : weight_sig +=0 * ({}) *({})'.format(weight_mu0id, weight_muid)

      hist_signal = self.tools.createHisto(tree_sig, self.quantity, hist_name=hist_signal_name, branchname=self.branchname, selection=selection_signal, weight=weight_sig)
      hist_signal.Sumw2()

      if self.add_overflow:
        overflow_signal_tot = hist_signal.GetBinContent(hist_signal.GetNbinsX()) + hist_signal.GetBinContent(hist_signal.GetNbinsX()+1)
        error_overflow_signal_tot = math.sqrt(math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()), 2) + math.pow(hist_signal.GetBinError(hist_signal.GetNbinsX()+1), 2)) 
        hist_signal.SetBinContent(hist_signal.GetNbinsX(), overflow_signal_tot)
        hist_signal.SetBinError(hist_signal.GetNbinsX(), error_overflow_signal_tot)
        hist_signal.SetBinContent(hist_signal.GetNbinsX()+1, 0)
        hist_signal.SetBinError(hist_signal.GetNbinsX()+1, 0)

      if self.do_shape: 
        int_signal = hist_signal.Integral()
        if int_signal != 0: hist_signal.Scale(1/int_signal)

      leg.AddEntry(hist_signal,'Signal')

      hist_signal.SetLineWidth(3)
      hist_signal.SetLineColor(ROOT.kOrange+1)
      hist_signal.SetFillColorAlpha(0, 0)
      signal_hists.append(hist_signal)

    # create frame  
    frame = hist_data_tot.Clone('frame')

    frame.SetTitle('')
    frame.GetXaxis().SetTitle(quantity.title)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleSize(0.047)
    frame.GetXaxis().SetTitleOffset(1.1)

    frame.GetYaxis().SetTitle('Normalised to unity')
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleSize(0.047)
    frame.GetYaxis().SetTitleOffset(1.5)

    if not self.quantity.do_log: frame.GetYaxis().SetRangeUser(1e-4, self.getMaxRangeY(signal_hists, hist_data_tot, self.quantity.do_log, use_sig=True)+0.2*self.getMaxRangeY(signal_hists, hist_data_tot, self.quantity.do_log, use_sig=True))
    else: frame.GetYaxis().SetRangeUser(1e-4, self.getMaxRangeY(signal_hists, hist_data_tot, self.quantity.do_log, use_sig=True)+50*self.getMaxRangeY(signal_hists, hist_data_tot, self.quantity.do_log, use_sig=True))

    ROOT.gStyle.SetOptStat(0)

    # draw the distributions
    frame.Draw()
    hist_data_tot.Draw('histo same')
    hist_data_tot.Draw('PE same')

    for hist_sig in signal_hists:
      hist_sig.Draw('histo same')
      hist_sig.Draw('PE1 same')

    # draw the legend
    leg.Draw('same')

    pad_up.cd()
    pad_up.RedrawAxis()

    canv.cd()
    canv.RedrawAxis()

    # add labels
    #if do_tdrstyle and not do_shape: self.tools.printLatexBox(0.36, 0.85,0Dimuon channel', size=0.038, pos='left', font=42)
    #if do_tdrstyle and do_shape: self.tools.printLatexBox(0.43, 0.88,0charge(#it{#mu}^{#pm}#it{#pi}^{#mp}) #neq 0', size=0.038, pos='left', font=42)
    #if not do_shape and do_tdrstyle: self.tools.printLatexBox(0.405, 0.575,0(#it{r}_{e}, #it{r}_{#mu}, #it{r}_{#tau}) = (0, 1, 0), Majorana', size=0.038, pos='left', font=42)
    #if add_CMSlabel and not do_tdrstyle: self.tools.printCMSTag(pad_up, CMS_tag, size=0.55 if plot_ratio else 0.43)
    #if add_CMSlabel and do_tdrstyle: self.tools.printInnerCMSTag(pad_up, CMS_tag, True, x_pos=0.17, y_pos=0.835, size=0.55)
    self.tools.printInnerCMSTag(pad_up, 'Preliminary', True, x_pos=0.17, y_pos=0.835, size=0.55)
    #if not do_shape: self.tools.printLumiTag(pad_up, 41.6, size=0.5, offset=0.507)

    save_name = '{}/{}'.format(self.outputdir, self.quantity.label)
    if quantity.do_log: save_name += '_log'
    canv.SaveAs('{}.png'.format(save_name))
    canv.SaveAs('{}.pdf'.format(save_name))
    #canv.SaveAs('{}/{}.C'.format(outputdir, self.quantity.label))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  #outdirlabel = 'study_Bs_selection_v0'
  outdirlabel = 'study_track_trgmu'

  signal_filenames = [
    #'/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_phi_selected_morevars.root',
    #'/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/flat_bparknano.root',
    '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_study_phi_ct.root'
    ]

  data_filenames = [
    #'/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/data/study_preselection_v0/ParkingBPH1_Run2018A/Chunk0_n419/merged/bparknano_phi_selected_morevars.root',
    #'/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/data/study_preselection_v0/ParkingBPH1_Run2018A/merged/bparknano_phi_selected_morevars.root',
    #'/eos/user/a/anlyon/CPVGen/V00_D1/ParkingBPH1/crab_20250225_225643/250225_215706/0000/merged/flat_bparknano.root',
    '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/data/V02/2018/D1/merged/small_bparknano_study_phi_ct.root',
    ]

  quantities = quantities['track_study']
  isnano = True

  for quantity in quantities:
    plotter = Plotter(
        quantity = quantity, 
        data_filenames = data_filenames,
        signal_filenames = signal_filenames,
        outdirlabel = outdirlabel,
        isnano = isnano,
        )
          
    plotter.plot()

    #quantities = [
    #   #Quantity(name_flat='bs_mass_corr', title='bs_mass_corr', label='bs_mass_corr', nbins=50, bin_min=5.05, bin_max=5.7),

    #   ##Quantity(name_flat='bs_beta', title='bs_beta', label='bs_beta', nbins=50, bin_min=0, bin_max=1),
    #   ##Quantity(name_flat='bs_beta', title='bs_beta', label='bs_beta', nbins=50, bin_min=0, bin_max=1, do_log=True),
    #   ##Quantity(name_flat='bs_eta', title='bs_eta', label='bs_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_charge ', title='bs_charge ', label='bs_charge ', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='bs_cos2d', title='bs_cos2d', label='bs_cos2d', nbins=50, bin_min=0.9, bin_max=1, do_log=True),

    #   ##Quantity(name_flat='bs_cxx', title='bs_cxx', label='bs_cxx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_cyx', title='bs_cyx', label='bs_cyx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_cyy', title='bs_cyy', label='bs_cyy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_czx', title='bs_czx', label='bs_czx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_czy', title='bs_czy', label='bs_czy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_czz', title='bs_czz', label='bs_czz', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='bs_lxysig', title='bs_lxysig', label='bs_lxysig', nbins=50, bin_min=0, bin_max=50, do_log=True),

    #   ##Quantity(name_flat='bs_mass', title='bs_mass', label='bs_mass', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_mass_err', title='bs_mass_err', label='bs_mass_err', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_mass_corr', title='bs_mass_corr', label='bs_mass_corr', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_phi', title='bs_phi', label='bs_phi', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='bs_pt', title='bs_pt', label='bs_pt', nbins=50, bin_min=0, bin_max=50, do_log=True),

    #   ##Quantity(name_flat='bs_sv_prob', title='bs_sv_prob', label='bs_sv_prob', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_sv_chi2', title='bs_sv_chi2', label='bs_sv_chi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_sv_ndof', title='bs_sv_ndof', label='bs_sv_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_vx', title='bs_vx', label='bs_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_vy', title='bs_vy', label='bs_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='bs_vz', title='bs_vz', label='bs_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='cos_theta_star_phi1', title='cos_theta_star_phi1', label='cos_theta_star_phi1', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='cos_theta_star_phi2', title='cos_theta_star_phi2', label='cos_theta_star_phi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='deltar_k1k2', title='deltar_k1k2', label='deltar_k1k2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='deltar_k1k3', title='deltar_k1k3', label='deltar_k1k3', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='deltar_k1k4', title='deltar_k1k4', label='deltar_k1k4', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='deltar_k2k3', title='deltar_k2k3', label='deltar_k2k3', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='deltar_k2k4', title='deltar_k2k4', label='deltar_k2k4', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='deltar_k3k4', title='deltar_k3k4', label='deltar_k3k4', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='deltar_min', title='deltar_min', label='deltar_min', nbins=50, bin_min=0, bin_max=0.2),
    #   #Quantity(name_flat='deltar_max', title='deltar_max', label='deltar_max', nbins=50, bin_min=0, bin_max=3),
    #   #Quantity(name_flat='deltar_phi1phi2', title='deltar_phi1phi2', label='deltar_phi1phi2', nbins=50, bin_min=0, bin_max=3),

    #   ##Quantity(name_flat='k1k3_mass', title='k1k3_mass', label='k1k3_mass', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1k3_pt', title='k1k3_pt', label='k1k3_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1k4_mass', title='k1k4_mass', label='k1k4_mass', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1k4_pt', title='k1k4_pt', label='k1k4_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2k3_mass', title='k2k3_mass', label='k2k3_mass', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2k3_pt', title='k2k3_pt', label='k2k3_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2k4_mass', title='k2k4_mass', label='k2k4_mass', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2k4_pt', title='k2k4_pt', label='k2k4_pt', nbins=50, bin_min=, bin_max=),

    #   ##Quantity(name_flat='phi1_eta', title='phi1_eta', label='phi1_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_phi', title='phi1_phi', label='phi1_phi', nbins=50, bin_min=, bin_max=),

    #   Quantity(name_flat='phi1_pt', title='phi1_pt', label='phi1_pt', nbins=50, bin_min=0, bin_max=30, do_log=True),

    #   #Quantity(name_flat='phi1_mass', title='phi1_mass', label='phi1_mass', nbins=50, bin_min=1.01, bin_max=1.03),
    #   #Quantity(name_flat='phi1_cos_theta_star_k1', title='phi1_cos_theta_star_k1', label='phi1_cos_theta_star_k1', nbins=50, bin_min=0, bin_max=1),
    #   ##Quantity(name_flat='phi1_cos_theta_star_k2', title='phi1_cos_theta_star_k2', label='phi1_cos_theta_star_k2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_deltar', title='phi1_deltar', label='phi1_deltar', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='phi1_cos2d', title='phi1_cos2d', label='phi1_cos2d', nbins=50, bin_min=-1, bin_max=1, do_log=True),

    #   ##Quantity(name_flat='phi1_cxx', title='phi1_cxx', label='phi1_cxx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_cyx', title='phi1_cyx', label='phi1_cyx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_cyy', title='phi1_cyy', label='phi1_cyy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_czx', title='phi1_czx', label='phi1_czx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_czy', title='phi1_czy', label='phi1_czy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_czz', title='phi1_czz', label='phi1_czz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_sv_chi2', title='phi1_sv_chi2', label='phi1_sv_chi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_sv_ndof', title='phi1_sv_ndof', label='phi1_sv_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_sv_prob', title='phi1_sv_prob', label='phi1_sv_prob', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_lxy', title='phi1_lxy', label='phi1_lxy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_lxysig', title='phi1_lxysig', label='phi1_lxysig', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='phi1_mass_err', title='phi1_mass_err', label='phi1_mass_err', nbins=50, bin_min=0, bin_max=0.15, do_log=True),

    #   ##Quantity(name_flat='phi1_vx', title='phi1_vx', label='phi1_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_vy', title='phi1_vy', label='phi1_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_vz', title='phi1_vz', label='phi1_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_eta', title='phi2_eta', label='phi2_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_phi', title='phi2_phi', label='phi2_phi', nbins=50, bin_min=, bin_max=),

    #   #Quantity(name_flat='phi2_pt', title='phi2_pt', label='phi2_pt', nbins=50, bin_min=0, bin_max=20, do_log=True),

    #   ##Quantity(name_flat='phi2_mass', title='phi2_mass', label='phi2_mass', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_cos_theta_star_k1', title='phi2_cos_theta_star_k1', label='phi2_cos_theta_star_k1', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_cos_theta_star_k2', title='phi2_cos_theta_star_k2', label='phi2_cos_theta_star_k2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_deltar', title='phi2_deltar', label='phi2_deltar', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_cos2d', title='phi2_cos2d', label='phi2_cos2d', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_cxx', title='phi2_cxx', label='phi2_cxx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_cyx', title='phi2_cyx', label='phi2_cyx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_cyy', title='phi2_cyy', label='phi2_cyy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_czx', title='phi2_czx', label='phi2_czx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_czy', title='phi2_czy', label='phi2_czy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_czz', title='phi2_czz', label='phi2_czz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_sv_chi2', title='phi2_sv_chi2', label='phi2_sv_chi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_sv_ndof', title='phi2_sv_ndof', label='phi2_sv_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_sv_prob', title='phi2_sv_prob', label='phi2_sv_prob', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_lxy', title='phi2_lxy', label='phi2_lxy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_lxysig', title='phi2_lxysig', label='phi2_lxysig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_mass_err', title='phi2_mass_err', label='phi2_mass_err', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_vx', title='phi2_vx', label='phi2_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_vy', title='phi2_vy', label='phi2_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi2_vz', title='phi2_vz', label='phi2_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='phi1_pt_times_phi2_pt', title='phi1_pt_times_phi2_pt', label='phi1_pt_times_phi2_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_eta', title='k1_eta', label='k1_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_phi', title='k1_phi', label='k1_phi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_pt', title='k1_pt', label='k1_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_dcasig', title='k1_dcasig', label='k1_dcasig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_charge', title='k1_charge', label='k1_charge', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_covlamlam', title='k1_covlamlam', label='k1_covlamlam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_covlamphi', title='k1_covlamphi', label='k1_covlamphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_covphiphi', title='k1_covphiphi', label='k1_covphiphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_covqoplam', title='k1_covqoplam', label='k1_covqoplam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_covqopphi', title='k1_covqopphi', label='k1_covqopphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_covqopqop', title='k1_covqopqop', label='k1_covqopqop', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_dxy', title='k1_dxy', label='k1_dxy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_dxysig', title='k1_dxysig', label='k1_dxysig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_dz', title='k1_dz', label='k1_dz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_dzsig', title='k1_dzsig', label='k1_dzsig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_normalisedchi2', title='k1_normalisedchi2', label='k1_normalisedchi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_pt_err', title='k1_pt_err', label='k1_pt_err', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_validfraction', title='k1_validfraction', label='k1_validfraction', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_vx', title='k1_vx', label='k1_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_vy', title='k1_vy', label='k1_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_vz', title='k1_vz', label='k1_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_highpurityflag', title='k1_highpurityflag', label='k1_highpurityflag', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_islost', title='k1_islost', label='k1_islost', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_ispacked', title='k1_ispacked', label='k1_ispacked', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_ndof', title='k1_ndof', label='k1_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_numberlosthits', title='k1_numberlosthits', label='k1_numberlosthits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_numberpixelhits', title='k1_numberpixelhits', label='k1_numberpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_numberpixellayers', title='k1_numberpixellayers', label='k1_numberpixellayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_numbertrackerlayers', title='k1_numbertrackerlayers', label='k1_numbertrackerlayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_numberofvalidhits', title='k1_numberofvalidhits', label='k1_numberofvalidhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_numberofvalidpixelhits', title='k1_numberofvalidpixelhits', label='k1_numberofvalidpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k1_qualityindex', title='k1_qualityindex', label='k1_qualityindex', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_eta', title='k2_eta', label='k2_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_phi', title='k2_phi', label='k2_phi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_pt', title='k2_pt', label='k2_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_dcasig', title='k2_dcasig', label='k2_dcasig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_charge', title='k2_charge', label='k2_charge', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_covlamlam', title='k2_covlamlam', label='k2_covlamlam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_covlamphi', title='k2_covlamphi', label='k2_covlamphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_covphiphi', title='k2_covphiphi', label='k2_covphiphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_covqoplam', title='k2_covqoplam', label='k2_covqoplam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_covqopphi', title='k2_covqopphi', label='k2_covqopphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_covqopqop', title='k2_covqopqop', label='k2_covqopqop', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_dxy', title='k2_dxy', label='k2_dxy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_dxysig', title='k2_dxysig', label='k2_dxysig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_dz', title='k2_dz', label='k2_dz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_dzsig', title='k2_dzsig', label='k2_dzsig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_normalisedchi2', title='k2_normalisedchi2', label='k2_normalisedchi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_pt_err', title='k2_pt_err', label='k2_pt_err', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_validfraction', title='k2_validfraction', label='k2_validfraction', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_vx', title='k2_vx', label='k2_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_vy', title='k2_vy', label='k2_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_vz', title='k2_vz', label='k2_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_highpurityflag', title='k2_highpurityflag', label='k2_highpurityflag', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_islost', title='k2_islost', label='k2_islost', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_ispacked', title='k2_ispacked', label='k2_ispacked', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_ndof', title='k2_ndof', label='k2_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_numberlosthits', title='k2_numberlosthits', label='k2_numberlosthits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_numberpixelhits', title='k2_numberpixelhits', label='k2_numberpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_numberpixellayers', title='k2_numberpixellayers', label='k2_numberpixellayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_numbertrackerlayers', title='k2_numbertrackerlayers', label='k2_numbertrackerlayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_numberofvalidhits', title='k2_numberofvalidhits', label='k2_numberofvalidhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_numberofvalidpixelhits', title='k2_numberofvalidpixelhits', label='k2_numberofvalidpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k2_qualityindex', title='k2_qualityindex', label='k2_qualityindex', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_eta', title='k3_eta', label='k3_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_phi', title='k3_phi', label='k3_phi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_pt', title='k3_pt', label='k3_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_dcasig', title='k3_dcasig', label='k3_dcasig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_charge', title='k3_charge', label='k3_charge', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_covlamlam', title='k3_covlamlam', label='k3_covlamlam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_covlamphi', title='k3_covlamphi', label='k3_covlamphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_covphiphi', title='k3_covphiphi', label='k3_covphiphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_covqoplam', title='k3_covqoplam', label='k3_covqoplam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_covqopphi', title='k3_covqopphi', label='k3_covqopphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_covqopqop', title='k3_covqopqop', label='k3_covqopqop', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_dxy', title='k3_dxy', label='k3_dxy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_dxysig', title='k3_dxysig', label='k3_dxysig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_dz', title='k3_dz', label='k3_dz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_dzsig', title='k3_dzsig', label='k3_dzsig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_normalisedchi2', title='k3_normalisedchi2', label='k3_normalisedchi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_pt_err', title='k3_pt_err', label='k3_pt_err', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_validfraction', title='k3_validfraction', label='k3_validfraction', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_vx', title='k3_vx', label='k3_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_vy', title='k3_vy', label='k3_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_vz', title='k3_vz', label='k3_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_highpurityflag', title='k3_highpurityflag', label='k3_highpurityflag', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_islost', title='k3_islost', label='k3_islost', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_ispacked', title='k3_ispacked', label='k3_ispacked', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_ndof', title='k3_ndof', label='k3_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_numberlosthits', title='k3_numberlosthits', label='k3_numberlosthits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_numberpixelhits', title='k3_numberpixelhits', label='k3_numberpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_numberpixellayers', title='k3_numberpixellayers', label='k3_numberpixellayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_numbertrackerlayers', title='k3_numbertrackerlayers', label='k3_numbertrackerlayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_numberofvalidhits', title='k3_numberofvalidhits', label='k3_numberofvalidhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_numberofvalidpixelhits', title='k3_numberofvalidpixelhits', label='k3_numberofvalidpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k3_qualityindex', title='k3_qualityindex', label='k3_qualityindex', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_eta', title='k4_eta', label='k4_eta', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_phi', title='k4_phi', label='k4_phi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_pt', title='k4_pt', label='k4_pt', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_dcasig', title='k4_dcasig', label='k4_dcasig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_charge', title='k4_charge', label='k4_charge', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_covlamlam', title='k4_covlamlam', label='k4_covlamlam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_covlamphi', title='k4_covlamphi', label='k4_covlamphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_covphiphi', title='k4_covphiphi', label='k4_covphiphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_covqoplam', title='k4_covqoplam', label='k4_covqoplam', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_covqopphi', title='k4_covqopphi', label='k4_covqopphi', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_covqopqop', title='k4_covqopqop', label='k4_covqopqop', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_dxy', title='k4_dxy', label='k4_dxy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_dxysig', title='k4_dxysig', label='k4_dxysig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_dz', title='k4_dz', label='k4_dz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_dzsig', title='k4_dzsig', label='k4_dzsig', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_normalisedchi2', title='k4_normalisedchi2', label='k4_normalisedchi2', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_pt_err', title='k4_pt_err', label='k4_pt_err', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_validfraction', title='k4_validfraction', label='k4_validfraction', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_vx', title='k4_vx', label='k4_vx', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_vy', title='k4_vy', label='k4_vy', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_vz', title='k4_vz', label='k4_vz', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_highpurityflag', title='k4_highpurityflag', label='k4_highpurityflag', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_islost', title='k4_islost', label='k4_islost', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_ispacked', title='k4_ispacked', label='k4_ispacked', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_ndof', title='k4_ndof', label='k4_ndof', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_numberlosthits', title='k4_numberlosthits', label='k4_numberlosthits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_numberpixelhits', title='k4_numberpixelhits', label='k4_numberpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_numberpixellayers', title='k4_numberpixellayers', label='k4_numberpixellayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_numbertrackerlayers', title='k4_numbertrackerlayers', label='k4_numbertrackerlayers', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_numberofvalidhits', title='k4_numberofvalidhits', label='k4_numberofvalidhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_numberofvalidpixelhits', title='k4_numberofvalidpixelhits', label='k4_numberofvalidpixelhits', nbins=50, bin_min=, bin_max=),
    #   ##Quantity(name_flat='k4_qualityindex', title='k4_qualityindex', label='k4_qualityindex', nbins=50, bin_min=, bin_max=),


    #]


    #for quantity in quantities:
    #  plotter = Plotter(
    #      quantity = quantity, 
    #      data_filenames = data_filenames,
    #      signal_filenames = signal_filenames,
    #      outdirlabel = outdirlabel,
    #      builder = 'Bs',
    #      )
    #      
    #  plotter.plot()

