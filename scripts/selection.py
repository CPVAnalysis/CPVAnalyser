import ROOT
import numpy as np
import os
from os import path
from tools import Tools
import sys
import math
from plotter import Quantity as Qte
from glob import glob

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to study the selection', add_help=True)
  parser.add_argument('--dirlabel'             , type=str, dest='dirlabel'            , help='name of the outdir'  , default=None)
  parser.add_argument('--mass'                 , type=str, dest='mass'                , help='mass'                , default=None)
  parser.add_argument('--category'             , type=str, dest='category'            , help='category'            , default=None)
  parser.add_argument('--quantity'             , type=str, dest='quantity'            , help='quantity'            , default=None)
  parser.add_argument('--action'               , type=str, dest='action'              , help='action'              , default=None)
  parser.add_argument('--cut_b_mass'           , type=str, dest='cut_b_mass'          , help='cut_b_mass'          , default=None)
  parser.add_argument('--cut_pi_pt'            , type=str, dest='cut_pi_pt'           , help='cut_pi_pt'           , default=None)
  parser.add_argument('--cut_sv_lxysig'        , type=str, dest='cut_sv_lxysig'       , help='cut_sv_lxysig'       , default=None)
  parser.add_argument('--cut_mu_dxysig'        , type=str, dest='cut_mu_dxysig'       , help='cut_mu_dxysig'       , default=None)
  parser.add_argument('--cut_pi_dxysig'        , type=str, dest='cut_pi_dxysig'       , help='cut_pi_dxysig'       , default=None)
  parser.add_argument('--cut_max_mu_pi_dxysig' , type=str, dest='cut_max_mu_pi_dxysig', help='cut_max_mu_pi_dxysig', default=None)
  parser.add_argument('--cut_min_mu_pi_dxysig' , type=str, dest='cut_min_mu_pi_dxysig', help='cut_min_mu_pi_dxysig', default=None)
  parser.add_argument('--cut_hnl_cos2d'        , type=str, dest='cut_hnl_cos2d'       , help='cut_hnl_cos2d'       , default=None)
  return parser.parse_args()


class Quantity(object):
  def __init__(self, name_flat='', label='', title='', logic='', units='', binMin=0., binMax=0.):
    self.name_flat = name_flat
    self.label = label
    self.title = title
    self.logic = logic
    self.units = units
    self.binMin = binMin
    self.binMax = binMax


class SelectedQuantity(object):
  '''
  this class will allow to impose a selection when studying 
  the impact of the cut on a given parameter
  '''
  def __init__(self, quantity, chosen_cut):
    self.quantity   = quantity
    self.chosen_cut = chosen_cut

#TODO apply weights on signal, move dxysig definition to BS (DCASig for pions), add selection on the b mass difference with mB

class Selection(object):
  def __init__(self, signal_files='', data_files='', baseline_selection='', quantity='', dirlabel='', preexisting_selection=None, proposed_cut=None):
    self.tools                 = Tools()
    self.signal_files          = signal_files
    self.data_files            = data_files
    self.baseline_selection    = baseline_selection
    self.quantity              = quantity
    self.dirlabel              = dirlabel
    self.preexisting_selection = preexisting_selection
    self.proposed_cut          = proposed_cut

    #self.weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysigbs_max5e6_v2_smalltable_v2'
    #print 'weight_hlt: {}'.format(self.weight_hlt)
    #self.weight_pusig = 'weight_pu_sig_D'
    #self.weight_mu0id = 'weight_mu0_softid'
    #self.weight_muid = 'weight_mu_looseid'

    #self.resolution_p0 = 0.0002747
    #self.resolution_p1 = 0.008302


  def createOutDir(self, outputdir):
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))


  def getLabel(self, str_):
    new_str_ = str_
    if 'abs(' in str_: new_str_ = new_str_.replace('abs(', '')
    if ')'in str_: new_str_ = new_str_.replace(')', '')
    return new_str_


  def getSelectionString(self):
    '''
    function to write the already-applied cuts in a string
    '''
    preselection_str = []
    for item, _ in enumerate(self.preexisting_selection):
      name_variable = self.preexisting_selection[item].quantity.name_flat
      preselection_str.append('{}{}{}'.format(name_variable,self.preexisting_selection[item].quantity.logic,self.preexisting_selection[item].chosen_cut))
    return ' && '.join(preselection_str)


  def getEff(self, entries_selected, entries_initial):
    eff = entries_selected / entries_initial if entries_initial!=0 else 0.
    return eff



  def getSignificanceDiffScan(self, npoints=20):
    '''
    plot the significance gain/loss as a function of the cut applied on the quantity of interest
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, npoints) 
      
    canv = self.tools.createTCanvas('canv', 900, 800)

    ROOT.gStyle.SetPadLeftMargin(0.8)
    ROOT.gStyle.SetPadRightMargin(0.8)
    
    gs_sig = []

    significance_diff_max = -1
    significance_diff_min = -1

    for signal_file in self.signal_files:

      g_sig = ROOT.TGraph()
      self.bs_mass_corr = Qte(name_flat='bs_mass_corr', label='bs_mass_corr', nbins=80, bin_min=5.05, bin_max=5.7)

      # compute the initial significance
      f_sig = ROOT.TFile.Open(signal_file, 'READ')
      tree_sig = self.tools.getTree(f_sig, 'signal_tree')
      signal_selection_ini = 'ismatched == 1' 
      if self.preexisting_selection !=None: signal_selection_ini += ' && {}'.format(self.getSelectionString())
      hist_sig_name = 'hist_sig_ini_{}'.format(self.quantity)
      hist_sig = self.tools.createHisto(tree_sig, self.bs_mass_corr, hist_name=hist_sig_name, branchname='flat', selection=signal_selection_ini)
      signal_yields_ini = hist_sig.Integral()

      background_yields_ini = 0
      for data_file in self.data_files:
         f_data = self.tools.getRootFile(data_file)
         tree_data = self.tools.getTree(f_data, 'signal_tree')
         background_selection_ini = 'fabs(bs_mass_corr - 5.367) > 0.2'
         if self.preexisting_selection !=None: background_selection_ini += ' && {}'.format(self.getSelectionString())
         hist_data_name = 'hist_data_ini_{}'.format(self.quantity)
         hist_data = self.tools.createHisto(tree_data, self.bs_mass_corr, hist_name=hist_data_name, branchname='flat', selection=background_selection_ini)
         background_yields_ini += hist_data.Integral()

      if background_yields_ini != 0:
        significance_ini = signal_yields_ini / math.sqrt(background_yields_ini)
      else:
        significance_ini = -99.

      for idx, cut in enumerate(points):
        print cut
        # compute signal yields
        signal_selection_sel = 'ismatched == 1' 
        if self.preexisting_selection !=None: signal_selection_sel += ' && {}'.format(self.getSelectionString())
        signal_selection_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
        f_sig = ROOT.TFile.Open(signal_file, 'READ')
        tree_sig = self.tools.getTree(f_sig, 'signal_tree')
        hist_sig_name = 'hist_sig_sel_{}'.format(self.quantity)
        hist_sig = self.tools.createHisto(tree_sig, self.bs_mass_corr, hist_name=hist_sig_name, branchname='flat', selection=signal_selection_sel)
        signal_yields_sel = hist_sig.Integral()

        # get background yields from unblinded data
        background_yields_sel = 0
        for data_file in self.data_files:
           f_data = self.tools.getRootFile(data_file)
           tree_data = self.tools.getTree(f_data, 'signal_tree')
           background_selection_sel = 'fabs(bs_mass_corr - 5.367) > 0.2'
           if self.preexisting_selection !=None: background_selection_sel += ' && {}'.format(self.getSelectionString())
           background_selection_sel += ' && {}{}{}'.format(self.quantity.name_flat, self.quantity.logic, cut)
           hist_data_name = 'hist_data_sel_{}_{}'.format(self.quantity, cut)
           hist_data = self.tools.createHisto(tree_data, self.bs_mass_corr, hist_name=hist_data_name, branchname='flat', selection=background_selection_sel)
           background_yields_sel += hist_data.Integral()

        # compute significance
        if background_yields_sel != 0:
          significance_sel = signal_yields_sel / math.sqrt(background_yields_sel)
        else:
          significance_sel = -99.

        # compute relative difference in significance
        significance_diff = ((significance_sel / significance_ini) - 1) * 100

        if significance_diff_max == -1 or significance_diff > significance_diff_max:
          significance_diff_max = significance_diff

        if significance_diff_min == -1 or significance_diff < significance_diff_min:
          significance_diff_min = significance_diff

        if significance_diff_min < -50: significance_diff_min = -50

        g_sig.SetPoint(idx, cut, significance_diff)
        g_sig.SetLineWidth(2)
        g_sig.SetLineStyle(9)
        g_sig.SetLineColor(ROOT.kOrange+1)
        g_sig.SetMarkerColor(ROOT.kOrange+1)
        g_sig.SetMarkerStyle(20)

      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on {}'.format(self.quantity.title))
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetTitle('Relative difference in significance [%]')
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetTitleSize(0.045)
      gs_sig[0].GetYaxis().SetTitleOffset(1.15)
      gs_sig[0].GetYaxis().SetRangeUser(significance_diff_min-0.15*significance_diff_min, significance_diff_max+0.15*significance_diff_max)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          g_sig.Draw('APL')
        else:
          g_sig.Draw('PL')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, significance_diff_min-0.15*significance_diff_min, self.proposed_cut, significance_diff_max+0.15*significance_diff_max)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      line.Draw('same')
    
    #legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.25, xmax=0.82, ymax=0.4, size=0.04)
    #for ifile, signal_file in enumerate(self.signal_files):
    #  legend.AddEntry(gs_sig[ifile], '{}GeV, {}mm'.format(signal_file.mass, signal_file.ctau))
    #legend.Draw()

    self.createOutDir('myPlots/selection/{}'.format(self.dirlabel))

    canv.cd()
    canv.SaveAs('myPlots/selection/{}/significance_diff_scan_{}.png'.format(self.dirlabel, self.getLabel(self.quantity.label)))



if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

    
  phi1_pt = Quantity('phi1_pt', 'phi1_pt', 'pt(phi1)', '>', 'GeV', 0, 10) 
  phi2_pt = Quantity('phi2_pt', 'phi2_pt', 'pt(phi2)', '>', 'GeV', 0, 10) 
  deltar_min = Quantity('deltar_min', 'deltar_min', 'min(deltaR)', '<', 'GeV', 0.01, 0.1) 
  deltar_max = Quantity('deltar_max', 'deltar_max', 'max(deltaR)', '<', 'GeV', 0.4, 1) 
  bs_lxysig = Quantity('bs_lxysig', 'Bs_lxy_sig', 'lxy_sig(Bs)', '>', 'GeV', 0, 50) 
  bs_cos2d = Quantity('bs_cos2d', 'bs_cos2d', 'cos2D(Bs)', '>', 'GeV', 0.95, 1) 
  sv_prob = Quantity('bs_sv_prob', 'bs_sv_prob', 'sv_prob(Bs)', '>', 'GeV', 0, 0.5) 

  #hnl_charge = Quantity('hnl_charge', 'hnl_charge', '#mu#pi charge', '==', '', -2, 2)
  ##hnl_cos2D = Quantity('hnl_cos2d', 'hnl_cos2D', 'SV cos2D', '>', '', 0.9993, 1) 
  #pi_pt = Quantity('pi_pt', 'pi_pt', 'displaced pion pT', '>', 'GeV', 0, 8)
  #mu_dxy = Quantity('abs(mu_dxy)', 'mu_dxy', 'displaced muon |IP| on xy', '>', 'cm', 0, 0.3)
  #pi_dxy = Quantity('abs(pi_dxy)', 'pi_dxy', 'displaced pion |IP| on xy', '>', 'cm', 0, 1)
  #max_mu_pi_dxysig = Quantity('max(abs(pi_dxysig), abs(mu_dxysig))', 'max_mu_pi_dxysig', 'max(muon, pion) |IP| significance on xy', '>', 'cm', 0, 80)
  #min_mu_pi_dxysig = Quantity('min(abs(pi_dcasig), abs(mu_dxysig_bs))', 'min_mu_pi_dxysig', 'min(muon, pion) |IP| significance on xy (BS)', '>', 'cm', 0, 200)
  #mu_dxysig = Quantity('abs(mu_dxysig)', 'mu_dxysig', 'displaced muon |IP| significance on xy', '>', 'cm', 0, 50)
  #pi_dxysig = Quantity('abs(pi_dxysig)', 'pi_dxysig', 'displaced pion |IP| significance on xy', '>', 'cm', 0, 100)
  #pi_dcasig = Quantity('abs(pi_dcasig)', 'pi_dcasig', 'pion DCA significance', '>', 'cm', 0, 50)
  #cos_theta_star_pion = Quantity('fabs(cos_theta_star_pion)', 'cos_theta_star_pion', '|cos(#theta_{#pi}*)|', '<', '', 0.5, 1) 
  #mu0_softid = Quantity('mu0_softid', 'mu0_softid', 'trigger muon soft ID', '==', '', 0, 1)
  #mu_looseid = Quantity('mu_looseid', 'mu_looseid', 'muon loose ID', '==', '', 0, 1)
  ##mu_customid = Quantity('mu_intimemuon==1 && mu_trackerhighpurityflag==1 && ((mu_isglobalmuon==1 && mu_numberofstations>0 && mu_numberoftrackerlayers<18) || (mu_isglobalmuon!=1 && mu_calocompatibility>0.05 && mu_numberoftrackerlayers>6 && mu_numberoftrackerlayers<16 && mu_numberofvalidpixelhits<6))', 'mu_customid', 'muon custom ID', '==', '', 0, 1)
  #mu_customid = Quantity('mu_customisedid', 'mu_customisedid', 'muon custom ID', '==', '', 0, 1)
  #mu_whnlid = Quantity('mu_whnlid', 'mu_whnlid', 'muon whnl ID', '==', '', 0, 1)
  #b_mass = Quantity('abs(b_mass-5.3)', 'b_mass', '#mu#mu#pi mass', '<', 'GeV', 0, 1.5) 
  #sv_lxy_sig = Quantity('sv_lxysig', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 500)
  #mu_numberofpixellayers = Quantity('mu_numberofpixellayers', 'mu_numberofpixellayers', 'mu_numberofpixellayers', '<', '', 0, 6)
  #mu_numberoftrackerlayers = Quantity('mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', 'mu_numberoftrackerlayers', '<', '', 0, 18)
  #deltar_trgmu_pi = Quantity('abs(deltar_trgmu_pi)', 'abs(deltar_trgmu_pi)', '|#DeltaR(trg#mu, #pi)|', '<', '', 0.3, 1.) 
  #deltar_trgmu_hnl = Quantity('abs(deltar_trgmu_hnl)', 'abs(deltar_trgmu_hnl)', '|#DeltaR(trg#mu, hnl)|', '<', '', 0., 0.5) 
  #mu0_pfiso03_rel = Quantity('mu0_pfiso03_rel', 'mu0_pfiso03_rel', 'primary muon relative PF iso03', '<', '', 0.3, 3)
  #mu0_pfiso04_rel = Quantity('mu0_pfiso04_rel', 'mu0_pfiso04_rel', 'primary muon relative PF iso04', '<', '', 0.3, 3)
  #mu_pfiso03_rel = Quantity('mu_pfiso03_rel', 'mu_pfiso03_rel', 'displaced muon relative PF iso03', '<', '', 0., 0.5)
  #mu_pfiso04_rel = Quantity('mu_pfiso04_rel', 'mu_pfiso04_rel', 'displaced muon relative PF iso04', '<', '', 0.3, 3)
  
  signal_files = [
     '/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/flat_bparknano.root',
  ]
   
  data_files = [f for f in glob('/eos/user/a/anlyon/CPVGen/V00_*/*/*/*/*/flat/flat_bparknano.root')]
  for f in data_files:
    print f


  opt = getOptions()

  do_scan_significance_diff = True

  dirlabel = 'study_sig_diff_scan'


  selection_sequential = [] 

  cut_phi1_pt = 6.5
  selection = Selection(signal_files=signal_files, data_files=data_files, quantity=phi1_pt, dirlabel=dirlabel, proposed_cut=cut_phi1_pt)
  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  selection_sequential.append(SelectedQuantity(phi1_pt, cut_phi1_pt))

  cut_phi2_pt = 5
  selection = Selection(signal_files=signal_files, data_files=data_files, quantity=phi2_pt, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_phi2_pt)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  selection_sequential.append(SelectedQuantity(phi2_pt, cut_phi2_pt))

  #cut_deltar_min = 0.05
  #selection = Selection(signal_files=signal_files, data_files=data_files, quantity=deltar_min, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_deltar_min)
  ##if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  #selection_sequential.append(SelectedQuantity(deltar_min, cut_deltar_min))

  #cut_deltar_max = 0.75
  #selection = Selection(signal_files=signal_files, data_files=data_files, quantity=deltar_max, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_deltar_max)
  ##if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  #selection_sequential.append(SelectedQuantity(deltar_max, cut_deltar_max))

  cut_bs_lxysig = 15
  selection = Selection(signal_files=signal_files, data_files=data_files, quantity=bs_lxysig, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_bs_lxysig)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  selection_sequential.append(SelectedQuantity(bs_lxysig, cut_bs_lxysig))

  cut_bs_cos2d = 0.995
  selection = Selection(signal_files=signal_files, data_files=data_files, quantity=bs_cos2d, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_bs_cos2d)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  selection_sequential.append(SelectedQuantity(bs_cos2d, cut_bs_cos2d))

  cut_sv_prob = 0.4
  selection = Selection(signal_files=signal_files, data_files=data_files, quantity=sv_prob, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_sv_prob)
  #if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=20)
  selection_sequential.append(SelectedQuantity(sv_prob, cut_sv_prob))

 #cut_mu_looseid = 1
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_looseid, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_looseid)
 ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
 #if do_print_cutflow: selection.printCutflowLine()
 #selection_sequential.append(SelectedQuantity(mu_looseid, cut_mu_looseid))

 #cut_mu_whnlid = 1
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_whnlid, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_whnlid)
 ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
 #if do_print_cutflow: selection.printCutflowLine()
 ###selection_sequential.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

 #cut_mu_customid = 1
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_customid, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_customid)
 ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
 #if do_print_cutflow: selection.printCutflowLine()
 ###selection_sequential.append(SelectedQuantity(mu_customid, cut_mu_customid)) # for the moment we don't use it

 #cut_hnl_charge = 0
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_charge, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_hnl_charge)
 ##if do_scan_efficiency: selection.getScanGraph(npoints=2)
 #if do_print_cutflow: selection.printCutflowLine(printNum=False)
 #selection_sequential.append(SelectedQuantity(hnl_charge, cut_hnl_charge))

 #cut_b_mass = float(opt.cut_b_mass)
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=b_mass, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_b_mass)
 #if do_b_mass:
 #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 #if do_print_cutflow: selection.printCutflowLine()
 #selection_sequential.append(SelectedQuantity(b_mass, cut_b_mass))

 #cut_pi_pt = float(opt.cut_pi_pt)
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_pt, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_pi_pt)
 #if do_pi_pt:
 #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 #if do_print_cutflow: selection.printCutflowLine()
 #selection_sequential.append(SelectedQuantity(pi_pt, cut_pi_pt))

 #cut_sv_lxy_sig = float(opt.cut_sv_lxysig)
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=sv_lxy_sig, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_sv_lxy_sig)
 #if do_sv_lxysig:
 #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 #if do_print_cutflow: selection.printCutflowLine()
 #selection_sequential.append(SelectedQuantity(sv_lxy_sig, cut_sv_lxy_sig))

 ##cut_mu_dxysig = float(opt.cut_mu_dxysig)
 ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=mu_dxysig, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_mu_dxysig)
 ##if do_mu_dxysig:
 ##  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 ##  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 ##  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 ##if do_print_cutflow: selection.printCutflowLine()
 ###selection_sequential.append(SelectedQuantity(mu_dxysig, cut_mu_dxysig))

 ##cut_pi_dxysig = float(opt.cut_pi_dxysig)
 ##selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dxysig, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_pi_dxysig)
 ##if do_pi_dxysig:
 ##  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 ##  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 ##  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 ##if do_print_cutflow: selection.printCutflowLine(printNum=False)
 ###selection_sequential.append(SelectedQuantity(pi_dxysig, cut_pi_dxysig))

 #cut_min_mu_pi_dxysig = float(opt.cut_min_mu_pi_dxysig)
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=min_mu_pi_dxysig, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_min_mu_pi_dxysig)
 #if do_min_mu_pi_dxysig:
 #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 #if do_print_cutflow: selection.printCutflowLine(printNum=False)
 #selection_sequential.append(SelectedQuantity(min_mu_pi_dxysig, cut_min_mu_pi_dxysig))

 #cut_hnl_cos2D = float(opt.cut_hnl_cos2d)
 #selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=hnl_cos2D, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_hnl_cos2D)
 #if do_hnl_cos2d:
 #  if do_scan_efficiency: selection.getScanGraph(npoints=30)
 #  if do_scan_significance: selection.getSignificanceScan(npoints=30)
 #  if do_scan_significance_diff: selection.getSignificanceDiffScan(npoints=30)
 #if do_print_cutflow: selection.printCutflowLine()
 #selection_sequential.append(SelectedQuantity(hnl_cos2D, cut_hnl_cos2D))

 ###cut_deltar_trgmu_pi = 1.5
 ###selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=deltar_trgmu_pi, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_deltar_trgmu_pi)
 ###if do_scan_efficiency: selection.getScanGraph(npoints=30)
 ###selection.printCutflowLine(printNum=False)
 ###selection_sequential.append(SelectedQuantity(deltar_trgmu_pi, cut_deltar_trgmu_pi))

 ###cut_pi_dcasig = 20
 ###selection = Selection(signal_label=signal_label, signal_files=signal_files, data_file=data_file, qcd_files=qcd_files, baseline_selection=baseline_selection, white_list=white_list, quantity=pi_dcasig, category=category, dirlabel=dirlabel, preexisting_selection=selection_sequential, proposed_cut=cut_pi_dcasig)
 ###if do_scan_efficiency: selection.getScanGraph(npoints=30)
 ####selection.printCutflowLine(printNum=False)
 ###selection_sequential.append(SelectedQuantity(pi_dcasig, cut_pi_dcasig))



