import numpy as np
import ROOT
import os
import os.path
from os import path
from tools import Tools


#TODO select bkg in sidebands only

class Quantity(object):
  def __init__(self, name_nano='', name_flat='', label='', title='', logic='', units='', binMin=0., binMax=0.):
    self.name_nano = name_nano
    self.name_flat = name_flat
    self.label = label
    self.title = title
    self.logic = logic
    self.units = units
    self.binMin = binMin
    self.binMax = binMax
                              

class Selection(Tools):
  def __init__(self, files, quantity, preexisting_selection=None, npoints=50, sample_type='nano', write_cut_analysis=False, proposed_cut=None, builder=''):
    self.tools                 = Tools()
    self.files                 = files
    self.quantity              = quantity
    self.preexisting_selection = preexisting_selection
    self.npoints               = npoints
    self.sample_type           = sample_type
    self.write_cut_analysis    = write_cut_analysis
    self.proposed_cut          = proposed_cut

    self.builder = 'Bs' #FIXME hardcoded
    #self.builder = builder
    if self.builder not in ['phi', 'Bs']:
      raise RuntimeError("Unknown builder. Please choose among ['phi', 'Bs']")
    print 'builder: {}'.format(self.builder)



    # efficiency will be assessed based on the mass distributions
    if self.builder == 'phi':
      self.mass = Quantity('PhiToKK_phi_mass', '', 'phi_mass', '', binMin=0, binMax=10000)
    elif self.builder == 'Bs':
      self.mass = Quantity('BsToPhiPhiTo4K_Bs_mass', '', 'Bs_mass', '', binMin=0, binMax=10000)
    
    # baseline selection
    #self.baseline_selection = ' && '.join(['{}==0'.format('hnl_charge' if str_=='sig' else 'b_hnl_charge')])
    #                                  , 'b_mass<6.35']) 


  def createOutDir(self, outputdir):
    if not path.exists(outputdir):
      os.system('mkdir {}'.format(outputdir))


  def getLabel(self, str_):
    new_str_ = str_
    if 'abs(' in str_: new_str_ = new_str_.replace('abs(', '')
    if ')'in str_: new_str_ = new_str_.replace(')', '')
    return new_str_


  def getPreselectionString(self, str_):
    '''
    function to write the already-applied cuts in a string
    str_ is to make the distinction between nano and matched samples,
    where variables are named differently
    '''
    preselection_str = []
    for item, _ in enumerate(self.preexisting_selection):
      name_variable = self.preexisting_selection[item].quantity.name_nano if self.sample_type == 'nano' else self.preexisting_selection[item].quantity.name_flat
      preselection_str.append('{}{}{}'.format(name_variable,self.preexisting_selection[item].quantity.logic,self.preexisting_selection[item].chosen_cut))
    return ' && '.join(preselection_str)


  def createHisto(self, file_, str_, with_preexisting_selection, with_extra_selection, cut=0):
    '''
    str_ makes the difference between signal and background
    '''

    filename = file_.sample_name 
    tree_name = 'Events' if self.sample_type == 'nano' else 'signal_tree'

    cut_variable = self.quantity.name_nano if self.sample_type == 'nano' else self.quantity.name_flat

    if self.sample_type == 'nano':
      if self.builder == 'phi':
        #baseline_selection = 'PhiToKK_isMatched==1' if str_=='sig' else 'fabs(PhiToKK_phi_mass - 1.02) > 0.007'
        baseline_selection = 'PhiToKK_isMatched==1' if str_=='sig' else 'fabs(PhiToKK_phi_mass - 1.02) > 0.015'
      elif self.builder == 'Bs':
        baseline_selection = 'BsToPhiPhiTo4K_isMatched==1' if str_=='sig' else 'fabs(BsToPhiPhiTo4K_Bs_mass - 5.367) > 0.2'
    else:
      baseline_selection = 'ismatched==1' if str_=='sig' else 'ismatched>-99' #FIXME
    
    c = ROOT.TCanvas()
    f = ROOT.TFile.Open(filename, 'READ')
    tree = f.Get(tree_name)
    
    # elements to fill the histogram
    distr_todraw = self.mass.name_nano if self.sample_type == 'nano' else self.mass.name_flat

    preselection = baseline_selection 
    if with_preexisting_selection and self.preexisting_selection != None: 
      preselection += ' && ' + self.getPreselectionString(str_)

    if with_extra_selection:
      preselection += ' && {}{}{}'.format(cut_variable, self.quantity.logic, cut)

    #print preselection
    
    hist = ROOT.TH1D('hist', 'hist', 1500, self.mass.binMin, self.mass.binMax)
    tree.Draw('{}>>hist'.format(distr_todraw), preselection)

    hist.SetDirectory(0)
    return hist


  def getEff(self, entries_selected, entries_initial):
    eff = entries_selected / entries_initial if entries_initial!=0 else 0.
    return eff


  def getScanGraph(self):

    '''
    plot the signal efficiency and background rejection as a function of the cut applied on the quantity of interest
    possibility to draw a line at the chosen cut
    '''

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, self.npoints) 

    canv = self.tools.createTCanvas('canv', 900, 800)
    canv.SetGrid()
    
    pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1)
    pad_up.SetBottomMargin(0.1)
    pad_up.Draw()
    canv.cd()
    pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
    pad_down.SetBottomMargin(0.15)
    pad_down.Draw()
    canv.cd()
    pad_leg = ROOT.TPad("pad_leg","pad_leg",0.5,0,1,0.25)
    pad_leg.SetBottomMargin(0.15)

    if self.write_cut_analysis:
      pad_down.cd()

    gs_sig = []
    gs_bkg = []

    # signal efficiency
    for ifile, file_ in enumerate(self.files):
      if file_.process != 'signal': continue
      g_sig = ROOT.TGraph()

      initial_sig_entries = self.createHisto(file_, 'sig', True, False).GetEntries()

      for idx, cut in enumerate(points):
        selected_sig_entries = self.createHisto(file_, 'sig', True, True, cut).GetEntries()
        g_sig.SetPoint(idx, cut, self.getEff(selected_sig_entries, initial_sig_entries))
        g_sig.SetLineWidth(0)
        g_sig.SetMarkerColor(ROOT.kOrange+1)
        g_sig.SetMarkerStyle(20+ifile)
  
      gs_sig.append(g_sig)

      gs_sig[0].GetXaxis().SetTitle('Cut on {}'.format(self.quantity.title))
      gs_sig[0].GetXaxis().SetLabelSize(0.045)
      gs_sig[0].GetXaxis().SetTitleSize(0.045)
      gs_sig[0].GetXaxis().SetTitleOffset(1.1)
      gs_sig[0].GetYaxis().SetLabelSize(0.045)
      gs_sig[0].GetYaxis().SetRangeUser(0, 1)
     
      for ig, g_sig in enumerate(gs_sig):
        if ig == 0:
          g_sig.Draw('AP')
        else:
          g_sig.Draw('P')

    # draw background rejection
    for ifile, file_ in enumerate(self.files):
      if file_.process != 'background': continue
      g_bkg = ROOT.TGraph()
      initial_bkg_entries = self.createHisto(self.files[0], 'bkg', True, False).GetEntries()

      for idx, cut in enumerate(points):
        selected_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, cut).GetEntries()
        g_bkg.SetPoint(idx, cut, 1-self.getEff(selected_bkg_entries, initial_bkg_entries))
        g_bkg.SetLineWidth(0)
        g_bkg.SetMarkerColor(ROOT.kBlue+2)
        g_bkg.SetMarkerStyle(20)

      g_bkg.Draw('P')

    if self.proposed_cut != None:
      line = ROOT.TLine(self.proposed_cut, 0, self.proposed_cut, 1)
      line.SetLineColor(2)
      line.SetLineWidth(3)
      line.Draw('same')
          
    legend = self.tools.getRootTLegend(xmin=0.5, ymin=0.37, xmax=0.82, ymax=0.61, size=0.03)
    for ifile, file_ in enumerate(self.files):
      if file_.process == 'signal':
        # only consistent if the background sample is fed before the signal samples
        legend.AddEntry(gs_sig[ifile-1], 'sig efficiency ({}GeV, {}mm)'.format(file_.signal_mass, file_.signal_ctau))
      else:
        legend.AddEntry(g_bkg, 'bkg rejection')
    legend.Draw()

    # write cutflow information
    if self.write_cut_analysis:
      pad_down.cd()
      proposed_sig_entries = self.createHisto(self.files[0], 'sig', True, True, self.proposed_cut).GetEntries()
      proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, self.proposed_cut).GetEntries()

      box1 = self.tools.getTextBox("NDC", 0.05, 0.7, 0.4, 0.98, 'Additional proposed cut: {}{}{}'.format(self.quantity.label, self.quantity.logic, self.proposed_cut), ROOT.kRed)
      box2 = self.tools.getTextBox("brNDC", 0.02, 0.3, 0.08, 0.4, '{}GeV'.format(self.files[0].signal_mass), ROOT.kBlack)
      box3 = self.tools.getTextBox("brNDC", 0.08, 0.45, 0.4, 0.65, 'N_matched ini: {}'.format(initial_sig_entries), ROOT.kOrange+1)
      box4 = self.tools.getTextBox("brNDC", 0.08, 0.25, 0.4, 0.55, 'N_matched new: {}'.format(proposed_sig_entries), ROOT.kOrange+1)
      box5 = self.tools.getTextBox("brNDC", 0.35, 0.35, 0.5, 0.65, '==> -{}%'.format(round((1 - proposed_sig_entries / initial_sig_entries)*100, 2)), ROOT.kOrange+1)
      box6 = self.tools.getTextBox("brNDC", 0.08, 0.15, 0.4, 0.3, 'N_nano ini: {}'.format(initial_bkg_entries), ROOT.kBlue+2)
      box7 = self.tools.getTextBox("brNDC", 0.08, 0., 0.4, 0.15, 'N_nano new: {}'.format(proposed_bkg_entries), ROOT.kBlue+2)
      box8 = self.tools.getTextBox("brNDC", 0.35, 0., 0.5, 0.3, '==> -{}%'.format(round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2)), ROOT.kBlue+2)

      box1.Draw('same')
      box2.Draw('same')
      box3.Draw('same')
      box4.Draw('same')
      box5.Draw('same')
      box6.Draw('same')
      box7.Draw('same')
      box8.Draw('same')

    canv.cd()
    self.createOutDir('myPlots/preselection')
    canv.SaveAs('myPlots/preselection/scan_{}.png'.format(self.getLabel(self.quantity.label)))


  def getROCGraph(self):

    points = np.linspace(self.quantity.binMin, self.quantity.binMax, self.npoints)

    canv = self.tools.createTCanvas('canv', 'canv', 1200, 1000)
    g = ROOT.TGraph2D()
    g.SetTitle('Cut on {} ({})'.format(self.quantity.title, self.quantity.logic))

    # do it for first file only
    file_ = self.files[0]
    
    initial_sig_entries = self.createHisto(file_, 'sig', True, False).GetEntries()
    initial_bkg_entries = self.createHisto(file_, 'bkg', True, False).GetEntries()

    for idx, cut in enumerate(points):
      selected_sig_entries = self.createHisto(file_, 'sig', True, True, cut).GetEntries()
      selected_bkg_entries = self.createHisto(file_, 'bkg', True, True, cut).GetEntries()
      g.SetPoint(idx, 1-self.getEff(selected_bkg_entries, initial_bkg_entries), self.getEff(selected_sig_entries, initial_sig_entries), cut)
    
    g.GetXaxis().SetTitle('background rejection')
    g.GetXaxis().SetLabelSize(0.038)
    g.GetXaxis().SetTitleSize(0.04)
    g.GetXaxis().SetTitleOffset(1.4)
    g.GetYaxis().SetTitle('signal efficiency')
    g.GetYaxis().SetLabelSize(0.038)
    g.GetYaxis().SetTitleSize(0.04)
    g.GetYaxis().SetTitleOffset(1.4)
    
    
    #ROOT.gStyle.SetPadRightMargin(5.5) 
    #g.GetZaxis().SetTitle('Cut on {}'.format(self.quantity.title))
    g.GetZaxis().SetTitle('Cut')
    #g.GetZaxis().SetLabelSize(0.038)
    g.GetZaxis().SetTitleSize(0.14)
    #g.GetZaxis().SetTitleOffset(1.2)

    ROOT.gStyle.SetPalette(53)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.7)
    g.Draw("PCOLZ")
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    ROOT.gPad.GetView().TopView()

    #ROOT.gPad.Update()
    #ROOT.gPad.Modified()
    ROOT.gStyle.SetPadRightMargin(1.5) 
   
    # draw diagonal line
    line = ROOT.TLine()
    line.SetLineColor(1)
    line.SetLineWidth(3)
    line.SetLineStyle(9)
    line.DrawLineNDC(0.1, 0.9, 0.9, 0.1)
    

    self.createOutDir('myPlots/preselection')
    canv.SaveAs('myPlots/preselection/roc_{}.png'.format(self.getLabel(self.quantity.label)))


  def printCutflowLine(self):
    initial_bkg_entries = self.createHisto(self.files[0], 'bkg', True, False).GetEntries()
    initial_sig_entries = self.createHisto(self.files[1], 'sig', True, False).GetEntries()
  
    proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, self.proposed_cut).GetEntries()
    proposed_sig_entries = self.createHisto(self.files[1], 'sig', True, True, self.proposed_cut).GetEntries()

    #print 'sig1 {} {}'.format(int(initial_sig_entries), int(initial_bkg_entries))
    #print 'sig1 {} {}'.format(int(proposed_sig_entries), int(proposed_bkg_entries))

    if len(self.files)==2:
      cutflow_line = '{qte} {log} {cut} & -{sig_per}\% & -{bkg_per}\% \\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2), 
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2),
      )

    elif len(self.files)==3:
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()

      ##print 'sig2 {}'.format(int(initial_sig1_entries))
      ##print 'sig2 {}'.format(int(proposed_sig1_entries))

      cutflow_line = '{qte} {log} {cut} & -{sig0_per}\% & -{sig1_per}\% & -{bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig0_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2), 
          sig1_per = round((1 - proposed_sig1_entries / initial_sig1_entries)*100, 2), 
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2),
      )

    else:
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()
      initial_sig2_entries = self.createHisto(self.files[3], 'sig', True, False).GetEntries()
      proposed_sig2_entries = self.createHisto(self.files[3], 'sig', True, True, self.proposed_cut).GetEntries()

      #print 'sig2 {}'.format(int(initial_sig1_entries))
      #print 'sig2 {}'.format(int(proposed_sig1_entries))
      #print 'sig3 {}'.format(int(initial_sig2_entries))
      #print 'sig3 {}'.format(int(proposed_sig2_entries))

      cutflow_line = '{qte} {log} {cut} {unit} & -{sig0_per}\% & -{sig1_per}\% & -{sig2_per}\% & -{bkg_per}\% \\\ '.format(
      #cutflow_line = '{qte} {log} {cut} {unit} & -{sig0_per}\% & -{sig1_per}\% & -{bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          unit = self.quantity.units,
          sig0_per = round((1 - proposed_sig_entries / initial_sig_entries)*100, 2) if initial_sig_entries!=0. else '-', 
          sig1_per = round((1 - proposed_sig1_entries / initial_sig1_entries)*100, 2) if initial_sig1_entries!=0. else '-',
          sig2_per = round((1 - proposed_sig2_entries / initial_sig2_entries)*100, 2) if initial_sig2_entries!=0. else '-',
          bkg_per = round((1 - proposed_bkg_entries / initial_bkg_entries)*100, 2) if initial_bkg_entries!=0. else '-',
      )

    cutflow_line = cutflow_line.replace('#eta', '$\eta$')
    cutflow_line = cutflow_line.replace('#mu#mu#pi', '$\mu\mu\pi$')
    cutflow_line = cutflow_line.replace('#mu#pi', '$\mu\pi$')
    cutflow_line = cutflow_line.replace('#Delta', '$\Delta$')
    cutflow_line = cutflow_line.replace('#pi', '$\pi$')
    cutflow_line = cutflow_line.replace('#Phi', '$\Phi$')
    print cutflow_line


  def printEfficiencyLine(self):
    # get the entries without cuts at all. Used for the cumulative efficiency
    nocut_bkg_entries = self.createHisto(self.files[0], 'bkg', False, False).GetEntries()
    nocut_sig_entries = self.createHisto(self.files[1], 'sig', False, False).GetEntries()

    initial_bkg_entries = self.createHisto(self.files[0], 'bkg', True, False).GetEntries()
    initial_sig_entries = self.createHisto(self.files[1], 'sig', True, False).GetEntries()
  
    proposed_bkg_entries = self.createHisto(self.files[0], 'bkg', True, True, self.proposed_cut).GetEntries()
    proposed_sig_entries = self.createHisto(self.files[1], 'sig', True, True, self.proposed_cut).GetEntries()

    if len(self.files)==1:
      cutflow_line = '{qte} {log} {cut} & {sig_per}\% & {bkg_per}\% \\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig_per = round((proposed_sig_entries / initial_sig_entries)*100, 1), 
          bkg_per = round((proposed_bkg_entries / initial_bkg_entries)*100, 1),
      )

    elif len(self.files)==2:
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()

      cutflow_line = '{qte} {log} {cut} & {sig0_per}\% & {sig1_per}\% & {bkg_per}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          sig0_per = round((proposed_sig_entries / initial_sig_entries)*100, 1), 
          sig1_per = round((proposed_sig1_entries / initial_sig1_entries)*100, 1), 
          bkg_per = round((proposed_bkg_entries / initial_bkg_entries)*100, 1),
      )

    else:
      nocut_sig1_entries = self.createHisto(self.files[2], 'sig', False, False).GetEntries()
      initial_sig1_entries = self.createHisto(self.files[2], 'sig', True, False).GetEntries()
      proposed_sig1_entries = self.createHisto(self.files[2], 'sig', True, True, self.proposed_cut).GetEntries()
      nocut_sig2_entries = self.createHisto(self.files[3], 'sig', False, False).GetEntries()
      initial_sig2_entries = self.createHisto(self.files[3], 'sig', True, False).GetEntries()
      proposed_sig2_entries = self.createHisto(self.files[3], 'sig', True, True, self.proposed_cut).GetEntries()

      cutflow_line = '{qte} {log} {cut} {unit} & {sig0_per}\% & {sig0_cumul}\% & {sig1_per}\% & {sig1_cumul}\% & {sig2_per}\% & {sig2_cumul}\% & {bkg_per}\% & {bkg_cumul}\% \\\ '.format(
          qte = self.quantity.title, 
          log = self.quantity.logic, 
          cut = self.proposed_cut, 
          unit = self.quantity.units,
          sig0_per = round((proposed_sig_entries / initial_sig_entries)*100, 1) if initial_sig_entries!=0. else '-', 
          sig1_per = round((proposed_sig1_entries / initial_sig1_entries)*100, 1) if initial_sig1_entries!=0. else '-',
          sig2_per = round((proposed_sig2_entries / initial_sig2_entries)*100, 1) if initial_sig2_entries!=0. else '-',
          bkg_per = round((proposed_bkg_entries / initial_bkg_entries)*100, 1) if initial_bkg_entries!=0. else '-',
          sig0_cumul = round((proposed_sig_entries / nocut_sig_entries)*100, 1),
          sig1_cumul = round((proposed_sig1_entries / nocut_sig1_entries)*100, 1),
          sig2_cumul = round((proposed_sig2_entries / nocut_sig2_entries)*100, 1),
          bkg_cumul = round((proposed_bkg_entries / nocut_bkg_entries)*100, 1),
      )

    cutflow_line = cutflow_line.replace('#eta', '$\eta$')
    cutflow_line = cutflow_line.replace('#mu#mu#pi', '$\mu\mu\pi$')
    cutflow_line = cutflow_line.replace('#mu#pi', '$\mu\pi$')
    cutflow_line = cutflow_line.replace('#Delta', '$\Delta$')
    cutflow_line = cutflow_line.replace('#pi', '$\pi$')
    cutflow_line = cutflow_line.replace('#Phi', '$\Phi$')
    print cutflow_line



class FileCollection(object): 
  '''
  this class allows to associate to the studied files the corresponding 
  mass and coupling points. Useful when scaning through e.g masses
  '''
  def __init__(self, sample_name='', process='', signal_mass='', signal_ctau=''):
    self.sample_name = sample_name
    self.process = process
    self.signal_mass = signal_mass
    self.signal_ctau = signal_ctau
    # add label for reweighted?
    if self.process not in ['signal', 'background']:
      raise RuntimeError("Unknown process. Please choose among ['signal', 'background']")


class PreselectedQuantity(object):
  '''
  this class will allow to impose a preselection when studying 
  the impact of the cut on a given parameter
  '''
  def __init__(self, quantity, chosen_cut):
    self.quantity   = quantity
    self.chosen_cut = chosen_cut


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  # builder
  do_phi = False
  do_Bs = True

  # what to study
  printCutflow = True
  printEfficiency = False
  printScan = False

  doPreselection = True
  doBaselineSelection = False
  doLooseSelection = False

  if do_Bs:

    file_background = FileCollection(
        sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/data/study_preselection_v0/ParkingBPH1_Run2018A/Chunk0_n419/merged/bparknano_phi_selected_morevars.root',
        process = 'background',
        )

    file_signal = FileCollection(
        sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_phi_selected_morevars.root',
        process = 'signal',
        signal_mass = '',
        signal_ctau = '',
        )
    
    files = []
    files.append(file_background)
    files.append(file_signal)

    phi1_mass = Quantity('abs(BsToPhiPhiTo4K_phi1_mass - 1.0195)', '', 'phi1_mass', '|mass(phi1) - 1.02|', '<', 'GeV', 0, 0.017) 
    phi2_mass = Quantity('abs(BsToPhiPhiTo4K_phi2_mass - 1.0195)', '', 'phi2_mass', '|mass(phi2) - 1.02|', '<', 'GeV', 0, 0.017) 
    phi1_pt = Quantity('BsToPhiPhiTo4K_phi1_pt', '', 'phi1_pt', 'pt(phi1)', '>', 'GeV', 0, 10) 
    phi1_eta = Quantity('fabs(BsToPhiPhiTo4K_phi1_eta)', '', 'phi1_eta', '|eta(phi1)|', '<', 'GeV', 0, 3) 
    phi2_pt = Quantity('BsToPhiPhiTo4K_phi2_pt', '', 'phi2_pt', 'pt(phi2)', '>', 'GeV', 0, 10) 
    phi1_pt_x_phi2_pt = Quantity('BsToPhiPhiTo4K_phi1_pt*BsToPhiPhiTo4K_phi2_pt', '', 'phi1_pt_x_phi2_pt', 'pt(phi1) * pt(phi2)', '>', 'GeV', 0, 10) 
    k1_pt = Quantity('BsToPhiPhiTo4K_k1_pt', '', 'k1_pt', 'pt(k1)', '>', 'GeV', 0, 10) 
    k2_pt = Quantity('BsToPhiPhiTo4K_k2_pt', '', 'k2_pt', 'pt(k2)', '>', 'GeV', 0, 10) 
    k3_pt = Quantity('BsToPhiPhiTo4K_k3_pt', '', 'k3_pt', 'pt(k3)', '>', 'GeV', 0, 10) 
    k4_pt = Quantity('BsToPhiPhiTo4K_k4_pt', '', 'k4_pt', 'pt(k4)', '>', 'GeV', 0, 10) 
    Bs_pt = Quantity('BsToPhiPhiTo4K_Bs_pt', '', 'Bs_pt', 'pt(Bs)', '>', 'GeV', 0, 10) 
    Bs_eta = Quantity('fabs(BsToPhiPhiTo4K_Bs_eta)', '', 'Bs_eta', '|eta(Bs)|', '<', 'GeV', 0, 3) 
    deltar_min = Quantity('BsToPhiPhiTo4K_deltaR_min', '', 'deltar_min', 'min(deltaR)', '<', 'GeV', 0, 0.5) 
    deltar_max = Quantity('BsToPhiPhiTo4K_deltaR_max', '', 'deltar_max', 'max(deltaR)', '<', 'GeV', 0, 3) 
    Bs_lxy_sig = Quantity('BsToPhiPhiTo4K_Bs_lxy_sig', '', 'Bs_lxy_sig', 'lxy_sig(Bs)', '>', 'GeV', 0, 10) 
    Bs_sv_prob = Quantity('BsToPhiPhiTo4K_Bs_sv_prob', '', 'Bs_sv_prob', 'sv_prob(Bs)', '>', 'GeV', 0, 10) 
    Bs_cos2D = Quantity('BsToPhiPhiTo4K_Bs_cos2D', '', 'Bs_cos2D', 'cos2D(Bs)', '>', 'GeV', 0, 1) 

    #k1_pt_prefit = Quantity('ProbeTracks_pt[PhiToKK_k1_idx]', '', 'k1_pt_prefit', 'p_{T}(k_{1}) (GeV) (prefit)', '>', 'GeV', 0.4, 3)
    #k2_pt_prefit = Quantity('ProbeTracks_pt[PhiToKK_k2_idx]', '', 'k2_pt_prefit', 'p_{T}(k_{1}) (GeV) (prefit)', '>', 'GeV', 0.4, 3)
    #k2_pt = Quantity('PhiToKK_phi_k2_pt', '', 'k2_pt', 'p_{T}(k_{2}) (GeV)', '>', 'GeV', 0.4, 3)
    #k1_pt = Quantity('PhiToKK_phi_k1_pt', '', 'k1_pt', 'p_{T}(k_{1}) (GeV)', '>', 'GeV', 0.4, 3)
    #k2_pt = Quantity('PhiToKK_phi_k2_pt', '', 'k2_pt', 'p_{T}(k_{2}) (GeV)', '>', 'GeV', 0.4, 3)
    #k1_eta_prefit = Quantity('abs(ProbeTracks_eta[PhiToKK_k1_idx])', '', 'k1_eta_prefit', '|#eta|(k_{1}) (prefit)', '<', '', 0, 2.5)
    #k2_eta_prefit = Quantity('abs(ProbeTracks_eta[PhiToKK_k2_idx])', '', 'k2_eta_prefit', '|#eta|(k_{1}) (prefit)', '<', '', 0, 2.5)
    #k1_eta = Quantity('abs(PhiToKK_phi_k1_eta)', '', 'k1_eta', '|#eta|(k_{1})', '<', '', 0, 2.5)
    #k2_eta = Quantity('abs(PhiToKK_phi_k2_eta)', '', 'k2_eta', '|#eta|(k_{2})', '<', '', 0, 2.5)

    #sv_prob = Quantity('PhiToKK_phi_sv_prob', '', 'sv_prob', 'SV probability', '>', '', 0, 0.05)
    #phi_cos2D = Quantity('PhiToKK_phi_cos2D', '', 'phi_cos2D', 'SV cos2D', '>', '', 0.93, 1) 
    #phi_mass = Quantity('abs(PhiToKK_phi_mass - 1.0195)', '', 'phi_mass', '|mass - 1.02|', '<', 'GeV', 0, 0.1) 
    #deltar = Quantity('PhiToKK_deltaR_postfit', '' , 'deltaR_postfit', '#Delta R (postfit)', '<', '', 0, 1)
    #lxysig = Quantity('PhiToKK_phi_lxy_sig', '', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 100)


    if doPreselection:

      preselection = [] 

      cut_phi1_mass = 0.012
      if printScan: Selection(files, phi1_mass, npoints=30, write_cut_analysis=False, proposed_cut=cut_phi1_mass).getScanGraph()
      if printCutflow: Selection(files, phi1_mass, proposed_cut=cut_phi1_mass).printCutflowLine()
      if printEfficiency: Selection(files, phi1_mass, proposed_cut=cut_phi1_mass).printEfficiencyLine()
      preselection.append(PreselectedQuantity(phi1_mass, cut_phi1_mass))

      cut_phi2_mass = 0.012
      if printScan: Selection(files, phi2_mass, npoints=30, write_cut_analysis=False, proposed_cut=cut_phi2_mass).getScanGraph()
      if printCutflow: Selection(files, phi2_mass, preexisting_selection=preselection, proposed_cut=cut_phi2_mass).printCutflowLine()
      if printEfficiency: Selection(files, phi2_mass, preexisting_selection=preselection, proposed_cut=cut_phi2_mass).printEfficiencyLine()
      preselection.append(PreselectedQuantity(phi2_mass, cut_phi2_mass))

      cut_phi1_pt = 2.5
      if printScan: Selection(files, phi1_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_phi1_pt).getScanGraph()
      if printCutflow: Selection(files, phi1_pt, preexisting_selection=preselection, proposed_cut=cut_phi1_pt).printCutflowLine()
      if printEfficiency: Selection(files, phi1_pt, preexisting_selection=preselection, proposed_cut=cut_phi1_pt).printEfficiencyLine()
      preselection.append(PreselectedQuantity(phi1_pt, cut_phi1_pt))

      #cut_phi1_eta = 2
      #if printScan: Selection(files, phi1_eta, npoints=30, write_cut_analysis=False, proposed_cut=cut_phi1_eta).getScanGraph()
      #if printCutflow: Selection(files, phi1_eta, preexisting_selection=preselection, proposed_cut=cut_phi1_eta).printCutflowLine()
      #if printEfficiency: Selection(files, phi1_eta, preexisting_selection=preselection, proposed_cut=cut_phi1_eta).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(phi1_eta, cut_phi1_eta))

      cut_phi2_pt = 1.8
      if printScan: Selection(files, phi2_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_phi2_pt).getScanGraph()
      if printCutflow: Selection(files, phi2_pt, preexisting_selection=preselection, proposed_cut=cut_phi2_pt).printCutflowLine()
      if printEfficiency: Selection(files, phi2_pt, preexisting_selection=preselection, proposed_cut=cut_phi2_pt).printEfficiencyLine()
      preselection.append(PreselectedQuantity(phi2_pt, cut_phi2_pt))

      cut_phi1_pt_x_phi2_pt = 6
      if printScan: Selection(files, phi1_pt_x_phi2_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_phi1_pt_x_phi2_pt).getScanGraph()
      if printCutflow: Selection(files, phi1_pt_x_phi2_pt, preexisting_selection=preselection, proposed_cut=cut_phi1_pt_x_phi2_pt).printCutflowLine()
      if printEfficiency: Selection(files, phi1_pt_x_phi2_pt, preexisting_selection=preselection, proposed_cut=cut_phi1_pt_x_phi2_pt).printEfficiencyLine()
      preselection.append(PreselectedQuantity(phi1_pt_x_phi2_pt, cut_phi1_pt_x_phi2_pt))

      #cut_k1_pt = 1.6
      #if printScan: Selection(files, k1_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_k1_pt).getScanGraph()
      ##if printCutflow: Selection(files, k1_pt, preexisting_selection=preselection, proposed_cut=cut_k1_pt).printCutflowLine()
      #if printEfficiency: Selection(files, k1_pt, preexisting_selection=preselection, proposed_cut=cut_k1_pt).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k1_pt, cut_k1_pt))

      #cut_k2_pt = 1.1
      #if printScan: Selection(files, k2_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_k2_pt).getScanGraph()
      ##if printCutflow: Selection(files, k2_pt, preexisting_selection=preselection, proposed_cut=cut_k2_pt).printCutflowLine()
      #if printEfficiency: Selection(files, k2_pt, preexisting_selection=preselection, proposed_cut=cut_k2_pt).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k1_pt, cut_k1_pt))

      #cut_k3_pt = 1.
      #if printScan: Selection(files, k3_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_k3_pt).getScanGraph()
      ##if printCutflow: Selection(files, k3_pt, preexisting_selection=preselection, proposed_cut=cut_k3_pt).printCutflowLine()
      #if printEfficiency: Selection(files, k3_pt, preexisting_selection=preselection, proposed_cut=cut_k3_pt).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k3_pt, cut_k3_pt))

      #cut_k4_pt = 0.75
      #if printScan: Selection(files, k4_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_k4_pt).getScanGraph()
      #if printCutflow: Selection(files, k4_pt, preexisting_selection=preselection, proposed_cut=cut_k4_pt).printCutflowLine()
      #if printEfficiency: Selection(files, k4_pt, preexisting_selection=preselection, proposed_cut=cut_k4_pt).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k4_pt, cut_k4_pt))

      cut_Bs_pt = 1.5
      if printScan: Selection(files, Bs_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_Bs_pt).getScanGraph()
      if printCutflow: Selection(files, Bs_pt, preexisting_selection=preselection, proposed_cut=cut_Bs_pt).printCutflowLine()
      if printEfficiency: Selection(files, Bs_pt, preexisting_selection=preselection, proposed_cut=cut_Bs_pt).printEfficiencyLine()
      preselection.append(PreselectedQuantity(Bs_pt, cut_Bs_pt))

      cut_Bs_eta = 2.5
      if printScan: Selection(files, Bs_eta, npoints=30, write_cut_analysis=False, proposed_cut=cut_Bs_eta).getScanGraph()
      if printCutflow: Selection(files, Bs_eta, preexisting_selection=preselection, proposed_cut=cut_Bs_eta).printCutflowLine()
      if printEfficiency: Selection(files, Bs_eta, preexisting_selection=preselection, proposed_cut=cut_Bs_eta).printEfficiencyLine()
      preselection.append(PreselectedQuantity(Bs_eta, cut_Bs_eta))

      cut_deltar_min = 0.15
      if printScan: Selection(files, deltar_min, npoints=30, write_cut_analysis=False, proposed_cut=cut_deltar_min).getScanGraph()
      if printCutflow: Selection(files, deltar_min, preexisting_selection=preselection, proposed_cut=cut_deltar_min).printCutflowLine()
      if printEfficiency: Selection(files, deltar_min, preexisting_selection=preselection, proposed_cut=cut_deltar_min).printEfficiencyLine()
      preselection.append(PreselectedQuantity(deltar_min, cut_deltar_min))

      cut_deltar_max = 2.5
      if printScan: Selection(files, deltar_max, npoints=30, write_cut_analysis=False, proposed_cut=cut_deltar_max).getScanGraph()
      if printCutflow: Selection(files, deltar_max, preexisting_selection=preselection, proposed_cut=cut_deltar_max).printCutflowLine()
      if printEfficiency: Selection(files, deltar_max, preexisting_selection=preselection, proposed_cut=cut_deltar_max).printEfficiencyLine()
      preselection.append(PreselectedQuantity(deltar_max, cut_deltar_max))

      cut_Bs_lxy_sig = 1.
      if printScan: Selection(files, Bs_lxy_sig, npoints=30, write_cut_analysis=False, proposed_cut=cut_Bs_lxy_sig).getScanGraph()
      if printCutflow: Selection(files, Bs_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_Bs_lxy_sig).printCutflowLine()
      if printEfficiency: Selection(files, Bs_lxy_sig, preexisting_selection=preselection, proposed_cut=cut_Bs_lxy_sig).printEfficiencyLine()
      preselection.append(PreselectedQuantity(Bs_lxy_sig, cut_Bs_lxy_sig))

      cut_Bs_sv_prob = 0.001
      if printScan: Selection(files, Bs_sv_prob, npoints=30, write_cut_analysis=False, proposed_cut=cut_Bs_sv_prob).getScanGraph()
      if printCutflow: Selection(files, Bs_sv_prob, preexisting_selection=preselection, proposed_cut=cut_Bs_sv_prob).printCutflowLine()
      if printEfficiency: Selection(files, Bs_sv_prob, preexisting_selection=preselection, proposed_cut=cut_Bs_sv_prob).printEfficiencyLine()
      preselection.append(PreselectedQuantity(Bs_sv_prob, cut_Bs_sv_prob))

      cut_Bs_cos2D = 0.9
      if printScan: Selection(files, Bs_cos2D, npoints=30, write_cut_analysis=False, proposed_cut=cut_Bs_cos2D).getScanGraph()
      if printCutflow: Selection(files, Bs_cos2D, preexisting_selection=preselection, proposed_cut=cut_Bs_cos2D).printCutflowLine()
      if printEfficiency: Selection(files, Bs_cos2D, preexisting_selection=preselection, proposed_cut=cut_Bs_cos2D).printEfficiencyLine()
      preselection.append(PreselectedQuantity(Bs_cos2D, cut_Bs_cos2D))

      print Selection(files, phi1_mass, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=0).getPreselectionString('bkg')




  if do_phi:

    file_background = FileCollection(
        #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/data/test_v0/ParkingBPH1_Run2018A/Chunk0_n500/merged/bparknano_morevars.root',
        sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/data/study_preselection_v0/ParkingBPH1_Run2018A/Chunk0_n3/merged/bparknano_loose.root',
        process = 'background',
        )

    file_signal = FileCollection(
        #sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_v0/BsToPhiPhiTo4K/nanoFiles/merged/bparknano.root',
        sample_name = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_loose.root',
        process = 'signal',
        signal_mass = '',
        signal_ctau = '',
        )
    
    files = []
    files.append(file_background)
    files.append(file_signal)

    k1_pt_prefit = Quantity('ProbeTracks_pt[PhiToKK_k1_idx]', '', 'k1_pt_prefit', 'p_{T}(k_{1}) (GeV) (prefit)', '>', 'GeV', 0.4, 3)
    k2_pt_prefit = Quantity('ProbeTracks_pt[PhiToKK_k2_idx]', '', 'k2_pt_prefit', 'p_{T}(k_{1}) (GeV) (prefit)', '>', 'GeV', 0.4, 3)
    k2_pt = Quantity('PhiToKK_phi_k2_pt', '', 'k2_pt', 'p_{T}(k_{2}) (GeV)', '>', 'GeV', 0.4, 3)
    k1_pt = Quantity('PhiToKK_phi_k1_pt', '', 'k1_pt', 'p_{T}(k_{1}) (GeV)', '>', 'GeV', 0.4, 3)
    k2_pt = Quantity('PhiToKK_phi_k2_pt', '', 'k2_pt', 'p_{T}(k_{2}) (GeV)', '>', 'GeV', 0.4, 3)
    k1_eta_prefit = Quantity('abs(ProbeTracks_eta[PhiToKK_k1_idx])', '', 'k1_eta_prefit', '|#eta|(k_{1}) (prefit)', '<', '', 0, 2.5)
    k2_eta_prefit = Quantity('abs(ProbeTracks_eta[PhiToKK_k2_idx])', '', 'k2_eta_prefit', '|#eta|(k_{1}) (prefit)', '<', '', 0, 2.5)
    k1_eta = Quantity('abs(PhiToKK_phi_k1_eta)', '', 'k1_eta', '|#eta|(k_{1})', '<', '', 0, 2.5)
    k2_eta = Quantity('abs(PhiToKK_phi_k2_eta)', '', 'k2_eta', '|#eta|(k_{2})', '<', '', 0, 2.5)

    sv_prob = Quantity('PhiToKK_phi_sv_prob', '', 'sv_prob', 'SV probability', '>', '', 0, 0.05)
    phi_cos2D = Quantity('PhiToKK_phi_cos2D', '', 'phi_cos2D', 'SV cos2D', '>', '', 0.93, 1) 
    phi_mass = Quantity('abs(PhiToKK_phi_mass - 1.0195)', '', 'phi_mass', '|mass - 1.02|', '<', 'GeV', 0, 0.1) 
    deltar = Quantity('PhiToKK_deltaR_postfit', '' , 'deltaR_postfit', '#Delta R (postfit)', '<', '', 0, 1)
    lxysig = Quantity('PhiToKK_phi_lxy_sig', '', 'sv_lxysig', 'significance of the SV displacement', '>', '', 0, 100)


    if doPreselection:

      preselection = [] 

      #cut_k1_pt_prefit = 0.75
      #if printScan: Selection(files, k1_pt_prefit, npoints=30, write_cut_analysis=False, proposed_cut=cut_k1_pt_prefit).getScanGraph()
      #if printCutflow: Selection(files, k1_pt_prefit, proposed_cut=cut_k1_pt_prefit).printCutflowLine()
      #if printEfficiency: Selection(files, k1_pt_prefit, proposed_cut=cut_k1_pt_prefit).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k1_pt_prefit, cut_k1_pt_prefit))

      #cut_k2_pt_prefit = 0.75
      #if printScan: Selection(files, k2_pt_prefit, npoints=30, write_cut_analysis=False, proposed_cut=cut_k2_pt_prefit).getScanGraph()
      #if printCutflow: Selection(files, k2_pt_prefit, proposed_cut=cut_k2_pt_prefit).printCutflowLine()
      #if printEfficiency: Selection(files, k2_pt_prefit, proposed_cut=cut_k2_pt_prefit).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k2_pt_prefit, cut_k2_pt_prefit))

      #cut_k1_eta_prefit = 2.5
      #if printScan: Selection(files, k1_eta_prefit, npoints=30, write_cut_analysis=False, proposed_cut=cut_k1_eta_prefit).getScanGraph()
      #if printCutflow: Selection(files, k1_eta_prefit, proposed_cut=cut_k1_eta_prefit).printCutflowLine()
      #if printEfficiency: Selection(files, k1_eta_prefit, proposed_cut=cut_k1_eta_prefit).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k1_eta_prefit, cut_k1_eta_prefit))

      #cut_k2_eta_prefit = 2.5
      #if printScan: Selection(files, k2_eta_prefit, npoints=30, write_cut_analysis=False, proposed_cut=cut_k2_eta_prefit).getScanGraph()
      #if printCutflow: Selection(files, k2_eta_prefit, proposed_cut=cut_k2_eta_prefit).printCutflowLine()
      #if printEfficiency: Selection(files, k2_eta_prefit, proposed_cut=cut_k2_eta_prefit).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k2_eta_prefit, cut_k2_eta_prefit))

      cut_phi_mass = 0.015
      if printScan: Selection(files, phi_mass, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_phi_mass).getScanGraph()
      if printCutflow: Selection(files, phi_mass, proposed_cut=cut_phi_mass).printCutflowLine()
      if printEfficiency: Selection(files, phi_mass, proposed_cut=cut_phi_mass).printEfficiencyLine()
      preselection.append(PreselectedQuantity(phi_mass, cut_phi_mass))

      cut_k1_pt = 0.8
      if printScan: Selection(files, k1_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_k1_pt).getScanGraph()
      if printCutflow: Selection(files, k1_pt, preexisting_selection=preselection, proposed_cut=cut_k1_pt).printCutflowLine()
      if printEfficiency: Selection(files, k1_pt, preexisting_selection=preselection, proposed_cut=cut_k1_pt).printEfficiencyLine()
      preselection.append(PreselectedQuantity(k1_pt, cut_k1_pt))

      #cut_k1_eta = 2.5
      #if printScan: Selection(files, k1_eta, npoints=30, write_cut_analysis=False, proposed_cut=cut_k1_eta).getScanGraph()
      #if printCutflow: Selection(files, k1_eta, preexisting_selection=preselection, proposed_cut=cut_k1_eta).printCutflowLine()
      #if printEfficiency: Selection(files, k1_eta, preexisting_selection=preselection, proposed_cut=cut_k1_eta).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k1_eta, cut_k1_eta))

      cut_k2_pt = 0.7
      if printScan: Selection(files, k2_pt, npoints=30, write_cut_analysis=False, proposed_cut=cut_k2_pt).getScanGraph()
      if printCutflow: Selection(files, k2_pt, preexisting_selection=preselection, proposed_cut=cut_k2_pt).printCutflowLine()
      if printEfficiency: Selection(files, k2_pt, preexisting_selection=preselection, proposed_cut=cut_k2_pt).printEfficiencyLine()
      preselection.append(PreselectedQuantity(k2_pt, cut_k2_pt))

      #cut_k2_eta = 2.5
      #if printScan: Selection(files, k2_eta, npoints=30, write_cut_analysis=False, proposed_cut=cut_k2_eta).getScanGraph()
      #if printCutflow: Selection(files, k2_eta, preexisting_selection=preselection,  proposed_cut=cut_k2_eta).printCutflowLine()
      #if printEfficiency: Selection(files, k2_eta, preexisting_selection=preselection, proposed_cut=cut_k2_eta).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(k2_eta, cut_k2_eta))

      cut_deltar = 0.25
      if printScan: Selection(files, deltar, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_deltar).getScanGraph()
      if printCutflow: Selection(files, deltar, preexisting_selection=preselection, proposed_cut=cut_deltar).printCutflowLine()
      if printEfficiency: Selection(files, deltar, preexisting_selection=preselection, proposed_cut=cut_deltar).printEfficiencyLine()
      preselection.append(PreselectedQuantity(deltar, cut_deltar))

      cut_sv_prob = 0.01
      if printScan: Selection(files, sv_prob, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_sv_prob).getScanGraph()
      if printCutflow: Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printCutflowLine()
      if printEfficiency: Selection(files, sv_prob, preexisting_selection=preselection, proposed_cut=cut_sv_prob).printEfficiencyLine()
      preselection.append(PreselectedQuantity(sv_prob, cut_sv_prob))

      #cut_lxysig = 0.5
      #if printScan: Selection(files, lxysig, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_lxysig).getScanGraph()
      #if printCutflow: Selection(files, lxysig, preexisting_selection=preselection, proposed_cut=cut_lxysig).printCutflowLine()
      #if printEfficiency: Selection(files, lxysig, preexisting_selection=preselection, proposed_cut=cut_lxysig).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(lxysig, cut_lxysig))

      #cut_phi_cos2D = 0.9
      #if printScan: Selection(files, phi_cos2D, npoints=30, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=cut_phi_cos2D).getScanGraph()
      #if printCutflow: Selection(files, phi_cos2D, preexisting_selection=preselection, proposed_cut=cut_phi_cos2D).printCutflowLine()
      #if printEfficiency: Selection(files, phi_cos2D, preexisting_selection=preselection, proposed_cut=cut_phi_cos2D).printEfficiencyLine()
      #preselection.append(PreselectedQuantity(phi_cos2D, cut_phi_cos2D))


      print Selection(files, phi_mass, npoints=2, write_cut_analysis=False, preexisting_selection=preselection, proposed_cut=0).getPreselectionString('bkg')

  #if doBaselineSelection:

  #  baseline_selection = []

  #  cut_hnl_charge = 0
  #  if printCutflow: Selection(files, hnl_charge, proposed_cut=cut_hnl_charge).printCutflowLine()
  #  if printEfficiency: Selection(files, hnl_charge, proposed_cut=cut_hnl_charge).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(hnl_charge, cut_hnl_charge))

  #  cut_mu0_softid = 1
  #  if printCutflow: Selection(files, mu0_softid, preexisting_selection=baseline_selection, proposed_cut=cut_mu0_softid).printCutflowLine()
  #  if printEfficiency: Selection(files, mu0_softid, preexisting_selection=baseline_selection, proposed_cut=cut_mu0_softid).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(mu0_softid, cut_mu0_softid))

  #  cut_mu_looseid = 1
  #  if printCutflow: Selection(files, mu_looseid, preexisting_selection=baseline_selection, proposed_cut=cut_mu_looseid).printCutflowLine()
  #  if printEfficiency: Selection(files, mu_looseid, preexisting_selection=baseline_selection, proposed_cut=cut_mu_looseid).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(mu_looseid, cut_mu_looseid))

  #  cut_pi_packedcandhashighpurity = 1
  #  if printCutflow: Selection(files, pi_packedcandhashighpurity, preexisting_selection=baseline_selection, proposed_cut=cut_pi_packedcandhashighpurity).printCutflowLine()
  #  if printEfficiency: Selection(files, pi_packedcandhashighpurity, preexisting_selection=baseline_selection, proposed_cut=cut_pi_packedcandhashighpurity).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(pi_packedcandhashighpurity, cut_pi_packedcandhashighpurity))

  #  cut_mu0mu_veto = 1
  #  if printCutflow: Selection(files, mu0mu_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0mu_veto).printCutflowLine()
  #  if printEfficiency: Selection(files, mu0mu_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0mu_veto).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(mu0mu_veto, cut_mu0mu_veto))

  #  cut_mu0pi_veto = 1
  #  if printCutflow: Selection(files, mu0pi_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0pi_veto).printCutflowLine()
  #  if printEfficiency: Selection(files, mu0pi_veto, preexisting_selection=baseline_selection, proposed_cut=cut_mu0pi_veto).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(mu0pi_veto, cut_mu0pi_veto))

  #  cut_sv_lxyz = 100
  #  if printCutflow: Selection(files, sv_lxyz, preexisting_selection=baseline_selection, proposed_cut=cut_sv_lxyz).printCutflowLine()
  #  if printEfficiency: Selection(files, sv_lxyz, preexisting_selection=baseline_selection, proposed_cut=cut_sv_lxyz).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(sv_lxyz, cut_sv_lxyz))

  #  cut_sv_err = 0
  #  if printCutflow: Selection(files, sv_err, preexisting_selection=baseline_selection, proposed_cut=cut_sv_err).printCutflowLine()
  #  if printEfficiency: Selection(files, sv_err, preexisting_selection=baseline_selection, proposed_cut=cut_sv_err).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(sv_err, cut_sv_err))

  #  cut_dr_mu0_pi = 0.01
  #  if printCutflow: Selection(files, dr_mu0_pi, preexisting_selection=baseline_selection, proposed_cut=cut_dr_mu0_pi).printCutflowLine()
  #  if printEfficiency: Selection(files, dr_mu0_pi, preexisting_selection=baseline_selection, proposed_cut=cut_dr_mu0_pi).printEfficiencyLine()
  #  baseline_selection.append(PreselectedQuantity(dr_mu0_pi, cut_dr_mu0_pi))


  #if doLooseSelection:

  #  loose_selection = []

  #  cut_trg_mu_pt = 0.3
  #  if printCutflow: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printCutflowLine()
  #  if printEfficiency: Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
  #  Selection(files, trg_mu_pt, proposed_cut=cut_trg_mu_pt).printEfficiencyLine()
  #  loose_selection.append(PreselectedQuantity(trg_mu_pt, cut_trg_mu_pt))

  #  cut_trg_mu_eta = 2.8
  #  if printScan: Selection(files, trg_mu_eta, npoints=30, write_cut_analysis=False, preexisting_selection=loose_selection, proposed_cut=cut_trg_mu_eta).getScanGraph()
  #  if printCutflow: Selection(files, trg_mu_eta, preexisting_selection=loose_selection, proposed_cut=cut_trg_mu_eta).printCutflowLine()
  #  if printEfficiency: Selection(files, trg_mu_eta, preexisting_selection=loose_selection, proposed_cut=cut_trg_mu_eta).printEfficiencyLine()
  #  loose_selection.append(PreselectedQuantity(trg_mu_eta, cut_trg_mu_eta))



