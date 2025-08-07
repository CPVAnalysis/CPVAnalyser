import os
import os.path
import sys
from os import path
import ROOT
from ROOT import gROOT, gStyle
import math
import numpy as np

from tools import Tools


class EfficiencyAnalyser(Tools):
  def __init__(self, filename='', cuts=None, eta_bins=None, pt_bins=None, title='', outdirlabel='default'):
    self.tools = Tools()
    self.filename = filename
    self.eta_bins = eta_bins
    self.pt_bins = pt_bins
    self.cuts = cuts
    self.title = title
    self.outdirlabel = outdirlabel
    self.outputdir = self.tools.getOutDir('./myPlots/efficiency', self.outdirlabel)

    self.treename = 'Events'

    if path.exists('{}/numbers.txt'.format(self.outputdir)):
      os.system('rm {}/numbers.txt'.format(self.outputdir))


  def getEfficiencyFromPOGFile(self, eta, pt):
     # for efficiency from tracking POG
     if abs(eta) < 0.9:
       file_eff = open('./data/track_reconstruction_efficiency_barrel.txt', 'read')
     elif abs(eta) > 0.9 and abs(eta) < 1.4:
       file_eff = open('./data/track_reconstruction_efficiency_transition.txt', 'read')
     elif abs(eta) > 1.4:
       file_eff = open('./data/track_reconstruction_efficiency_endcap.txt', 'read')

     lines = file_eff.readlines()
     list_pt = []
     list_eff = []

     for line in lines:
      idx = line.find(' ')
      line_pt = line[0:idx-1]
      line_eff = line[idx+1:len(line)-1]

      list_pt.append(float(line_pt))
      list_eff.append(float(line_eff))

     for ipt, line_pt in enumerate(list_pt):
       if pt > list_pt[ipt] and pt < list_pt[ipt+1]:
         eff_avg = (list_eff[ipt] + list_eff[ipt+1]) / 2.
         #print '{} < {} < {}: {}'.format(list_pt[ipt], entry.GenPart_pt[part_idx], list_pt[ipt+1], eff_avg)

     return eff_avg


  def getEfficiency(self, eta_bins=None, pt_bins=None):
    f = self.tools.getRootFile(self.filename, with_ext=False) 
    tree = self.tools.getTree(f, self.treename)

    if eta_bins != None and pt_bins != None:
      n_num = np.zeros((len(eta_bins), len(pt_bins)))
      n_deno = np.zeros((len(eta_bins), len(pt_bins)))
      efficiency = np.zeros((len(eta_bins), len(pt_bins))) # efficiency retrieved from nanoAOD builder
      efficiency_trkPOG = np.zeros((len(eta_bins), len(pt_bins))) # efficiency retrieved from tracking POG
      error = np.zeros((len(eta_bins), len(pt_bins))) 
    else:
      n_num = np.zeros(len(self.cuts))
      n_deno = np.zeros(len(self.cuts))
      efficiency = np.zeros(len(self.cuts))
      error = np.zeros(len(self.cuts))

    eff_POG = 0
    count_gen = 0

    for ientry, entry in enumerate(tree):
      if ientry%1000 == 0:
        print '{}% events processed'.format(round(float(ientry)/float(tree.GetEntries())*100, 1))
      
      cand_ismatched = 0
      for icand in range(0, entry.nBsToPhiPhiTo4K):
        if entry.BsToPhiPhiTo4K_isMatched[icand] == 1: cand_ismatched = 1

      matching_cond = {}
      matching_cond['candidate'] = cand_ismatched==1 

      has_trgmu = False

      # searching for trigger muon
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) == 13 and entry.GenPart_pt[icand] > 7.5 and abs(entry.GenPart_eta[icand]) < 1.5:
        #if abs(entry.GenPart_pdgId[icand]) == 13 and entry.GenPart_pt[icand] > 10.5 and abs(entry.GenPart_eta[icand]) < 1.5:
          has_trgmu = True

      # skip events without trigger muon
      #if not has_trgmu: continue

      Bs_idx = -1
      phi1_idx = -1
      phi2_idx = -1
      k1_idx = -1
      k2_idx = -1
      k3_idx = -1
      k4_idx = -1

      phi_idx = []
      count_phi = 0
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) == 333:
          mother_idx = entry.GenPart_genPartIdxMother[icand]
          mother_pdgid = entry.GenPart_pdgId[mother_idx]
          if abs(mother_pdgid) == 531:
            count_phi += 1
            phi_idx.append(icand)
            Bs_idx = mother_idx

      if count_phi != 2:
        print 'WARNING - did not find exactly two phi daughters'
        continue

      phi1_idx = phi_idx[0]
      phi2_idx = phi_idx[1]

      count_phi1_daughters = 0
      count_phi2_daughters = 0

      k1_idx_list = []
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) == 321 and entry.GenPart_genPartIdxMother[icand] == phi1_idx:
          count_phi1_daughters += 1
          k1_idx_list.append(icand)

      if count_phi1_daughters != 2:
        print 'WARNING - did not find exactly two K daughters of phi1'
        continue

      k1_idx = k1_idx_list[0]
      k2_idx = k1_idx_list[1]

      k2_idx_list = []
      for icand in range(0, entry.nGenPart):
        if abs(entry.GenPart_pdgId[icand]) == 321 and entry.GenPart_genPartIdxMother[icand] == phi2_idx:
          count_phi2_daughters += 1
          k2_idx_list.append(icand)

      if count_phi2_daughters != 2:
        print 'WARNING - did not find exactly two K daughters of phi2'
        continue

      k3_idx = k2_idx_list[0]
      k4_idx = k2_idx_list[1]

      if eta_bins != None and pt_bins != None:
        # fetch n_deno and n_num with acceptance cuts
        if entry.GenPart_pt[k1_idx] > 0.5 and abs(entry.GenPart_eta[k1_idx]) < 2.5 \
        and entry.GenPart_pt[k2_idx] > 0.5 and abs(entry.GenPart_eta[k2_idx]) < 2.5 \
        and entry.GenPart_pt[k3_idx] > 0.5 and abs(entry.GenPart_eta[k3_idx]) < 2.5 \
        and entry.GenPart_pt[k4_idx] > 0.5 and abs(entry.GenPart_eta[k4_idx]) < 2.5 \
        and has_trgmu:
             count_gen += 1

             obj_pt = entry.GenPart_pt[Bs_idx]
             obj_eta = abs(entry.GenPart_eta[Bs_idx])

             for ibin_eta, eta_bin in enumerate(eta_bins):
               for ibin_pt, pt_bin in enumerate(pt_bins):
                 bin_min_eta, bin_max_eta = eta_bin
                 bin_min_pt, bin_max_pt = pt_bin

                 if obj_eta > bin_min_eta and obj_eta < bin_max_eta and obj_pt > bin_min_pt and obj_pt < bin_max_pt: 
                   # for nanoAOD efficiency
                   n_deno[ibin_eta][ibin_pt] = n_deno[ibin_eta][ibin_pt] + 1
                   if matching_cond['candidate']: n_num[ibin_eta][ibin_pt] = n_num[ibin_eta][ibin_pt] + 1

                   # for efficiency from tracking POG
                   k1_eff = self.getEfficiencyFromPOGFile(eta = abs(entry.GenPart_eta[k1_idx]), pt = entry.GenPart_pt[k1_idx])
                   k2_eff = self.getEfficiencyFromPOGFile(eta = abs(entry.GenPart_eta[k2_idx]), pt = entry.GenPart_pt[k2_idx])
                   k3_eff = self.getEfficiencyFromPOGFile(eta = abs(entry.GenPart_eta[k3_idx]), pt = entry.GenPart_pt[k3_idx])
                   k4_eff = self.getEfficiencyFromPOGFile(eta = abs(entry.GenPart_eta[k4_idx]), pt = entry.GenPart_pt[k4_idx])
                   eff_POG += (k1_eff * k2_eff * k3_eff * k4_eff)

                   
    efficiency_POG = eff_POG / count_gen
    print 'efficiency_POG = {}'.format(efficiency_POG)

    # compute efficiency
    f = open('{}/numbers.txt'.format(self.outputdir), 'a')
    f.write('\n')
    if eta_bins != None and pt_bins != None:
        for ibin_eta, eta_bin in enumerate(eta_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            efficiency[ibin_eta][ibin_pt] = float(n_num[ibin_eta][ibin_pt]) / float(n_deno[ibin_eta][ibin_pt]) if float(n_deno[ibin_eta][ibin_pt]) != 0. else 0.
            if n_num[ibin_eta][ibin_pt] == 0: n_num[ibin_eta][ibin_pt] = 1e-11
            if float(n_num[ibin_eta][ibin_pt]) != 0 and float(n_deno[ibin_eta][ibin_pt]) != 0:
              error[ibin_eta][ibin_pt] = efficiency[ibin_eta][ibin_pt] * ( math.sqrt(float(n_num[ibin_eta][ibin_pt]))/float(n_num[ibin_eta][ibin_pt])  + math.sqrt(float(n_deno[ibin_eta][ibin_pt]))/float(n_deno[ibin_eta][ibin_pt]) )     
            else:
              error[ibin_eta][ibin_pt] = 0
            # for aesthetics
            if efficiency[ibin_eta][ibin_pt] == 0.: efficiency[ibin_eta][ibin_pt] = 1e-9

        for ibin_eta, eta_bin in enumerate(eta_bins):
          for ibin_pt, pt_bin in enumerate(pt_bins):
            f.write('\n {} {} {} {} {}+-{}'.format(eta_bin, pt_bin, n_deno[ibin_eta][ibin_pt], n_num[ibin_eta][ibin_pt], efficiency[ibin_eta][ibin_pt], error[ibin_eta][ibin_pt]))
            print '{} {} {} {} {}+-{}'.format(eta_bin, pt_bin, n_deno[ibin_eta][ibin_pt], n_num[ibin_eta][ibin_pt], efficiency[ibin_eta][ibin_pt], error[ibin_eta][ibin_pt])
    f.close()



    return efficiency, error


  def plot2DEfficiency(self):
    gStyle.SetPadRightMargin(0.16)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetOptStat(0)
    gStyle.SetPaintTextFormat(".2f")

    efficiency = -99.
    error = -99.
    efficiency, error = self.getEfficiency(eta_bins=self.eta_bins, pt_bins=self.pt_bins)

    canv_name = '2d_canv_{}'.format(self.outdirlabel)
    canv = self.tools.createTCanvas(canv_name, 900, 800)

    bin_min_eta, a = self.eta_bins[0]
    b ,bin_max_eta = self.eta_bins[len(self.eta_bins)-1]

    hist_name = 'hist{}'.format(self.outdirlabel)
    hist = ROOT.TH2D(hist_name, hist_name, len(self.eta_bins), 1, 100, len(self.pt_bins), 1, 100)
    for ibin_eta, eta_bin in enumerate(self.eta_bins):
      for ibin_pt, pt_bin in enumerate(self.pt_bins):
        hist.Fill(str(eta_bin), str(pt_bin), efficiency[ibin_eta][ibin_pt])

    hist.Draw('colz')
    hist.SetMarkerSize(2)
    hist.Draw('text' +'same')

    hist.SetTitle(self.title)
    hist.GetXaxis().SetTitle('B_{s} |#eta|')
    hist.GetXaxis().SetLabelSize(0.037)
    hist.GetXaxis().SetTitleSize(0.042)
    hist.GetXaxis().SetTitleOffset(1.1)
    ylabel = 'B_{s} #it{p}_{T} (GeV)'
    hist.GetYaxis().SetTitle(ylabel)
    hist.GetYaxis().SetLabelSize(0.037)
    hist.GetYaxis().SetTitleSize(0.042)
    hist.GetYaxis().SetTitleOffset(1.8)
    hist.GetZaxis().SetTitle('Efficiency')
    hist.GetZaxis().SetLabelSize(0.037)
    hist.GetZaxis().SetTitleSize(0.042)
    hist.GetZaxis().SetTitleOffset(1.2)
    hist.GetZaxis().SetRangeUser(-1e-9, 1)

    ROOT.gPad.Modified()
    ROOT.gPad.Update()

    canv.cd()
    canv.SaveAs('{}/2d_pt_eta.png'.format(self.outputdir))
    canv.SaveAs('{}/2d_pt_eta.pdf'.format(self.outputdir))

    
  def plot1DEfficiency(self, binning='eta'):
    if binning not in ['eta', 'pt']:
      raise RuntimeError("Please choose in which quantity to bin ('eta', 'bin')")

    if binning == 'eta':
      bins = self.eta_bins
      eta_bins = bins
      pt_bins = [(0, 1000)]
    elif binning == 'pt':
      bins = self.pt_bins
      eta_bins = [(0, 1000)]
      pt_bins = bins

    efficiency, error = self.getEfficiency(eta_bins=eta_bins, pt_bins=pt_bins)

    canv_name = 'canv_{}_{}'.format(binning, self.outdirlabel)
    canv = self.tools.createTCanvas(canv_name, 900, 800)

    pad = ROOT.TPad('pad', 'pad', 0, 0, 1, 1)
    pad.Draw()
    pad.SetGrid()
    pad.cd()

    graph = ROOT.TGraphAsymmErrors()
    for ibin, bin_ in enumerate(bins):
      bin_min, bin_max = bin_
      x = (float(bin_min) + float(bin_max))/2.
      y = efficiency[ibin][0] if binning == 'eta' else efficiency[0][ibin]
      err = error[ibin][0] if binning == 'eta' else error[0][ibin]
      point = graph.GetN()
      graph.SetPoint(point, x, y)
      graph.SetPointError(point, (bin_max - bin_min)/2., (bin_max - bin_min)/2., err, err)
      
    graph.Draw('AP')  

    graph.SetTitle(self.title)
    graph.SetLineColor(ROOT.kBlue+2)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue+2)
    if binning == 'eta': xlabel = 'B_{s} |#eta|'
    elif binning == 'pt': xlabel = 'B_{s} #it{p}_{T} (GeV)'
    graph.GetXaxis().SetTitle(xlabel)
    graph.GetXaxis().SetLabelSize(0.037)
    graph.GetXaxis().SetTitleSize(0.042)
    graph.GetXaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetTitle('Efficiency')
    graph.GetYaxis().SetLabelSize(0.037)
    graph.GetYaxis().SetTitleSize(0.042)
    graph.GetYaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetRangeUser(0, 1.)

    #label = ROOT.TPaveText(0.6,0.7,0.85,0.8,"brNDC")
    #text = 'inclusive pT' if binning == 'eta' else 'inclusive eta'
    #label.AddText(text)
    #label.SetBorderSize(0)
    #label.SetFillColorAlpha(0, 0)
    #label.SetTextSize(0.035)
    #label.SetTextFont(42)
    #label.SetTextAlign(11)
    #label.Draw()

    canv.cd()
    canv.SaveAs('{}/1d_{}.png'.format(self.outputdir, binning))
    canv.SaveAs('{}/1d_{}.pdf'.format(self.outputdir, binning))




if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  do_studyRecoEfficiency = True
  do_studyMuonIDEfficiency = False

  #TODO implement check that all gen entries are in nano

  if do_studyRecoEfficiency:
    filename = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_loose.root' #TODO take fullgen?

    eta_bins = [(0, 10)]
    pt_bins = [(0, 1000)] 
    #eta_bins = [(0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 2.5)]
    #pt_bins = [(0.5, 5), (5, 10), (10, 15), (15, 30), (30, 100)] 

    outdirlabel = '102X_crab_trgmu_filter_loose'

    title = ''

    analyser = EfficiencyAnalyser(filename=filename, eta_bins=eta_bins, pt_bins=pt_bins, title=title, outdirlabel=outdirlabel)

    analyser.plot2DEfficiency()
    #analyser.plot1DEfficiency(binning='eta')
    #analyser.plot1DEfficiency(binning='pt')
