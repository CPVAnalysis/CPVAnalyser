import os
import os.path
from os import path
import ROOT


class Tools(object):
  def getRootFile(self, filename, with_ext=False):
    if not with_ext:
      f = ROOT.TFile.Open(filename, 'READ')
    else:
      f = ROOT.TFile.Open('root://t3dcachedb.psi.ch:1094/'+filename, 'READ')
    return f


  def getTree(self, rootfile, tree_name):
    tree = rootfile.Get(tree_name)
    if not tree:
      raise RuntimeError('Tree of name "{}" was not found in {}'.format(tree_name, rootfile.GetName()))
    return tree


  def createHisto(self, tree, quantity, hist_name='hist', branchname='flat', weight=-99, selection=''):
    ROOT.TH1.SetDefaultSumw2()
    hist = ROOT.TH1D(hist_name, hist_name, quantity.nbins, quantity.bin_min, quantity.bin_max)

    if selection == '' and weight == -99: selection_string = selection
    elif selection != '' and weight == -99: selection_string = '({})'.format(selection)
    elif selection == ''  and weight != -99: selection_string = '({})'.format(weight)
    else: selection_string = '({sel}) * ({wght})'.format(sel=selection, wght=weight)

    print(selection_string)

    qte = quantity.name_flat if branchname != 'nano' else quantity.name_nano
    tree.Project(hist_name, qte, selection_string)

    hist.SetDirectory(0)
    return hist


  def getRatioHistogram(self, hist1, hist2): 
    hist_ratio = hist1.Clone('hist_ratio')
    hist_ratio.Divide(hist2)
    return hist_ratio


  def createTCanvas(self, name, dimx=800, dimy=700):
    canv = ROOT.TCanvas(name, name, dimx, dimy)
    ROOT.SetOwnership(canv, False)
    return canv


  def getRootTLegend(self, xmin=0.65, ymin=0.7, xmax=0.85, ymax=0.9, size=0.03, do_alpha=True):
    legend = ROOT.TLegend(xmin, ymin, xmax, ymax)
    legend.SetTextSize(size)
    legend.SetLineColor(0)
    if do_alpha: legend.SetFillColorAlpha(0, 0)
    else: legend.SetFillColor(0)
    legend.SetBorderSize(0)
    return legend
    

  def getTextBox(self, style, xmin, ymin, xmax, ymax, text, size=0.11, colour=ROOT.kBlack):
    box = ROOT.TPaveText(xmin, ymin, xmax, ymax, style)
    box.AddText(text)
    box.SetBorderSize(0)
    box.SetFillColorAlpha(0, 0)
    box.SetTextColor(colour)
    box.SetTextSize(size)
    box.SetTextFont(42)
    box.SetTextAlign(11)
    return box


  def printLatexBox(self, x, y, text, size=0.04, pos='center', font=42):
    box = ROOT.TLatex()
    box.SetNDC()
    box.SetTextFont(font)
    if pos == 'center': box.SetTextAlign(22) 
    else: box.SetTextAlign(12) 
    box.SetTextSize(size)    
    box.DrawLatex(x, y, text)
  
  
  def getRootXAxis(self, hist, title='', label_size=0.037, title_size=0.042, offset=1.1, xmin=-99, xmax=-99): 
    ROOT.SetOwnership(hist, True)
    hist.GetXaxis().SetTitle(title)
    hist.GetXaxis().SetLabelSize(label_size)
    hist.GetXaxis().SetTitleSize(title_size)
    hist.GetXaxis().SetTitleOffset(offset)
    if xmin != -99 and xmax != -99:
      hist.GetXaxis().SetRangeUser(xmin, xmax)
    return hist


  def getRootYAxis(self, hist, title='', label_size=0.037, title_size=0.042, offset=1.1, ymin=-99, ymax=-99): 
    ROOT.SetOwnership(hist, False)
    hist.GetYaxis().SetTitle(title)
    hist.GetYaxis().SetLabelSize(label_size)
    hist.GetYaxis().SetTitleSize(title_size)
    hist.GetYaxis().SetTitleOffset(offset)
    if ymin != -99 and ymax != -99:
      hist.GetXaxis().SetRangeUser(ymin, ymax)
    return hist


  def printCMSTag(self, pad, cms_tag, size=0.55, offset=0.08, font=52):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    # print CMS
    tag.SetTextFont(61)
    tag.SetTextAlign(11) 
    tag.SetTextSize(size*pad.GetTopMargin())    
    tag.DrawLatex(pad.GetLeftMargin(), 1-pad.GetTopMargin()+0.2*pad.GetTopMargin(), 'CMS')
    # print CMS tag
    tag.SetTextFont(font)
    tag.SetTextSize(0.9*size*pad.GetTopMargin())
    tag.SetTextAlign(11)
    tag.DrawLatex(pad.GetLeftMargin()+offset, 1-pad.GetTopMargin()+0.2*pad.GetTopMargin(), cms_tag)      
    pad.Update()


  def printInnerCMSTag(self, pad, cms_tag, print_tag=False, x_pos=0.15, y_pos=0.83, size=0.55):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    # print CMS
    tag.SetTextFont(61)
    tag.SetTextAlign(11) 
    tag.SetTextSize(size*pad.GetTopMargin())    
    tag.DrawLatex(x_pos, y_pos, 'CMS')
    ## print CMS tag
    tag.SetTextFont(52)
    tag.SetTextSize(0.9*size*pad.GetTopMargin())
    tag.SetTextAlign(11)
    x_pos_tag = x_pos
    y_pos_tag = y_pos - 0.06
    tag.DrawLatex(x_pos_tag, y_pos_tag, cms_tag)      
    pad.Update()


  def printCMSTagInFrame(self, pad, cms_tag, size=0.55, offset=0.11):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    # print CMS
    tag.SetTextFont(61)
    tag.SetTextAlign(11) 
    tag.SetTextSize(size*pad.GetTopMargin())    
    tag.DrawLatex(pad.GetLeftMargin()+0.2*pad.GetLeftMargin(), 1-pad.GetTopMargin()-0.8*pad.GetTopMargin(), 'CMS')
    # print CMS tag
    tag.SetTextFont(52)
    tag.SetTextSize(0.9*size*pad.GetTopMargin())
    tag.SetTextAlign(11)
    tag.DrawLatex(pad.GetLeftMargin()+0.2*pad.GetLeftMargin()+offset, 1-pad.GetTopMargin()-0.8*pad.GetTopMargin(), cms_tag)      
    pad.Update()


  def printLumiTag(self, pad, lumi, size=0.43, offset=0.57):
    pad.cd()
    tag = ROOT.TLatex()
    tag.SetNDC()
    lumi_text = str(round(lumi, 2)) + ' fb^{-1} (13 TeV)'
    tag.SetTextFont(42)
    tag.SetTextAlign(11) 
    tag.SetTextSize(0.9*size*pad.GetTopMargin())    
    tag.DrawLatex(pad.GetLeftMargin()+offset, 1-pad.GetTopMargin()+0.2*pad.GetTopMargin(), lumi_text)
    pad.Update()


  def getOutDir(self, maindir, outdirlabel, do_shape=False, do_luminorm=False, do_stack=False, do_log=False, add_overflow=False):
    if not path.exists(maindir):
      os.system('mkdir -p {}'.format(maindir))
    os.system('cp ../data/index.php {}'.format(maindir))
    dirlabel = outdirlabel

    outputdir = '{}/{}'.format(maindir, dirlabel)
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    os.system('cp ../data/index.php {}'.format(outputdir))
    os.system('cp ../data/index.php {}/..'.format(outputdir))

    norm = None
    if do_shape: norm = 'shape'
    elif do_luminorm: norm = 'luminorm'

    #if not do_shape and not do_stack and not do_log: dirlabel += '/plain'
    if do_stack: 
      if norm==None and not do_log and not do_luminorm: dirlabel += '/stack'
      elif norm==None and do_log:  dirlabel += '/stack/log'
      elif norm!=None and not do_log and not add_overflow: dirlabel += '/stack_{}'.format(norm)
      elif norm!=None and not do_log and add_overflow: dirlabel += '/stack_{}_overflow'.format(norm)
      elif norm!=None and do_log and not add_overflow: dirlabel += '/stack_{}_log'.format(norm)
      else: dirlabel += '/stack_{}_log_overflow'.format(norm)
    else:
      if norm==None and do_log: dirlabel += '/log'
      elif norm==None and not do_log: dirlabel += '/lin'
      elif norm!=None and not do_log and not add_overflow: dirlabel += '/{}'.format(norm)
      elif norm!=None and do_log and not add_overflow: dirlabel += '/{}_log'.format(norm)
      elif norm!=None and not do_log and add_overflow: dirlabel += '/{}_overflow'.format(norm)
      else: dirlabel += '/{}_log_overflow'.format(norm)

    outputdir = '{}/{}'.format(maindir, dirlabel)
    
    if not path.exists(outputdir):
      os.system('mkdir -p {}'.format(outputdir))
    os.system('cp ../data/index.php {}'.format(outputdir))

    return outputdir


