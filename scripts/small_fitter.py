import ROOT
from ROOT import RooFit
import math

from tools import Tools

class Fitter(Tools):
  def __init__(self, filename='', nbins=250, title=' ', outdirlabel='testing', label=''):
    self.tools = Tools()
    self.filename = filename
    self.nbins = nbins
    self.title = title
    self.outdirlabel = outdirlabel
    self.label = label

  def performFit(self):
    # open the file and get the tree
    inputfile = ROOT.TFile.Open(self.filename)
    treename = 'signal_tree'
    #treename = 'tree'
    tree = self.tools.getTree(inputfile, treename)

    # get signal mass
    signal_mass = 5.367

    # we declare invariant mass as a RooRealVar (for the residual) and as a RooDataHist (for the fit):
    binMin = 5.1 #5.07 #signal_mass - 0.15*signal_mass
    binMax = 5.6 #5.67 #signal_mass + 0.15*signal_mass
    nbins = self.nbins

    signal_min = 5.25
    signal_max = 5.5

    hist = ROOT.TH1D("hist", "hist", nbins, binMin, binMax)
    c1 = self.tools.createTCanvas(name='c1', dimx=700, dimy=600)
    branch_name = 'bs_mass_corr'
    #tree.Draw("{}>>hist".format(branch_name))
    cond = 'score > 0.96'
    tree.Draw("{}>>hist".format(branch_name), cond)

    invmass = ROOT.RooRealVar("invmass","invmass", binMin, binMax)
    invmass.setBins(nbins)

    rdh = ROOT.RooDataHist("rdh", "rdh", ROOT.RooArgList(invmass), hist)

    mult = 1.
    nsig = ROOT.RooRealVar('nsig','number of signal events', 20000*mult,0,1E07)
    nbkg = ROOT.RooRealVar('nbkg','number of bkg events (combinatorial)', 1000*mult,0,1.0E07) 

    # Define the PDF to fit: 
    # Double sided crystal ball
    # we declare all the parameters needed for the fits 
    #mean_init = signal_mass - 0.01*signal_mass
    #mean_min = signal_mass - 0.05*signal_mass
    #mean_max = signal_mass + 0.05*signal_mass
    #mean  = ROOT.RooRealVar("mean","mean", mean_init, mean_min, mean_max)

    #sigma = ROOT.RooRealVar("sigma","sigma", 0.01, 0.005, 0.15)

    #alpha_1 = ROOT.RooRealVar("alpha_1", "alpha_1", -1, -3, 3)
    ##RooRealVar *n_1      = new RooRealVar("n_1", "n_1", 0, 15)
    #n_1 = ROOT.RooRealVar("n_1", "n_1", 0, 3)
    #alpha_2 = ROOT.RooRealVar("alpha_2", "alpha_2", 1, -3, 3)
    ##RooRealVar *n_2      = new RooRealVar("n_2", "n_2", 0, 15)
    #n_2 = ROOT.RooRealVar("n_2", "n_2", 0, 3)

    #CBpdf_1 = ROOT.RooCBShape("CBpdf_1", "CBpdf_1", invmass, mean, sigma, alpha_1, n_1)
    #CBpdf_2 = ROOT.RooCBShape("CBpdf_2", "CBpdf_2", invmass, mean, sigma, alpha_2, n_2)

    ## defines the relative importance of the two CBs
    #sigfrac = ROOT.RooRealVar("sigfrac","sigfrac", 0.5, 0.0 ,1.0)

    # we add the two CB pdfs together
    #fitmodel_signal = ROOT.RooAddPdf("fitmodel_signal", "fitmodel_signal", CBpdf_1, CBpdf_2, sigfrac, ROOT.RooArgList(nsig))

    voigt_mean  = ROOT.RooRealVar('voigt_mean', 'voigt_mean', signal_mass, 5.3, 5.4)
    #voigt_sigma = ROOT.RooRealVar('voigt_sigma', 'voigt_sigma', 0.02, 0.01, 0.03)
    #voigt_sigma = ROOT.RooRealVar('voigt_sigma', 'voigt_sigma', 0.025, 0.0245, 0.0255)
    voigt_sigma = ROOT.RooRealVar('voigt_sigma', 'voigt_sigma', 0.025)
    #voigt_width = ROOT.RooRealVar('voigt_width', 'voigt_width', 0.02, 0.01, 0.03)
    #voigt_width = ROOT.RooRealVar('voigt_width', 'voigt_width', 0.023, 0.02, 0.0235)
    voigt_width = ROOT.RooRealVar('voigt_width', 'voigt_width', 0.0228)
    voigt       = ROOT.RooVoigtian('voigt', 'voigt', invmass, voigt_mean, voigt_sigma, voigt_width)
    
    # choose signal model
    fitmodel_signal = ROOT.RooAddPdf('fitmodel_signal','signal',ROOT.RooArgList(voigt),ROOT.RooArgList(nsig))


    # background
    exp_a = ROOT.RooRealVar('exp_a', 'exp_a', -0.0001, -5.,0.)
    exp = ROOT.RooExponential('exp', 'Exponential', invmass, exp_a)
    fitmodel_bkg = ROOT.RooAddPdf('fitmodel_bkg', 'Combinatorial bkg', ROOT.RooArgList(exp), ROOT.RooArgList(nbkg))
        
    # model
    model = ROOT.RooAddPdf('model',
                           'signal + bkg(comb)',
                           ROOT.RooArgList(fitmodel_signal,fitmodel_bkg),
                           ROOT.RooArgList(nsig,nbkg)) 

    # we define the frame where to plot
    canv = self.tools.createTCanvas(name="canv", dimx=900, dimy=800)

    # and the two pads
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.99, 0.99)
    pad1.SetLeftMargin(0.15)
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.2)
    pad2.SetLeftMargin(0.15)
    pad1.Draw()
    pad2.Draw()

    frame = invmass.frame(ROOT.RooFit.Title(self.title))

    # plot the data
    rdh.plotOn(frame, ROOT.RooFit.Name("data"))

    # fit the PDF to the data
    #RooFitResult *result
    result = model.fitTo(rdh, ROOT.RooFit.Extended(True))
    #fit_range_min = signal_mass - 0.1*signal_mass
    #fit_range_max = signal_mass + 0.1*signal_mass
    #invmass.setRange("peak", fit_range_min, fit_range_max)
    #result = model.fitTo(rdh, ROOT.RooFit.Extended(True), ROOT.RooFit.Range("peak"))

    # plot the fit    
    model.plotOn(frame, ROOT.RooFit.Components('fitmodel_signal'),ROOT.RooFit.LineColor(ROOT.kOrange+7), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name('signal'))
    model.plotOn(frame, ROOT.RooFit.Components('fitmodel_bkg'),ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name('bkg_comb'))
    model.plotOn(frame, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("model"), ROOT.RooFit.Components("model"))

    # and write the fit parameters
    model.paramOn(frame,   
         ROOT.RooFit.Layout(0.2, 0.4, 0.8),
         ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))
         )

    frame.getAttText().SetTextSize(0.03)
    frame.getAttLine().SetLineColorAlpha(0, 0)
    frame.getAttFill().SetFillColorAlpha(0, 0)

    # we compute the chisquare
    chisquare = frame.chiSquare("model","data")

    # and print it
    label1 = ROOT.TPaveText(0.62,0.65,0.72,0.8,"brNDC")
    label1.SetBorderSize(0)
    label1.SetFillColor(ROOT.kWhite)
    label1.SetTextSize(0.03)
    label1.SetTextFont(42)
    label1.SetTextAlign(11)
    qte = '#chi^{2}/ndof'
    #label1.AddText('Double sided Crystal Ball')
    #label1.AddText('{} = {}'.format(qte, round(chisquare, 2)))
    #label1.AddText(self.label)
    print "chisquare = {}".format(chisquare)

    n_sig = nsig.getVal()
    n_bkg = nbkg.getVal() * (signal_max - signal_min) / (binMax - binMin) # number of background events in signal window
    significance = n_sig / (n_sig + n_bkg) #math.sqrt(n_sig + n_bkg)

    lumi = 33.6
    lumi_proj = 41.6 * 2 + 112

    n_sig_proj = n_sig * lumi_proj / lumi
    n_bkg_proj = n_bkg * lumi_proj / lumi 
    significance_proj = n_sig_proj / (n_sig_proj + n_bkg_proj) # math.sqrt(n_sig_proj + n_bkg_proj)

    print '\n significance: '
    print '{} & {} & {} & {} & {} & {}'.format(int(round(n_sig, 0)), int(round(n_bkg, 0)), round(significance, 2), int(round(n_sig_proj, 0)), int(round(n_bkg_proj, 0)), round(significance_proj, 2))
    print '\n'

    # We define and plot the residuals    
    # construct a histogram with the pulls of the data w.r.t the curve
    hpull = frame.pullHist()
    for i in range(0, frame.GetNbinsX()):
       hpull.SetPointError(i,0,0,0,0)

    # create a new frame to draw the pull distribution and add the distribution to the frame
    frame2 = invmass.frame(ROOT.RooFit.Title(" "))
    frame2.addPlotable(hpull,"P")#,"E3")

    # plot of the curve and the fit
    canv.cd()
    pad1.cd()

    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetTitle("#it{m}(K^{+}K^{-}K^{+}K^{-}) (GeV)")
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.Draw()
    label1.Draw()

    # plot of the residuals
    pad2.cd()
    ROOT.gPad.SetLeftMargin(0.15) 
    ROOT.gPad.SetPad(0.01,0.01,0.99,0.2)

    frame2.GetYaxis().SetNdivisions(3)
    frame2.GetYaxis().SetLabelSize(0.17)
    frame2.GetYaxis().SetTitleSize(0.17)
    frame2.GetYaxis().SetTitleOffset(0.24)
    frame2.GetYaxis().SetRangeUser(-5,5)  
    frame2.GetYaxis().SetTitle("Pulls") 
    frame2.GetXaxis().SetTitle("")  
    frame2.GetXaxis().SetLabelOffset(5) 
    frame2.Draw()

    line = ROOT.TLine()
    line.DrawLine(binMin,0,binMax,0)
    line.SetLineColor(2)
    line.DrawLine(binMin,-3,binMax,-3)
    line.DrawLine(binMin,3,binMax,3)

    # save output
    canv.cd()
    outputdir = self.tools.getOutDir('./myPlots/fits', self.outdirlabel)
    canv.SaveAs("{}/fit_{}.png".format(outputdir, label))
    canv.SaveAs("{}/fit_{}.pdf".format(outputdir, label))

    # additionally, get the pull histogram
    canv_pull = self.tools.createTCanvas(name="canv_pull", dimx=700, dimy=600)

    hist_pull = ROOT.TH1D("hist_pull", "hist_pull", 120, -5, 5)

    for i in range(0, hpull.GetN()):

      x = ROOT.Double()
      point = ROOT.Double()
      hpull.GetPoint(i,x,point) 
      if x<binMin or x>binMax: continue

      hist_pull.Fill(point)

    hist_pull.SetTitle("")
    hist_pull.SetLineColor(4)
    hist_pull.SetLineWidth(2)
    hist_pull.Draw()

    Xaxis = hist_pull.GetXaxis()
    Yaxis = hist_pull.GetYaxis()
    Xaxis.SetTitle("pulls")
    Xaxis.SetTitleSize(0.045)
    Xaxis.SetLabelSize(0.045)
    Xaxis.SetTitleOffset(1.1)
    Yaxis.SetTitleSize(0.045)
    Yaxis.SetLabelSize(0.045)
    Yaxis.SetTitleOffset(1.26)
    ROOT.gStyle.SetOptStat(0)
    hist_pull.Draw()

    # and fit it
    fgauss= ROOT.TF1("fgauss", "gaus", -5, 5)
    fgauss.SetLineColor(2)
    hist_pull.Fit("fgauss")
    fgauss.Draw("same")
    ROOT.gStyle.SetOptFit(0011)

    #canv_pull.SaveAs("test.png")  
    canv_pull.SaveAs("{}/pulls_{}.png".format(outputdir, label))
    canv_pull.SaveAs("{}/pulls_{}.pdf".format(outputdir, label))


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)

  outdirlabel = './myPlots/selection/' 

  #label = 'signal'
  #filename = '/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/flat_bparknano.root'
  #fitter = Fitter(filename=filename, nbins=100, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #label = 'cutbased1' 
  #filename = './myPlots/selection/tree_cutbased1.root'
  #fitter = Fitter(filename=filename, nbins=100, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  #label = 'tree_GK2' 
  #filename = './myPlots/selection/tree_GK2.root'
  #fitter = Fitter(filename=filename, nbins=100, outdirlabel=outdirlabel, label=label)
  #fitter.performFit()

  label = 'mva_test_2025Mar11_16h33m47s' 
  filename = './outputs/test_2025Mar11_16h33m47s.root'
  fitter = Fitter(filename=filename, nbins=100, outdirlabel=outdirlabel, label=label)
  fitter.performFit()


