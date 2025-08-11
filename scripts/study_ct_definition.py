import os
from os import path
import ROOT

from tools import Tools

class Quantity(object):
    def __init__(self, name, title, legend, label, colour):
        self.name = name
        self.title = title
        self.legend = legend
        self.label = label
        self.colour = colour


class CtAnalyser(object):
    def __init__(self, filename, branch_gen, quantities_reco):
        self.tools = Tools()
        self.filename = filename
        self.branch_gen = branch_gen
        self.quantities_reco = quantities_reco

        self.tree_name = 'Events'

        self.outdir = './myPlots/ct_definition'
        if not path.exists(self.outdir):
            os.system('mkdir -p {}'.format(self.outdir))


    def plot_ct(self, quantity_reco, extra_selection=None):
        f = self.tools.getRootFile(self.filename)
        tree = self.tools.getTree(f, self.tree_name)

        canv_name = 'canv' + quantity_reco.label
        if extra_selection != None: canv_name += '_selected'
        canv = self.tools.createTCanvas(canv_name)
        canv.SetLogy()

        bin_min = -0.1
        bin_max = 0.5
        nbins = 100

        hist_name = 'hist'
        hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
        hist.SetLineColor(quantity_reco.colour)
        hist.SetLineWidth(2)

        selection = 'BsToPhiPhiTo4K_isMatched == 1' 
        if extra_selection != None:
            selection += ' && {}'.format(extra_selection)

        tree.Draw('{} >> hist'.format(quantity_reco.name), selection)

        frame = hist.Clone('frame')
        frame.SetTitle('')

        frame.GetXaxis().SetTitle(quantity_reco.title)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetTitleOffset(1.1)

        frame.GetYaxis().SetTitle('Entries')
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleOffset(1.)

        frame.Draw()

        ROOT.gStyle.SetOptStat(0)

        figname = quantity_reco.label
        if extra_selection != None: figname += '_selected'
        canv.SaveAs('{}/{}.png'.format(self.outdir, figname))
        canv.SaveAs('{}/{}.pdf'.format(self.outdir, figname))


    def plot_ct_all(self, do_extra_selection=False):
        f = self.tools.getRootFile(self.filename)
        tree = self.tools.getTree(f, self.tree_name)

        canv_name = 'canv_all'
        if do_extra_selection: canv_name += '_selected'
        canv = self.tools.createTCanvas(canv_name)
        canv.SetLogy()

        leg = self.tools.getRootTLegend(0.6, 0.5, 0.9, 0.8, 0.035)

        bin_min = -0.1
        bin_max = 0.5
        nbins = 100

        hists = []
        for quantity_reco in self.quantities_reco:
            hist_name = 'hist' + quantity_reco.name
            hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
            hist.SetLineColor(quantity_reco.colour)
            hist.SetLineWidth(2)

            selection = 'BsToPhiPhiTo4K_isMatched == 1' 
            if do_extra_selection:
                extra_selection = 'abs({}/{} -1) < 1'.format(quantity_reco.name, self.branch_gen)
                selection += ' && {}'.format(extra_selection)

            tree.Draw('{} >> {}'.format(quantity_reco.name, hist_name), selection)
            hist.SetDirectory(0)
            hists.append(hist)
            leg.AddEntry(hist, quantity_reco.legend)

            #print('{} -- Entries: {} Mean: {} StdDev: {}'.format(quantity_reco.label, hist.GetEntries(), hist.GetMean(), hist.GetStdDev()))

        frame = hists[0].Clone('frame')
        frame.SetTitle('')

        frame.GetXaxis().SetTitle('B_{s} ct [cm]')
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetTitleOffset(1.)

        frame.GetYaxis().SetTitle('Entries')
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleOffset(1.)

        frame.Draw()
        for hist in hists:
            hist.Draw('same')

        leg.Draw()

        ROOT.gStyle.SetOptStat(0)

        figname = 'ct_all'
        if do_extra_selection: figname += '_selected'
        canv.SaveAs('{}/{}.png'.format(self.outdir, figname))
        canv.SaveAs('{}/{}.pdf'.format(self.outdir, figname))


    def plot_2D_reco_vs_gen(self, quantity_reco, extra_selection=None):
        f = self.tools.getRootFile(self.filename)
        tree = self.tools.getTree(f, self.tree_name)

        canv_name = 'canv' + quantity_reco.label
        if extra_selection != None: canv_name += '_selected'
        canv = self.tools.createTCanvas(canv_name)

        bin_min = -0.1
        bin_max = 0.5
        nbins = 100

        hist_name = 'hist'
        hist = ROOT.TH2D(hist_name, hist_name, nbins, bin_min, bin_max, nbins, bin_min, bin_max)

        selection = 'BsToPhiPhiTo4K_isMatched == 1' 
        if extra_selection != None:
            selection += ' && {}'.format(extra_selection)

        tree.Draw('{}:{} >> hist'.format(self.branch_gen, quantity_reco.name), selection, 'colz')

        frame = hist.Clone('frame')
        frame.SetTitle('')

        frame.GetXaxis().SetTitle(quantity_reco.title)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetTitleOffset(1.1)

        frame.GetYaxis().SetTitle('gen B_{s} ct [cm]')
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleOffset(1.)

        frame.Draw()

        ROOT.gStyle.SetOptStat(0)

        figname = '2d_gen_vs_{}'.format(quantity_reco.label)
        if extra_selection != None: figname += '_selected'
        canv.SaveAs('{}/{}.png'.format(self.outdir, figname))
        canv.SaveAs('{}/{}.pdf'.format(self.outdir, figname))


    def plot_1D_reco_over_gen(self, quantity_reco, extra_selection=None):
        f = self.tools.getRootFile(self.filename)
        tree = self.tools.getTree(f, self.tree_name)

        canv_name = 'canv' + quantity_reco.label
        if extra_selection != None: canv_name += '_selected'
        canv = self.tools.createTCanvas(canv_name)
        canv.SetLogy()

        if extra_selection == None:
            bin_min = -100
            bin_max = 100
        else:
            bin_min = -10
            bin_max = 10
        nbins = 100

        hist_name = 'hist'
        hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)

        selection = 'BsToPhiPhiTo4K_isMatched == 1' 
        if extra_selection != None:
            selection += ' && {}'.format(extra_selection)

        tree.Draw('{}/{} >> hist'.format(quantity_reco.name, self.branch_gen), selection)

        frame = hist.Clone('frame')
        frame.SetTitle('')

        frame.GetXaxis().SetTitle(quantity_reco.title + ' / gen B_{s} ct'.replace('[cm]', ''))
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetTitleOffset(1.1)

        frame.GetYaxis().SetTitle('Entries')
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleOffset(1.)

        frame.Draw()

        ROOT.gStyle.SetOptStat(0)

        figname = '1d_{}_over_gen'.format(quantity_reco.label)
        if extra_selection != None: figname += '_selected'
        canv.SaveAs('{}/{}.png'.format(self.outdir, figname))
        canv.SaveAs('{}/{}.pdf'.format(self.outdir, figname))


    def plot_1D_reco_over_gen_all(self, do_log, extra_selection_cut):
        f = self.tools.getRootFile(self.filename)
        tree = self.tools.getTree(f, self.tree_name)

        canv_name = 'canv_all'
        canv = self.tools.createTCanvas(canv_name)
        if do_log:
            canv.SetLogy()

        canv_selected_name = 'canv_selected_all'
        canv_selected = self.tools.createTCanvas(canv_selected_name)
        if do_log:
            canv_selected.SetLogy()

        leg = self.tools.getRootTLegend(0.6, 0.5, 0.9, 0.8, 0.035)
        leg_selected = self.tools.getRootTLegend(0.6, 0.5, 0.9, 0.8, 0.035)

        hists = []
        hists_selected = []
        for quantity_reco in self.quantities_reco:
            bin_min = -100
            bin_max = 100
            nbins = 100
            hist_name = 'hist' + quantity_reco.name
            hist = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
            hist.SetLineColor(quantity_reco.colour)
            hist.SetLineWidth(2)

            selection = 'BsToPhiPhiTo4K_isMatched == 1' 

            tree.Draw('{}/{} >> {}'.format(quantity_reco.name, self.branch_gen, hist_name), selection)
            hist.SetDirectory(0)
            hists.append(hist)
            leg.AddEntry(hist, quantity_reco.legend)

            bin_min = -1 * extra_selection_cut + 1
            bin_max = extra_selection_cut + 1
            nbins = 100
            hist_selected_name = 'hist_selected' + quantity_reco.name
            hist_selected = ROOT.TH1D(hist_selected_name, hist_selected_name, nbins, bin_min, bin_max)
            hist_selected.SetLineColor(quantity_reco.colour)
            hist_selected.SetLineWidth(2)

            selection = 'BsToPhiPhiTo4K_isMatched == 1' 
            extra_selection = 'abs({}/{} -1) < {}'.format(quantity_reco.name, self.branch_gen, extra_selection_cut)
            selection += ' && {}'.format(extra_selection)

            tree.Draw('{}/{} >> {}'.format(quantity_reco.name, self.branch_gen, hist_selected_name), selection)
            hist_selected.SetDirectory(0)
            hists_selected.append(hist_selected)
            leg_selected.AddEntry(hist_selected, quantity_reco.legend)

            bin_min = -100
            bin_max = 100
            nbins = 100
            hist_zero_name = 'hist_zero' + quantity_reco.name
            hist_zero = ROOT.TH1D(hist_zero_name, hist_zero_name, nbins, bin_min, bin_max)

            selection_zero = 'BsToPhiPhiTo4K_isMatched == 1 && {}<0'.format(quantity_reco.name) 

            tree.Draw('{} >> {}'.format(quantity_reco.name, hist_zero_name), selection_zero)

            bin_min = -1 * extra_selection_cut + 1
            bin_max = extra_selection_cut + 1
            nbins = 100
            hist_zero_selected_name = 'hist_zero_selected' + quantity_reco.name
            hist_zero_selected = ROOT.TH1D(hist_zero_selected_name, hist_zero_selected_name, nbins, bin_min, bin_max)

            selection_zero_selected = 'BsToPhiPhiTo4K_isMatched == 1 && {}<0'.format(quantity_reco.name) 
            extra_selection = 'abs({}/{} -1) < {}'.format(quantity_reco.name, self.branch_gen, extra_selection_cut)
            selection_zero_selected += ' && {}'.format(extra_selection)

            tree.Draw('{} >> {}'.format(quantity_reco.name, hist_zero_selected_name), selection_zero_selected)

            print('{} -- not selected -- Entries: {} Mean: {} StdDev: {} Zero: {} ({}%)'.format(quantity_reco.label, hist.GetEntries(), round(hist.GetMean(), 3), round(hist.GetStdDev(), 3), hist_zero.GetEntries(),  round( float(hist_zero.GetEntries()) / float(hist.GetEntries())  , 2)*100    ))
            print('{} -- selected     -- Entries: {} Mean: {} StdDev: {} Efficiency: {} Zero: {} ({}%)'.format(quantity_reco.label, hist_selected.GetEntries(), round(hist_selected.GetMean(), 3), round(hist_selected.GetStdDev(), 3), round(float(hist_selected.GetEntries()) / float(hist.GetEntries()), 2), hist_zero_selected.GetEntries(),  round( float(hist_zero_selected.GetEntries()) / float(hist_selected.GetEntries()), 2)*100  ))
            print('\n')

        canv.cd()
        frame = hists[0].Clone('frame')
        frame.SetTitle('')

        frame.GetXaxis().SetTitle('reco B_{s} ct / gen B_{s} ct')
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetTitleOffset(0.9)

        frame.GetYaxis().SetTitle('Entries')
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleOffset(1.)

        frame.Draw()
        for hist in hists:
            hist.Draw('same')

        leg.Draw()

        ROOT.gStyle.SetOptStat(0)

        figname = '1d_reco_over_gen_ct_all'
        if do_log: figname += '_log'
        canv.SaveAs('{}/{}.png'.format(self.outdir, figname))
        canv.SaveAs('{}/{}.pdf'.format(self.outdir, figname))

        canv_selected.cd()
        frame = hists_selected[0].Clone('frame')
        frame.SetTitle('')

        frame.GetXaxis().SetTitle('reco B_{s} ct / gen B_{s} ct')
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleSize(0.047)
        frame.GetXaxis().SetTitleOffset(0.9)

        frame.GetYaxis().SetTitle('Entries')
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleSize(0.047)
        frame.GetYaxis().SetTitleOffset(1.)

        frame.Draw()
        for hist in hists_selected:
            hist.Draw('same')

        leg.Draw()

        ROOT.gStyle.SetOptStat(0)

        figname = '1d_reco_over_gen_ct_all_selected'
        if do_log: figname += '_log'
        canv_selected.SaveAs('{}/{}.png'.format(self.outdir, figname))
        canv_selected.SaveAs('{}/{}.pdf'.format(self.outdir, figname))


    def process(self):
        for quantity_reco in quantities_reco:
            extra_selection = 'abs({}/{} -1) < 1'.format(quantity_reco.name, self.branch_gen)

            self.plot_ct(quantity_reco, extra_selection=None)
            self.plot_ct(quantity_reco, extra_selection=extra_selection)
            self.plot_2D_reco_vs_gen(quantity_reco, extra_selection=None)
            self.plot_2D_reco_vs_gen(quantity_reco, extra_selection=extra_selection)
            self.plot_1D_reco_over_gen(quantity_reco, extra_selection=None)
            self.plot_1D_reco_over_gen(quantity_reco, extra_selection=extra_selection)

        self.plot_ct_all(do_extra_selection=False)
        self.plot_ct_all(do_extra_selection=True)
        self.plot_1D_reco_over_gen_all(do_log=False, extra_selection_cut=1)
        self.plot_1D_reco_over_gen_all(do_log=True, extra_selection_cut=1)



if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)

    filename = '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_study_ct.root'
    branch_gen = 'BsToPhiPhiTo4K_gen_Bs_ct'

    quantities_reco = [
        Quantity('BsToPhiPhiTo4K_Bs_ct_2D_cm', 'B_{s} ct (computed wrt BS) [cm]', 'BS', 'ct_2D_cm', ROOT.kBlue-3),
        Quantity('BsToPhiPhiTo4K_Bs_ct_2D_cm_posbsz', 'B_{s} ct (computed wrt BS(z(B_{s})) [cm]', 'BS(z(B_{s}))', 'ct_2D_cm_posbsz', ROOT.kGreen+3),
        Quantity('BsToPhiPhiTo4K_Bs_ct_2D_cm_posbspv', 'B_{s} ct (computed wrt BS(z(PV)) [cm]', 'BS(z(PV))', 'ct_2D_cm_posbspv', ROOT.kRed+1),
        Quantity('BsToPhiPhiTo4K_Bs_ct_2D_cm_posthepv', 'B_{s} ct (computed wrt best PV [cm]', 'best PV', 'ct_2D_cm_posthepv', ROOT.kMagenta+1),
    ]

    analyser = CtAnalyser(
            filename = filename,
            branch_gen = branch_gen,
            quantities_reco = quantities_reco,
            )

    analyser.process()
    
