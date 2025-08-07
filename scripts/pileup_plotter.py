import os
from glob import glob
import ROOT

import sys
sys.path.append('..')
from tools import Tools



class Plotter(object):
    def __init__(self, year):
        self.tools = Tools()
        self.year = year

        if self.year == '2024':
            self.parts = ['B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        else:
            self.parts = ['C', 'D', 'E', 'F', 'G']

        self.colours = [ROOT.kBlue+4, ROOT.kBlue+1, ROOT.kBlue-1, ROOT.kBlue-10, ROOT.kOrange-9, ROOT.kOrange+1, ROOT.kOrange+9, ROOT.kRed]

        ROOT.gStyle.SetOptStat(0)


    def get_label(self, name):
        idx1 = name.find(self.year) + 5
        idx2 = name.find('.root')
        label = name[idx1:idx2]

        return label


    def plot_2018(self):

        quantity = 'PV_npvs'
        #quantity = 'PV_npvsGood'

        f_profile_data = ROOT.TFile.Open('./myPlots/pileup/2018/pileup_2018.root', 'READ')
        profile_data = f_profile_data.Get('pileup_2018total')
        profile_data.Scale(1./profile_data.Integral())
        profile_data.SetFillStyle(3001)
        profile_data.SetFillColorAlpha(ROOT.kAzure+8, 0.5)

        profile_data.GetXaxis().SetRangeUser(0, 90)
        profile_data.GetXaxis().SetTitle(' ')
        profile_data.GetYaxis().SetTitle(' ')
        profile_data.SetTitle(' ')

        bin_min = 0
        bin_max = 200
        nbins = 200
        hist_name = 'hist_tot'
        hist_tot = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)

        hist_tot.SetLineColor(1)
        hist_tot.SetLineWidth(2)

        profile_filenames = [f for f in glob('myPlots/pileup/profile_2018*.root')]
        #profile_filenames = [f for f in glob('myPlots/pileup/profile*D1*.root')]

        for profile_filename in profile_filenames:
            print(profile_filename)
        
            f = ROOT.TFile.Open(profile_filename, 'READ')

            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')

            hist_tot.Add(profile)
            f.Close()

        factor = 1.2
        hist_scaled = ROOT.TH1D('hist_scaled', 'hist_scaled', nbins, bin_min*factor, bin_max*factor)
        
        for i in range(1, nbins+1):
          hist_scaled.SetBinContent(i, hist_tot.GetBinContent(i))
        
        hist_scaled.SetLineColor(2)
        hist_scaled.SetLineWidth(2)

        # normalise distributions
        hist_scaled.Scale(1./hist_scaled.Integral())
        hist_tot.Scale(1./hist_tot.Integral())

        maximum = max([hist.GetMaximum() for hist in [profile_data, hist_tot, hist_scaled]])
        profile_data.GetYaxis().SetRangeUser(0, maximum+maximum*0.1)

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.55, xmax=0.8, ymax=0.75, size=0.035)
        leg.AddEntry(hist_tot, '{}'.format(quantity))
        leg.AddEntry(hist_scaled, '{} x {}'.format(quantity, factor))
        leg.AddEntry(profile_data, 'PU profile in data')

        canv_name = 'canv'
        canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
        canv.cd()

        profile_data.Draw('hist')
        hist_tot.Draw('hist same')
        hist_scaled.Draw('hist same')

        leg.Draw()

        figname = 'comparison_pileup_2018_{}'.format(quantity)
        canv.SaveAs('myPlots/pileup/{}.png'.format(figname))
        canv.SaveAs('myPlots/pileup/{}.pdf'.format(figname))

        


    def plot_profile_per_part(self):
        
        profile_filenames = [f for f in glob('myPlots/pileup/profile_{}_*.root'.format(self.year))]

        canv_name = 'canv'
        canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
        canv.cd()

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.45, xmax=0.8, ymax=0.75, size=0.035)

        hists = []
        int_tot = 0
        for ipart, part in enumerate(self.parts):
            print('\n' + part)

            bin_min = 0
            bin_max = 200
            nbins = 200
            hist_name = 'hist_tot_part_{}'.format(part)
            print(hist_name)
            hist_tot_part = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
            hist_tot_part.SetLineColor(self.colours[ipart])
            hist_tot_part.SetLineWidth(2)

            leg.AddEntry(hist_tot_part, part)

            for profile_filename in profile_filenames:
                if part not in profile_filename: continue

                label = self.get_label(profile_filename)
                f = ROOT.TFile.Open(profile_filename, 'READ')
                profile = f.Get('hist_{}'.format(label))

                int_tot += profile.Integral()

                hist_tot_part.Add(profile)
                f.Close()

            hists.append(hist_tot_part)


        maximum = max([hist.GetMaximum() for hist in hists])
        for ihist, hist in enumerate(hists):
            #print(hist.Integral()/float(int_tot))
            #hist.Scale(hist.Integral()/float(int_tot))
            ##hist.Scale(1./hist.Integral())
            #hist.Scale((hist.Integral()/float(int_tot))/hist.Integral())
            #print((hist.Integral()/float(int_tot))/hist.Integral())

            #print( 1./hist.Integral() * hist.Integral()/float(int_tot) )

            if ihist == 0:
                hist.GetXaxis().SetRangeUser(0, 90)
                hist.GetYaxis().SetRangeUser(0, maximum+maximum*0.1)
                hist.SetTitle(year)
                hist.Draw('hist')
            else:
                hist.Draw('hist same')

            leg.Draw()


        figname = 'pileup_{}'.format(self.year)
        canv.SaveAs('myPlots/pileup/{}.png'.format(figname))
        canv.SaveAs('myPlots/pileup/{}.pdf'.format(figname))


    def plot_profile_per_year(self):
        quantity = 'PV_npvs'
        
        profile_filenames = [f for f in glob('myPlots/pileup/profile_{}_*.root'.format(self.year))]

        canv_name = 'canv'
        canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
        canv.cd()

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.45, xmax=0.8, ymax=0.75, size=0.035)

        bin_min = 0
        bin_max = 200
        nbins = 200
        hist_name = 'hist_tot_{}'.format(self.year)
        hist_tot = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
        hist_tot.SetLineColor(1)
        hist_tot.SetLineWidth(2)

        #hists = []
        for profile_filename in profile_filenames:
            print(profile_filename)

            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')

            hist_tot.Add(profile)
            f.Close()

        factor = 1.2
        hist_scaled = ROOT.TH1D('hist_scaled', 'hist_scaled', nbins, bin_min*factor, bin_max*factor)
        
        for i in range(1, nbins+1):
          hist_scaled.SetBinContent(i, hist_tot.GetBinContent(i))
        
        hist_scaled.SetLineColor(2)
        hist_scaled.SetLineWidth(2)

        # normalise distributions
        hist_scaled.Scale(1./hist_scaled.Integral())
        hist_tot.Scale(1./hist_tot.Integral())

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.55, xmax=0.8, ymax=0.75, size=0.035)
        leg.AddEntry(hist_tot, '{}'.format(quantity))
        leg.AddEntry(hist_scaled, '{} x {}'.format(quantity, factor))

        hist_tot.GetXaxis().SetRangeUser(0, 100)
        hist_tot.SetTitle(self.year)

        hist_tot.Draw('hist')
        hist_scaled.Draw('hist same')

        leg.Draw()

        figname = '{}_{}'.format(quantity, self.year)
        canv.SaveAs('myPlots/pileup/{}.png'.format(figname))
        canv.SaveAs('myPlots/pileup/{}.pdf'.format(figname))


    def plot_profile_all_years(self):
        quantity = 'PV_npvs'
        
        profile_filenames_2018 = [f for f in glob('myPlots/pileup/profile_2018_*.root')]
        profile_filenames_2022 = [f for f in glob('myPlots/pileup/profile_2022_*.root')]
        profile_filenames_2024 = [f for f in glob('myPlots/pileup/profile_2024_*.root')]

        canv_name = 'canv'
        canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
        canv.cd()

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.45, xmax=0.8, ymax=0.75, size=0.035)

        bin_min = 0
        bin_max = 200
        nbins = 200

        hist_name_2018 = 'hist_tot_2018'
        hist_tot_2018 = ROOT.TH1D(hist_name_2018, hist_name_2018, nbins, bin_min, bin_max)

        hist_name_2022 = 'hist_tot_2022'
        hist_tot_2022 = ROOT.TH1D(hist_name_2022, hist_name_2022, nbins, bin_min, bin_max)

        hist_name_2024 = 'hist_tot_2024'
        hist_tot_2024 = ROOT.TH1D(hist_name_2024, hist_name_2024, nbins, bin_min, bin_max)

        for profile_filename in profile_filenames_2018:
            print(profile_filename)
            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')
            hist_tot_2018.Add(profile)
            f.Close()

        for profile_filename in profile_filenames_2022:
            print(profile_filename)
            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')
            hist_tot_2022.Add(profile)
            f.Close()

        for profile_filename in profile_filenames_2024:
            print(profile_filename)
            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')
            hist_tot_2024.Add(profile)
            f.Close()

        factor = 1.2

        hist_scaled_2018 = ROOT.TH1D('hist_scaled_2018', 'hist_scaled_2018', nbins, bin_min*factor, bin_max*factor)
        for i in range(1, nbins+1):
          hist_scaled_2018.SetBinContent(i, hist_tot_2018.GetBinContent(i))
        hist_scaled_2018.SetLineColor(1)
        hist_scaled_2018.SetLineWidth(2)

        hist_scaled_2022 = ROOT.TH1D('hist_scaled_2022', 'hist_scaled_2022', nbins, bin_min*factor, bin_max*factor)
        for i in range(1, nbins+1):
          hist_scaled_2022.SetBinContent(i, hist_tot_2022.GetBinContent(i))
        hist_scaled_2022.SetLineColor(2)
        hist_scaled_2022.SetLineWidth(2)

        hist_scaled_2024 = ROOT.TH1D('hist_scaled_2024', 'hist_scaled_2024', nbins, bin_min*factor, bin_max*factor)
        for i in range(1, nbins+1):
          hist_scaled_2024.SetBinContent(i, hist_tot_2024.GetBinContent(i))
        hist_scaled_2024.SetLineColor(4)
        hist_scaled_2024.SetLineWidth(2)

        # normalise distributions
        hist_scaled_2018.Scale(1./hist_scaled_2018.Integral())
        hist_scaled_2022.Scale(1./hist_scaled_2022.Integral())
        hist_scaled_2024.Scale(1./hist_scaled_2024.Integral())

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.55, xmax=0.8, ymax=0.75, size=0.035)
        leg.AddEntry(hist_scaled_2018, '2018')
        leg.AddEntry(hist_scaled_2022, '2022')
        leg.AddEntry(hist_scaled_2024, '2024')

        maximum = max([hist.GetMaximum() for hist in [hist_scaled_2018, hist_scaled_2022, hist_scaled_2024]])
        hist_scaled_2018.GetXaxis().SetRangeUser(0, 100)
        hist_scaled_2018.GetYaxis().SetRangeUser(0, maximum+maximum*0.1)
        hist_scaled_2018.GetXaxis().SetTitle('PV_npvs (x {})'.format(factor))
        hist_scaled_2018.SetTitle(' ')

        hist_scaled_2018.Draw('hist')
        hist_scaled_2022.Draw('hist same')
        hist_scaled_2024.Draw('hist same')

        leg.Draw()

        figname = '{}_all_years'.format(quantity)
        canv.SaveAs('myPlots/pileup/{}.png'.format(figname))
        canv.SaveAs('myPlots/pileup/{}.pdf'.format(figname))


    def compare_central_to_custom(self):
        quantity = 'PV_npvs'
        
        profile_filenames = [f for f in glob('myPlots/pileup/profile_{}_*.root'.format(self.year))]

        canv_name = 'canv'
        canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
        canv.cd()

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.45, xmax=0.8, ymax=0.75, size=0.035)

        bin_min = 0
        bin_max = 200
        nbins = 200

        hist_name_central = 'hist_central_{}'.format(self.year)
        hist_central = ROOT.TH1D(hist_name_central, hist_name_central, nbins, bin_min, bin_max)
        hist_central.SetLineColor(2)
        hist_central.SetLineWidth(2)
        
        if self.year == '2022':
            values_central = [7.075550618391933e-8, 1.8432226484975646e-7, 4.6156514471969593e-7, 0.0000011111611991838491, 
                 0.0000025719752161798103, 0.000005724865812608344, 0.000012255841383374045, 0.000025239403069596116, 
                 0.00005001054998201597, 0.00009536530158990567, 0.00017505633393457624, 0.00030942214916825035, 
                 0.0005268123536229287, 0.0008642843968521786, 0.0013669182280399903, 0.0020851167548246985, 
                 0.0030695148409245446, 0.004363635945105083, 0.005995143197404548, 0.007967247822222358, 
                 0.010252302872826594, 0.01278957659177177, 0.015488544412469806, 0.01823784978331645, 
                 0.020918669702105028, 0.023420019399650906, 0.025652949149203495, 0.027560835627835043, 
                 0.02912397347687914, 0.030358091266301533, 0.03130778480604892, 0.03203676872496023, 
                 0.0326170853351521, 0.03311902652393314, 0.033602777248239, 0.0341120235754556, 
                 0.03466927947785801, 0.03527261707506484, 0.035893786618889145, 0.03647817900850185, 
                 0.036947435730750315, 0.03720550450678737, 0.037148460727673235, 0.03667753703450604, 
                 0.03571377296329832, 0.034211859754226276, 0.032170439241889726, 0.029636506070368274, 
                 0.02670262519076345, 0.023497154911314072, 0.020169158697337236, 0.016870783471647905, 
                 0.013740289679427057, 0.010888563843704815, 0.008390977574442656, 0.006285186751143873, 
                 0.004574246293656772, 0.003233538335807419, 0.002219622271900557, 0.0014792038980537092, 
                 0.0009568560481315006, 0.0006007171037926386, 0.00036596934105178995, 0.0002163349104153549, 
                 0.00012407362512604619, 0.0000690356949524181, 0.000037263645547231494, 0.00001951170588910065, 
                 0.000009910336118978026, 0.0000048826244075428666, 0.0000023333596885075797, 
                 0.0000010816029570543702, 4.863048449289416e-7, 2.1208148308081624e-7, 8.97121135679932e-8, 
                 3.6809172420519874e-8, 1.4649459937201982e-8, 5.655267024863598e-9, 2.117664468591336e-9, 
                 7.692038404370259e-10, 2.7102837405697987e-10, 9.263749466613295e-11, 3.071624552355945e-11, 
                 9.880298997379985e-12, 3.0832214331312204e-12, 9.33436314183754e-13, 2.7417209623761203e-13, 
                 7.813293248960901e-14, 2.1603865264197903e-14, 5.796018523167997e-15, 1.5088422256459697e-15, 
                 3.811436255838504e-16, 9.342850737730402e-17, 2.2224464483477953e-17, 5.130498608124184e-18, 
                 1.1494216669980747e-18, 2.499227229379666e-19, 5.2741621866055994e-20, 1.080281961755894e-20, 
                 2.1476863811171814e-21 ]
            
        elif self.year == '2024':
            values_central = [1.0599126204703409e-05, 3.951260716440464e-05, 5.582896642311492e-05, 6.453438928373494e-05, 7.94442211109479e-05,
                7.937281404906837e-05, 8.907245068638336e-05, 9.415642420373598e-05, 9.502623228017968e-05, 9.199478729752352e-05,
                0.00013938048564508576, 0.00017010543636561607, 0.0001134056824651737, 0.00014296639146898143, 0.0002581456179543262,
                0.0005678611952417843, 0.0008332593678487858, 0.0009099520114358141, 0.0009172050801512289, 0.0009505506810280174,
                0.0010644229070947433, 0.001339110596250398, 0.0018834003581107888, 0.0029865186505389474, 0.004897011268804792,
                0.007488712649153006, 0.010316155472531027, 0.01280202341832599, 0.014747096050878332, 0.016326943326136174,
                0.017545204330173675, 0.018522496650041513, 0.019313535713508023, 0.01993990862588472, 0.020383348878016044,
                0.020721612374020128, 0.02102766065973508, 0.02150180376692927, 0.022073985047278605, 0.0226770721698269,
                0.023284487595937403, 0.02396015609623849, 0.024817192902298803, 0.026070521311442062, 0.02789493736924733,
                0.03041415774056745, 0.03431433930850965, 0.03926963065106093, 0.04413912495173111, 0.0482701694536395,
                0.051618664470990054, 0.05375451405860295, 0.052590289193478454, 0.049365338159627375, 0.04486833402546383,
                0.03819745898668014, 0.030826335590683835, 0.023458382461983215, 0.017085594960378436, 0.01195868965716555,
                0.0078471509432766, 0.005033866678064158, 0.0031381464673017096, 0.001932392575354725, 0.0012357872700131897,
                0.0007126320704765079, 0.0003612978990027777, 0.000182623168268646, 8.818852384472863e-05, 3.151355162852159e-05,
                1.0187289280910454e-05, 4.7968992292803826e-06, 1.8577163576511699e-06, 5.729305524684033e-07, 2.602401818463374e-07,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 9.949892938900103e-09,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0]


        for x in list(range(100)):
            hist_central.SetBinContent(x, values_central[x])


        hist_name = 'hist_tot_{}'.format(self.year)
        hist_tot = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)

        for profile_filename in profile_filenames:
            print(profile_filename)

            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')

            hist_tot.Add(profile)
            f.Close()

        factor = 1.2
        hist_scaled = ROOT.TH1D('hist_scaled', 'hist_scaled', nbins, bin_min*factor, bin_max*factor)
        
        for i in range(1, nbins+1):
          hist_scaled.SetBinContent(i, hist_tot.GetBinContent(i))
        
        hist_scaled.SetLineColor(1)
        hist_scaled.SetLineWidth(2)

        # normalise distribution
        hist_scaled.Scale(1./hist_scaled.Integral())

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.55, xmax=0.8, ymax=0.75, size=0.035)
        leg.AddEntry(hist_central, 'central profile')
        leg.AddEntry(hist_scaled, 'custom profile')

        hist_scaled.GetXaxis().SetRangeUser(0, 100)
        maximum = max([hist.GetMaximum() for hist in [hist_central, hist_scaled]])
        hist_scaled.GetYaxis().SetRangeUser(0, maximum + maximum*0.1)
        hist_scaled.SetTitle(self.year)

        hist_scaled.Draw('hist')
        hist_central.Draw('hist same')

        leg.Draw()

        figname = 'central_vs_custom_{}'.format(self.year)
        canv.SaveAs('myPlots/pileup/{}.png'.format(figname))
        canv.SaveAs('myPlots/pileup/{}.pdf'.format(figname))


    def compare_pv_mc_to_data(self):
        quantity = 'PV_npvs'
        factor = 1.2

        bin_min = 0
        bin_max = 200
        nbins = 200

        canv_name = 'canv'
        canv = self.tools.createTCanvas(name=canv_name, dimx=800, dimy=700)
        canv.cd()

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.45, xmax=0.8, ymax=0.75, size=0.035)

        ## pv profile in data
        profile_filenames = [f for f in glob('myPlots/pileup/profile_{}_*.root'.format(self.year))]
        hist_name = 'hist_tot_{}'.format(self.year)
        hist_tot = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)

        for profile_filename in profile_filenames:
            print(profile_filename)

            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')

            print(profile.GetEntries())
            hist_tot.Add(profile)
            f.Close()

        hist_scaled = ROOT.TH1D('hist_scaled', 'hist_scaled', nbins, bin_min*factor, bin_max*factor)
        
        for i in range(1, nbins+1):
          hist_scaled.SetBinContent(i, hist_tot.GetBinContent(i))
        
        hist_scaled.SetLineColor(1)
        hist_scaled.SetLineWidth(2)

        # pv profile in simulation
        mc_file = '/eos/cms/store/group/phys_bphys/anlyon/CPVGen/sin2b/2024/bparknano_tot.root'
        f = ROOT.TFile.Open(mc_file, 'READ')
        tree = f.Get('Events')

        hist_name = 'hist_mc'
        hist_mc = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
        
        tree.Draw('{} >> hist_mc'.format(quantity))

        hist_scaled_mc = ROOT.TH1D('hist_scaled_mc', 'hist_scaled_mc', nbins, bin_min*factor, bin_max*factor)
        
        for i in range(1, nbins+1):
          hist_scaled_mc.SetBinContent(i, hist_mc.GetBinContent(i))
        
        hist_scaled_mc.SetLineColor(4)
        hist_scaled_mc.SetLineWidth(2)


        # normalise distributions
        hist_scaled_mc.Scale(1./hist_scaled_mc.Integral())
        hist_scaled.Scale(1./hist_scaled.Integral())

        leg = self.tools.getRootTLegend(xmin=0.6, ymin=0.55, xmax=0.8, ymax=0.75, size=0.035)
        leg.AddEntry(hist_scaled_mc, 'MC')
        leg.AddEntry(hist_scaled, 'data')

        hist_scaled.GetXaxis().SetRangeUser(0, 100)
        hist_scaled.GetXaxis().SetTitle('number of PV (scaled by factor of 1.2)')
        maximum = max([hist.GetMaximum() for hist in [hist_scaled_mc, hist_scaled]])
        hist_scaled.GetYaxis().SetRangeUser(0, maximum + maximum*0.1)
        hist_scaled.SetTitle(self.year)

        hist_scaled.Draw('hist')
        hist_scaled_mc.Draw('hist same')

        leg.Draw()

        figname = 'pv_data_vs_mc_{}'.format(self.year)
        canv.SaveAs('myPlots/pileup/{}.png'.format(figname))
        canv.SaveAs('myPlots/pileup/{}.pdf'.format(figname))



    def create_profile_file(self):
        quantity = 'PV_npvs'
        
        profile_filenames = [f for f in glob('myPlots/pileup/profile_{}_*.root'.format(self.year))]


        bin_min = 0
        bin_max = 200
        nbins = 200
        hist_name = 'hist_tot_{}'.format(self.year)
        hist_tot = ROOT.TH1D(hist_name, hist_name, nbins, bin_min, bin_max)
        hist_tot.SetLineColor(1)
        hist_tot.SetLineWidth(2)

        for profile_filename in profile_filenames:
            print(profile_filename)

            f = ROOT.TFile.Open(profile_filename, 'READ')
            if quantity == 'PV_npvs':
                profile = f.Get('hist')
            elif quantity == 'PV_npvsGood':
                profile = f.Get('hist_good')

            hist_tot.Add(profile)
            f.Close()

        outfile = ROOT.TFile.Open('myPlots/pileup/pileup_{}.root'.format(self.year), 'RECREATE')
        factor = 1.2
        pileup = ROOT.TH1D('pileup', 'pileup', nbins, bin_min*factor, bin_max*factor)
        
        for i in range(1, nbins+1):
          pileup.SetBinContent(i, hist_tot.GetBinContent(i))
        
        pileup.SetLineColor(1)
        pileup.SetLineWidth(2)
        pileup.GetXaxis().SetRangeUser(0, 200)

        # normalise distributions
        pileup.Scale(1./pileup.Integral())

        outfile.Write()
        outfile.Close()
        print('--> myPlots/pileup/pileup_{}.root created'.format(self.year))


    def create_string(self):
        infile = ROOT.TFile.Open('myPlots/pileup/pileup_{}.root'.format(self.year), 'READ')
        hist = infile.Get('pileup')

        bins = []
        probValues = []

        nbins = hist.GetNbinsX()
        for i in range(0, nbins):
            bins.append(i)
            probValues.append(hist.GetBinContent(i))

        print('\n bins:')
        print(bins)
        print('\n values:')
        print(probValues)
        print(max(probValues))



if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)

    year = '2024'
    plotter = Plotter(year = year)
    #plotter.process()
    #plotter.plot_2018()
    #plotter.plot_profile_per_year()
    #plotter.plot_profile_all_years()
    #plotter.create_profile_file()
    #plotter.create_string()
    #plotter.compare_central_to_custom()
    plotter.compare_pv_mc_to_data()
