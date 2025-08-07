import os
import ROOT

#TODO test sv_prob (like sin2b)

infilename = '/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/nanoFiles/merged/bparknano_Bs_selected.root'

f = ROOT.TFile.Open(infilename, 'READ')
tree = f.Get('Events')

count_tot = 0
count_tot_multcand = 0
count_matched = 0
count_matched_multcand = 0

for ientry, entry in enumerate(tree):
   if ientry%1000 == 0:
     print '{}% events processed'.format(round(float(ientry)/float(tree.GetEntries())*100, 1))

   ncand = entry.nBsToPhiPhiTo4K

   has_at_least_one_matched = False
   for icand in range(0, ncand):
      if entry.BsToPhiPhiTo4K_isMatched[icand]: has_at_least_one_matched = True

   if not has_at_least_one_matched: continue

   #NOTE reverse=True: decreasing order

   #sorted_idx = sorted([ii for ii in range(0, ncand)], key = lambda x : entry.BsToPhiPhiTo4K_Bs_pt[x], reverse=True)
   #sorted_idx = sorted([ii for ii in range(0, ncand)], key = lambda x : entry.BsToPhiPhiTo4K_Bs_cos2D[x], reverse=True)
   sorted_idx = sorted([ii for ii in range(0, ncand)], key = lambda x : entry.BsToPhiPhiTo4K_phi1_pt[x] * entry.BsToPhiPhiTo4K_phi2_pt[x], reverse=True)

   cand_idx = sorted_idx[0]

   count_tot += 1 
   if entry.BsToPhiPhiTo4K_isMatched[cand_idx]: count_matched += 1

   if ncand < 2: continue
   count_tot_multcand += 1 
   if entry.BsToPhiPhiTo4K_isMatched[cand_idx]: count_matched_multcand += 1

   #print '\n'
   #for icand in range(0, ncand):
   #   print 'pt: {}'.format(entry.BsToPhiPhiTo4K_Bs_pt[icand])

   #print 'cand_idx: {}'.format(cand_idx)

   #for icand in range(0, ncand):
   #   print 'isMatched: {}'.format(entry.BsToPhiPhiTo4K_isMatched[icand])



print 'count_matched / count_tot = {} / {} = {}'.format(count_matched, count_tot, count_matched/float(count_tot))
print 'count_matched_multcand / count_tot_multcand = {} / {} = {}'.format(count_matched_multcand, count_tot_multcand, count_matched_multcand/float(count_tot_multcand))


