import ROOT

class Quantity(object):
  def __init__(self, name_nano='', name_flat='', title='', label='', nbins=0, bin_min=0., bin_max=0., do_log=False):
    self.name_nano = name_nano
    self.name_flat = name_flat
    self.title = title
    self.label = label
    self.nbins = nbins
    self.bin_min = bin_min
    self.bin_max = bin_max
    self.do_log = do_log

quantities = {}
quantities['track_study'] = [
  Quantity(name_nano='PhiToKK_phi_k1_drTrg[BsToPhiPhiTo4K_phi1_idx]', label='k1_drtrg', title='#DeltaR(k_{1}, trg#mu) ', nbins=80, bin_min=-0, bin_max=5),
  Quantity(name_nano='PhiToKK_phi_k2_drTrg[BsToPhiPhiTo4K_phi1_idx]', label='k2_drtrg', title='#DeltaR(k_{2}, trg#mu) ', nbins=80, bin_min=-0, bin_max=5),
  Quantity(name_nano='PhiToKK_phi_k1_drTrg[BsToPhiPhiTo4K_phi2_idx]', label='k3_drtrg', title='#DeltaR(k_{3}, trg#mu) ', nbins=80, bin_min=-0, bin_max=5),
  Quantity(name_nano='PhiToKK_phi_k2_drTrg[BsToPhiPhiTo4K_phi2_idx]', label='k4_drtrg', title='#DeltaR(k_{4}, trg#mu) ', nbins=80, bin_min=-0, bin_max=5),
  Quantity(name_nano='PhiToKK_phi_k1_dzTrg[BsToPhiPhiTo4K_phi1_idx]', label='k1_dztrg', title='#Deltaz(k_{1}, trg#mu) ', nbins=80, bin_min=-10, bin_max=10, do_log=True),
  Quantity(name_nano='PhiToKK_phi_k2_dzTrg[BsToPhiPhiTo4K_phi1_idx]', label='k2_dztrg', title='#Deltaz(k_{2}, trg#mu) ', nbins=80, bin_min=-10, bin_max=10, do_log=True),
  Quantity(name_nano='PhiToKK_phi_k1_dzTrg[BsToPhiPhiTo4K_phi2_idx]', label='k3_dztrg', title='#Deltaz(k_{3}, trg#mu) ', nbins=80, bin_min=-10, bin_max=10, do_log=True),
  Quantity(name_nano='PhiToKK_phi_k2_dzTrg[BsToPhiPhiTo4K_phi2_idx]', label='k4_dztrg', title='#Deltaz(k_{4}, trg#mu) ', nbins=80, bin_min=-10, bin_max=10, do_log=True),
 ]

