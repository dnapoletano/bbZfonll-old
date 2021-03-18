import sys
sys.path.append('../../lib')

import _bbHFONLL as fb
import numpy as np
import pandas as pd
# import reader

class A0_coeff :
    def __init__(self,scale,pdf):
        self.scale=scale
        self.pdf=pdf

    def sigma_m_lim(self,order,muR,muF,mH):
        fb.Scales_i_fonll_set(self.scale,order)
        fb.Scales_i_or_set(self.scale,-1)
        if(muF <=  fb.Scales_mbmb_get(self.scale)):
            return 0.
        else:
            fb.Scales_m_H_set(self.scale,mH)
            fb.Scales_mu_R_set(self.scale,muR)
            fb.Scales_mu_F_set(self.scale,muF)
            return -fb.xsecDiff(self.pdf,self.scale)
