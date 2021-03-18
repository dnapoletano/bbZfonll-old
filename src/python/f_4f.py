import sys
sys.path.append('../../lib')

import _bbHFONLL as fb
import numpy as np
import pandas as pd
# import reader

class B_coeff :
    def __init__(self,NLO_file,LO_file,gg_file,pdf,ref_scale):
        self.nlo=pd.DataFrame.from_csv(NLO_file, sep='\t')
        self.lo=pd.DataFrame.from_csv(LO_file, sep='\t')
        self.gg=pd.DataFrame.from_csv(gg_file, sep='\t')
        self.pdf=pdf
        self.mur=ref_scale[0]
        self.muf=ref_scale[1]
        self.m2b=np.square(ref_scale[2])
        
    def xs4FNLO(self,m2r,m2f):
        '''
        Description : Calculate the B coefficient at a given scale
        Arguments   : m2r = squared renormalisation scale
                      m2f = squared factorisation   scale
        '''
        alphar = fb.LUMI_ass(self.pdf,m2r)*2*np.pi
        af     = fb.LUMI_ass(self.pdf,m2f)
        NLO,LO,gg=self.get_4f(np.sqrt(m2r),np.sqrt(m2f))
        xs     = NLO + alphar*(-np.log(m2r/self.m2b)*LO +
                               np.log(m2f/self.m2b)*gg)/(3.*np.pi)
        return xs

    def xs4FNLO_mh(self,m2r,m2f,mh):
        '''
        Description : Calculate the B coefficient at a given scale and higgs mass
        Arguments   : m2r = squared renormalisation scale
                      m2f = squared factorisation   scale
                      mh  = mass of the higgs
        '''
        NLO,LO,gg   = self.get_4f_mh(np.sqrt(m2r),np.sqrt(m2f),mh)
        alphar = fb.LUMI_ass(self.pdf,m2r)*2*np.pi
        af     = fb.LUMI_ass(self.pdf,m2f)
        xs     = NLO + alphar*(-np.log(m2r/self.m2b)*LO +
                               np.log(m2f/self.m2b)*gg)/(3.*np.pi)
        return xs

    def get_4f(self,mur,muf,mode=0):
        MUR=mur/self.mur
        MUF=muf/self.muf
        xsnlo  = float(self.nlo[self.nlo.muR==MUR][self.nlo[self.nlo.muR==MUR].muF==MUF].XSNLO) 
        if(mode==1):
            xsnlo = xsnlo + float(self.nlo[self.nlo.muR==MUR][self.nlo[self.nlo.muR==MUR].muF==MUF].XSint)
            xslo   = float(self.lo[self.lo.muR==MUR][self.lo[self.lo.muR==MUR].muF==MUF].XSNLO)
            xsgg   = float(self.gg[self.gg.muR==MUR][self.gg[self.gg.muR==MUR].muF==MUF].XSNLO)
        return xsnlo,xslo,xsgg

    def get_4f_mh(self,mur,muf,mh,mode=0):
        MUR=mur/self.mur
        MUF=muf/self.muf
        xsnlo  = float(self.nlo[self.nlo.muR==MUR][self.nlo[self.nlo.muR==MUR].muF==MUF].XSNLO[int(mh)]) 
        if(mode==1):
            xsnlo = xsnlo + float(self.nlo[self.nlo.muR==MUR][self.nlo[self.nlo.muR==MUR].muF==MUF].XSint[int(mh)])
        xslo   = float(self.lo[self.lo.muR==MUR][self.lo[self.lo.muR==MUR].muF==MUF].XSNLO[int(mh)])
        xsgg   = float(self.gg[self.gg.muR==MUR][self.gg[self.gg.muR==MUR].muF==MUF].XSNLO[int(mh)])
        return xsnlo,xslo,xsgg
