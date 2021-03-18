import sys
sys.path.append('../../lib')

import _bbHFONLL as fb
import numpy as np
import pandas as pd
import messages
import reader
import f_4f
import f_5f
import massless_lim

class fonll :
    def __init__(self,infile='../../examples/test.info'):
        self.msg=messages.messages()
        self.msg.print_intro()
        self.reader=reader.reader(infile)
        
        pdf_set=self.reader.get_str_variable('pdf_set','NNPDF30_nnlo_as_0118')
        pdf_mem=self.reader.get_i_variable('pdf_mem',0)
        order_5f=self.reader.get_i_variable('order_5f',2)
        order_fonll=self.reader.get_i_variable('order_fonll',2)
        mb=self.reader.get_d_variable('mb',4.58)
        mur=self.reader.get_d_variable('mur',125.)
        muf=self.reader.get_d_variable('muf',125.)
        mh=self.reader.get_d_variable('mh',125.)
        sqrts=self.reader.get_d_variable('sqrts',13000)
        iter=self.reader.get_i_variable('iter',100000)
        a_min=self.reader.get_d_variable('a_min',0.7e-2)
        proc=self.reader.get_str_variable('proc','bbz') 
        self.scale=fb.new_Scales(2,2,mb,mur,muf,mh,sqrts,iter,pdf_mem,a_min,proc)
        self.calc=self.reader.get_str_variable('calc_type','var_mass')              
        self.pdf=fb.new_LUMI(pdf_set, pdf_mem, muf, mb, False);
        nlo_file=self.reader.get_str_variable('4f_nlo_file','4f/4F_NLO.txt')
        lo_file=self.reader.get_str_variable('4f_lo_file','4f/4F_LO.txt')
        gg_file=self.reader.get_str_variable('4f_gg_file','4f/4F_gg.txt')
        
        self.B=f_4f.B_coeff(nlo_file,lo_file,gg_file,self.pdf,[mur,muf,mb])
        self.A=f_5f.A_coeff(self.scale,pdf_set)
        self.A0=massless_lim.A0_coeff(self.scale,pdf_set)
        
        self.mh=mh
        self.mur=mur
        self.muf=muf
        self.order_5f=order_5f
        self.order_fonll=order_fonll
        
        
    def calc_fonll(self,mur,muf,mh):
        # mur=self.mur
        # muf=self.muf
        # mh=self.mh

        fonll=0.
        if(self.calc=='fixed_mass'):
            fonll=self.B.xs4FNLO(np.square(mur),np.square(muf)) + \
                   self.A.sigma5F(self.order_5f,mur,muf,mh) - \
                   self.A0.sigma_m_lim(order_fonll,mur,muf,mh)
        elif(self.calc=='var_mass'):
            fonll=self.B.xs4FNLO_mh(np.square(mur),np.square(muf),mh) + \
                   self.A.sigma5F(self.order_5f,mur,muf,mh) - \
                   self.A0.sigma_m_lim(self.order_fonll,mur,muf,mh)
        text_res='----------------------------\n'
        text_res+='  XS FONLL :      ' + str(fonll) +'\n'
        text_res+='----------------------------\n'
        self.msg.infog(text_res)
        return fonll
            
    def set_pdf_mem(mem):
        fb.Scales_i_member_set(self.scale,mem)

    def __del__(self):
        fb.delete_LUMI(self.pdf)
        fb.delete_Scales(self.scale)
        self.msg.print_outro()
