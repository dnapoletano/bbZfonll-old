import sys
sys.path.append('./../lib')
sys.path.append('../src/python')

import fonll

fnll=fonll.fonll('test.info')

fnll.calc_fonll(125.,125.,125.)
