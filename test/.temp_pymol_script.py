

import sys
sys.path.append('..')
import src.bioiain as bi


pdb = cmd.load('/cri4/iain/scripts/bioiain/test/data/test_list/8A9G.pdb', 'original',)
print(pdb,)
cmd.disable('(all)',)
bi.log('header', "I'm a log",)
