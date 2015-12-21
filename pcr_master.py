#!/usr/bin/python

"""
tested with Python 2.7
dependencies: numpy, pandas, matplotlib, scipy, statsmodules
sample group information in template/mapping1
"""
import os
path1 = '/home/joe/p/template'
os.chdir(path1)
import pandas as pd
import mapping1, pcr1, pcr2, pwttest1, anova1, barplot1, boxplot1, barplot2

class RtPcr():
    """

    """
    def __init__(self, endo='det_ctrl', cali='treat1 condi1', mapping=mapping1.mapping_1):
        self.endo = endo
        self.cali = cali
        self.mapping = mapping
        
    def main(self, file1='pcrraw.csv', path=''):
        
        p1 = pcr1.Pcr1(self.endo, self.cali, self.mapping)
        p2 = pcr2.Pcr2(self.cali)
        gg = p2.pcr_2(p1.pcr_1(file1)[0])[0]

        genes = pd.Series([group1 for group1, group2, group3 in gg.index]).unique()
        genes.sort()

        groups = list(pd.Series([group2 for group1, group2, group3 in gg.index]).unique())
        groups = [groups[i] for i in [1,3,0,2]]


        ## student t-test
        t1 = pwttest1.Pwttest1(['condi1', 'condi2'], ['treat1', 'treat2']).pwttest_1(gg,cat=0)

        # two-way ANOVA
        a1 = anova1.Anova1().anova_1(gg)

        # figure 1
        f1 = barplot1.Barplot1(genes, groups, path).barplot_1(gg)

        # figure 2
        f2 = boxplot1.Boxplot1(genes, groups, path).boxplot_1(gg)

        # figure 3
        f3 = barplot2.Barplot2(genes, groups, path).barplot_2(gg)

        return gg, t1, a1, f1, f2, f3


RtPcr().main()
        
if __name__ == "__main__":
    try:
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("endo", help="endogeneous control")
        parser.add_argument("cali", help="group name of calibrator samples")
        parser.add_argument("file1", help="path and file name of raw data")
        args = parser.parse_args()
        RtPcr(endo=args.endo, cali=args.cali).main(file1=args.file1)
    except IOError:
        import sys
        RtPcr(endo=sys.argv[1], cali=sys.argv[2]).main(file1=sys.argv[3])
    except:
        RtPcr(endo='det_ctrl', cali='treat1 condi1').main(file1='pcrraw.csv')
    
