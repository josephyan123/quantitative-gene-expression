import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Barplot1():
    """

    """
    def __init__(self, genes, groups, path):
        self.genes = genes
        self.groups = groups
        self.path = path

    def barplot_1(self, df):
        
        def detgrp(x):
            lst1 = [range(value+1)[1:] for value in x['Rel_exp'].count(level='Grp')]
            lst2 = [item for sublist in lst1 for item in sublist]
            x['repeats'] = lst2
            x.index = x.index.droplevel(level='Sample')
            x = x.set_index('repeats', append=True)
            x = x.unstack()
            return x

        newformat = df.groupby(level='Detector').apply(detgrp)
        newformat.index = newformat.index.droplevel(0)
        newformat = newformat.sortlevel(ascending=False)

        # draw figure 1
        fig, axes = plt.subplots(len(self.genes))
        plt.subplots_adjust(wspace=0.4, hspace=0.3)
        inde = 0
        for name, group in newformat.groupby(level='Detector'):
            mi = pd.MultiIndex(levels=[list(self.genes), self.groups], labels=[[inde]*len(self.groups), range(len(self.groups))], names=[u'Detector', u'Grp'])
            group = group.reindex(mi)
            yerr=np.array([group['Neg_err'].values, group['Pos_err'].values]).swapaxes(0,1).swapaxes(0,2)
            group['Rel_exp'].plot(kind='bar', ax=axes[inde], color = 'b', legend=False, title=name, yerr=yerr, grid=False)
            axes[inde].set_xticklabels(self.groups, rotation=0)
            axes[inde].set_ylabel('RNA Level')
            axes[inde].set_xlabel('')
            inde += 1

        plt.savefig(self.path+'bar_chart.png', dpi=400, bbox_inches='tight')  ### change file name
        plt.close()
