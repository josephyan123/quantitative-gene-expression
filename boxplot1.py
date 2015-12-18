import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Boxplot1():
    """
    draw a boxplot with whiskers
    """
    def __init__(self, genes, groups, path):
        self.genes = genes
        self.groups = groups
        self.path = path

    def boxplot_1(self, df):
        
        gg0 = df.reset_index('Sample')

        def arbsort(x):
            df=pd.DataFrame()
            for grp_nam in self.groups:
                df = df.append(x.xs(grp_nam, level='Grp', drop_level=False))
            return df
        gg0 = gg0.groupby(level='Detector', group_keys=False).apply(arbsort)

        fig, axes = plt.subplots(len(self.genes), figsize=(8, 12))
        plt.subplots_adjust(wspace=0.4, hspace=0.3)

        ind3 = 0
        for nam_det, grp_det in gg0.groupby(level='Detector'):
            grp_det.columns = [u'Sample', u'Value', u'Poserr', u'Negerr']
            ind4 = 0
            for group2 in self.groups:
                grp_grp = grp_det.xs(group2, level='Grp')
                grp_grp['x'] = np.random.randn(len(grp_grp))/10+2*ind4
                grp_grp.index = np.arange(len(grp_grp))
                yerr2 = np.array([grp_grp.Poserr, grp_grp.Negerr])
                grp_grp.plot(kind='scatter', x='x', y='Value', yerr=yerr2, ax=axes[ind3])
                bp = grp_grp[['Value']].boxplot(sym='', widths=1.5, positions=[ind4*2.0], whis=[5,95], ax=axes[ind3])
                plt.setp(bp['boxes'], linewidth=2)
                plt.setp(bp['whiskers'], linewidth=2)
                plt.setp(bp['caps'], color='blue', linewidth=3)
                for j in np.arange(len(grp_grp)):
                    axes[ind3].text(grp_grp.x[j]+0.1, grp_grp.Value[j], grp_grp.Sample[j], verticalalignment='center')
                ind4 += 1
            axes[ind3].set_xlim([-1, len(self.groups)*2-1])
            axes[ind3].set_xticks(range(0, len(self.groups)*2, 2))
            axes[ind3].set_xticklabels(self.groups)
            axes[ind3].set_xlabel('')
            axes[ind3].set_ylabel('RNA Level')
            axes[ind3].set_title(nam_det)
            ind3 += 1

        # save figure to a .png file
        plt.savefig(self.path+'boxplot.png', dpi=400, bbox_inches='tight')
        plt.close()
