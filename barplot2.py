import numpy as np
import matplotlib.pyplot as plt

class Barplot2():
    """

    """
    def __init__(self, genes, groups, path):
        self.genes = genes
        self.groups = groups
        self.path = path

    def barplot_2(self, df):
	hh = df.groupby(level=['Detector', 'Grp'])['Rel_exp'].agg(['mean', 'std', 'count'])
        def barfigfmt(x, col):
            x1 = x[col]
            x1_1 = x1.unstack()
            try:
                del x.columns.name
            except:
                pass
            x1_1 = x1_1[self.groups]
            x1_1 = x1_1.reindex(self.genes)
            return x1_1

        hh1_1 = barfigfmt(hh, 'mean')
        hh1_1e = barfigfmt(hh, 'std')/np.sqrt(barfigfmt(hh, 'count'))

        # plot the bars with data from hh1_1 and error bar from hh1_1e
        #plt.rcParams.update({'font.size': 10})
        fig = plt.figure(); ax = fig.add_subplot(1, 1, 1)
        hh1_1.plot(kind='bar', ax=ax, yerr=hh1_1e, ylim=[0, hh1_1.max().max()*1.6], grid=False)
        ax.set_xticklabels(hh1_1.index, rotation=0)
        ax.set_ylabel('RNA Level')
        ax.set_xlabel('')
        ax.legend(loc='best', frameon=False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        # label the bars with values
        for i, each in enumerate(hh1_1.index):
            for j, col in enumerate(hh1_1.columns):
                y = hh1_1.ix[each][col]
                ye = hh1_1e.ix[each][col]
                ax.text(i+(0.12*j-0.18), y+ye+0.02, "{0:.2g}".format(y), ha='center', va='bottom')

            
        # # label p-values from pairwise student t-test
        # xmin, xmax = ax.get_xlim(); xlen = xmax - xmin
        # dfy = hh1_1 + hh1_1e
        # #psign1 = psign[['condi1', 'treat1', 'condi2', 'treat2']]
        # colo = ['red', 'blue', 'red', 'blue']
        # for i in range(len(hh1_1)):
        #     for h, j in enumerate([1, 2, 3, 3]):
        #         ax.text((-0.2+0.135*j+i)/len(self.genes)*xlen,
        #                 dfy.iloc[i,j]+(0.2 if h==3 else 0.1),
        #                 psign1.ix[i,h],
        #                 fontsize=8, color=colo[h], weight='bold', ha='center', va='bottom')

        # add notes:
        ax.text(-0.5, -0.4, '*, **, ***, ****, indicate p-value < 0.05, 0.01, 0.001, or 0.0001 respectively.',
                fontsize=10, ha='left', va='bottom')
        ax.text(-0.5, -0.5, 'red color letters show knockout vs. floxed; ',
                fontsize=8, color='red', ha='left', va='bottom')
        ax.text(0.28, -0.5, 'blue color letters show high-fat diet vs. low-fat diet.',
                fontsize=8, color='blue', ha='left', va='bottom')

        # label with p-values from ANOVA analysis
        # indf = 0
        # for name, group in anova2.groupby(level='Detector'):
        #     group.columns=['p-value']
        #     height = (hh1_1.max(axis=1)[indf])* 1.6
        #     for m in range(len(group)):
        #         ax.text(indf, height-0.1*m, group.index.get_level_values(1)[m] + str('%f' % (group['p-value'][m]))[:5],
        #                 ha='center', va='bottom', fontsize=6, fontweight='bold')
        #     indf += 1

        # save the figure
        plt.savefig(self.path+'figure.png', dpi=400, bbox_inches='tight')
        plt.close()
