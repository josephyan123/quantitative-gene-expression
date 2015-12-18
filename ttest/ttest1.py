import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from scipy import stats

class Ttest1():
    """
    pairwise t-test without any corrections
    first argument is the dataframe with MultiIndex and values to test
    'gbg': name of groupby column
    'cpg': name of column to contain groups to be compared with its first group value
    'value': name of column containing the values to be compared
    'ctrl': group value in cpg used as control group when doing pairwise comparison
    return 'pvalue', 'psign'
        'pvalue': dataframe containing p-value from pairwise t-test
        'psign': convert 'pvalue' dataframe to a new dataframe with '*'
            '': >= 0.05, '*': <0.05>=0.01, '**': <0.01>=0.001, '***': <0.001>=0.0001, '****': <0.0001
    """
    
    def __init__(self, df, value):
        self.df = df
        self.value = value

    def onetotwo(self, df):
        df1 = df.reset_index(level=['Grp', 'Sample'])
        df1['Grp1'] = df1['Grp'].apply(lambda x: x.split()[1])
        df1['Grp2'] = df1['Grp'].apply(lambda x: x.split()[0])
        df1 = df1.drop('Grp', axis=1)
        df1 = df1.set_index(['Grp1', 'Grp2', 'Sample'], append=True)
        return df1

    def ttest_1(self, gbg='Condi', cpg='Treat', ctrl='control'):
        df1 = self.onetotwo(self.df).reset_index([gbg, cpg])
        group1 = list(set(df1[gbg]))
        group2 = list(set(df1[cpg]))
        tt2 = np.ones((len(group1), len(group2)))
        pv2 = np.ones((len(group1), len(group2)))
        for i, g1 in enumerate(group1):
                group = df1[(df1[gbg]==g1)]
                for j, g2 in enumerate(group2):
                        tt2[i,j], pv2[i, j] = stats.ttest_ind(group[self.value][group[cpg]==ctrl], group[self.value][group[cpg]==g2])
                        
        pv3 = pd.DataFrame(pv2, index=group1, columns=group2)
        
        ps = range(len(pv3))

        for m in list(np.arange(len(pv3))):
                ps[m] = Series(pd.cut(pv3.ix[m, :], [0, 0.0001, 0.001, 0.01, 0.05, 1], right=False, labels=['****', '***', '**', '*', '']))

        psign = pd.DataFrame(ps)
        psign.index = group1
        psign.columns = group2

        pvalue = pv3.copy()

        return pvalue, psign
