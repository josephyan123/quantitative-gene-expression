
import pandas as pd
import numpy as np
from pandas import Series, DataFrame

class Pcr1():
    """
    file, the .csv file containing raw data
    det_ctrl, internal control, like 'cyclophilin B' or 'Cyclophilin A'
    sam_ctrl, group label of samples used as control, its average value should be 1 in the delta-delta-Ct method
    mapping, dictionary for each sameple belong to a specific group
    return two dataframe, first is the value for each sample, second is the average for each group
    """
    def __init__(self, det_ctrl, sam_ctrl, mapping):
        self.det_ctrl = det_ctrl
        self.sam_ctrl = sam_ctrl
        self.mapping = mapping

    def pcr_1(self, file):

        pcr = pd.read_csv(file)
        pcr1 = pcr[['Sample', 'Detector', 'Ct']].set_index(['Sample', 'Detector'])
        pcr1['ctmedian'] = pcr1.Ct.median(level=['Sample', 'Detector'])
        pcr1['ind1'] = pcr1.ctmedian * 0.02
        pcr1['ctdfm'] = abs(pcr1.Ct - pcr1.ctmedian)
        pcr1['ctdfmmax'] = pcr1.ctdfm.max(level=['Sample', 'Detector'])
        pcr1['ctdfmmedian'] = pcr1.ctdfm.median(level=['Sample', 'Detector'])
        pcr1['ind2'] = pcr1.ctdfmmax / pcr1.ctdfmmedian
        pcr1.Ct[(pcr1.ctdfm > pcr1.ind1) & (pcr1.ind2 > 2.5)] = np.nan
        pcr1 = pcr1[['Ct']]

        me = pcr1.mean(level=['Sample', 'Detector'])
        me.columns = ['ctmean']
        st = pcr1.std(level=['Sample', 'Detector'])
        st.columns = ['ctstd']
        pcr2 = me.join(st)
        pcr2 = pcr2.reset_index(level='Detector')
        pcr2['Grp'] = Series(self.mapping)
        pcr2 = pcr2.set_index(['Detector', 'Grp'], append=True)

        c = pcr2.xs(self.det_ctrl, level='Detector'); c.columns = ['c_mean', 'c_std']
        t = pcr2[[group2 not in [self.det_ctrl] for group1, group2, group3 in pcr2.index]]; t.columns = ['t_mean', 't_std']
        pcrm = c.join(t.reset_index(level='Detector')).set_index(['Detector'], append=True).reorder_levels(['Detector', 'Grp', 'Sample']).sortlevel()

        def f2(x):
            deltact = x.t_mean - x.c_mean
            pldstd = (x.c_std ** 2 + x.t_std ** 2) ** 0.5
            deltact_c = deltact.xs(self.sam_ctrl, level='Grp').mean()
            ddct = deltact - deltact_c
            ddct_p = ddct + pldstd
            ddct_m = ddct - pldstd
            reex_r = 2 ** (-ddct)
            reex_c = reex_r.xs(self.sam_ctrl, level='Grp').mean()
            reex = reex_r/reex_c
            reex_p = 2 ** (-ddct_m)/reex_c
            reex_m = 2 ** (-ddct_p)/reex_c
            eb_p = reex_p - reex
            eb_n = reex - reex_m
            pcrind = DataFrame({'Rel_exp': reex, 'Pos_err': eb_p, 'Neg_err': eb_n}, columns=['Rel_exp', 'Pos_err', 'Neg_err'])
            return pcrind

        pcri = pcrm.groupby(level='Detector').apply(f2)

        def f4(x):
            return Series([x.mean(), x.std()], index=['mean', 'std'])

        pcrr = pcri.groupby(level=['Detector', 'Grp'])['Rel_exp'].apply(f4).unstack()

        return pcri, pcrr



if __name__ == '__main__':
    import sys
    fpcr1(sys.argv[1:])
