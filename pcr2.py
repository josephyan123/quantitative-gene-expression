import numpy as np
class Pcr2():

    """
    remove outliers automatically
    """

    def __init__(self, sam_ctrl):
        self.sam_ctrl = sam_ctrl

    def pcr_2(self, df):
        def demean(arr):
            return arr - arr.mean()
        grpd = df.groupby(level=['Detector', 'Grp'])['Rel_exp']
        dm1 = grpd.transform(demean)
        std1 = grpd.transform(np.std)
        use = df[abs(dm1)<=2*std1]
        use = use/use['Rel_exp'].xs(self.sam_ctrl, level='Grp').mean()
        discard = df[abs(dm1)>2*std1]
        discard = discard/use['Rel_exp'].xs(self.sam_ctrl, level='Grp').mean()
        return use, discard
