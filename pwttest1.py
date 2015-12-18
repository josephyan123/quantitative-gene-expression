import pandas as pd
from ttest import ttest1

class Pwttest1():
    """
    object initiation:
        lst1 & lst2: split the 'Grp' into 2 variables, eg. 'treat1 condi1' -> 'treat1', 'condi1'
            lst1: list of unique value of rear part, in above example [condi1, cond2, ...] with
                first element is the control value
                the name of this level of index is 'Grp1'
            lst2: list of unique value of front part, in above example [treat1, treat2, ...] with
                first element is the control value
                the name of this level of index is 'Grp2'
    method arguments:
        df: input dataframe
        cat: =0 (default value) return p-value; =1 return symbol (*, <0.05 etc)
    """
    def __init__(self, lst1, lst2):
        self.lst1 = lst1
        self.lst2 = lst2

    def pwttest_1(self, df, cat):
        
        def pwt(x):
            z = ttest1.Ttest1(x, 'Rel_exp')
            s1 = z.ttest_1(gbg='Grp1', cpg='Grp2', ctrl=self.lst2[0])[cat][self.lst2[1]]
            s2 = z.ttest_1(gbg='Grp2', cpg='Grp1', ctrl=self.lst1[0])[cat][self.lst1[1]]
            s = pd.concat([s1, s2])
            return s

        return df.groupby(level='Detector').apply(pwt)
