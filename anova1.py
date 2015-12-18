from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

class Anova1():
    """

    """
    def __init__(self):
        pass
    
    def anova_1(self, df):
        anv = df[['Rel_exp']].reset_index('Grp')
        fmt1 = lambda x: x.split()[0]
        fmt2 = lambda x: x.split()[1]
        anv['Treat'] = anv['Grp'].map(fmt1)
        anv['Condi'] = anv['Grp'].map(fmt2)
        anv = anv.rename(columns = {'Rel_exp':'Exp'})

        def f6(x):
            formula = 'Exp ~ C(Treat) + C(Condi) + C(Treat):C(Condi)'
            lm = ols(formula, x).fit()
            return anova_lm(lm)

        grouped= anv.groupby(level='Detector')

        anova1 = grouped.apply(f6)

        anova2 = anova1[['PR(>F)']]
        anova2 = anova2[[group2 not in ['Residual'] for group1, group2 in anova2.index]]

        return anova2
