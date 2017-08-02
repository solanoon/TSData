from analyzer.Trait import Trait
import os, sys
import math
#
# do advanced analyze for TSExp
# (2nd experiment modules)
#
# INFO: sub-experiment's path must be relative 
#       with parent's workpath.
#



#
# @description
# Do Trait analysis for each (gene) clusters
#
# input: (gene cluster)
# output: (result table)
class TSExpTrait(object):
    def __init__(self, exp, traitname, workdir=None):
        self.exp = exp
        self.traitname = traitname
        self.params = {}

        # work directory is automatically set
        # (under parent's workdir)
        # it's relative to parent's workdir.
        name = "%s_%s" % (exp.name, self.traitname)
        if (workdir is None):
            self.workdir = name
        else:
            self.workdir = workdir
        self.workdir = os.path.join(self.exp.workdir, self.workdir)
        print self.workdir
        if (not os.path.exists(self.workdir)):
            os.mkdir(self.workdir)

    def SetSpecies(self, species):
        self.params['species'] = species

    def run(self):
        # get param from tsd_out
        self.exp.params.update(self.params)
        if ('species' not in self.params):
            raise Exception("'species' param is required")
        species = self.params['species']
        # initialize GOTerm object
        go = Trait()
        go.load("%s-%s" % (self.traitname, species))
        # add column for GOTerm analysis
        self.exp.AddColumn(self.traitname, 'value', 'goterm')
        self.exp.AddColumn('%s-file' % self.traitname, 'file', 'goterm')
        self.exp.AddColumn('%s-pvalue' % self.traitname, 'pvalue', 'goterm')
        # run GOTerm analysis for each cluster (gene)
        # ... and save each table.
        for idx,row in self.exp.GetTable().iterrows():
            name = row.index.values.astype(str)[0]
            genelist_fn = row['cluster']
            with open(os.path.join(self.exp.workdir,genelist_fn),"r") as f:
                genelist = f.readlines()
            result_table = go.GOanalysis(genelist)
            if (not result_table.empty):
                fn = os.path.join( self.workdir, '%s.csv' % name )
                result_table.to_csv( fn , index=False)
                desc = str(result_table['GOTerm'].iloc[0])
                pval = -math.log10( float(result_table['pvalue'].iloc[0]) ) # convert to log-scale
            else:
                fn = ''
                desc = '-'
                pval = 0    # 'zero' as log10 value
            self.exp.df.loc[idx][self.traitname] = desc
            self.exp.df.loc[idx]['%s-file' % self.traitname] = fn
            self.exp.df.loc[idx]['%s-pvalue' % self.traitname] = pval
        # append result table into TSExp ...
        #self.exp.tables = result_table_exp
        #self.exp.parent = self.exp_in.name



def TSExpGOTerm(exp, workdir=None):
    return TSExpTrait(exp, 'GOTerm', workdir)

def TSExpPathway(exp, workdir=None):
    return TSExpTrait(exp, 'Pathway', workdir)

def TSExpMotif(exp, workdir=None):
    return TSExpTrait(exp, 'Motif', workdir)


