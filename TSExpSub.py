from analyzer.GOTerm import GOTerm
import os, sys
#
# do advanced analyze for TSExp
# (2nd experiment modules)
#
# INFO: sub-experiment's path must be relative 
#       with parent's workpath.
#



#
# @description
# Do GOTerm analysis for each (gene) clusters
#
# input: (gene cluster)
# output: (result table)
class TSExpGOTerm(object):
    def __init__(self, exp, workdir=None):
        self.exp = exp
        self.params = {}

        # work directory is automatically set
        # (under parent's workdir)
        # it's relative to parent's workdir.
        exp.name = "%s_goterm" % exp_in.name
        if (workdir is None):
            self.exp.workdir = exp.name
        else:
            self.exp.workdir = workdir
        self.workdir = os.path.join(self.exp_in.workdir, self.exp.workdir)
        if (not os.path.exists(self.workdir)):
            os.mkdir(self.workdir)

    def SetSpecies(species):
        self.params['species'] = species

    def run(self):
        # get param from tsd_out
        self.exp.params.update(self.params)
        if ('species' not in self.params):
            raise Exception("'species' param is required")
        species = self.params['species']
        # initialize GOTerm object
        go = GOTerm()
        go.load(species)
        # add column for GOTerm analysis
        self.exp.AddColumn('goterm', 'value')
        self.exp.AddColumn('goterm-file', 'file')
        self.exp.AddColumn('goterm-pvalue', 'value')
        # run GOTerm analysis for each cluster (gene)
        # ... and save each table.
        result_table_exp = []
        for idx,row in self.exp_in.GetTable().iterrows():
            name = row.index
            genelist_fn = row['cluster']
            with open(os.path.join(self.exp_in.workdir,genelist_fn),"r") as f:
                genelist = f.readlines()
            result_table = go.GOanalysis(genelist)
            if (not result_table.empty):
                fn = os.path.join( self.exp.workdir, '%s.csv' % name )
                result_table.to_csv(os.path.join( self.exp_in.workdir, fn ), index=False)
                desc = str(result_table['GOTerm'].iloc[0])
                pval = float(result_table['GOTerm'].iloc[1])
            else:
                fn = ''
                desc = '(no GOTerm found)'
                pval = 1
            row['goterm'] = desc
            row['goterm-file'] = fn
            row['goterm-pvalue'] = pval
        # append result table into TSExp ...
        self.exp.tables = result_table_exp
        self.exp.parent = self.exp_in.name


#
# @description
# Do KEGG Pathway analysis for each (gene) clusters
#
# input: (gene cluster)
# output: (result table)
class TSExpPathway(object):
    def __init__(self):
        pass



#
# @description
# Do Promotor analysis for each (gene) clusters
#
# input: (gene cluster)
# output: (result table)
class TSExpPromotor(object):
    def __init__(self):
        pass


