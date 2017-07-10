from analyzer.GOTerm import GOTerm
import os, sys
#
# do advanced analyze for TSExp
# (2nd experiment modules)
#



#
# @description
# Do GOTerm analysis for each (gene) clusters
#
# input: (gene cluster)
# output: (result table)
class TSExpGOTerm(object):
    def __init__(self, exp_in, exp, workdir=None):
        self.exp_in = exp_in
        self.exp = exp
        # work directory is automatically set
        # (under parent's workdir)
        # it's relative to parent's workdir.
        exp.name = "%s_goterm" % exp_in.name
        if (workdir is None):
            self.exp.workdir = os.path.join(exp_in.workdir, exp.name)
        else:
            self.exp.workdir = workdir
        os.mkdir(self.exp.workdir)
    def run(self):
        # get param from tsd_out
        species = self.exp.params['species']
        # initialize GOTerm object
        go = GOTerm()
        go.load(species)
        # run GOTerm analysis for each cluster (gene)
        # ... and save each table.
        result_table_exp = []
        for cluster in self.exp_in.clusters:
            name = cluster['name']
            genelist = cluster['cluster']
            result_table = go.GOanalysis(genelist)
            fn = os.path.join( self.exp.workdir, '%s.csv' % name )
            result_table.to_csv(fn, index=False)
            result_table_exp.append({
                'name': 'GOTerm of cluster %s' % name, 
                'desc': '',
                'path':fn
            })
        # append result table into TSExp ...
        self.exp.tables = result_table_exp


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


