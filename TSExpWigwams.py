#import TSData
import os
from wigwams.scripts import wigwams_wrapper
#import unittest.mock
#from mock import patch
import sys

# 
# @description
# save & process TSData file for wigwams processable format
#

class TSExpWigwams(object):
    def __init__(self, tsexp, tsd, workdir):
        self.exp = tsexp
        self.tsd = tsd  # TODO: copy
        self.exp.params['replication'] = "flatten"
        self.exp.name = "wigwams"
        self.exp.desc = "All time replication are averaged as wigwams doesn't supports replication."
        self.workdir = workdir

    # "flatten", "rescale", "none"
    def SetReplicationProcess(self, replication):
        self.exp.params['replication'] = replication

    # summarize result and save
    def Summarize(self):
        wigwams_out = 'exported_modules.tsv'
        # load cluster info
        try:
            with open(wigwams_out,'r') as f:
                clusters = []
                for l in f.readlines()[1:]:
                    cluster_idx, cluster_groups, cluster_gn = l.strip().split('\t')
                    cluster_idx = int(cluster_idx)
                    while (len(clusters) < cluster_idx):
                        clusters.append({
                            'cluster': [],
                            'name': '%d cluster' % (cluster_idx+1),
                            'desc': cluster_groups,
                        })
                    clusters[cluster_idx-1]['cluster'].append(cluster_gn)
        except Exception as e:
            error_msg = str(e)
            self.exp.SetError(error_msg)
            clusters = []
            raise e
        # store cluster info
        self.exp.clusters = clusters
        self.exp.graphs = []
        for i in range(len(clusters)):
            self.exp.graphs.append({
                'path': 'plots/Module%d.eps' % (i+1),
                'desc': 'eps plot file'
            })

    def run(self):
        if (self.tsd == None):
            raise Exception('No TSData to experiment')
        if (self.workdir == None):
            raise Exception('No workdir specified to process')
        if (self.exp == None):
            raise Exception('No TSExp to experiment')
        # TODO: use params
        # tidy tsd data into wigwams input format (workdir changed)
        wigwams_input_path = 'wigwams_input.csv'
        # must process replication
        rep_type = self.exp.params['replication']
        if (rep_type == "flatten"):
            self.tsd.flatten_replication()
        elif (rep_type == "rescale"):
            self.tsd.rescale_replication()
        elif (rep_type == "none"):
            pass
        else:
            errormsg = "Unknown replication process command: %s" % rep_type
            self.exp.SetError(errormsg)
        # must convert timedata into float format & save
        self.tsd.convert_timedata_float()
        with open(wigwams_input_path, 'w') as f:
            # drop a row - SampleID (header)
            df_meta_2 = self.tsd.df_meta
            f.write(df_meta_2.to_csv(header=None))
            f.write(self.tsd.df.to_csv(header=None))
            f.close()

        # feed parameter & execute wigwams
        #with patch('sys.argv', [
        #    '--Expression', wigwams_input_path]):
        old_sys_argv = sys.argv
        old_workdir = os.getcwd()
        wigwams_workdir = os.path.join( os.getcwd(), self.workdir )
        os.chdir(wigwams_workdir)
        sys.argv = [wigwams_workdir] + ('--Expression %s' % wigwams_input_path).split()
        wigwams_wrapper.main()
        sys.argv = old_sys_argv
        os.chdir(old_workdir)

        # summarize output data and finish
        self.Summarize()
        self.exp.SetFinish()
