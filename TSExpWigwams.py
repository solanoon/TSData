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
        # set default params
        self.params = {}
        self.params['replication'] = "flatten"
        self.exp = tsexp
        self.tsd = tsd  # TODO: copy
        self.exp.name = "wigwams"
        self.exp.desc = "All time replication are averaged as wigwams doesn't supports replication."
        self.exp.workdir = self.workdir = workdir

    # "flatten", "rescale", "none"
    def SetReplicationProcess(self, replication):
        self.params['replication'] = replication

    # summarize result and save
    def Summarize(self):
        self.exp.AddColumn('group', 'value')
        self.exp.AddColumn('cluster', 'file')
        self.exp.AddColumn('epsplot', 'file')
        self.exp.AddColumn('image', 'image')
        self.exp.AddColumn('pvalue', 'value')

        # load cluster info
        # - legacy code reads exported_modules
        # - new code reads swept-modules
        _legacy = False
        if (_legacy):
            wigwams_out = os.path.join(self.workdir, 'exported_modules.tsv')
            with open(wigwams_out,'r') as f:
                clusters = []
                for l in f.readlines()[1:]:
                    cluster_idx, cluster_groups, cluster_gn = l.strip().split('\t')
                    cluster_idx = int(cluster_idx)
                    while (len(clusters) < cluster_idx):
                        clusters.append({'group':cluster_groups, 'cluster':[]})
                    clusters[cluster_idx-1]['cluster'].append(cluster_gn)
        else:
            wigwams_out = os.path.join(self.workdir, 'intermediate-module-structures', 'swept_modules.tsv')
            with open(wigwams_out,'r') as f:
                clusters = []
                for l in f.readlines():
                    clustername,seed,setsize,pvalue,genelist = l.strip().split('\t')
                    pvalue = float(pvalue) * -1 # convert to log-scale value to positive
                    clusters.append({
                        'group': clustername,
                        'cluster': genelist.split(','),
                        'pvalue': pvalue
                    })
        # store cluster info
        self.exp.clusters = clusters
        self.exp.graphs = []
        for i in range(len(clusters)):
            cluster_path = "clusters/%03d.txt" % (i+1)
            plot_path = "plots/Module%03d.eps" % (i+1)
            image_path = "plots/Module%03d.png" % (i+1)
            abs_cluster_folder = os.path.join(self.workdir, "clusters")
            if (not os.path.exists(abs_cluster_folder)):
                os.mkdir(abs_cluster_folder)
            abs_cluster_path = os.path.join(self.workdir, cluster_path)
            with open(abs_cluster_path,"w") as f:
                f.write('\n'.join(clusters[i]['cluster']))
            self.exp.AddRow('%03d' % (i+1), [
                # cluster, image, eps data, pvalue
                clusters[i]['group'],
                cluster_path,
                plot_path,
                image_path,
                clusters[i]['pvalue']
            ])

    def conv_eps2png(self):
        # start converting
        # requires: ghostscript
        for idx,row in self.exp.GetTable().iterrows():
            image_path = row['image']
            abs_image_path = os.path.join(self.workdir, image_path)
            data_path = row['epsplot']
            abs_data_path = os.path.join(self.workdir, data_path)
            if (os.path.exists(abs_data_path)):
                cmd = "gs -o %s -sDEVICE=pngalpha %s" % (abs_data_path, abs_image_path)
                os.system(cmd)
            else:
                # image unavailable
                row['image'] = ''
                row['epsplot'] = ''

    def run(self):
        if (self.tsd == None):
            raise Exception('No TSData to experiment')
        if (self.workdir == None):
            raise Exception('No workdir specified to process')
        if (self.exp == None):
            raise Exception('No TSExp to experiment')
        # tidy tsd data into wigwams input format (workdir changed)
        wigwams_input_path = os.path.abspath(os.path.join(self.workdir, 'wigwams_input.csv'))
        # TODO: use params
        self.exp.params = self.params
        # must process replication
        rep_type = self.params['replication']
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
        wigwams_workdir = os.path.abspath(self.workdir)
        os.chdir(wigwams_workdir)
        sys.argv = [wigwams_workdir] + ('--Expression %s' % wigwams_input_path).split()
        wigwams_wrapper.main()
        sys.argv = old_sys_argv
        os.chdir(old_workdir)

        # summarize output data and finish
        self.Summarize()
        self.conv_eps2png()
