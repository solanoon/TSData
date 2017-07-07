import TSData
import os
from wigwams.scripts import wigwams_wrapper
#import unittest.mock
#from mock import patch
import sys

# 
# @description
# save & process TSData file for wigwams processable format
#

class TSExpWigwams(TSData.TSExp):
    def __init__(self, tsd, workdir):
        self.tsd = tsd  # TODO: copy
        self.replication = "flatten"
        TSData.TSExp.expname = "wigwams"
        TSData.TSExp.desc = "All time replication are averaged as wigwams doesn't supports replication."
        TSData.TSExp.workdir = workdir

    # "flatten", "rescale", "none"
    def SetReplicationProcess(replication):
        self.replication = replication

    def run(self):
        # tidy tsd data into wigwams input format
        wigwams_input_path = os.path.join(TSData.TSExp.workdir, 'wigwams_input.csv')
        # must process replication
        if (self.replication == "flatten"):
            self.tsd.flatten_replication()
        elif (self.replication == "rescale"):
            self.tsd.rescale_replication()
        elif (self.replication == "none"):
            pass
        else:
            raise Exception("Unknown replication process command: %s" % self.replication)
        # must convert timedata into float format
        self.tsd.convert_timedata_float()
        with open(wigwams_input_path, 'w') as f:
            # drop a row - SampleID (header)
            df_meta_2 = self.tsd.df_meta
            f.write(df_meta_2.to_csv(header=None))
            f.write(self.tsd.df.to_csv(header=None))
            f.close()
        # feed parameter
        #with patch('sys.argv', [
        #    '--Expression', wigwams_input_path]):
        old_sys_argv = sys.argv
        sys.argv = [sys.argv[0]] + ('--Expression %s' % wigwams_input_path).split()

        # execute wigwams
        wigwams_wrapper.main()

        sys.argv = old_sys_argv

        # TODO summarize output data
        return

