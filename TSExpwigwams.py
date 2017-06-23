import TSData
import os
from wigwams.scripts import wigwams_wrapper
#import unittest.mock
from mock import patch

# 
# @description
# save & process TSData file for wigwams processable format
#

class TSExpWigwams(TSData.TSExp):
    def __init__(self, tsd, workdir):
        self.tsd = tsd
        TSData.TSExp.expname = "wigwams"
        TSData.TSExp.desc = "All time replication are averaged as wigwams doesn't supports replication."
        TSData.TSExp.workdir = workdir

    def run(self):
        # tidy tsd data into wigwams input format
        wigwams_input_path = os.path.join(TSData.TSExp.workdir, 'wigwams_input.csv')
        with open(wigwams_input_path, 'w') as f:
            f.write(self.tsd.df_meta.to_csv())
            f.write(self.tsd.df.to_csv(header=None))
            f.close()
        # feed parameter
        with patch('sys.argv', [
            '--Expression', wigwams_input_path,
            '--SizeThresholds', 50]):
            # execute wigwams
            wigwams_wrapper.main()

        # TODO summarize output data
        return

