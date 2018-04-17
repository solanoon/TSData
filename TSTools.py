import TSData
import pandas as pd
import numpy as np
from datetime import datetime
import os

# TODO
# 1. loadjson / savejson implement
# 2. readmatrix implement
# 3. convertRawData implement.

# convert to TSD files
# without including microarray data
def ConvertCSVtoTSDs(fp, dest_dir=''):
    df_csv = pd.read_csv(fp)
    tsds = []

    # essential tests
    keys_to_record = ['SeriesID', 'SampleID', 'Time', 'Species', 'Date', 'Desc', 'Tissue', 'Genotype', 'Ecotype', 'Filepath', 'Filetype']
    if 'SeriesID' not in df_csv:
        raise Exception('Column "SeriesID" is must required.')
    if 'Date' not in df_csv:
        df_csv['Date'] = str(datetime.now())
    for k in keys_to_record[2:]:
        if  k not in df_csv:
            df_csv[k] = np.nan

    for (name, df) in df_csv.groupby('SeriesID'):
        # Name / h / stress / species / tissue / genotype / ecotype / filename / filetype
        dat = {}
        for k in {'SeriesID','Date','Desc'}:
            dat[k.lower()] = df.iloc[0][k]
        dat['name']=dat['seriesid']
        # save to metadata / conditions
        tsd = TSData()
        tsd.metadata.set(dat)
        df_new = df[keys_to_record].set_index('SampleID').transpose()
        tsd.appendsample(df_new)
        #print df_new
        print tsd.df_meta
        tsds.append(tsd)
        # save TSData (without metadata)
        tsd.save(dest_dir+name+'.tsd')
    return tsds

def ConvertRawData(data, datatype):
    pass


#
# loads Timeseries folder and manages them
#
class TSLoader(list):
    def __init__(self):
        super(TSLoader, self).__init__()
        # available filter condition:
        # - Stress
        # - Sample
        # - Species
        # - MinRepCnt (replication count min)
        # - ExpExist (expression data should be exist)
        self.filters = {}
        # filter mode: 'drop' / 'modify'
        self.filter_mode = 'drop'
    def loadpath(self, path='data'):
        # read Timeseries files from given path(directory)
        fps = [f[:-4] for f in os.listdir(path) if f[-4:] == '.tsd']
        for fp in fps:
            tsd = TSData.load(os.path.join(path, fp+'.tsd'))
            # check out filter
            logics = self._filter_check_logic(tsd)
            if (not logics.all()):
                if (self.filter_mode == 'modify'):
                    tsd.filter_by_logic(logics)
                elif (self.filter_mode == 'drop'):
                    print 'dropped: %s (%s)' % (fp, tsd.df_meta.loc['Species'].tolist()[0])
                    continue
            self.append(tsd)
    def addfilter(self,k,v):
        self.filters[k]=v
    def clearfilter(self):
        self.filters = {}
    def filter(self):
        # iterate list and drop TS in case of not matching to filter
        if (self.filter_mode == 'modify'):
            self = [x.filter_by_logic(self._filter_check_logic(x)) for x in self]
        elif (self.filter_mode == 'drop'):
            self = filter(lambda t: self._filter_check_logic(t).all(), self)
    def export_merged(self):
        # return a big merged single TSData
        tsd = TSData.TSData()
        for t in self:
            tsd.append(t)
        return tsd
    def export_tfrecord(self):
        raise Exception('NotImplemented')
    # (internal function)
    # check is TSData suitable to filter
    def _filter_check_logic(self,tsd):
        # only filter for metadata
        filter_metadata = dict(self.filters)
        if ('MinRepCnt' in filter_metadata):
            del filter_metadata['MinRepCnt']
        if ('ExpExist' in filter_metadata):
            del filter_metadata['ExpExist']
        # first check metadata valid
        logics = []
        for k,v in filter_metadata.items():
            logics.append(tsd.df_meta.loc[k].str.match(v, case=False))
        logic = np.logical_and.reduce( logics )
        # then check replication test
        if ('MinRepCnt' in self.filters):
            minrepcnt = int(filter_metadata['MinRepCnt'])
            df_rep = tsd.getSeries()
            logic_minrep = (df_rep.loc['repcnt'] < minrepcnt)
            logic = np.logical_and(logic, logic_minrep)
        if ('ExpExist' in self.filters):
            logic = logic & tsd.df.empty
        return logic
