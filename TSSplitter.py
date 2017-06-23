#
# used to split microarray array into TSDatas
#

import pandas as pd
import numpy as np
import TSData
import os

class TSSplitter:
    def __init__(self):
        self.df = None

    def LoadMatrix(self, path):
        self.df = pd.read_csv(path,header=0,index_col=0)

    # @description
    # Extract microarray file from big-microarray file using TS metadata.
    # @argument
    # if exactname false, then it'll no-case-sensitive, startswith character.
    # if includedf true, then TS will be resaved. else, extract microarray file will be saved.
    def Extract(self, TSpath, exactname=False, includedf=False):
        if (self.df == None):
            raise Exception("Must LoadMatrix(...) first")

        # read TSdata (should include only header file)
        tsd = TSData()
        tsd.load(TSpath)
        """
        if (tsd.metadata['dfpath'] == None):
            raise Exception("Should only use TSdata without microarray-less data")
        """

        df_out = None
        if (exactname):
            df_out = self.df[TSpath.df.index]
        else:
            col_valid = []
            col_target = [s.upper() for s in list(TSpath.df.index)]
            col_valid_cnt = 0
            for col in tsd.df.index:
                col_up = col.upper()
                col_res = False
                for col_comp in col_target:
                    if (col_up.startswith(col_comp)):
                        col_res = True
                        col_valid_cnt += 1
                        break
                col_valid.append(col_res)
            if (col_valid_cnt < col_target):
                raise Exception("Not all columns included in matrix file")
            elif (col_valid_cnt > col_target):
                raise Exception("Too many columns selected; maybe invalid filter")
            df_out = self.df.loc[:,col_valid]

        # in case of includedf
        if (includedf):
            tsd.df = df_out
            tsd.save()
        else:
            fn = os.path.basename(TSpath)
            df_out.to_csv(
                os.path.dirname(TSpath) + "/" + os.path.splitext(fn)[0] + '.txt',
                sep='\t'
                )

    # creates multiple csv file, maybe
    # @argument
    # include: include microarray data into TS file. if not, generate separate microarray file
    def ExtractFromCSV(self, outdir="./", included=False):
        csvs = CSVtoTS(csvpath)
        for csv in csvs:
            print 'processing series: %s' % csv.metadata['name']
            csvpath = outdir + '/' + csv.metadata['name']
            csv.save(csvpath)
            self.Extract(csvpath)


#
# utilities start
#

# Requires [CID, SampleID, Time]
# Automatically will be ordered by time
def CSVtoTS(csvpath):
    r = []
    csv = pd.read_csv(csvpath, header=0)
    CIDs = set(csv['CID'])
    for CID in list(CIDs):
        df = csv.loc[csv['CID'] == CID]
        # generate Timeval
        timeval = [] #TODO
        df['Timeval'] = timeval
        # order replication first, then time
        df.sort(['Timeval', 'Rep'], inplace=True)
        # generate df_meta
        TSDate = None
        TSSource = None
        TSType = 'CEL'
        TSDesc = ''
        if ('Date' in df):
            TSDate = df['Date'][0]
        if ('Source' in df):
            TSSource = df['Source'][0]
        if ('Type' in df):
            TSSource = df['Type'][0]
        if ('Desc' in df):
            TSSource = df['Desc'][0]
        df_meta = pd.DataFrame(index=['SampleID','CID','Time'],columns=df['SampleID'])
        for ind, row in df.iterrows():
            df_meta[row['SampleID']] = [row['CID'], row['Time']]
        # fill metadata
        tsdat = TSData()
        tsdat.df_meta = df_meta
        tsdat.metadata['name'] = CID
        tsdat.metadata['date'] = TSDate
        tsdat.metadata['source'] = TSSource
        tsdat.metadata['type'] = TSType
        tsdat.metadata['desc'] = TSDesc
        # append and finish
        r.append(tsdat)
    return r

#
# requires [CID, Date, Type, Author, Source, Desc]
#
def FillTSMetadata(csvpath, tsdats):
    pass
