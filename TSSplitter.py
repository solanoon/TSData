#
# used to split microarray array into TSDatas
#

import pandas as pd
import numpy as np
import TSData
import os

class TSSplitter:
    def __init__(self):
        self.df = pd.DataFrame()
        self.groups = {}

    def LoadMatrix(self, path, sep="\t"):
        self.df = pd.read_csv(path,header=0,index_col=0,sep=sep)

    def LoadGroup_csv(self, path):
        csv_df = pd.read_csv(path, index_col=False, header=0)
        for idx, row in csv_df.iterrows():
            # generate df_meta
            TSDate = None
            TSSource = None
            TSType = 'CEL'
            TSDesc = ''
            if ('Date' in row):
                TSDate = row['Date'][0]
            if ('Source' in row):
                TSSource = row['Source'][0]
            if ('Type' in row):
                TSType = row['Type'][0]
            if ('Desc' in row):
                TSDesc = row['Desc'][0]
            self.groups[ row[[0]][0] ] = {
                'Species': row[[1]][0],
                'Stress': row[[2]][0],
                'Tissue': row[[6]][0],
                'Genotype': row[[8]][0],
                'Expid': row[[10]][0],
                'Date': TSDate,
                'Source': TSSource,
                'Type': TSType,
                'Desc': TSDesc,
                'geneinfo': [],     # genename, time, rep
                }

    # @description for depreciated LoadGroup method.
    def LoadSample_csv(self, path):
        def gene_sorter(e):
            ts = str(e[1])
            t_ = (12-len(ts))*'0' + ts
            return t_ + '_' + str(e[2])

        csv_df = pd.read_csv(path, index_col=False, header=0) # 1: SampleID, 3: CID
        CIDs = list(set(csv_df.ix[:,2].tolist()))
        for CID in CIDs:
            if (CID not in self.groups):
                print '%s CID not exists, ignore.' % CID
                continue
            group = self.groups[CID]
            b_CID = csv_df.ix[:,2] == CID
            csv_cond = csv_df[b_CID]
            genelist = csv_cond.ix[:,0].tolist()
            rep = csv_cond.ix[:,5].fillna(0).tolist()   # fill NA for use
            time = csv_cond.ix[:,10].fillna(0).tolist()
            ignore_list = (csv_cond.ix[:,5].isnull())
            # sort with rep/time
            geneinfo = zip(genelist, time, rep, ignore_list)
            geneinfo = filter(lambda x: not x[3], geneinfo) # filter invalid genes
            geneinfo.sort(key=gene_sorter)
            group['geneinfo'] = geneinfo

    # @description
    # Extract microarray file from big-microarray file using TS metadata.
    # @argument
    # if exactname false, then it'll no-case-sensitive, startswith character.
    # if includedf true, then TS will be resaved. else, extract microarray file will be saved.
    def Extract(self, tsd, exactname=False, includedf=False):
        if (self.df.empty):
            raise Exception("Must LoadMatrix(...) first")
        if (tsd.df_meta.empty):
            raise Exception("Not allow empty dataframe (GENENAME info necessary)")

        # read TSdata (should include only header file)
        #tsd = TSData.TSData()
        #tsd.load(TSpath)
        """
        if (tsd.metadata['dfpath'] == None):
            raise Exception("Should only use TSdata without microarray-less data")
        """

        df_out = None
        if (exactname):
            df_out = self.df[tsd.df_meta.columns]
        else:
            col_valid = []
            col_target = [s.upper() for s in list(tsd.df_meta.columns)]
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
            #tsd.save()
        else:
            fn = os.path.basename(TSpath)
            df_out.to_csv(
                os.path.dirname(TSpath) + "/" + os.path.splitext(fn)[0] + '.txt',
                sep='\t'
                )

    # creates multiple csv file, maybe
    # @argument
    # include: include microarray data into TS file. if not, generate separate microarray file
    def ExtractAll(self, outdir="./", included=False):
        for CID,group in self.groups.items():
            # generate df_meta
            genenames = map(lambda x: x[0], group['geneinfo'])
            if (len(genenames) == 0):
                print "%s CID has no genes, ignored" % CID
                continue
            df_meta = pd.DataFrame(index=['CID','Time'],columns=genenames)
            df_meta.index.name = 'SampleID'
            for genename,time,rep,valid in group['geneinfo']:
                df_meta[genename] = [CID, time]
            # make TSData
            tsdat = TSData.TSData()
            tsdat.df_meta = df_meta
            tsdat.metadata['name'] = CID
            tsdat.metadata['date'] = group['Date']
            tsdat.metadata['source'] = group['Source']
            tsdat.metadata['type'] = group['Type']
            tsdat.metadata['desc'] = group['Desc']
            # extract df
            self.Extract(tsdat)
            print df_meta.columns
            print tsdat.df.columns
            tsdat.save()


