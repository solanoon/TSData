import numpy as np
import pandas as pd
import json
import os
from io import StringIO

##
# description about TSData format
# (first line: JSON formatted metadatas)
# (from second line: tabled data of metadata)
#
# Timeseries Metadata: Time,SampleID,CID,Offset
#
# Series: Package of multiple(or single) condition(s).
#
# 1. all conditions should have same index
# 2. time column not have be the same one (between conditions).
# 3. one file can have multiple 'series' ( can see series names in GetSeries() )
# 4. meta / event data should be in 'series' metadata.
# 5. replication event - just need to be in same 'time'.
#

##
# general utility
#

# @description: Get time data in hours
def GetTSFromString(s):
    if (s[-1] == 'm'):
        return float(s[:-1])
    elif (s[-1] == 'h'):
        return float(s[:-1]) * 60
    elif (s[-1] == 'D'):
        return float(s[:-1]) * 3600
    elif (s[-1] == 'M'):
        return float(s[:-1]) * 3600 * 30
    else:
        raise Exception("Cannot extract time info from string: %s" % s)

def GetRepFromString(s):
    raise Exception("Cannot extract replication from string: %s" % s)


# Condition desc (kind of metadata, not must required)
# used in: TSData.conditions
class TSCondition():
    def __init__(self):
        # string desc (Should be filled)
        self.species = None     # Arabidopsis, Oryza sativa, ...
        self.datatype = 'CEL'
        self.tissue = None      # Leaf, Root, ...
        self.genotype = None    # WT(None) or (modified genename)
        self.age = None         # organ age (None if not provided)
        # None or valid desc
        self.temp = None        # temperature (ex: Hot_24, Cold_4, Cold, ...)
        self.water = None       # watering level (Drought, ...)
        self.nutrition = None   # food/nutrition level

    def __repr__(self):
        #return self.__dict__
        return {
            'species': self.species,
            'tissue': self.tissue,
            'datatype': self.datatype,
            'genotype': self.genotype,
            'age': self.age,
            'temp': self.temp,
            'water': self.water,
            'nutrition': self.nutrition,
            }

    def Set(self, key, value):
        if (key == 'species'):
            self.species = value
        elif (key == 'datatype'):
            self.datatype = value
        elif (key == 'tissue'):
            self.tissue = value
        elif (key == 'genotype'):
            self.genotype = value
        elif (key == 'age'):
            self.age = value
        elif (key == 'temp'):
            self.temp = value
        elif (key == 'water'):
            self.water = value
        elif (key == 'nutrition'):
            self.nutrition = value


##
# TSData
# @description: main class that loads / saves / calculates Timeseries differential
#
class TSData:
    def metadata_init(self):
        # metadata
        self.metadata = {
            'name':None,        # name of this series data
            'date':None,        # generated or measured date
            'desc':'',          # description of data
            'author':None,      # generator of data
            'source':None,      # source of data
            'dfpath':None,      # dataframe(microarray) path. if None, then it is included, if 'default', then it is same-named txt file.
            'version': 0.1,
            }

        self.conditions = {}    # key: CID (condition ID)

    def __init__(self):
        self.metadata_init()

        # matrix includes: [SampleID, CID, Time]
        # column is ordered in this way:
        # 1h_1r, 1h_2r, 1h_3r, 2h_1r, 2h_2r, 3h_1r, ...
        self.df_meta = pd.DataFrame()
        # a well-known microarray
        self.df = pd.DataFrame()

        # options for save
        # (not inherited, not saved)
        self.cur_path = None
        self.sep = ','
        self.save_df_header_prop = None # if none, all of them are saved
        self.save_SampleID = None       # only save with columns (SampleID)


    def __str__(self):
        # tell how many samples / conditions are existing
        return "Metadata Info:\n"+\
                str(self.metadata)+"\n"+\
                "Series count: %d\n"+\
                "Gene count: %d" % (len(self.df_meta.columns), len(self.df.index))


    # @description load for general TS file
    def load(self, path):
        # internal function
        def save_section(sectionname, mat_data):
            if (sectionname == "###TSJsonData"):
                self.metadata = json.loads(mat_data)
            elif (sectionname == "###TSDataCondition"):
                datas = {}
                for l in mat_data.split('\n'):
                    datas[l[0].lower()] = l[1:]
                if ('cid' not in datas):
                    print "###TSDataCondition does not have CID, ignored."
                else:
                    for cid in datas['cid']:
                        tscond = TSCondition
                        for k in datas:
                            if (k == 'cid'):
                                continue
                            if (len(datas[k]) < idx):
                                tscond.Set(k, datas[k][idx])
                        self.conditions[cid] = tscond
                        idx += 1
            elif (sectionname == "###TSDataHeader"):
                # TODO merge matrix with previous one
                self.df_meta = pd.read_csv(StringIO(mat_data), sep=self.sep)
            elif (sectionname == "###TSDataMatrix"):
                # TODO merge matrix with previous one
                self.df = pd.read_csv(StringIO(mat_data), sep=self.sep)
            elif (sectionname == "###TSDataEvent"):
                print("###TSDataEvent section is not currently supported, sorry.")


        self.cur_path = path
        if (self.cur_path[-4:] == '.csv'):
            self.sep = ','
        else:
            self.sep = '\t'
        with open(path, "r") as f:
            l = f.readline()
            # check signature
            if (l[:9] != "###TSData"):
                raise Exception("Invalid TSData Format")
            param = l.split(',')
            if (len(param) > 1):
                self.metadata['version'] = float(param[1])
            self.metadata_init()
            # read all section to end of the file
            cmd = "###TSJsonData"
            mat_data = ""
            for l in f.readline():
                if (l[:9] == "###TSData"):
                    if (mat_data == ""):
                        print "Section %s is empty, ignored." % cmd
                    else:
                        # new section starts, save previous section
                        save_section(cmd, mat_data)
                        cmd = l
                        mat_data = ""
                else:
                    mat_data += l
            # if previous section remains, then save.
            if (mat_data != ""):
                save_section(cmd, mat_data)
            # fin
            f.close()

        # if dfpath exists, then load it.
        if (self.metadata['dfpath'] is not None):
            dfpath = self.metadata['dfpath']
            if (dfpath == 'default'):
                print("dfpath: default parameter found.")
                dfpath = path[:-4] + '.csv'
            if (dfpath == path):
                raise Exception('dfpath cannot be same with original file path')
            self.df = pd.read_csv(dfpath)

        # sanity check
        if (len(self.df.columns) != len(self.df_meta.columns)):
            raise Exception('DataHeader and DataMatrix column count is different!')


    # load data from CSV file (Header: SampleID, CID, Time)
    def load_csv(self, path):
        self.cur_path = path
        self.sep=','
        with open(path, "r") as f:
            l = f.readline() + f.readline() + f.readline()
            self.df_meta = pd.read_csv(StringIO(unicode(l)), sep=self.sep, index_col=0)
            self.df = pd.read_csv(StringIO(unicode(f.read())), sep=self.sep, index_col=0, header=None)
            self.df.columns = self.df_meta.columns  # fit column for same


    # @description save as TS file format
    def save(self, path=None):
        if (path == None):
            path = self.cur_path
            if (path == None):
                raise Exception("Should once open a file if None-path specified")
        self.cur_path = path
        
        with open(path, "w") as f:
            f.write('###TSData,0.1')
            f.write( json.dumps(self.metadata ) )
            f.write('###TSDataCondition')
            # TODO write conditions
            f.write('###TSDataHeader')
            f.write( self.df_head.to_csv() )
            if (self.metadata['dfpath'] == None):
                f.write('###TSDataMatrix')
                f.write( self.df.to_csv() ) 
            f.close()

    def save_df(self, path=None):
        if (path == None):
            path = self.metadata['dfpath']
        if (path == None):
            print("path and dfpath is all null, Aborted.")
            return False
        if (path == 'default'):
            if (self.cur_path == None):
                print("current file path isn't specified.")
            path = self.cur_path[:-4] + '.csv'
        self.df.to_csv(path, sep='\t')
        return True

    # @description clear myself.
    def clear(self):
        self.__init__()

    # @description metadata related (TODO)
    def setName(self, name):
        self.name = name
    def setType(self, type):
        self.type = type

    # @description get CIDs
    def GetCIDs(self):
        return list(set(self.df_meta.loc['CID']))
    # @description get SampleIDs
    def GetSampleIDs(self):
        return list(self.df_meta.index)

    # -------------------------
    # modifiers
    # -------------------------

    # @description Merge other TS object into current one.
    def merge(self, ts):
        # raise error if row index(genename) set is different
        if (self.df.index != ts.df.index):
            raise Exception("row index(genename) set is different")
        self.df[ts.df.columns] = ts.df
        self.df_meta[ts.df.columns] = ts.df_meta

    # @description flatten replicated TS column into single one
    # (in case of replication is unsupported, ex: wigwams)
    # @argument func: first / avg(default) / max / min
    def flatten_replication(self, func='avg'):
        def removeduplicated(a):
            r = []
            for e in a:
                if (e not in r):
                    r.append(e)
            return r
        df = self.df
        dfh = self.df_meta
        # newly generating data
        n_dfh = pd.DataFrame(index=dfh.index)
        n_df = pd.DataFrame(index=df.index)
        # filter from CID and Time
        CID_set = removeduplicated(dfh.loc['CID'])  # to be in order
        i=0
        n_df_cols = []
        for CID in list(CID_set):
            arr_bool = dfh.loc['CID'] == CID
            dfh_cond = dfh.loc[:, arr_bool]
            df_cond = df.loc[:, arr_bool]
            time_set = removeduplicated(dfh_cond.loc['Time'])
            for t in list(time_set):
                # set data column name automatically, using time & CID
                col_name = '%s_%s' % (CID, t)
                df_cond_time = df_cond.loc[:,dfh_cond.loc['Time']==t]
                if (func == 'avg'):
                    n_df[col_name] = df_cond_time.mean(axis=1)
                elif (func == 'first'):
                    raise Exception('First: NotImplemented')
                elif (func == 'max'):
                    raise Exception('Max: NotImplemented')
                elif (func == 'min'):
                    raise Exception('Min: NotImplemented')
                else:
                    raise Exception('%s: NotSupported' % func)
                n_df_cols.append(col_name)
                n_dfh[col_name] = [CID, t]
        # use result
        self.df = n_df
        self.df_meta = n_dfh

    def rescale_replication(self):
        raise Exception("replication rescaling is not supported now!")

    # @description split this timeseries from CID
    # in case of program requires each condition as separated file (ex: EDISA)
    def SplitByCID(self):
        r = []
        df = self.df
        dfh = self.df_meta
        CID_set = removeduplicated(dfh.loc['CID'])  # to be in order
        # newly generating data
        for CID in list(CID_set):
            arr_bool = dfh.loc['CID'] == CID
            dfh_cond = dfh.loc[:, arr_bool]
            df_cond = df.loc[:, arr_bool]
            tsdat = TSData()
            tsdat.df = df_cond
            tsdat.df_meta = dfh_cond
            tsdat.metadata['dfpath'] = None
            r.append(tsdat)
        return r

    def CutByCID(self, cid_group):
        # TODO bugcheck (index order may be changed)
        # create new TSData with loss of conditional data
        tsd = TSData()
        # calculate cid_intersection
        cid_available = list(self.df_meta.loc['CID'])
        cid_group = list(set(cid_group) & set(cid_available))
        if (len(cid_group) == 0):
            print "CutByCID warning - No cid group set found."
            return None
        tsd.df_meta = self.df_meta[cid_group]
        tsd.df = self.df[cid_group]
        return tsd
    
    def CutBySampleID(self, sid_group):
        # TODO bugcheck (index order may be changed)
        # create new TSData with loss of conditional data
        tsd = TSData()
        # calculate cid_intersection
        sid_group = list(set(sid_group) & set(self.df_meta.columns))
        if (len(sid_group) == 0):
            print "CutByCID warning - No cid group set found."
            return None
        tsd.df_meta = self.df_meta[sid_group]
        tsd.df = self.df[sid_group]
        return tsd

    def RemoveByCID(self, cid_group):
        raise Exception('not supported yet')

    def RemoveBySampleID(self, sid_group):
        raise Exception('not supported yet')

    #
    # depreciated
    #

    # @description load genomic data(microarray data) for time data
    def loadmatrix(arr_files, discard_cols=None):
        df_r = None
        # load files
        # merge files into a big dataframe
        for arr_file in arr_files:
            df = pd.read_table(arr_file)
            df.set_index(df.columns[0], inplace=True, drop=True)
            if (df_r == None):
                df_r = df
            else:
                df_r = df_r.merge(df, left_index=True, right_index=True, how='outer')
        # replace all NaN to Zero
        df_r.fillna(0)

        if (df_r == None):
            raise Exception("should load files")

        # discard some columns
        if (discard_cols):
            for colname in discard_cols:
                if colname in df_r.columns:
                    del df_r

        self.TS_matrix = df_r
        self.TS_available = False

    def setmatrixinfo(arr_TS_col, arr_TS_time=None, arr_TS_rep=None):
        if (self.TS_matrix == None):
            raise Exception("load matrix first")
        df = self.TS_matrix
        col_size = len(df.columns)
        if (len(arr_TS_col) != col_size):
            raise Exception("timeseries data doesn't match with orignial TS matrix")
        ts_available = True

        # if timeseries information is null, then fill it
        if arr_TS_time == None:
            arr_TS_time = []
            for col in arr_TS_col:
                arr_TS_time.append(GetTSFromString(col))
        elif (type(arr_TS_time) is list):
            if (len(arr_TS_time) == 1):
                print('[WARNING] timeseries data doesnt fit with origin file dimension (%s)', str(self.metadata['name']))
                ts_available = False

        # if replication info is null, then fill it
        if (type(arr_TS_rep) is int):
            arr_TS_rep = [arr_TS_rep,]*col_size
        if (type(arr_TS_rep) is list):
            if (len(arr_TS_rep) == 1):
                arr_TS_rep = [arr_TS_rep[0],]*col_size
            elif (len(arr_TS_rep) != col_size):
                raise Exception('[ERROR] timeseries replication data doesnt fit with origin file dimension')
        if (arr_TS_rep == None):
            arr_TS_rep = []
            for col in arr_TS_col:
                arr_TS_rep.append(GetRepFromString(col))

        self.TS_name = arr_TS_col
        self.TS_time = arr_TS_time
        self.TS_rep = arr_TS_rep

        # resort column:
        # sort first with time, then replication
        sortdata = []
        for i in range(col_size):
            sortdata.append([
                df.columns[i],
                self.TS_time[i] * 1000 + arr_TS_rep
                ])
        sorted(sortdata, key=lambda x: x[1])
        sortcols = []
        for l in sortdata:
            sortcols.append(l[0])
        df.reindex(columns=sortcols)

        # so timeseries now available
        self.TS_available = ts_available
    


# implementations for extensions with TSData
class TSExp:
    def __init__(self):
        self.expname = "TSExp_base"
        self.workdir = ""
        self.params = {}
        self.values = []
        self.clusters = []
        self.graphs = []
        # 0: finished
        # 1: pending
        # 2: processing
        # 3: (error or something else)
        self.status = 0
        self.progress = 0
        self.desc = ""

    def __str__(self):
        return "Expname: %s\n"\
                "WorkDir: %s\n"\
                "ExpParam: %s\n"\
                "Exp Result: Value %d, Clusters %d, Graphs %d"\
                % (self.expname, self.workdir, json.dumps(self.params),
                        len(self.values), len(self.clusters), len(self.graphs))

    # @description
    # load TSExp result from path
    def LoadExp(self, path):
        raise Exception("NotImplemented")

    def GetStatus(self):
        return self.status
    def SetError(self, code=3):
        self.status = code
    def SetFinish(self):
        self.status = 0
    def SetPending(self):
        self.status = 1
    def SetProcessing(self):
        self.status = 2

    def GetProgress(self):
        return self.progress
    def SetProgress(self, v):
        self.progress = v

    def GetDesc(self):
        return self.desc

    def GetLog(self):
        return ""
