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
    s = s.strip()
    if (s[-1] == 'm'):
        return float(s[:-1]) / 60.0
    elif (s[-1] == 'h'):
        return float(s[:-1])
    elif (s[-1] == 'D'):
        return float(s[:-1]) * 60
    elif (s[-1] == 'M'):
        return float(s[:-1]) * 60 * 30
    else:
        return float(s)
        #raise Exception("Cannot extract time info from string: %s" % s)

def GetRepFromString(s):
    raise Exception("Cannot extract replication from string: %s" % s)


# Condition desc (kind of metadata, not must required)
# used in: TSData.conditions
class TSCondition(object):
    def __init__(self):
        # string desc (Should be filled)
        self.species = None     # Arabidopsis, Oryza sativa, ...
        self.tissue = None
        self.datatype = 'CEL'
        self.genotype = 'WT'
        self.ecotype = None
        self.age = None
        self.stress = []        # kind of 'tag'

    def __repr__(self):
        #return self.__dict__
        return json.dumps(Get(self))

    # get dict array
    def Get(self):
        return {
            'age': self.age,
            'species': self.species,
            'tissue': self.tissue,
            'datatype': self.datatype,
            'genotype': self.genotype,
            'ecotype': self.ecotype,
            'stress': self.stress,
            }

    # load from json dict array
    def Set(self, d):
        self.age = d['age']
        self.species = d['species']
        self.tissue = d['tissue']
        self.datatype = d['datatype']
        self.genotype = d['genotype']
        self.ecotype = d['ecotype']
        self.stress = d['stress']

    def SetStress(self, value):
        self.stress = value.split(',')


##
# TSData
# @description: main class that loads / saves / calculates Timeseries differential
#
class TSData(object):
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
        return ("Metadata Info:\n"+\
                str(self.metadata)+"\n"+\
                "Series(Condition) count: %d\n"+\
                "Sample count: %d\n"+\
                "Gene count: %d") % (self.GetConditionCount(), len(self.df_meta.columns), len(self.df.index))

    def getCondition(self, key):
        if (key not in self.conditions):
            self.conditions[key] = TSCondition()
        return self.conditions[key]

    def removeCondition(self, key):
        if (key in self.conditions):
            del self.conditions[key]

    # fit condition metadata with df_meta.columns
    # used when load / save (to keep integrity)
    def fitCondition(self):
        cols = list(set(self.df_meta.ix['CID',:].tolist()))
        condition_valid = {}
        for k in self.conditions:
            if (k in cols):
                condition_valid[k] = self.conditions[k]
            else:
                print 'remove invalid condition %s' % k
        self.conditions = condition_valid
        for c in cols:
            if (c not in self.conditions):
                self.conditions[c] = TSCondition()


    # @description load for general TS file
    def load(self, path):
        # internal function
        def save_section(sectionname, mat_data):
            if (sectionname == "###TSJsonData"):
                self.metadata = json.loads(mat_data)
            elif (sectionname == "###TSDataCondition"):
                conds = json.loads(mat_data)
                for cid,dat in conds.items():
                    tscond = TSCondition()
                    tscond.Set(dat)
                    self.conditions[cid] = tscond
            elif (sectionname == "###TSDataHeader"):
                # TODO merge matrix with previous one
                self.df_meta = pd.read_csv(StringIO(unicode(mat_data)), sep=self.sep, index_col=0)
            elif (sectionname == "###TSDataMatrix"):
                # TODO merge matrix with previous one
                self.df = pd.read_csv(StringIO(unicode(mat_data)), sep=self.sep, index_col=0)
            elif (sectionname == "###TSDataEvent"):
                print("###TSDataEvent section is not currently supported, sorry.")


        self.cur_path = path
        if (self.cur_path[-4:] == '.txt'):
            self.sep = '\t'
        else:
            self.sep = ','
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
            for l in f.readlines():
                l = l.strip()
                if (l[:9] == "###TSData"):
                    if (mat_data == ""):
                        print "Section %s is empty, ignored." % cmd
                    else:
                        # new section starts, save previous section
                        save_section(cmd, mat_data)
                        cmd = l
                        mat_data = ""
                else:
                    mat_data += l + '\n'
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
        self.fitCondition()


    # load data from CSV file (Header: SampleID, CID, Time)
    # No condition metadata, in this case.
    def load_csv(self, path):
        self.cur_path = path
        self.sep=','
        with open(path, "r") as f:
            l = f.readline() + f.readline() + f.readline()
            self.df_meta = pd.read_csv(StringIO(unicode(l)), sep=self.sep, index_col=0)
            self.df = pd.read_csv(StringIO(unicode(f.read())), sep=self.sep, index_col=0, header=None)
            self.df.columns = self.df_meta.columns  # fit column for same
        self.fitCondition()


    # @description save as TS file format
    def save(self, path=None):
        self.fitCondition() # condition integrity

        if (path == None):
            path = self.cur_path
            if (path == None):
                raise Exception("Should once open a file if None-path specified")
        self.cur_path = path
        
        with open(path, "w") as f:
            f.write('###TSData,0.1\n')
            f.write( json.dumps(self.metadata )+'\n' )
            f.write('###TSDataCondition\n')
            d_cond = {}
            for k in self.conditions:
                d_cond[k] = self.conditions[k].Get()
            f.write( json.dumps(d_cond) )
            f.write('\n')
            f.write('###TSDataHeader\n')
            f.write( self.df_meta.to_csv() )
            if (self.metadata['dfpath'] == None):
                f.write('###TSDataMatrix\n')
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
    def SetName(self, name):
        self.metadata['name'] = name
    # @description get CIDs
    def GetCIDs(self):
        return list(set(self.df_meta.loc['CID']))
    # @description get SampleIDs
    def GetSampleIDs(self):
        return list(self.df_meta.index)
    def GetGeneCount(self):
        return len(self.df.index)
    def GetConditionCount(self):
        return len(self.conditions)

    # -------------------------
    # modifiers
    # -------------------------

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
        n_dfh = pd.DataFrame(index=['CID','Time'])
        n_dfh.index.name = 'SampleID'
        n_df = pd.DataFrame(index=df.index)
        # filter from CID and Time
        CID_set = removeduplicated(dfh.loc['CID'])  # to be in order
        i=0
        n_df_cols = []
        for CID in list(CID_set):
            arr_bool = dfh.loc['CID'] == CID
            dfh_cond = dfh.loc[:, list(arr_bool)]
            df_cond = df.loc[:, list(arr_bool)]
            time_set = removeduplicated(dfh_cond.loc['Time'])
            for t in list(time_set):
                # set data column name automatically, using time & CID
                col_name = '%s_%s' % (CID, t)
                df_cond_time = df_cond.loc[:,list(dfh_cond.loc['Time']==t)]
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
        print self.df_meta

    def rescale_replication(self):
        raise Exception("replication rescaling is not supported now!")

    #
    # @description
    # converts timedata into float(hour) format
    # raise Exception if inavailable format.
    #
    def convert_timedata_float(self):
        arr_t = self.df_meta.loc['Time'].tolist()
        arr_t = map(lambda x: GetTSFromString(x), arr_t)
        self.df_meta.loc['Time'] = arr_t

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
    


#
# ------ General Functions ------
#

# @description Merge TS objects and creates new one
def merge(arr_ts):
    if (len(arr_ts) == 0):
        raise Exception("merging TS is zero; give valid TS set.")

    ts_r = TSData()

    # check validation
    arr_idx = None
    for ts in arr_ts:
        if arr_idx is None:
            arr_idx = list(ts.df.index)
        else:
            if (arr_idx != list(ts.df.index)):
                raise Exception('Cannot merge series, as Genename index is different.')

    # sum'em all
    for ts in arr_ts:
        ts_r.conditions.update(ts.conditions)
    arr_df = map(lambda x: x.df, arr_ts)
    arr_df_meta = map(lambda x: x.df_meta, arr_ts)
    ts_r.df = pd.concat(arr_df, axis=1)
    ts_r.df_meta = pd.concat(arr_df_meta, axis=1)

    return ts_r

