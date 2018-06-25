import numpy as np     #numpy是提供矩阵运算的库
import pandas as pd    #基于numpy，内含dataframe和series两种数据
import json            #json是一种轻量级的数据交换格式，易于人阅读和编写
import os              #os模块包含普遍的操作系统功能
from io import StringIO   #io模块是用来处理各种类型的I/O操作流
import collections        #collections提供了许多有用的集合类

##
# description about TSData format
# (first line: JSON formatted metadatas)         第一行 json编码的元数据
# (from second line: tabled data of metadata)    第二行起 表格数据
#
# Timeseries Metadata: Time,SampleID,CID,Offset   时间序列元数据-时间，样本id，CID，弥补(?)
#
# Series: Package of multiple(or single) condition(s).               series: 单一/多个条件
#
# 1. all conditions should have same index                           所有的条件 指数一致
# 2. time column not have be the same one (between conditions).      时间行不同一
# 3. one file can have multiple 'series' ( can see series names in GetSeries() )    一个文件可以有很多series
# 4. meta / event data should be in 'series' metadata.               元/事件 数据需在‘series’元数据里
# 5. replication event - just need to be in same 'time'.             republication事件需‘时间’相同
#

##
# general utility    一般效用
#

# @description: Get time data in hours     获取时间信息
def GetTSFromString(s):                    #定义函数s
    s = s.strip()                          #用于移除字符串头尾指定的字符
    if (s[-1] == 's'):                     #结尾是s/负数
        s = s[:-1]  # days, months...
    if (s[-1] == 'm'):
        return float(s[:-1]) / 60.0
    elif (s[-3:] == 'min'):                #结尾的三个字符是min->变换成小时
        return float(s[:-3]) / 60.0
    elif (s[-1] == 'h'):
        return float(s[:-1])
    elif (s[-4:] == 'hour'):               #单位是小时hour则保持不变
        return float(s[:-4])
    elif (s[-1] == 'D' or s[-1] == 'd'):
        return float(s[:-1]) * 60
    elif (s[-3:] == 'day'):
        return float(s[:-3]) * 60
    elif (s[-1] == 'M'):
        return float(s[:-1]) * 60 * 30     #这三个elif不太理解 难道不是x24x30吗
    else:
        return float(s)
        #raise Exception("Cannot extract time info from string: %s" % s)

def GetRepFromString(s):
    raise Exception("Cannot extract replication from string: %s" % s)    #raise引发exception异常

def convertTime2Int(l):
    return [GetTSFromString(e) for e in l]    #将时间数据的格式转换成int形式

def convertTime2Str(i):
    # per hour
    r = []
    if (i > 24):
        r.append(str(int(i/24))+'d')     #输入的数据若大于24，则转换为天数并在后面加d
        i = i % 24
    if (i >= 1):
        r.append(str(int(i))+'h')        #输入的数据在1-24之间则在后面加h
        i = i % 1
    if (i > 0):
        r.append(str(int(i*60))+'m')     #输入的数据在0-1之间则乘60(?)+m
    return ' '.join(r)                   #将空格和r字符串相连


#
# gene matrix refiner
# ( internally used in TSData::readmatrix() )
#
class GeneMatrix(pd.DataFrame):
    def __init__(self):
        super(GeneMatrix, self).__init__()
        self._files_read = []
        self._refine_columns = False
        self._refine_index = False

    def load_from_path(self, path, sep=','):
        if (path not in self._files_read):
            df = pd.read_csv(path, index_col=0, sep=sep)
            self.load_from_df(df)
            self._files_read.append(path)

    def load_from_df(self, df):
        if (self._refine_columns):
            l = [x.replace('.CEL.gz','') for x in df.columns.tolist()]
            for i in range(len(l)):
                if (l[i][:3] == "GSM"):
                    l[i] = l[i].split('_',1)[0]
            df.columns = l
        if (self._refine_index):
            df.index = [x.split('_',1)[0] for x in df.index.tolist()]
        print df.columns
        self[df.columns] = df
        #self = pd.concat( (self,df), axis=1 )

    def set_refine_columns(self, v=True):
        self._refine_columns = v

    def set_refine_index(self, v=True):
        self._refine_index = v

    def reindex_from_listfile(self, fp):
        with open(fp,'r') as f:
            l = f.read().strip().split('\n')
        self = self.reindex(l)

    def reindex_from_df(self, df):
        self = self.reindex(df.index)

    # it's no-inplace function
    def drop_duplicated_index(self):
        return self[~self.index.duplicated(keep='first')]


# (DEPRECIATED method)
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
        self.exception = False  # is it 'true' timeseries, or just control/stress?

    def __repr__(self):
        return json.dumps(self.Get())

    # get dict array (just for DEPRECIATED case)
    def get(self):
        return {
            'Age': self.age,
            'Species': self.species,
            'Tissue': self.tissue,
            'Filepath': '',
            'Filetype': self.datatype,
            'Genotype': self.genotype,
            'Ecotype': self.ecotype,
            'Stress': self.stress,
            'Exception': None,
            }

    # load from json dict array
    def load(self, d):
        self.age = d['age']
        self.species = d['species']
        self.tissue = d['tissue']
        self.datatype = d['datatype']
        self.genotype = d['genotype']
        self.ecotype = d['ecotype']
        self.stress = d['stress']
        if ('exception' in d):
            self.exception = d['exception']

    def SetStress(self, value):
        self.stress = value.split(',')

    # use only in case of None value isn't supported
    def FillNaN(self):
        if (self.age == None):
            self.age = ''
        if (self.species == None):
            self.species = ''
        if (self.tissue == None):
            self.tissue = ''
        if (self.datatype == None):
            self.datatype = ''
        if (self.genotype == None):
            self.genotype = ''
        if (self.ecotype == None):
            self.ecotype = ''
        if (self.stress == None):
            self.stress = []
        if (self.exception == None):
            self.exception = ''

class TSDict(dict):
    def set(self, d):
        keys = self.keys()
        for k in keys:
            if (k in d):
                self[k] = d[k]
    def fillNaN(self):
        keys = self.keys()
        for k in keys:
            if (self[k] == None):
                self[k] = ''

##
# TSData
# @description: main class that loads / saves / calculates Timeseries differential
#
class TSData(object):
    def metadata_init(self):
        # metadata
        self.metadata = TSDict({
            'name':None,        # name of this series data
            'date':None,        # generated or measured date
            'desc':'',          # description of data
            'version': 0.2,
            })

    def __init__(self):
        self.metadata_init()

        #
        # - matrix includes: [
        #   SampleID
        #   SeriesID
        #   Replication
        #   Time            (Control type should be set to 0)
        #   Stress          (comma separation accepted)
        #   Stressdetail    (detail information on stress; Temperature, Recovering ...)
        #   Tissue
        #   Genotype
        #   Ecotype
        #   Exception       (timeseries: None, else string; control-exp?)
        # ]
        #
        # - column ordering: (series name, time point, replication point)
        # - column name: SampleID
        #
        self.df_meta = pd.DataFrame(index=[
            'SeriesID', 'Replication', 'Time',
            'Stress', 'Stressdetail',
            'Tissue', 'Genotype', 'Ecotype', 'Age', 'Species',
            'Exception', 'Filepath', 'Filetype',
            'Desc'
            ])
        self.df_meta.index.name = 'SampleID'
        # a well-known microarray
        self.df = pd.DataFrame()

        # options for save
        # (not inherited, not saved)
        self.cur_path = None
        self.workdir = None
        self.sep = ','


    def __str__(self):
        # tell how many samples / conditions are existing
        print self.getReplication()
        return (("Metadata Info:\n"+\
                str(self.metadata)+"\n"+\
                "Series Info:\n%s") % (str(self.df_meta))
                + ('\nTimePoints:\n%s' % str(self.getTimePoints()))
                + ('\nReplications: %d' % len(self.getReplication().columns.tolist()))
                )

    def getSampleNames(self):
        return self.df_meta.iloc[0]

    # if force_convert=='integer', all item is converted into sequential integer number.
    #   ex) 0, 0, later, later, last, last --> 0, 0, 1, 1, 2, 1
    # if force_convert=='string', all item is converted into readable string.
    #   ex) 0, 60, 120 --> 0, 1h, 2h
    def getTimePoints(self, force_convert=None):
        tps = self.df_meta.loc['Time']
        if (force_convert == 'integer'):
            i = 0
            tpnow = tp[0]
            r = []
            for tp in tps:
                if tp != tpnow:
                    tpnow = tp
                    i += 1
                r.append(i)
            return r
        elif (force_convert == 'string'):
            tpnow = tp[0]
            conv = lambda x: float(x) if x.replace('.','').isdigit() else x
            t = conv(tpnow)
            r = []
            for tp in tps:
                if tp != tpnow:
                    tpnow = tp
                    t = conv(tpnow)
                    i += 1
                r.append(t)
            return r
        return tps

    def getConditionNames(self):
        #return self.conditions.keys()
        return self.df_meta.columns.tolist()

    def getExpression(self, gn):
        return self.df.loc()[gn]

    #
    # get information about replication (means same series & timepoint)
    # and about genemic expression (if provided)
    #
    # returns: DataFrame
    # - columns: SeriesID
    # - index: count, time, avg, min, max, std, samples(joined with comma)
    #
    # gn: specific gene-name to get avg/min/max/std
    #
    def getSeries(self, gn=None):
        df_rep = pd.DataFrame(index=['count','time','avg','min','max','std','samples'])
        df_rep.index.name = 'SeriesID'
        # group by 'SeriesID' and 'Time'
        groups = list(self.df_meta.transpose().groupby(['SeriesID','Time']))
        # sort once more (as timepoint is wrong recognized as string)
        def sort_group_time(a,b):
            r = float(a[0][1]) - float(b[0][1])
            if (r > 0):
                return 1
            elif (r == 0):
                return 0
            else:
                return -1
        groups.sort(sort_group_time)
        for k,df_rep_s in groups:
            samples = df_rep_s.transpose().columns
            cnt = len(samples)
            sample_str = ','.join( samples.tolist() )
            time = k[1]
            # split samples (row)
            df_expr_split = self.df
            if (gn is not None):
                df_expr_split = self.getExpression(gn)
            # split samples (column)
            df_expr_split = np.array(self.df[samples])
            # get statistics from split
            _avg = np.mean(df_expr_split)
            _min = np.min(df_expr_split)
            _max = np.max(df_expr_split)
            _std = 0
            if (df_expr_split.shape[1] > 1):
                _std = np.std(df_expr_split, ddof=1)
            # add new record
            df_rep[k[0]+'_'+k[1]] = [
                cnt, time, _avg, _min, _max, _std, sample_str
                ]
        return df_rep
    # @description (DEPRECIATED)
    def getReplication(self, gn=None, force_convert=None):
        return self.getSeries(gn)

    def fix(self):
        # sort series metadata : by seriesname / time / replication
        def sort_manually(a,b):
            al = (a[1],a[2],a[3]) 
            bl = (b[1],b[2],b[3])
            if (al > bl):
                return 1
            else:
                return 0
        col_sort_data = zip(
            self.df_meta.columns.tolist(),
            self.df_meta['SeriesID'].tolist(),
            self.df_meta['Time'].tolist(),
            self.df_meta['Replication'].tolist())
        col_sort.sort(sort_manually)
        cols = [x[0] for x in col_sort]
        self.df_meta = self.df_meta[cols]
        self.df = self.df[ self.df_meta.columns ]
        # last: check validation of datatable top-column
        if (self.df.index.name):
            self.df.index.name = self.df.index.name.replace('#','_')
        idx=self.df.index.tolist()
        for i in range(len(idx)):
            idx[i] = idx[i].replace('#','_')
        self.df.index = idx

    def appendsample(self, df_meta_new):
        self.df_meta = pd.concat([self.df_meta, df_meta_new], axis=1, join_axes=[self.df_meta.index])

    # @description load for general TS file
    def load(self, path):
        # internal function
        _cond = {}
        def save_section_v01(sectionname, mat_data):
            if (sectionname == "###TSJsonData"):
                self.metadata = json.loads(mat_data)
            elif (sectionname == "###TSDataCondition"):
                # DEPRECIATED section, need converting to dataframe
                conds = json.loads(mat_data)
                for cid,dat in conds.items():
                    tscond = TSCondition()
                    tscond.load(dat)
                    tscond.stress = ','.join(tscond.stress)
                    _cond[cid] = tscond
            elif (sectionname == "###TSDataHeader"):
                # add column if not exists
                _uni = unicode(mat_data.decode('utf-8'))
                df_meta_old = pd.read_csv(StringIO(_uni), sep=self.sep, index_col=0)
                # need to change some row name
                df_meta_old = df_meta_old.rename(index={'CID':'SeriesID', 'Title':'Desc'})
                # merge df_meta, with adding column & left-innerjoin row
                #self.df_meta = pd.concat([self.df_meta, df_meta_old], axis=1, join_axes=[self.df_meta.index])
                self.appendsample(df_meta_old)
            elif (sectionname == "###TSDataMatrix"):
                _uni = unicode(mat_data.decode('utf-8'))
                self.df = pd.read_csv(StringIO(_uni), sep=self.sep, index_col=0)
            elif (sectionname == "###TSDataEvent"):
                print("###TSDataEvent section is not currently supported, sorry.")

        def save_section_v02(sectionname, mat_data):
            if (sectionname == "###TSJsonData"):
                self.metadata = json.loads(mat_data)
            elif (sectionname == "###TSDataHeader"):
                _uni = unicode(mat_data.decode('utf-8'))
                self.df_meta = pd.read_csv(StringIO(_uni), sep=self.sep, index_col=0)
            elif (sectionname == "###TSDataMatrix"):
                _uni = unicode(mat_data.decode('utf-8'))
                self.df = pd.read_csv(StringIO(_uni), sep=self.sep, index_col=0)
            elif (sectionname == "###TSDataEvent"):
                print("###TSDataEvent section is not currently supported, sorry.")

        self.cur_path = path
        self.workdir = os.path.dirname(path)
        save_section = save_section_v02
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
                if (self.metadata['version'] < 0.2):
                    save_section = save_section_v01
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

        # COMPATIBILITY WORK
        if (self.metadata['version'] == 0.1):
            # make temporary dict for (sample - series) pair
            _samplepair = {}
            for sampleid in self.df_meta.columns:
                _samplepair[sampleid] = self.df_meta[sampleid]['SeriesID']
            # fill df_meta (using series & sample name)
            for sampleid in self.df.columns:
                _seriesid = _samplepair[sampleid]
                _tscond = _cond[_seriesid]
                _metadata = _tscond.get()
                self.df_meta[sampleid]['SeriesID'] = _seriesid
                for k,d in _metadata.items():
                    self.df_meta[sampleid][k] = d
            self.metadata['version'] = 0.2

        # sanity check
        if (not self.df.empty):
            if (len(self.df.columns) != len(self.df_meta.columns)):
                raise Exception('DataHeader and DataMatrix column count is different!')
        # check and rematch column order
        #self.df = self.df[ self.df_meta.columns ]
        # fix order / wrong character
        #self.fix()

    # (DEPRECIATED METHOD)
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


    # @description save as TS file format
    def save(self, path=None):
        if (path == None):
            path = self.cur_path
            if (path == None):
                raise Exception("Should once open a file if None-path specified")
        self.cur_path = path
        
        with open(path, "w") as f:
            f.write('###TSData,0.2\n')
            f.write( json.dumps(self.metadata )+'\n' )
            f.write('###TSDataHeader\n')
            f.write( self.df_meta.to_csv(encoding="utf-8") )
            if (not self.df.empty):   # only add matrix data if dataframe exists
                f.write('###TSDataMatrix\n')
                f.write( self.df.to_csv(encoding="utf-8") ) 
            f.close()
        return True

    def readmatrix(self, fix_genename=True):
        df_g = GeneMatrix()

        # TODO: before reading, preprocess & read matrix files ...?
        for c in self.df_meta.columns:
            fp = self.df_meta[c]['Filepath']
            if (fp == None or pd.isnull(fp)):
                print '[WARNING] %s filepath is NaN. Canceled.' % c
                return
            df_g.load_from_path(fp)
        self.df = df_g[self.df_meta.columns]

    def readmatrix_from_matrix(self, df_mat):
        # set df with given matrix (samplename x genename)
        self.df = df_mat[ self.df_meta.columns ]

    # @description clear myself.
    def clear(self):
        self.__init__()

    # @descript filter logic to df / df_meta dataframe.
    def filter_by_logic(self, logic):
        self.filter(self.df.columns[logic])
    # @descript filter by SampleID s.
    def filter(self, names):
        self.df = self.df[names]
        self.df_meta = self.df_meta[names]

    # @description get CIDs
    def GetSeries(self):
        return list(set(self.df_meta.loc['SeriesID']))
    # @description get SampleIDs
    def GetSampleIDs(self):
        return list(self.df_meta.columns)
    def get_index(self):
        return self.df.index.tolist()
    def GetGeneCount(self):
        return len(self.df.index)
    def IsTSExists(self, tsname):
        return tsname in self.conditions
    # returns timepoint, replicates
    def GetTSTimepoint(self, tsname):
        arr_bool = list(self.df_meta.loc['SeriesID',:] == tsname)
        df_extracted = self.df_meta.loc[:,arr_bool]
        timepoints = []
        replicates = {}
        for i in list(df_extracted.loc['Time',:]):
            if i not in timepoints:
                timepoints.append(i)
                replicates[i] = 0
            replicates[i] += 1
        replicates_arr = map(lambda x: replicates[x], timepoints)
        return (timepoints, replicates_arr)
    # get all TS description(including timepoints/replicates)
    # returns: [ { dict, id, timepoints:[], replicates:[] } ]
    def GetAllTSDescription(self):
        r = []
        for k in self.conditions:
            d = self.conditions[k].Get()
            timepoints, replicates = self.GetTSTimepoint(k)
            d['id'] = k
            d['timepoints'] = timepoints
            d['replicates'] = replicates
            r.append(d)
        return r
    # get Gene expression profile for time using specific gene markers
    # (expression value is averaged, x pos calculated with: (index)/(colsize) )
    # returns: [ (x:0~1 float, y:avg. gene expr.) ]
    def GetProfile(self, gnames):
        r = []
        df_extracted = self.df.loc[gnames,:]
        df_mean = df_extracted.mean(axis=0)
        y = list(df_mean)
        x = [float(i)/len(y) for i in range(len(y))]
        return zip(x,y)

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

def load(fpath):
    tsd = TSData()
    tsd.load(fpath)
    return tsd





if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process TSData.')
    parser.add_argument('command', type=str, default='open',
            help='open / fix')
    parser.add_argument('file', type=str, help='file for command')
    args = parser.parse_args()
    if (args.command == 'open'):
        tsd = load(args.file)
        print tsd
    elif (args.command == 'fix'):
        tsd = load(args.file)
        print 'fixing time order / name / etc ...'
        tsd.fix()
        tsd.save()
    else:
        parser.print_help()
