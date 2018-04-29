import TSData
import pandas as pd
import numpy as np
from datetime import datetime
import os

def activateR():
    global r
    global pandas2ri
    global importr
    from rpy2.robjects import r, pandas2ri
    from rpy2.robjects.packages import importr
    # initialize pandas2ri
    pandas2ri.activate()



# TODO
# 1. loadjson / savejson implement
# 2. readmatrix implement
# 3. convertRawData implement.

# convert to TSD files
# without including microarray data
def ConvertCSVtoTSDs(fp, dest_dir='', save=True):
    df_csv = pd.read_csv(fp,encoding='utf-8')
    tsds = []

    # essential tests
    keys_to_record = ['SeriesID', 'SampleID', 'Time', 'Species', 'Stress', 'Date', 'Desc', 'Tissue', 'Genotype', 'Ecotype', 'Filepath', 'Filetype']
    if 'SeriesID' not in df_csv:
        raise Exception('Column "SeriesID" is must required.')
    if 'Date' not in df_csv:
        df_csv['Date'] = str(datetime.now())
    for k in keys_to_record[2:]:
        if  k not in df_csv:
            df_csv[k] = np.nan

    for (name, df) in df_csv.groupby('SeriesID'):
        # Name / h / stress / species / tissue / genotype / ecotype / filename / filetype
        # parse json metadata first
        dat = {}
        for k in {'SeriesID','Date','Desc'}:
            dat[k.lower()] = df.iloc[0][k]
        # VALIDATION CHECK: no comma allowed
        if (',' in dat['seriesid']):
            raise Exception("Comma in name isn't allowed!")
        dat['name']=dat['seriesid']
        # parse series metadata frame
        tsd = TSData.TSData()
        tsd.metadata.set(dat)
        df_new = df[keys_to_record].set_index('SampleID').transpose()
        tsd.appendsample(df_new)
        #print df_new
        print tsd.df_meta
        tsds.append(tsd)
        # save TSData (without metadata)
        if (save):
            tsd.save(dest_dir+name+'.tsd')
    return tsds

# read tsd file's gene matrix (at once)
def FillGenematrix(tsl, sep=','):
    df_g = TSData.GeneMatrix()
    df_g.set_refine_columns()
    df_g.set_refine_index()
    for tsd in tsl:
        for fp in tsd.df_meta.loc['Filepath']:
            if (fp == None or pd.isnull(fp)):
                continue
            df_g.load_from_path(fp, sep)
    df_g = df_g.drop_duplicated_index()      # drop duplicated genename
    for tsd in tsl:
        print 'saving %s' % tsd.metadata['name']
        tsd.readmatrix_from_matrix(df_g)
        tsd.save()

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
    def merge(self):
        # return a big merged single TSData
        tsd = TSData.TSData()
        for t in self:
            tsd.append(t)
        return tsd
    # column: timeseries, row: [labelname]
    # if mixed, returns 'mixed'
    def get_labels(self, labelname='Stress'):
        r = []
        series = []
        for tsd in self:
            g = tsd.df_meta.transpose().groupby(by='SeriesID')
            for k,df in g:
                series.append(k)
                r.append(df[labelname].tolist()[0])
        return pd.DataFrame(r, columns=[labelname], index=series)
    # column: timeseries X time, row: [labelname]
    def get_labels_per_time(self, labelname='Stress'):
        r = []
        for tsd in self:
            r.append( tsd.df_meta.loc[labelname] )
        return pd.concat(r, axis=1)
    def get_timepoints(self):
        return self.get_labels_per_time('Time')
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
            # COMMENT: must match name with getSeries() dataframe column
            minrepcnt = int(self.filters['MinRepCnt'])
            df_rep = tsd.getSeries().loc['count']
            logic_minrep_group = (df_rep >= minrepcnt)
            logic_minrep = [ logic_minrep_group[x[0]+'_'+x[1]] for x in zip(tsd.df_meta.loc['SeriesID'].tolist(),tsd.df_meta.loc['Time'].tolist()) ]
            logic = np.logical_and(logic, logic_minrep)
        if ('ExpExist' in self.filters):
            logic = logic & tsd.df.empty
        return logic




#
# DEG(gene expression) processing part
#

def _do_DEG_per_timepoint(tsd, strict_rep):
    default_rep_std = 0.2
    default_rep_cnt = 2
    strict_rep = False

    reps = tsd.getReplication()
    use_bonf_adj=False
    if (use_bonf_adj):
        thres = thres/len(reps)
    # check for exception (single timepoint not allowed)
    if (len(reps) == 1):
        raise Exception('single timepoint is not allowed!')
    # gather tp/std/mean for each tp
    rep_means = []
    rep_stds = []
    rep_cnts = []
    rep_x = []
    for rep in reps:
        rep_tp = float(rep['x'])
        samples = tsd.df[rep['samples']]
        rep_cnt = samples.shape[1]
        rep_mean = np.mean(samples,axis=1)
        if (rep_cnt == 1):
            if (strict_rep):
                raise Exception('Replication must be over 2; single replication not allowed.')
            rep_std = np.full(rep_mean.shape, default_rep_std)
            rep_cnt = default_rep_cnt
        else:
            rep_std = np.std(samples,axis=1)
        rep_means.append(rep_mean)
        rep_stds.append(rep_std)
        rep_cnts.append(rep_cnt)
        rep_x.append(rep_tp)
    # calculate Ttest p-value for each timepoint
    _shape = rep_means[0].shape
    ts_pval = [] #[np.full(_shape,0),]  # default pvalue(non_deg)
    ts_tval = []                        # t value
    ts_tp = rep_x[1:] #[-1.0,]+rep_x[1:]          # default timepoint(zero)
    mean0 = rep_means[0]
    std0 = rep_stds[0]
    cnt0 = rep_cnts[0]
    for (mean,std,cnt) in zip(rep_means[1:],rep_stds[1:],rep_cnts[1:]):
        val_t, val_p = scipy.stats.ttest_ind_from_stats(
                mean0,std0,cnt0,
                mean,std,cnt)
        """
        if (thres > 0):
            val_p[val_p < thres] = 0
            val_p[val_p >= thres] = 1
        val_p = 1-val_p
        val_sign = (val_t>0)*2-1     # get sign - down DEG or up DEG?
        val_sign *= -1
        arr_p = val_p*val_sign
        """
        ts_pval.append(ts_pval)
        ts_tval.append(ts_tval)
    # regenerate as dataframe (colname: tsd.metadata['name'], rowname: tsd.get_index())
    df_pval = pd.DataFrame(np.array(ts_pval), index=tsd.get_index(), columns=[tsd.metadata['name'],])
    df_tval = pd.DataFrame(np.array(ts_tval), index=tsd.get_index(), columns=[tsd.metadata['name'],])
    return {'pvalue': df_pval, 'tvalue': df_tval}
def do_DEG_per_timepoint(tsl, strict_rep=False):
    if (type(tsl) is TSData):
        return _do_DEG_per_timepoint(tsl)
    else:
        r_pv = []
        r_tv = []
        for tsd in tsl:
            d = _do_DEG_per_timepoint(tsd)
            r_pv.append( d['pvalue'] )
            r_tv.append( d['tvalue'] )
        r_pv = pd.concat(r_pv, axis=1) # TODO: is this correct?
        r_tv = pd.concat(r_tv, axis=1)
        return {'pvalue':r_pv, 'tvalue':r_tv}


def _do_limma(tsd, pval_mode='adj.P.Val'):
    # calculate design matrix
    df_model = modelmatrix( tsd.df_meta.loc['Time'] )
    df_model[df_model.columns[0]] = 1       # REAL design matrix
    # run rscript
    importr('limma')
    r_df = pandas2ri.py2ri(tsd.df)
    r_df_design = pandas2ri.py2ri(df_model)
    r_fit = r['lmFit'](r_df, r_df_design)
    r_fit = r['eBayes'](r_fit)
    r_fit = r['topTable'](r_fit, coef=2, adjust="fdr", n=np.inf)

    # use columns: adj.P.Val, t
    r_idx = r_fit.rownames
    c_pval = r_fit.rx2(pval_mode)
    c_tval = r_fit.rx2('t')
    df_pval = pd.Series( c_pval, index=r_idx )
    df_tval = pd.Series( c_tval, index=r_idx )
    # reorder row index
    df_pval = df_pval.reindex(index = tsd.df.index)
    df_tval = df_tval.reindex(index = tsd.df.index)
    return {'pvalue':df_pval, 'tvalue':df_tval}
# adj.P.Val or P.Value
def do_limma(tsl, pval_mode='adj.P.Val'):
    # COMMENT:
    # when appending new Series to previous pandas dataframe,
    # 1. when dataframe is empty --> All data is saved properly
    # 2. when already filled dataframe --> left-joined (some data may be removed)
    if (type(tsl) is TSData):
        return _do_limma(tsl)
    else:
        df_pv = pd.DataFrame()
        df_tv = pd.DataFrame()
        for tsd in tsl:
            d = _do_limma(tsd, pval_mode)
            df_pv[tsd.metadata['name']] = d['pvalue']
            df_tv[tsd.metadata['name']] = d['tvalue']
        return {'pvalue':df_pv, 'tvalue':df_tv}



def do_DESeq(tsl):
    raise Exception("NotImplemented")


# DEG with foldchange with t-value (TODO: max - min ?)
def _do_DEG_FC(tsd, fc=0.6):
    # get all time series
    # and get average of these std.
    series = tsd.getSeries()
    samples_first = np.mean( tsd.df[series.loc["samples"][0].split(',')] ,axis=1)
    samples_last = np.mean( tsd.df[series.loc["samples"][-1].split(',')] ,axis=1)
    series_FC = samples_last - samples_first
    #print series_FC[:20]
    df_pv = pd.Series(fc / np.abs(series_FC) * 0.05, index=tsd.df.index)
    #print df_pv
    df_tv = pd.Series(series_FC, index=tsd.df.index)
    #samples_std = np.std(samples, axis=1, ddof=1)
    return {'pvalue': df_pv, 'tvalue': df_tv}
def do_DEG_FC(tsl, fc=0.6):
    if (type(tsl) is TSData):
        return _do_DEG_FC(tsl, fc)
    else:
        df_pv = pd.DataFrame()
        df_tv = pd.DataFrame()
        for tsd in tsl:
            d = _do_DEG_FC(tsd, fc)
            df_pv[tsd.metadata['name']] = d['pvalue']
            df_tv[tsd.metadata['name']] = d['tvalue']
        return {'pvalue':df_pv, 'tvalue':df_tv}


# creates design matrix
def factor(table):
    if (type(table) == pd.DataFrame):
        return pd.unique(table.values.flatten())
    else:
        return pd.unique(table)
# creates model matrix (same with R - model.matrix)
# @input single columns / index: (sample names)
# @output columns:factor_list, index:(sample names) 
def modelmatrix(table, axis=0, factor_list=None):
    if (factor_list is None):
        factor_list = factor(table)
    # do factor encoding
    r = []
    if (type(table) == pd.DataFrame or type(table) == pd.Series):
        idx = table.index
        table = table.values.flatten()
    else:
        raise Exception('Only DataFrame/Series format is supported')
    for x in table.tolist():
        r.append(factor_list == x)
    return pd.DataFrame(np.array(r)*1, columns=factor_list, index=idx)


def compile_tfrecord(tsl):
    raise Exception('NotImplemented')

def batch(batch_size, np_arrs, dim=0):
    # COMMENT: all numpy arrays should have same row count
    idx = np.random.randint(np_arrs[0].shape[dim], size=batch_size)
    if (dim == 0):
        return [x[idx,:] for x in np_arrs]
    else:
        return [x[:,idx] for x in np_arrs]
