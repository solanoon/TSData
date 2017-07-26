import numpy as np
import pandas as pd
import json
import os, sys
from itertools import groupby
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


# implementations for extensions with TSData
# file format
# - first line: metadata json
# - second line ~ end: result dataframe table
class TSExp(object):
    def __init__(self):
        self.name = "TSExp_base"
        self.desc = ""
        self.log = []       # experiment result log (string list)
        #self.parent = None  # (string or None) experiment parent of this exp   # (depreciated)
        self.params = {}    # (key, value) params used in experiments
        # data contains here
        # first 2 rows:
        # (type: file, string, number, image, stringwithfile, imagewithfile, stringhide)
        # (display: 0, 1) - show default?
        # (toolname: string) - analysis tool name
        # (desc: string) - some detailed explanation with data
        # *withfile: filename separated with last ':' character.
        self.df = pd.DataFrame(index=['_toolname','_type','_display','_desc'], columns=[])


        # non-saved statuses
        self.workdir = ""
        self.children = []

    def __str__(self):
        return "Expname: %s\n"\
                "WorkDir: %s\n"\
                "ExpParam: %s\n"\
                "Exp Result: %dx%d\n"\
                "Logs: %s"\
                % (self.name, self.workdir, json.dumps(self.params), self.df.shape[0], self.df.shape[1], str(self.log))

    def to_json(self):
        r = {}
        # each cell should contain metadata consistently
        r = json.loads( self.GetTable().to_json(orient='split') )
        _toolname = self.df.loc['_toolname']
        _display = self.df.loc['_display']
        _desc = self.df.loc['_desc']
        _type = self.df.loc['_type']    # append to metadata
        data_new = []
        for row,gname in zip(r['data'], r['index']):
            row_new = [{'value':gname, 'type':'string'}]
            for v,t in zip(row, _type):
                row_new.append({'value':v, 'type':t})
            data_new.append(row_new)
        r['data'] = data_new
        del r['index']
        r['columns'] = [{'name':x} for x in r['columns']]
        r['columns'].insert(0,{'name':'name'})
        # set column information
        idx = 0
        for display,desc in zip(_display, _desc):
            r['columns'][idx]['display'] = display
            r['columns'][idx]['desc'] = desc
            idx += 1
        # use toolname to make columns_top
        columns_top = [(i, len(list(_))) for i,_ in groupby(_toolname)]
        r['columns_top'] = columns_top
        return r

    # @description
    # load TSExp result from path
    def load(self, path):
        self.workdir = os.path.dirname(path)
        with open(path,'r') as f:
            d = json.loads(f.readline())
            self.df = pd.read_csv(StringIO(unicode(f.read())), index_col=0)
        self.name = d['name']
        self.desc = d['desc']
        self.log = d['log']
        #self.parent = d['parent']
        self.params = d['params']
    def save(self, path=None):
        if (path is None):
            path = os.path.join(self.workdir, self.name+'.json')
        with open(path,'w') as f:
            json.dump({
                'name': self.name,
                'desc': self.desc,
                'log': self.log,
                #'parent': self.parent,
                'params': self.params
            }, f)
            f.write('\n')
            f.write(self.df.to_csv())

    def SetParam(self, k, v):
        self.params[k] = v
    def SetDefaultParam(self, k, v):
        if (k not in self.params):
            self.SetParam(k, v)
    def GetDesc(self):
        return self.desc

    # @description
    # add column with types and etc ...
    def AddColumn(self, cname, ctype, display=1, desc=''):
        self.df[cname] = ''
        self.df[cname]['_toolname'] = self.name
        self.df[cname]['_type'] = ctype
        self.df[cname]['_display'] = display
        self.df[cname]['_desc'] = desc

    # @description
    # add generic data row
    def AddRow(self, name, row_array):
        self.df.loc[name] = row_array

    # @description
    # get addon-analysis-data from dataframe
    # returns empty dataframe if not exists
    def GetAnalysis(self, name):
        return self.df[self.df.loc['_toolname'] == name]

    def GetOriginalTable(self):
        return self.GetAnalysis(self.name)

    # @description
    # get table without metadata (not starts with '_')
    def GetTable(self):
        return self.df.loc[ [not str(x).startswith('_') for x in self.df.index] ]

    # @description
    # concat addon-analysis-data into dataframe
    def ConcatAnalysis(self, df):
        self.df = pd.concat([self.df, df], join='outer', axis=1)

    # @desciption
    # delete addon-analysis-data from dataframe
    def DeleteAnalysis(self, name):
        self.df.drop(self.df.loc['_toolname'] == name,inplace=True,axis=1)


# @description
# gather all TSExp from `directory`, and tidy them to process
# read-only.
class TSExpGroup(object):
    def __init__(self):
        self.list_exp = []  # listed exp for visualizing
        self.dict_exp = {} # (key: expname) total experiments, just read from path...
        self.path = None

    def __repr__(self):
        return 'path: %s\nexp total count:%d\nexp list count: %d' % (self.path, len(self.dict_exp), len(self.list_exp))

    # @description
    # read experiments from path directory
    # and tidy them...
    def read(self, path):
        fp_list = []
        for fp in os.listdir(path):
            fp_rel = os.path.join(path, fp)
            if (os.path.isfile(fp_rel) and fp.endswith(".exp")):
                fp_list.append(fp_rel)
        for fp in fp_list:
            exp = TSExp()
            exp.load(fp)
            if (exp.name in self.dict_exp):
                print 'WARNING: expname %s(%s) is already exists' % (exp.name, fp)
            self.dict_exp[exp.name] = exp

        # make hierarcical list from loaded experiment
        # (ex)
        # parent
        # - child
        # - child
        # parent
        # ...
        # -> depreciated
        for _, exp in self.dict_exp.items():
            """
            if (exp.parent is not None):
                parent_name = exp.parent
                if (parent_name not in self.dict_exp):
                    print 'WARNING: expname %s has parent %s, but parent not exists.'\
                        % (exp.name, parent_name)
                    continue
                self.dict_exp[parent_name].children.append(exp)
        for _, exp in self.dict_exp.items():
            if (exp.parent is None):
                self.list_exp.append(exp)
                for child in exp.children:
                    self.list_exp.append(child)
            """
            self.list_exp.append(exp)

        self.path = path

    def to_json(self):
        r = []
        for exp in self.list_exp:
            r.append( exp.to_json() )
        return r
