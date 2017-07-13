import numpy as np
import pandas as pd
import json
import os, sys


# implementations for extensions with TSData
class TSExp(object):
    def __init__(self):
        self.name = "TSExp_base"
        self.parent = None  # (string or None) experiment parent of this exp
        self.workdir = ""
        self.params = {}
        self.values = []    # (value: string, name, desc: string) single value
        self.clusters = []  # (cluster: [string], name, desc: string) clusters
        self.graphs = []    # (path: string, name, desc) cytoscape.js ?
        self.images = []    # (path: string, name, desc) only png, jpg, ...
        self.files = []     # (path: string, name, desc) other types
        self.tables = []    # (path: string, name, desc) mostly csv-formatted data

        # non-saved statuses

        # 0: finished
        # 1: pending
        # 2: processing
        # 3: (error or something else)
        self.status = 0
        self.stat_msg = ""
        self.progress = 0
        self.desc = ""
        self.children = []

    def __str__(self):
        return "Expname: %s\n"\
                "WorkDir: %s\n"\
                "ExpParam: %s\n"\
                "Exp Result: Value %d, Clusters %d, Graphs %d"\
                % (self.name, self.workdir, json.dumps(self.params),
                        len(self.values), len(self.clusters), len(self.graphs))

    # @description
    # load TSExp result from path
    def load(self, path):
        with open(path,'r') as f:
            d = json.load(f)
        self.values = d['values']
        self.clusters = d['clusters']
        self.graphs = d['graphs']
        self.images = d['images']
        self.files = d['files']
        self.tables = d['tables']
        self.params = {}

        status = d['status']
        self.name = status['name']
        self.parent = status['parent']
        self.status = status['status']
        self.stat_msg = status['stat_msg']
        self.progress = status['progress']
        self.desc = status['desc']
    def save(self, path=None):
        if (path is None):
            path = os.path.join(self.workdir, self.name+'.json')
        with open(path,'w') as f:
            json.dump({
                'status': {
                    'name': self.name,
                    'parent': self.parent,
                    'status': self.status,
                    'stat_msg': self.stat_msg,
                    'progress': self.progress,
                    'desc': self.desc
                },
                'values': self.values,
                'clusters': self.clusters,
                'graphs': self.graphs,
                'images': self.images,
                'files': self.files,
                'tables': self.tables
            }, f)

    # @description
    # Get as formatted one (this can be customized ...)
    def format(self):
        r = []
        for v in self.values:
            d = dict(v)
            d['type'] = 'value'
            r.append(d)
        for v in self.clusters:
            d = dict(v)
            d['type'] = 'clusters'
            r.append(d)
        for v in self.graphs:
            d = dict(v)
            d['type'] = 'graphs'
            r.append(d)
        for v in self.images:
            d = dict(v)
            d['type'] = 'images'
            r.append(d)
        for v in self.files:
            d = dict(v)
            d['type'] = 'files'
            r.append(d)
        for v in self.tables:
            d = dict(v)
            d['type'] = 'tables'
            r.append(d)
        return r

    def GetStatus(self):
        return self.status
    def SetError(self, stat_msg='', code=3):
        self.status = code
        self.stat_msg = stat_msg
    def SetFinish(self):
        self.status = 0
    def SetPending(self):
        self.status = 1
    def SetProcessing(self):
        self.status = 2
    def IsFinish(self):
        return self.status == 0
    def IsProcessing(self):
        return self.status == 2
    def IsError(self):
        return self.status >= 3

    def GetProgress(self):
        return self.progress
    def SetProgress(self, v):
        self.progress = v

    def SetParam(self, k, v):
        self.params[k] = v

    def GetDesc(self):
        return self.desc
    def GetMessage(self):
        return self.stat_msg


# @description
# gather all TSExp from `directory`, and tidy them to process
# read-only.
class TSExpGroup(object):
    def __init__(self):
        self.list_exp = []  # listed exp for visualizing
        self.total_exp = {} # (key: expname) total experiments, just read from path...
        self.path = None

    def __repr__(self):
        return 'path: %s\nexp count: %d' % (self.path, len(self.list_exp))

    # @description
    # read experiments from path directory
    # and tidy them...
    def read(self, path):
        fp_list = []
        for fp in fp_list:
            exp = TSExp()
            exp.load(fp)
            if (exp.name in self.total_exp):
                print 'WARNING: expname %s(%s) is already exists' % (exp.name, fp)
            self.total_exp[exp.name] = exp

        # make hierarcical list from loaded experiment
        for _, exp in self.total_exp.items():
            if (exp.parent is None):
                self.list_exp.append(exp)
            else:
                parent_name = self.parent
                if (parent_name not in self.total_exp):
                    print 'WARNING: expname %s has parent %s, but parent not exists.'\
                        % (exp.name, parent_name)
                    continue
                self.total_exp[parent_name].children.append(exp)
