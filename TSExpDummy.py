class TSExpDummy(object):
    def __init__(self, exp):
        self.exp = exp
        exp.name = "_dummy"
        exp.desc = "do nothing!"
    
    def run(self):
        # add dummy columns
        self.exp.AddColumn('abc', 'text')
        self.exp.AddColumn('def', 'text')
        # add dummy rows
        self.exp.AddRow('item-1', ['s','e'])
        self.exp.AddRow('item-2', ['s','e'])
        self.exp.AddRow('item-3', ['s','e'])
