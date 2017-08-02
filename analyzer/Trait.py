import scipy.stats as stats
import pandas as pd
import os

# default file lists to load
list_t2g_alias = {
    "GOTerm-rice": os.path.join(os.path.dirname(__file__), "data/GOBPname2gene.rice.txt"),
    "GOTerm-arabidopsis": os.path.join(os.path.dirname(__file__), "data/GOBPname2gene.arabidopsis.txt"),
    "GOTerm-human": os.path.join(os.path.dirname(__file__), "/data/GOBPname2gene.human.txt"),
    "GOTerm-mouse": os.path.join(os.path.dirname(__file__), "/data/GOBPname2gene.rice.txt"),
    "Pathway-arabidopsis": os.path.join(os.path.dirname(__file__), "/data/pathwayname2gene.ath.txt"),
    "Motif-arabidopsis": os.path.join(os.path.dirname(__file__), "/data/kmer2genes.k6.n1.txt")
}

# trait to genes file list
class Trait(object):
    def __init__(self):
        self.dic_trait2gene = {}
        self.dic_gene2trait = {}
        self.backgroundGenes = None
        self.dic_trait2ratio = {}   # gene-ratio of each trait
        self.method = 'fisher'      # method for statistical test

    def load(self, alias):
        if (alias not in list_t2g_alias):
            raise Exception('invalid alias')
        else:
            self.load_t2g(list_t2g_alias[alias])

    # load trait to genes file
    def load_t2g(self, fn):
        with open(fn, 'r') as f:
            for line in f:
                trait,genes = line.rstrip('\n').split('\t')
                lst_gene=genes.split(',')
                for gene in lst_gene:
                    if (self.backgroundGenes!= None and gene not in self.backgroundGenes):
                        continue
                    if trait not in self.dic_trait2gene:
                        self.dic_trait2gene[trait] = set()
                    if gene not in self.dic_gene2trait:
                        self.dic_gene2trait[gene] = set()
                    self.dic_gene2trait[gene].add(trait)
                    self.dic_trait2gene[trait].add(gene)
        for trait in self.dic_trait2gene.keys():
            self.dic_trait2ratio[trait] = float(len(self.dic_trait2gene[trait]))/len(self.dic_gene2trait)

    # trait, pvalue, (fisher 4 elements)
    def GOanalysis(self, lst_gene):
        return pd.DataFrame(
            self.GOanalysis_list(lst_gene),
            columns=('GOTerm', 'pvalue', 'occured_in_tested', 'total_tested', 'occured_in_background', 'total_background')
            )

    # trait, pvalue, (fisher 4 elements)
    def GOanalysis_list(self, lst_gene):
        dic_trait2count = {}
        set_tested_gene = set()
        for gid in lst_gene:
            if gid not in self.dic_gene2trait:
                # invalid gene name, no GOterm
                continue
            # count related traits with given genes
            for trait in self.dic_gene2trait[gid]:
                if not trait in dic_trait2count:
                    dic_trait2count[trait] = 0
                dic_trait2count[trait] += 1
            set_tested_gene.add(gid)
        lst_out = []
        for trait,count in dic_trait2count.items():
            occured_in_tested, total_tested, occured_in_background, total_background = count, len(set_tested_gene), len(self.dic_trait2gene[trait]), len(self.dic_gene2trait)
            if occured_in_tested == 0:
                pval = 1.0
            else:
                if self.method == 'binomial':
                    pval=1.0-stats.binom.cdf(occured_in_tested-1, total_tested, self.dic_trait2ratio[trait])     # ovccured_in_test-1 means p(X>=n) i.e. contain
                elif self.method == 'fisher':
                    oddratio,pval=stats.fisher_exact([[occured_in_tested, total_tested-occured_in_tested], [occured_in_background-occured_in_tested, total_background-total_tested-occured_in_background+occured_in_tested]], alternative='greater')        # 2X2 fisher's exact test
                lst_out.append([trait, pval, occured_in_tested, total_tested, occured_in_background, total_background])
        # sort with pvalue
        return sorted(lst_out, key=lambda x: x[1])

