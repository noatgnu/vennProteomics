import pandas as pd
from matplotlib_venn import venn2, venn3
from uniprotparser.betaparser import UniprotParser, UniprotSequence
from io import StringIO

class Analysis:
    def __init__(self):
        self.data = []
        self.labels = []
        pass

    def add_data(self, filename, accession_id_column, label="", gene_name_column="Gene names", has_gene_name=False):
        df = pd.read_csv(filename, sep="\t")
        if not has_gene_name:
            parser = UniprotParser()
            acc_set = set()
            for accession in df[accession_id_column]:
                for acc in accession.split(";"):
                    uni = UniprotSequence(acc, parse_acc=True)
                    acc_set.add(uni.accession)
            data = []
            for d in parser.parse(acc_set):
                data.append(pd.read_csv(StringIO(d.text), sep="\t"))
            if len(data) > 0:
                if len(data) > 1:
                    data = pd.concat(data, ignore_index=True)
                else:
                    data = data[0]
            self.data.append(list(data["Gene Names"].unique()))
        else:
            self.data.append(list(df[gene_name_column].unique()))
        if label == "":
            self.labels.append(f"Dataset #{len(self.data)}")
        else:
            self.labels.append(label)

    def create_venn_diagram(self):
        if 1 > len(self.data) < 4:
            if len(self.data) == 2:
                venn2(self.data, set_labels=self.labels)
            else:
                venn3(self.data, set_labels=self.labels)


if __name__ == "__main__":
    venn2()