from unittest import TestCase
from vennproteomics.analysis import Analysis
import matplotlib.pyplot as plt
import pandas as pd
class TestAnalysis(TestCase):
    def test_add_data(self):
        vp = Analysis()
        vp.add_data(r"C:\Users\Toan Phung\Downloads\GolgiIP_HEK293_Alessi Lab.txt", gene_name_column="Gene.names", has_gene_name=True, label="Experimental")
        vp.add_data(r"C:\Users\Toan Phung\Downloads\Golgi-IP_Wade Harper.txt", gene_name_column="FC>1", label="FC>1", has_gene_name=True)
        vp.add_data(r"C:\Users\Toan Phung\Downloads\Golgi-IP_Wade Harper.txt", gene_name_column="FC>1", label="FC<1",
                    has_gene_name=True)


    def test_create_venn_diagram(self):
        vp = Analysis()
        data = pd.read_csv(r"C:\Users\Toan Phung\Downloads\GolgiIP_HEK293_Alessi Lab.txt", sep="\t")
        print(data.columns)
        vp.add_data(data[(data["Significant GolgiTAG-IP/Control-IP"]=="+") & (data["Fold enrichment (Log2): GolgiTAG-IP/Control-IP"]>=2)], gene_name_column="Gene.names",
                    has_gene_name=True, label="Significant GolgiTAG-IP/Control-IP")
        vp.add_data(r"C:\Users\Toan Phung\Downloads\Golgi-IP_Wade Harper.txt", gene_name_column="FC>1", label="FC>1",
                    has_gene_name=True)
        vp.add_data(r"C:\Users\Toan Phung\Downloads\Golgi-IP_Wade Harper.txt", gene_name_column="FC<1", label="FC<1",
                    has_gene_name=True)
        fig, ax = plt.subplots()
        vp.create_venn_diagram(ax)
        plt.tight_layout()
        fig.savefig("result.svg")