from molvs.curator import CurationPipeline, MixturesFilter, Neutralize
import pandas as pd
import numpy as np


class ReachDatabaseTestUM:

    def setup(self):
        steps = [
            #MixturesFilter(),
            Neutralize(),
            #TautomerCheck(),
            #DuplicatesFilter(1)
             ]
        self.cc = CurationPipeline(steps=steps)


    def test_get_file(self):
        import csv
        df = pd.read_csv("tests/test_data/ReachChemicalIds_dropsmiles.csv",
                  sep=",", quotechar='"', escapechar='\\', encoding='utf-8',
                  quoting=csv.QUOTE_NONNUMERIC, doublequote=False)

        assert not df.empty


    def test_reach_pipeline(self):
        import csv
        df = self.cc.runFile("tests/test_data/ReachChemicalIds_dropsmiles.csv", smiles_col=8, index_col=0,
                         sep=",", quotechar='"', escapechar='\\', encoding='utf-8',
                         quoting=csv.QUOTE_NONNUMERIC, doublequote=False)
        print(df)
        assert not df.empty


