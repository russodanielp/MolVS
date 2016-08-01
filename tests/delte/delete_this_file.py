from molvs.curator import CurationPipeline, MixturesFilter, Neutralize, StandardizeStep, CreateRDKitMols
import pandas as pd
import numpy as np


class ReachDatabaseTestUM:

    def setup(self):
        steps = [
            CreateRDKitMols(),
            StandardizeStep(),
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


    # def test_reach_pipeline(self):
    #     import csv
    #     try:
    #         df = self.cc.runFile("tests/test_data/ReachChemicalIds_dropsmiles.csv", smiles_col=8, index_col=0,
    #                          sep=",", quotechar='"', escapechar='\\', encoding='utf-8',
    #                      quoting=csv.QUOTE_NONNUMERIC, doublequote=False)
    #     except Exception as e:
    #         print(str(e))
    #     #print(df)
    #     assert False

    def test_reach_pipeline(self):
        import csv
        df = pd.read_csv("tests/test_data/ReachChemicalIds_dropsmiles.csv",
                         sep=",", quotechar='"', escapechar='\\', encoding='utf-8',
                         quoting=csv.QUOTE_NONNUMERIC, doublequote=False)
        self.cc.run(df, 9)
        assert df.empty

