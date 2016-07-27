from molvs.curator import *
import pandas as pd
import numpy as np


class PipelineTestUM:

    def setup(self):
        self.cc = CurationPipeline(steps=[])
        self.Mols =  [
            ['[Na+].[Cl-]', 1], # NaCl - inorganic molecule
            ['[NH4+].[SH-]', 1],# Ammonium Bisulfide - Inorganic molecule
            ['[C-]#[O+].C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3.C1=CC=C(C=C1)'  # Tris(triphenylphosphine)rhodium carbonyl hydride - organometallic
             'P(C2=CC=CC=C2)C3=CC=CC=C3.C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3.[Rh]', 0],
            ['C1=C[C-]=CC=C[C-]=C1.C1=C[C-]=CC=C[C-]=C1.[U+4]', 0], # Uranocene - organometallic
            ['CC(=O)OC1=CC=CC=C1C(=O)O.CC(=O)NC1=CC=C(C=C1)O', 0],  # Aspirin/Acetaminophen - mixture
            ['C1CN(CCN1CCOCCO)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl.Cl.Cl', 1], # Hydroxyzine Dihydrochloride - salt
            ['CC(=O)OC1=CC=CC=C1C(=O)O', 1],# Aspirin - organic compound
            ['CC(=O)NC1=CC=C(C=C1)O', 0]# Acetaminophen - organic compound
        ]
        self.pre_neut_mols = ["c1cccc[nH+]1",
                              "C[N+](C)(C)C",
                              "c1ccccc1[NH3+]",
                              "CC(=O)[O-]",
                              "c1ccccc1[O-]",
                              "CCS",
                              "C[N-]S(=O)(=O)C",
                              "C[N-]C=C",
                              "C[N-]N=C",
                              "c1ccc[n-]1",
                              "CC[N-]C(=O)CC"]


        self.neut_mols = ["c1ncccc1",
                          "C[N+](C)(C)C",
                          "Nc1ccccc1",
                          "CC(O)=O",
                          "Oc1ccccc1",
                          "CCS",
                          "CNS(=O)(C)=O",
                          "CNC=C",
                          "CNN=C",
                          "c1[nH]ccc1",
                          "CCNC(CC)=O"]

        self.Mols = [[Chem.MolFromSmiles(row[0]), row[1]] for row in self.Mols]
        self.df = pd.DataFrame(self.Mols)

        self.pre_neutralized = pd.DataFrame([Chem.MolFromSmiles(mol) for mol in self.pre_neut_mols])
        self.neutralized = pd.DataFrame([Chem.MolFromSmiles(mol) for mol in self.neut_mols])

        self.inorganic_salts = self.df.iloc[:2, :]
        self.organometallics = self.df.iloc[2:4]
        self.mixtures =  self.df.iloc[4:6]
        self.controls = self.df.iloc[6:, :]



        self.mol_idx = 0
    # self.testFile = 'activity_file.txt'



    def test_mixtures_filter(self):
        """ Testing class MixturesFilter"""


        df = pd.concat([self.mixtures, self.controls])

        self.cc = self.cc.addStep(MixturesFilter())
        testMols = self.cc.run(df, self.mol_idx)

        # Hydroxyzine Dihydrochloride is a salt/mixture that will have a large organic component to be kepy
        # this is that component that needs to be added to the controls
        hydroxy_dichlor = pd.DataFrame([[Chem.MolFromSmiles(mol), 1] for mol in ['C1CN(CCN1CCOCCO)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl']])
        other_compound = pd.DataFrame([[Chem.MolFromSmiles(mol), 0] for mol in ['CC(=O)Oc1ccccc1C(=O)O']])
        self.controls = pd.concat([other_compound, hydroxy_dichlor, self.controls])


        testMols.iloc[:, self.mol_idx] = testMols.iloc[:, self.mol_idx].map(Chem.MolToSmiles)
        self.controls.iloc[:, self.mol_idx] = self.controls.iloc[:, self.mol_idx].map(Chem.MolToSmiles)
        print(testMols.values, self.controls.values)
        assert np.array_equal(testMols.values, self.controls.values)

    def test_duplicates_filter(self):
        """ Testing class DuplicatesFilter"""
        df = pd.concat([self.controls, self.controls])

        self.cc = self.cc.addStep(DuplicatesFilter(1))
        testMols = self.cc.run(df, 0)

        testMols.iloc[:, self.mol_idx] = testMols.iloc[:, self.mol_idx].map(Chem.MolToSmiles)
        self.controls.iloc[:, self.mol_idx] = self.controls.iloc[:, self.mol_idx].map(Chem.MolToSmiles)

        testMols.sort_index(inplace=True)
        self.controls.sort_index(inplace=True)

        assert np.array_equal(testMols.values, self.controls.values)

    def test_neutralize(self):
        """ Testing class Neutralize"""

        self.cc = self.cc.addStep(Neutralize())
        testMols = self.cc.run(self.pre_neutralized, 0)

        testMols.iloc[:, self.mol_idx] = testMols.iloc[:, self.mol_idx].apply(Chem.MolToSmiles)
        self.neutralized.iloc[:, self.mol_idx] = self.neutralized.iloc[:, self.mol_idx].apply(Chem.MolToSmiles)

        assert np.array_equal(testMols.values, self.neutralized.values)

    def full_test_from_dataframe(self):
        """ Testing full pipeline from DataFrame """
        steps = [
            MixturesFilter(),
            Neutralize(),
            TautomerCheck(),
            DuplicatesFilter(1)]
        pipeline = CurationPipeline(steps=steps)
        df = pipeline.run(self.df, 0)
        assert not df.empty

    def full_test_from_file(self):
        """ Testing full pipeline from file """
        steps = [
            MixturesFilter(),
            Neutralize(),
            TautomerCheck(),
            DuplicatesFilter(1)]
        pipeline = CurationPipeline(steps=steps)
        df = pipeline.runFile(r'tests/test_data/database_file.csv', smiles_col=0)
        assert not df.empty



