"""
Author: Daniel P. Russo
"""
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import Chem
from rdkit.Chem import AllChem
from molvs.tautomer import TautomerCanonicalizer
from molvs.charge import Uncharger
from molvs.fragment import LargestFragmentChooser
import pandas as pd
import re


class CurationStep:

    def runStep(self, df, cmp_index):
        raise NotImplementedError("Please implement this method")



class MixturesFilter(CurationStep):
    """class to remove mixtures"""
    def __init__(self, **kwargs):
        self.chooser = LargestFragmentChooser(**kwargs)

    def filterFunction(self, cmp):
        return self.chooser.choose(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) for mol in df[cmp_index]]
        return df.dropna(subset=[cmp_index])


class Neutralize(CurationStep):
    """ Class to neutralize compounds  """
    def __init__(self, **kwargs):
        self.neutralizer = Uncharger(**kwargs)

    def filterFunction(self, cmp):
        return self.neutralizer.uncharge(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) for mol in df[cmp_index]]
        return df

class TautomerCheck(CurationStep):
    """ class to standaridze tautomers"""
    def __init__(self, **kwargs):
        self.tautomerizer = TautomerCanonicalizer(**kwargs)

    def filterFunction(self, cmp):
        return self.tautomerizer.canonicalize(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) for mol in df[cmp_index]]
        return df

class DuplicatesFilter(CurationStep):
    """ Class to eliminate duplicates """
    def __init__(self, cmp_index, act_index, take='highest'):
        self._cmp_index = cmp_index
        self._act_index = act_index
        if take == 'highest':
            self._keep = 'last'
        elif take == 'lowest':
            self._keep = 'first'

    def filterFunction(self, df):
        df['smiles'] = [Chem.MolToSmiles(mol) for mol in df[self._cmp_index]]
        df.sort_values(['smiles', self._act_index], inplace=True)
        df.drop_duplicates(subset=['smiles'], keep=self._keep, inplace=True)
        df.drop('smiles', axis=1, inplace=True)
        return df

    def runStep(self, df):
        return self.filterFunction(df)


CURATION_STEPS=(
    MixturesFilter(),
    Neutralize(),
    TautomerCheck(),
    DuplicatesFilter()
)

class CurationPipeline:
    """class for creating a MolVS curation pipeline"""
    def __init__(self, steps=[]):

        """Initialize a CurationPipeline with a custom list of :class:`~molvs.curator.CurationStep`.

        :param steps: A list of list of :class:`~molvs.curator.CurationStep`.
        """
        self._steps = steps

    def addStep(self, curationStep):
        return CurationPipeline(self._steps + [curationStep])

    def run(self, df, cmp_index):
        new_df = df.copy()
        for step in self._steps:
            new_df = step.runStep(new_df, cmp_index)
        return new_df

    def runFile(self, filename, delimiter='\t'):
        """
        :param filename: name of activity file
        :param delimiter: delimiter of the file
        :return: None
        """
        df = pd.read_csv(filename, delimiter=delimiter, header=None, dtype=str)
        df.columns = ['smiles', 'activity']
        df['rdkit_mols'] = [Chem.MolFromSmiles(smi) for smi in df['smiles']]
        return self.run(df)
