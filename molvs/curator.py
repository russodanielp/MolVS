"""
Author: Daniel P. Russo
"""

from rdkit import Chem

# had to change imports for pychemsim project
from molvs.tautomer import TautomerCanonicalizer
from molvs.charge import Uncharger
from molvs.fragment import LargestFragmentChooser
from molvs.standardize import Standardizer
import pandas as pd
import logging


log = logging.getLogger(__name__)


# TODO Gracefully keep track of compounds lost in each step through a log explaining why each compound was eliminated

class CurationStep:

    def runStep(self, df, cmp_index):
        raise NotImplementedError("Please implement this method")


class CreateRDKitMols(CurationStep):
    """Instantiates RDKit Mols from Smiles for an entire Pandas DF column"""

    def __init__(self):
        log.debug("Initializing CreateRDKitMols")
        self._rdkitMols = []

    def stepFunction(self, cmp):
        log.debug("Running CreateRDKitMols on {0}".format(cmp))
        return Chem.MolFromSmiles(cmp)

    def runStep(self, df, cmp_index):
        # drop null values in column containing smiles
        df.dropna(subset=[df.columns[cmp_index]], inplace=True)
        df.iloc[:, cmp_index] = [self.stepFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        return df

class ReturnToSmiles(CurationStep):
    """Converts RDKit Mols to Smiles for an entire Pandas DF column"""

    def __init__(self):
        log.debug("Initializing ReturnToSmiles")
        self._rdkitMols = []

    def stepFunction(self, cmp):
        log.debug("Running ReturnToSmiles on {0}".format(cmp))
        return Chem.MolToSmiles(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.stepFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        return df


class StandardizeStep(CurationStep):
    """class to standadize all molecules before use in the pipeline"""
    def __init__(self, **kwargs):
        log.debug("Initializing StandardizeStep")
        self.standardizer = Standardizer(**kwargs)


    def filterFunction(self, cmp):
        try:
            return self.standardizer.standardize(cmp)
        except ValueError as e:
            log.error("{0} on {1}".format(e, cmp))


    def runStep(self, df, cmp_index):
        log.debug("Running StandardizeStep")
        df.iloc[:, cmp_index] = [self.filterFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        return df


class MixturesFilter(CurationStep):
    """class to remove mixtures"""
    def __init__(self, **kwargs):
        self.chooser = LargestFragmentChooser(**kwargs)

    def filterFunction(self, cmp):
        return self.chooser.choose(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        subset = df.columns[cmp_index]
        return df.dropna(subset=[subset])


class Neutralize(CurationStep):
    """ Class to neutralize compounds  """
    def __init__(self, **kwargs):
        self.neutralizer = Uncharger(**kwargs)

    def filterFunction(self, cmp):
        return self.neutralizer.uncharge(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        return df

class TautomerCheck(CurationStep):
    """ class to standaridze tautomers"""
    def __init__(self, **kwargs):
        self.tautomerizer = TautomerCanonicalizer(**kwargs)

    def filterFunction(self, cmp):
        return self.tautomerizer.canonicalize(cmp)

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        return df

class InorganicFilter(CurationStep):
    """ Filter removes any inorganic compounds, i.e., compounds not containingg Carbon"""

    def __init__(self):
        pass

    def filterFunction(self, cmp):
        C = Chem.MolFromSmarts('[#6]')
        return cmp if cmp.HasSubstructMatch(C) else None

    def runStep(self, df, cmp_index):
        df.iloc[:, cmp_index] = [self.filterFunction(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        return df

class DuplicatesFilter(CurationStep):
    """ Class to eliminate duplicates """
    def __init__(self, act_index=None, take='highest'):
        self._act_index = act_index
        if take == 'highest':
            self._keep = 'last'
        elif take == 'lowest':
            self._keep = 'first'

    def filterFunction(self, df, cmp_index):
        # TODO Fix this SettingWithCopyWarning
        df.loc[:, 'smiles'] = [Chem.MolToSmiles(mol) if mol else None for mol in df.iloc[:, cmp_index]]
        # if there is an activity column, sort the values so the appropriate activity is kept
        if self._act_index:
            df.sort_values(['smiles', self._act_index], inplace=True)
        df.drop_duplicates(subset=['smiles'], keep=self._keep, inplace=True)
        df.drop('smiles', axis=1, inplace=True)
        return df

    def runStep(self, df, cmp_index):
        return self.filterFunction(df, cmp_index)



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

    def runFile(self, filename, smiles_col=0, **kwargs):
        """
        Reads database file performs specified curation steps.

        :param filename: name of activity file
        :param smiles_col: column index of smiles

        :param **kwargs: ::ref:: Pandas.read_csv()

        :return: Curated database in Pandas DataFrame
        """
        log.debug("CurationPipeline.runFile()")
        df = pd.read_csv(filename, **kwargs)
        df.columns = [i for i in range(len(df.columns))]
        df[smiles_col] = [Chem.MolFromSmiles(smi) for smi in df[smiles_col]]
        return self.run(df, smiles_col)
