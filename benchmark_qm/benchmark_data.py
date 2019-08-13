from yaml import load, dump
from yaml import CLoader as Loader, CDumper as Dumper
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem
from warnings import warn
import numpy as np
import pybel as pyb
import openbabel as ob
import os
import rdkit


class BenchmarkEntry:
    """
    A class for storing info related to a molecule for QM benchmarking
    """
    def __init__(self, name=None, smiles=None, inchi=None, index=None, categories=None, multiplicity=None, charge=None,
                 rmg_symmetry=None, expt_sources=None, preferred_expt_source=None, ff_xyz=None, dft_xyzs=None,
                 dft_freqs=None, yaml_file=None, qm_files=None):
        """

        :param name: Name of the species
        :param smiles: SMILES string of the species
        :param inchi: InChi string of molecule
        :param index: Species index in the benchmarking set
        :param categories: list of categories (e.g. aromatics) that the species belongs to
        :param multiplicity: Spin multiplicity (2S+1)
        :param charge: The total net charge on the molecule
        :param rmg_symmetry: The symmetry number of the molecule as calculated by rmg
        :param expt_sources: A list of sources with experimental source objects
        :param preferred_expt_source: The chosen experimental data used for benchmarking this species
        :param ff_xyz: The xyz coordinates of the lowest energy conformer from force fields (as a string)
        :param dft_xyzs: A dictionary with the DFT method string as keys and xyz geometry as values
        :param dft_freqs: A dictionary with the DFT method string as keys and frequencies as values
        :param yaml_file: The relative location of the yaml file (used to store the data of this object)
        :param qm_files: A dictionary for the mapping {str(file_description):str(file_path)}
        """
        self.name = name
        self.smiles = smiles
        self.inchi = inchi
        self.index = index
        self.categories = categories or []
        self.multiplicity = multiplicity
        self.charge = charge
        self.rmg_symmetry = rmg_symmetry
        self.expt_sources = expt_sources or []
        self.preferred_expt_source = preferred_expt_source
        self.ff_xyz = ff_xyz
        self.dft_xyzs = dft_xyzs or {}
        self.dft_freqs = dft_freqs or {}
        self.yaml_file = yaml_file
        self.qm_files = qm_files or {}

    def save_to_yaml(self, path=None):
        """
        Save the benchmark entry to a YAML file
        :param path: The relative location of the YAML file (if not specified self.yaml_file is used)
        :return: None
        """
        if path is None:
            if self.yaml_file is None:
                raise ValueError('YAML path not specified for BenchmarkEntry {0}'.format(self.name))
            else:
                path = self.yaml_file

        # Make parent directories if they don't exist
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        with open(path, 'w') as f:
            # Dump the species name at the top of the file
            f.write('name: {}\n'.format(self.name))

            # Dump the remaining attributes
            remaining_attributes = deepcopy(self.__dict__)
            _ = remaining_attributes.pop('name')
            dump(remaining_attributes, f, Dumper=Dumper)

    def load_from_yaml(self, path=None):
        if path is None:
            if self.yaml_file is None:
                raise ValueError('YAML path not specified for BenchmarkEntry {0}'.format(self.name))
            else:
                path = self.yaml_file

        with open(path, 'r') as f:
            attributes = load(f, Loader=Loader)
            attributes['yaml_file'] = path

        self.__init__(**attributes)

    def generate_ff_xyz(self, nconf=None):
        #  Try using RDKit first
        try:
            mol = smiles_to_rdkit(self.smiles, nconf=nconf)
            conformer = mol.GetConformers()[0]
            coords = conformer.GetPositions()
            elements = [atom.GetSymbol() for atom in mol.GetAtoms()]

        except:
            print 'RDKit failed to generate 3D geometry for {0}, trying OpenBabel'.format(self.smiles)
            pybmol = smiles_to_openbabel(self.smiles, nconf=nconf)
            coords = [list(atom.coords) for atom in pybmol.atoms]
            elements = [atomic_symbol[atom.atomicnum] for atom in pybmol.atoms]

        xyz_body_lines = ['{0}  {1[0]: .10f}  {1[1]: .10f}  {1[2]: .10f}'.format(symbol, coord) for symbol, coord in
                          zip(elements, coords)]
        xyz_body = '\n'.join(xyz_body_lines)
        xyz_header = '{0}\n\n'.format(len(elements))
        xyz_string = ''.join([xyz_header, xyz_body, '\n'])

        self.ff_xyz = xyz_string

    def generate_qchem_file(self, template_file):
        mol_xyz = self.ff_xyz
        mol_xyz = mol_xyz.split('\n')[1:]  # Remove header with num_atoms
        mol_xyz[0] = '$molecule\n{0} {1}'.format(self.charge, self.multiplicity)
        mol_xyz = '\n'.join(mol_xyz)
        mol_xyz += '$end\n\n'

        with open(template_file, 'r') as f:
            file_text = f.readlines()

        file_text = mol_xyz + ''.join(file_text)

        method_name = os.path.split(template_file)[-1].split('.')[0]
        file_name = 'spc_{0}_{1}.inp'.format(self.index, method_name)

        with open(os.path.join(os.path.dirname(self.yaml_file), file_name), 'w') as f:
            f.write(file_text)


def smiles_to_rdkit(smi, gen_3d=True, nconf=None):
    """
    Convert smiles to RDKit molecule.
    Tries to generate the lowest-energy conformer.
    Original Author: Colin Grambow
    Modified by Mark Payne

    :param smi: Smiles string
    :param gen_3d: Generate 3D coordinates and perform a conformer search
    :param nconf: number of conformers to generate (default n_atoms*10)
    :return mol: RDkit molecule object
    """
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)

    if nconf is None:
        nconf = mol.GetNumAtoms()*10

    if gen_3d:
        cids = AllChem.EmbedMultipleConfs(mol, nconf, AllChem.ETKDG())

        AllChem.MMFFSanitizeMolecule(mol)
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)

        energies = []
        use_uff = False
        if AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=0) is None:
            #  Use UFF instead
            print 'MMFF not defined for {0}, switching to UFF instead'.format(smi)
            use_uff = True
            AllChem.UFFOptimizeMoleculeConfs(mol)

        for cid in cids:
            if use_uff:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            else:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
                ff.Minimize()

            energy = ff.CalcEnergy()
            energies.append(energy)

        energies = np.asarray(energies)
        min_energy_idx = np.argsort(energies)[0]

        new_mol = Chem.Mol(mol)
        new_mol.RemoveAllConformers()
        min_conf = mol.GetConformer(cids[min_energy_idx])
        new_mol.AddConformer(min_conf, assignId=True)
        mol = new_mol

    return mol


def smiles_to_openbabel(smi, nconf=None, geom_steps=500):
    """
    Uses OpenBabel to automatically return the geometry of the lowest energy conformer using a weighted
    rotor search with the mmff94s force field.
    Original Author(s): Colin Grambow, Matt Johnson, Ryan Gillis
    Modified by: Mark Payne

    :param smi: Smiles string
    :param nconf: The number of conformers to generate (default n_atoms*10)
    :param geom_steps: The number of optimization steps for each conformer
    :return pybmol: pybel mol object set to the lowest energy conformer
    """
    energies = []
    pybmol = pyb.readstring("smi", smi)
    pybmol.make3D()
    obmol = pybmol.OBMol

    if nconf is None:
        n_atoms = len(pybmol.atoms)
        nconf = n_atoms*10

    ff = ob.OBForceField.FindForceField("mmff94s")
    ff.Setup(obmol)
    ff.WeightedRotorSearch(nconf, geom_steps)
    ff.GetConformers(obmol)
    for n in xrange(obmol.NumConformers()):
        obmol.SetConformer(n)
        ff.Setup(obmol)
        energies.append(ff.Energy())

    lowest_energy_index = np.argmin(energies)
    obmol.SetConformer(lowest_energy_index)

    return pybmol


try:
    atomic_symbol = {i: Chem.Atom(i).GetSymbol() for i in range(1, 117)}
except RuntimeError:
    rdkit_version = rdkit.__version__
    if int(rdkit_version.split('.')[0]) < 2018:
        warn(Warning('RDKit version {0} detected. BenchmarkEntry objects require RDKit 2018 or later to generate force'
                     'field geometries'.format(rdkit_version)))


if __name__ == '__main__':
    pass
