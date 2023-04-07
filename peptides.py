#!/usr/bin/env python

"""
Data structures and utilities related to peptides
"""

import logging
from pyteomics import parser as pyparser
from pyteomics import mass as pyteomics_mass
from pyteomics import achrom
from pyteomics import electrochem
from pyteomics import mass as massUtil
from numpy import median
import copy

import re


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)

all_amino_acids = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
    'P', 'S', 'T', 'W', 'Y', 'V', 'U'
]

# mass of a hydrogen atom
HYDROGEN_MASS = 1.00794

# dictionary from amino acid symbols to positions in all_amino_acids
all_aa_pos_dict = dict()
for i in range(0, len(all_amino_acids)):
    all_aa_pos_dict[all_amino_acids[i]] = i

# constraints on tryptic peptide length for fragment ion / database search purposes.
# These are subject to debate. 7 is probably a pretty hard minimum, but analysis of Chip 7
# MS data indicated that we're less likely to see knottins that we make with no unique
# peptides length 8 or longer. Max could be anywhere in the 25-35 range, probably.
DEFAULT_MIN_PEPTIDE_LENGTH = 8
DEFAULT_MAX_PEPTIDE_LENGTH = 30

# constraints on tryptic peptide mass for purposes of seeing it in the mass spec
DEFAULT_MIN_PEPTIDE_MASS = 800
DEFAULT_MAX_PEPTIDE_MASS = 3100

# minimum calculated charge for a peptide we want to observe by MS. 0.0 means
# no minimum.
DEFAULT_MIN_PEPTIDE_CHARGE = 0.0

# maximum number of M residues in a tryptic peptide for MS observation
DEFAULT_MAX_PEPTIDE_METHIONINES = 1
# maximum number of W residues in a tryptic peptide for MS observation
DEFAULT_MAX_PEPTIDE_TRYPTOPHANS = 1

# From: Engelman, D.M., Steitz, T.A., and Goldman, A. (1986). 
# Identifying nonpolar transbilayer helices in amino acid sequences of 
# membrane proteins. Ann. Rev. Biophys. Biophys. Chem. 15, 321-353
AA_HYDROPHOBICITIES = {'F': 3.7, 'M': 3.4, 'I': 3.1, 'L': 2.8, 'V': 2.6,
                       'C': 2.0, 'W': 1.9, 'A': 1.6, 'T': 1.2, 'G': 1.0,
                       'S': 0.6, 'P': -0.2, 'Y': -0.7, 'H': -3.0, 'Q': -4.1,
                       'N': -4.8, 'E': -8.2, 'K': -8.8, 'D': -9.2, 'R': -12.3}

AA_UNMOD_MASSES = {
    'A': 71.03711,
    'C': 103.00919, # note no iodoacetamide
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406,
    'M': 131.04049,
    'N': 114.04293,
    'P': 97.05276,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333,
    'U': 150.95363	
}

# Table 4 Measured Molar Extinction Coefficients (214 nm) of protein and peptides from
# J. Agric Food Chem. Vol. 55. No. 14, 2007
EXTINCTION_COEFFICIENT = {'P': 2675, 'H': 5125, 'F': 5200, 'Y': 5375, 'W': 29050, 'M': 980, 'R': 102, 'N': 136,
                          'Q': 142, 'C': 225, 'G': 21, 'A': 32, 'S': 34, 'K': 41, 'T': 41, 'V': 43, 'I': 45,
                          'L': 45, 'D': 58, 'E': 78}


def make_masstable_with_mods(aa_modifications):
    """
    Make an amino acid mass table with the indicated modifications
    :param aa_modifications: 
    :return: 
    """
    result = AA_UNMOD_MASSES.copy()
    for aa_mod in aa_modifications:
        result[aa_mod.aa] += aa_mod.massdiff
    return result


def make_masstable_with_iodo_c():
    """
    Make an amino acid mass table with iodoacetamide on C
    :return: 
    """
    return make_masstable_with_mods([MODIFICATION_IODOACETAMIDE_STATIC])


def calc_peptide_rt(peptide_sequence):
    """
    Calculate the predicted retention time of a peptide  
    :param peptide_sequence: 
    :return: 
    """
    return achrom.calculate_RT(peptide_sequence, achrom.RCs_guo_ph7_0)


def calc_peptide_extinction_coeff(peptide_sequence):
    """Calculate the predicted extinction coefficient of a peptide"""
    #peptide bond molar extinction coefficient
    peptide_bond = 923

    # EQUATION for e = extinction coefficient :
    #e  = (923)*(n-1)+(frequency of AA i)*(e value of AA i)

    # where (n-1) is the number of peptide bonds (i.e. number of letters in the AA sequence string minus 1)
    # will need a check for N-terminal proline

    absorbances = []
    for index, AA in enumerate(peptide_sequence):
        if index == 0 and AA == 'P':
            absorbances.append(30)
        else:
            absorbances.append(EXTINCTION_COEFFICIENT[AA])
    return peptide_bond*(len(peptide_sequence)-1) + sum(absorbances)


def calc_mean_hydrophobicity(sequence):
    """Calculate the mean AA hydrophobicity for a peptide. This will
    allow direct comparison of peptides of the same length. Not sure yet
    how to compare peptides of different lengths"""
    sum_hydrophobicity = 0
    for aa in list(sequence):
        sum_hydrophobicity += AA_HYDROPHOBICITIES[aa]
    return sum_hydrophobicity / len(sequence)


def is_cleavagepep_good_for_ms(pepseq,
                               min_peptide_charge=DEFAULT_MIN_PEPTIDE_CHARGE,
                               max_cleavage_methionines=DEFAULT_MAX_PEPTIDE_METHIONINES,
                               max_cleavage_tryptophans=DEFAULT_MAX_PEPTIDE_TRYPTOPHANS):
    """Given a variant sequence and a base sequence, find all tryptic peptides
    of the variant that:
    1. satisfy minimum predicted charge constraints
    2. satisfy amino acid constraints
    Assumes mass and length have already been checked."""
    # check charge
    if min_peptide_charge > 0 and predict_charge(pepseq) < min_peptide_charge:
        return False
    # check # Ms
    if pepseq.count('M') > max_cleavage_methionines:
        return False
    # check # Ws
    if pepseq.count('W') > max_cleavage_tryptophans:
        return False
    return True


def calc_theoretical_peak_mzs(peptide_sequence, charges, position_massdiffs, min_mz, max_mz,
                              nterm_deltamass=0.0, cterm_deltamass=0.0,
                              ion_types=('b', 'y')):
    """
    Calculate the theoretical m/z values of all b- and y-ions for the given peptide sequence,
    in the fragment charge states specified, with the modifications specified, that occur between
    min_mz and max_mz
    Todo: make ion types a parameter
    Todo: support variable modifications. But that would mean multiple results
    Todo: make the main loop more efficient, if anyone cares.
    :param peptide_sequence: 
    :param charges: 
    :param position_massdiffs: a list of delta masses for each sequence position. Must be the same
    length as peptide_sequence.
    :param min_mz: 
    :param max_mz: 
    :param nterm_deltamass: delta mass for the N terminus
    :param cterm_deltamass: delta mass for the C terminus
    :param ion_types: a list of ion types to include. Currently only supports 'b' and 'y'
    :return: a list of fragment m/zs
    """
    if not ion_types:
        raise ValueError("No ion types supplied")
    if len(peptide_sequence) != len(position_massdiffs):
        raise ValueError("peptide_sequence has length %d; position_massdiffs has length %d. Must be same." %
                         (len(peptide_sequence), len(position_massdiffs)))
    for ion_type in ion_types:
        if ion_type not in ('b', 'y'):
            raise ValueError('calc_theoretical_peak_mzs got unknown ion type "%s". Only "b" and "y" supported.' %
                             ion_type)
    peak_mzs = []

    # loop on charges. This would be far more efficient if I looped on amino acids, instead, and
    # calculated the mass. I haven't done that because the first line of the for loop, which affects
    # all ions, uses charge. It wouldn't be that hard, but so far I don't need the performance boost. *shrug*
    for charge in charges:
        nterm_mass = HYDROGEN_MASS * charge + nterm_deltamass
        b_ion_mass = nterm_mass
        # walk forward through the sequence, adding mass to b_ion_mass for each amino acid
        for position in range(0, len(peptide_sequence) - 1):
            b_ion_mass += AA_UNMOD_MASSES[peptide_sequence[position]]
            b_ion_mass += position_massdiffs[position]
            b_ion_mz = b_ion_mass / charge
            if min_mz <= b_ion_mz <= max_mz and 'b' in ion_types:
                peak_mzs.append(b_ion_mz)
        # calculate the mass of the last amino acid, which isn't part of b_ion_mass
        last_aa_mass = AA_UNMOD_MASSES[peptide_sequence[-1]]
        last_aa_mass += position_massdiffs[len(peptide_sequence) - 1]
        # start y_ion_mass out with the full mass of the peptide.
        # take the highest b ion mass, add the mass of the C terminus, and subtract the mass of the N terminus
        cterm_mass =  2 * HYDROGEN_MASS + AminoacidModification.MOD_OXIDATION_MASSDIFF + cterm_deltamass
        y_ion_mass = b_ion_mass + last_aa_mass + cterm_mass
        # walk forward to calculate all the y-ion masses and mzs
        for position in range(0, len(peptide_sequence) - 1):
            y_ion_mass -= AA_UNMOD_MASSES[peptide_sequence[position]]
            y_ion_mass -= position_massdiffs[position]
            y_ion_mz = y_ion_mass / charge
            if min_mz <= y_ion_mz <= max_mz and 'y' in ion_types:
                peak_mzs.append(y_ion_mz)
    return peak_mzs


def predict_charge(sequence):
    """Predict the charge we're most likely to view this peptide in. """
    # This has to be extremely wacky, because algorithmic prediction works
    # terribly
    #todo: hone this based on more experimental results
    #todo: is pH 2?

    pep_mass = massUtil.fast_mass(sequence)
    most_likely_charge = electrochem.charge(sequence, pH=2)

    if most_likely_charge > 7:
        most_likely_charge = 7
    if pep_mass < 3000:
        most_likely_charge = min(most_likely_charge, 4)
    if pep_mass < 4800:
        most_likely_charge = min(most_likely_charge, 6)
    if pep_mass > 3000:
        most_likely_charge = max(most_likely_charge, 4)
    if pep_mass > 3500:
        most_likely_charge = max(most_likely_charge, 5)
    if pep_mass > 5000:
        most_likely_charge = max(most_likely_charge, 6)

    return most_likely_charge


def get_run_sequences(run, min_pprophet=None):
    """load the *set* of peptide identifications from a run,
    with optional PeptideProphet cutoff"""

    result = set()
    for peptide_id in run.peptide_ids:
        if not min_pprophet or peptide_id.pprophet >= min_pprophet:
            result.add(peptide_id.sequence)

    return result


def get_analysis_sequences(analysis, min_pprophet=None):
    result = set()
    for run in analysis.runs:
        result = result.update(get_run_sequences(run, min_pprophet=min_pprophet))
    return result


class MSMSPipelineAnalysis:
    """Stores what we need from a PepXML file"""

    def __init__(self):
        self.runs = list()


class MSRunSearchResult:
    """Stores what we need from a set of peptide search results"""

    def __init__(self, name=None, search_engine=None, fasta=None,
                 max_missed_cleavages=0, min_tryptic_termini=2,
                 modifications=None):
        self.peptide_ids = list()
        self.name = name
        self.search_engine = search_engine
        self.fasta = fasta
        self.max_missed_cleavages = max_missed_cleavages
        self.min_tryptic_termini = min_tryptic_termini
        if modifications is None:
            modifications = []
        self.modifications = modifications

    def get_peptide_sequences(self):
        """return the *set* of peptide sequences in this run"""
        result = set()
        for peptide_id in self.peptide_ids:
            result.add(peptide_id.sequence)
        return result


    def calc_peptide_median_times(self):
        """return a map from peptide sequences to median retention time"""
        seq_id_map = self.get_peptide_ids_map()
        result = dict()
        for peptide in seq_id_map:
            result[peptide] = median([myid.time for myid in seq_id_map[peptide]])
        return result

    def get_peptide_ids_map(self):
        """return a map from peptide sequences to all their IDs"""
        result = dict()
        for peptide_id in self.peptide_ids:
            if peptide_id.sequence not in result:
                result[peptide_id.sequence] = list()
            result[peptide_id.sequence].append(peptide_id)
        return result

    def get_peptide_proteins_map(self):
        """return a map from peptide sequences to sets of all proteins
        associated with them"""
        result = dict()
        for peptide_id in self.peptide_ids:
            if peptide_id.sequence not in result:
                result[peptide_id.sequence] = set()
            for protein in peptide_id.proteins:
                result[peptide_id.sequence].add(protein)
        return result

    def get_pepids_for_sequence(self, peptide_sequence):
        allidsmap = self.get_peptide_ids_map()
        if not peptide_sequence in allidsmap:
            return None
        return allidsmap[peptide_sequence]

    def get_peptide_spectralcount_map(self):
        """return a map from peptide sequences to spectral counts"""
        result = dict()
        for peptide_id in self.peptide_ids:
            if peptide_id.sequence not in result:
                result[peptide_id.sequence] = 0
            result[peptide_id.sequence] += 1
        return result

    def remove_nonunique_peptide_ids(self):
        new_peptide_ids = []
        for peptide_id in self.peptide_ids:
            if peptide_id.proteins and len(peptide_id.proteins) == 1:
                new_peptide_ids.append(peptide_id)
        self.peptide_ids = new_peptide_ids


class PeptideIdentification:
    """Stores what we need from a pepXML peptide identification.
          mass: theoretical mass
          observed mass: the mass actually seen by the MS
          delta_mass is observed - theoretical mass"""

    def __init__(self, scan, time, sequence, probability, mass, charge,
                 proteins, observed_mass=None, prev_aa=None,
                 next_aa=None, modifications=None, spectrum_name=None,
                 num_tol_term=None, ratio_heavy_light=None,
                 quant_heavy_area=None, quant_light_area=None,
                 quant_labelfree_peakintensity=None,
                 quant_labelfree_peakarea=None,
                 quant_labelfree_peak_rt_seconds=None,
                 quant_labelfree_start_rt_seconds=None,
                 quant_labelfree_end_rt_seconds=None):
        self.scan = int(scan)
        self.time = time
        self.sequence = sequence
        self.probability = probability
        self.mass = mass
        self.charge = charge
        self.proteins = proteins
        self.observed_mass = observed_mass
        self.prev_aa = prev_aa
        self.next_aa = next_aa
        if modifications is None:
            modifications = []
        self.modifications = modifications
        self.spectrum_name = spectrum_name
        self.num_tol_term = num_tol_term
        self.ratio_heavy_light = ratio_heavy_light
        self.quant_heavy_area = quant_heavy_area
        self.quant_light_area = quant_light_area
        self.quant_labelfree_peakintensity = quant_labelfree_peakintensity
        self.quant_labelfree_peakarea = quant_labelfree_peakarea
        self.quant_labelfree_peak_rt_seconds = quant_labelfree_peak_rt_seconds
        self.quant_labelfree_start_rt_seconds = quant_labelfree_start_rt_seconds
        self.quant_labelfree_end_rt_seconds = quant_labelfree_end_rt_seconds

    def get_modpeptide_string(self):
        return ModifiedPeptide(self.sequence, self.modifications).to_string()

    def get_delta_mass(self):
        if self.observed_mass:
            return self.observed_mass - self.mass
        return 0

    def calc_n_missed_cleavages(self):
        """calculate internal missed cleavages"""
        return calc_peptide_missed_cleavages(self.sequence)

    def has_labelfree_quant(self):
        return self.quant_labelfree_peakintensity is not None

    def __str__(self):
        return "PeptideIdentification: %s, prob=%f" % (self.sequence,
                                                       self.probability)


def calc_peptide_missed_cleavages(pepseq):
    """
    Count the missed cleavages in a peptide sequence
    :param pepseq:
    :return:
    """
    fully_tryptic_peps = calc_tryptic_peptides(pepseq, 0)
    return len(fully_tryptic_peps) - 1


def calc_tryptic_peptides(proteinseq, n_missed_cleavages):
    """
    Trypsinize, returning an unordered list of all peptide sequences with <= n_missed_cleavages
    :param proteinseq:
    :param n_missed_cleavages:
    :return:
    """
    return pyparser.cleave(proteinseq, pyparser.expasy_rules["trypsin"], missed_cleavages=n_missed_cleavages)


class AminoacidModification:
    """ stores information about an amino acid or n-terminal modification that
        may (or always does) happen to peptides.
        special character values for aa  n and c (lowercase) indicate
        n-terminal and c-terminal mods """

    NTERM_CHAR = 'n'
    CTERM_CHAR = 'c'

    MOD_IODOACETAMIDE_MASSDIFF = 57.021
    MOD_OXIDATION_MASSDIFF = 15.994915
    MOD_DISULFIDE_C_MASSDIFF = -1.00794

    def __init__(self, aa, massdiff, is_variable):
        self.aa = aa
        self.massdiff = massdiff
        self.is_variable = is_variable

        self.modified_mass = massdiff
        if not self.is_terminal():
            self.modified_mass = massdiff + pyteomics_mass.std_aa_mass[aa]

    def is_terminal(self):
        """Is this a terminal mod? (if not, aa mod) """
        return self.aa == AminoacidModification.NTERM_CHAR or \
               self.aa == AminoacidModification.CTERM_CHAR

    def get_variable_Y_N(self):
        """return is_variable mapped to Y or N"""
        if self.is_variable:
            return 'Y'
        return 'N'

    def to_string(self):
        var_string = ''
        if self.is_variable:
            var_string = 'V'
        return self.aa + str(self.massdiff) + var_string

    @staticmethod
    def parse_modlist_string(modlist_string):
        """parse a comma-separated string representing a list of mods"""
        result = []
        for chunk in modlist_string.split(","):
            result.append(AminoacidModification.parse_string(chunk))
        return result

    @staticmethod
    def parse_string(mod_string):
        """Parse a string defining an AminoacidModification. Format:
        XDIFF[V]
        where X is the amino acid symbol, DIFF is the mass difference, and
        V, if present, indicates variable.
        Examples: C57.021, M15.99V"""
        aa = mod_string[0]
        is_variable = False
        if mod_string.endswith("V"):
            is_variable = True
            mod_string = mod_string[:len(mod_string) - 1]
        massdiff = float(mod_string[1:])
        return AminoacidModification(aa, massdiff, is_variable)

    def create_modified_aa(self, position_0based):
        return ModifiedAminoacid(position_0based,
                                 self.modified_mass, self.massdiff)


class ModifiedAminoacid:
    """ 
    stores information about a specific instance of an amino acid
    modification (or n-term or cc-term mod) on a particular peptide ID
    """

    # special 'position' values that indicate N-terminal or C-terminal
    POSITION_NTERM = -1
    POSITION_CTERM = -2

    def __init__(self, position_0based, modified_mass, mass_diff):
        """position: 0-based.
           modified_mass: modified mass
           mass_diff: mass difference (for c or n term, same as modified_mass)"""
        self.position_0based = position_0based
        self.modified_mass = modified_mass
        self.mass_diff = mass_diff

    def get_position_1based(self):
        return self.position_0based + 1

    def is_terminal(self):
        return self.position_0based in (ModifiedAminoacid.POSITION_NTERM,
                                        ModifiedAminoacid.POSITION_CTERM)

    def mass_to_string(self):
        return "[" + str(int(round(self.modified_mass))) + "]"


def apply_modifications_to_sequence(sequence, aa_modifications):
    """
    Given a sequence and a list of AminoacidModifications, apply the mods to each position
    in the sequence, as well as the N and C termini 
    Does not handle variable modifications -- raises an exception.
    Does not handle multiple modifications on the same position -- raises an exception
    :param sequence: peptide aa string
    :param aa_modifications: list of AminoacidModifications
    :return: a tuple of:
    (list of mass differences with the same length as the sequence,
    N-terminal delta mass,
    C-terminal delta mass)
    """
    position_massdiff_list = [0] * len(sequence)
    nterm_deltamass = 0.0
    cterm_deltamass = 0.0
    aa_aamod_map = {}
    # build a map from amino acids to ModifiedAminoAcids, and also add n-terminal and c-terminal modifications
    for aa_modification in aa_modifications:
        if aa_modification.aa in aa_aamod_map:
            raise ValueError('make_modified_peptide got multiple modifications for the same position.')
        if aa_modification.is_variable:
            raise ValueError('make_modified_peptide got a variable modification')
        aa_aamod_map[aa_modification.aa] = aa_modification
        # handle N- or C-terminal modifications here
        if aa_modification.aa == AminoacidModification.NTERM_CHAR:
            nterm_deltamass = aa_modification.massdiff
        elif aa_modification.aa == AminoacidModification.CTERM_CHAR:
            cterm_deltamass = aa_modification.massdiff
    # step through the sequence and create the ModifiedAminoacids for each position
    for i in range(0, len(sequence)):
        if sequence[i] in aa_aamod_map:
            position_massdiff_list[i] = aa_aamod_map[sequence[i]].massdiff
    return position_massdiff_list, nterm_deltamass, cterm_deltamass


class ModifiedPeptide:
    """ stores information about a specific peptide with a specific set of modifications.
    Rather brittle. Multiple modifications on the same residue will mess this up. """

    def __init__(self, sequence, modified_aas=None, mass=None):
        self.sequence = sequence
        if modified_aas is None:
            modified_aas = []
        self.modified_aas = modified_aas
        self.mass = mass
        if not mass:
            self.mass = calc_mass_with_mods(sequence, modified_aas)

    def to_string(self):
        result = ""
        for modified_aa in self.modified_aas:
            if modified_aa.position_0based == ModifiedAminoacid.POSITION_NTERM:
                result = result + modified_aa.mass_to_string()
        for i in range(0, len(self.sequence)):
            aa = self.sequence[i]
            result = result + aa
            for modified_aa in self.modified_aas:
                if modified_aa.position_0based == i:
                    result = result + modified_aa.mass_to_string()
        for modified_aa in self.modified_aas:
            if modified_aa.position_0based == ModifiedAminoacid.POSITION_CTERM:
                result = result + modified_aa.mass_to_string()
        return result


def calc_mass_with_mods(peptide_sequence, modified_aas):
    result = calc_unmod_monoisotopic_mass(peptide_sequence)
    for mod_aa in modified_aas:
        result += mod_aa.mass_diff
    return result


def calc_uniquemass_modpeps(peptide_seq, peptide_mods):
    mods_thus_far = []
    remaining_varmods = []
    for mod in peptide_mods:
        if mod.is_variable:
            remaining_varmods.append(mod)
        else:
            if mod.is_terminal():
                position = ModifiedAminoacid.POSITION_NTERM
                if mod.aa == AminoacidModification.CTERM_CHAR:
                    position = ModifiedAminoacid.POSITION_CTERM
                mods_thus_far.append(mod.create_modified_aa(position))
            else:
                for position in [m.start() for m in re.finditer(mod.aa, peptide_seq)]:
                    mods_thus_far.append(mod.create_modified_aa(position))
    return add_varmods_recurse(peptide_seq, mods_thus_far, remaining_varmods)


def add_varmods_recurse(peptide_seq, mods_thus_far, remaining_varmods):
    """helper method for calc_uniquemass_modpeps. Does *not* create all possible
    combinations of modified and unmodified amino acids. Instead, creates one
    representative for each unique-mass-producing arrangement of modifications"""
    if not remaining_varmods:
        return [ModifiedPeptide(peptide_seq, mods_thus_far)]

    mod = remaining_varmods[0]
    new_remaining_varmods = []
    if len(remaining_varmods) > 1:
        new_remaining_varmods = remaining_varmods[1:]

    if mod.is_terminal():
        position = ModifiedAminoacid.POSITION_NTERM
        if mod.aa == AminoacidModification.CTERM_CHAR:
            position = ModifiedAminoacid.POSITION_CTERM
        result = add_varmods_recurse(peptide_seq, mods_thus_far, new_remaining_varmods)
        result.extend(
            add_varmods_recurse(peptide_seq, mods_thus_far + (mod.create_modified_aa(position)), new_remaining_varmods))
        return result
    else:
        aa_positions = [m.start() for m in re.finditer(mod.aa, peptide_seq)]
        thismod_accum = []
        result = []
        # this for loop will lopp from no occurrences of the modification
        # through all of them, adding from left to right in the sequence
        for i in range(0, len(aa_positions) + 1):
            result.extend(add_varmods_recurse(peptide_seq, mods_thus_far + thismod_accum, new_remaining_varmods))
            if i < len(aa_positions):
                thismod_accum.append(mod.create_modified_aa(aa_positions[i]))
        return result


#restriction: no more than one variable modification per residue/terminus
def calc_all_possible_masses(peptide_seq, peptide_mods):
    """calculate all possible monoisotopic masses of a peptide"""
    mass_no_var_mods = calc_unmod_monoisotopic_mass(peptide_seq)

    var_mods = list()
    for mod in peptide_mods:
        if mod.is_variable:
            var_mods.append(mod)
        else:
            if mod.aa in [AminoacidModification.NTERM_CHAR,
                          AminoacidModification.CTERM_CHAR]:
                mass_no_var_mods = mass_no_var_mods + mod.mass_diff
            else:
                mass_no_var_mods = (mass_no_var_mods +
                                    (peptide_seq.count(mod.aa) * mod.massdiff))

    return list(add_masses_for_mods_recurse(peptide_seq, var_mods, mass_no_var_mods))


def add_masses_for_mods_recurse(pep_sequence, remaining_varmods,
                                mass_without_remaining_mods):
    """recursively create a *set* of all masses for all modifications"""
    if not remaining_varmods:
        return set([mass_without_remaining_mods])
    mod = remaining_varmods[0]
    remaining_varmods.remove(mod)

    mod_locs_in_seq_count = 0
    if mod.aa in [AminoacidModification.NTERM_CHAR,
                  AminoacidModification.CTERM_CHAR]:
        mod_locs_in_seq_count = 1
    else:
        mod_locs_in_seq_count = pep_sequence.count(mod.aa)

    result = set()

    numbers_of_mods_to_add = range(0, mod_locs_in_seq_count + 1)

    for i in numbers_of_mods_to_add:
        result = result.union(result, add_masses_for_mods_recurse(pep_sequence,
                                                                  remaining_varmods,
                                                                  mass_without_remaining_mods + (i * mod.massdiff)))
    return result


def calc_peptide_mass_no_var_mods(pep, is_nterminal, nterm_conjugate_mass):
    """calculate the mass we'll see if there are no oxidized Methionines"""
    unmodified_mass = pyteomics_mass.fast_mass(pep)
    if is_nterminal:
        unmodified_mass += nterm_conjugate_mass
    return unmodified_mass + AminoacidModification.MOD_IODOACETAMIDE_MASSDIFF * pep.count("C")


def calc_peptide_masses(pep, is_nterminal, nterm_conjugate_mass,
                        should_use_ox_m):
    """Calculate all possible masses for this peptide, including modified and
    unmodified Methionines if should_use_ox_m"""
    unmodified_mass = pyteomics_mass.fast_mass(pep)

    unmodified_mass += AminoacidModification.MOD_IODOACETAMIDE_MASSDIFF * pep.count("C")
    if not should_use_ox_m:
        return [unmodified_mass]
    if is_nterminal:
        unmodified_mass += nterm_conjugate_mass
    methionineCount = pep.count("M")
    result = list()
    for i in range(0, methionineCount + 1):
        result.append(unmodified_mass + (AminoacidModification.MOD_OXIDATION_MASSDIFF * i))
    return result


def count_cysteines_in_pepseqs(pepseqs):
    result = []
    for pepseq in pepseqs:
        n_cys = pepseq.count('C')
        while n_cys > len(result) - 1:
            result.append(0)
        result[n_cys] += 1
    return result




# modifications corresponding to iodoacetamide and variable ox methionine

MODIFICATION_IODOACETAMIDE_STATIC = AminoacidModification("C", AminoacidModification.MOD_IODOACETAMIDE_MASSDIFF,
                                                          False)
MODIFICATION_IODOACETAMIDE_VAR = AminoacidModification("C", AminoacidModification.MOD_IODOACETAMIDE_MASSDIFF,
                                                       False)
MODIFICATION_OXM_VAR = AminoacidModification("M", AminoacidModification.MOD_OXIDATION_MASSDIFF,
                                                       False)

MODIFICATIONS_IODOACETAMIDE_OXMVAR = [AminoacidModification("C",
                                                            AminoacidModification.MOD_IODOACETAMIDE_MASSDIFF,
                                                            False),
                                      AminoacidModification("M",
                                                            AminoacidModification.MOD_OXIDATION_MASSDIFF,
                                                            True)]


def calc_pepmasses_partial_reduction_nooxm(peptide_sequence, is_alkylated=False):
    """Calculate all possible masses for peptide_sequence for all reduction states from fully-disulfide-bonded
    to fully-reduced. Either assume sample alkylation or don't"""
    disulfide_count = peptide_sequence.count('C') / 2
    unreduced_mass = calc_unmod_monoisotopic_mass(peptide_sequence) + peptide_sequence.count('C') * AminoacidModification.MOD_DISULFIDE_C_MASSDIFF
    disulfide_reduction_massdiff = 2 * (-AminoacidModification.MOD_DISULFIDE_C_MASSDIFF)
    if is_alkylated:
        disulfide_reduction_massdiff += 2 * AminoacidModification.MOD_IODOACETAMIDE_MASSDIFF
    result = list()
    for i in range(0, disulfide_count+1):
        result.append(unreduced_mass + i * disulfide_reduction_massdiff)
    return result

def calc_unmod_monoisotopic_mass(peptide_sequence):
    """return unmodified monoisotpic mass"""
    return calc_monoisotopic_mass(peptide_sequence, massUtil.std_aa_comp)


def calc_monoisotopic_mass(modified_pep_sequence, modified_mass_dict):
    """return monoisotopic mass of a modified peptide"""
    return massUtil.calculate_mass(modified_pep_sequence,
                               aa_comp=modified_mass_dict)


def get_aa_mass(aa):
    """convenience method to get the monoisotopic mass of an amino acid"""
    return pyteomics_mass.std_aa_mass[aa]


def calc_mz_from_neutralmass_charge(neutral_mass, charge):
    """
    Given a neutral mass and a charge, calculate mz
    :param neutral_mass:
    :param charge:
    :return:
    """
    return neutral_mass / charge + HYDROGEN_MASS


def calc_mz_from_mplush_charge(m_plus_h, charge):
    """
    Given an M+H and a charge, calculate mz
    :param m_plus_h:
    :param charge:
    :return:
    """
    return calc_mz_from_neutralmass_charge(m_plus_h - HYDROGEN_MASS, charge)


def calc_neutralmass_from_mz_charge(mz, charge):
    """
    Given an mz and a charge, calculate the neutral mass of the ion
    :param mz:
    :param charge:
    :return:
    """
    return (mz - HYDROGEN_MASS) * charge


def calc_mplush_from_mz_charge(mz, charge):
    """
    Given an mz and a charge, calculate the M+H mass of the ion
    :param mz:
    :param charge:
    :return:
    """
    return calc_neutralmass_from_mz_charge(mz, charge) + HYDROGEN_MASS



