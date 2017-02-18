""" Python module for using Probanno-Standalone

This module has been adapted from Probanno-Standalone for use as a python module. Probabilistic Annotation is an algorithm
used to assign probabilities to complexes and reactions given a genome for use in reconstructing genome-scale metabolic
models. Probabilistic Annotation was originally written by Mike Mundy, Matt Benedict, and a stand-alone version adapted
from this by Terry Farrah. This module is an adaptation using the stand-alone version to run Probabilistic Annotation
from within a python application.

Example: Downloading a genome, running probabilistic annotation, and exporting results to JSON and TSV files.

      >>> rxn_probs = probanno.generate_reaction_likelihoods(probanno.get_fasta_by_id('267377', 'my.fasta'), 'templates/GramNegative.json', genome_id='267377')
      >>> probanno.export_json(rxn_probs, 'myprobabilities.json')



"""

import json
import re
import requests
import cobra

from lib.ProbAnnotationWorker import ProbAnnotationWorker

# Constants
UNIPROT_BASE_URL = 'http://www.uniprot.org/uniprot/?format=fasta&query=organism:'


def get_fasta_by_id(proteome_id, output_file):
    """
    Downloads a FASTA file for the proteome by organism ID
    :param proteome_id: ID of the organism, e.g. 267377
    :param output_file: path to a destination for the file
    :return: the name of the file where the FASTA is saved (should be output_file)
    """
    taxid_pattern = re.compile('^\d{1,7}$')
    # if not taxid_pattern.match(proteome_id):  # fetch file from Uniprot
    #     raise ValueError(str(proteome_id) + ' is not a valid proteome identifier')
    url = UNIPROT_BASE_URL + proteome_id
    attempts = 0
    while attempts < 3:
        try:
            response = requests.get(url)
            if response.status_code > 399 or response.status_code < 200:
                raise requests.HTTPError(response.status_code + ': ' + response.content)
            content = response.content
            if len(content) < 10:
                raise FastaNotFoundError()
            with open(output_file, 'w') as f:
                f.write(content)
            break
        except requests.HTTPError as e:
            attempts += 1
            if attempts >= 3:
                raise FastaNotFoundError('Failed to download fasta: ' + response.status_code + ' response.content')
    return output_file


def generate_reaction_probabilities(fasta_file, template_model_file, genome_id=None):
    """
    A function for generating reaction likelihoods for a given genome according to the Probabilistic Annotation
    algorithm as
    :param fasta_file: file name of a proteome sequence for an organism in FASTA format
    :param template_model_file: filename of a template model. Some are included, e.g. templates/GramNegative.json
    :param genome_id: (optional) genome id for the organism. Used in naming intermediate files
    :return: ReactionProbabilities object

    """
    if genome_id is None:
        # Use fasta_file name minus extension. worker uses only for file names and logging
        genome_id = '.'.join(fasta_file.split('.')[0:-1])
    # Create a worker for running the algorithm.
    worker = ProbAnnotationWorker(genome_id)
    try:
        template_model = _load_template_file(template_model_file)

        # Run blast using the fasta file.
        blast_result_file = worker.runBlast(fasta_file)

        # Calculate roleset probabilities.
        rolestring_tuples = worker.rolesetProbabilitiesMarble(blast_result_file)

        # Calculate per-gene role probabilities.
        role_probs = worker.rolesetProbabilitiesToRoleProbabilities(rolestring_tuples)

        # Calculate whole cell role probabilities.
        total_role_probs = worker.totalRoleProbabilities(role_probs)

        # Calculate complex probabilities.
        complex_probs = worker.complexProbabilities(total_role_probs, complexesToRequiredRoles=_complex_to_roles_dict(template_model))

        # Calculate reaction probabilities.
        rxn_probs = worker.reactionProbabilities(complex_probs, rxnsToComplexes=_reactions_to_complexes_dict(template_model))

        # Store in dictionary for better serialization
        return ReactionProbabilities([{'reaction': r[0], 'probability': r[1], 'type': r[2], 'complexes': _deserialize_cplx(r[3], worker.config['separator']), 'gpr': r[4]} for r in rxn_probs])
    finally:
        worker.cleanup()  # worker creates lots of temporary and intermediate files. Allow it to clean up


def probabilistic_gapfill(model, universal_model, reaction_probabilities, default_penalties=None, dm_rxns=False, ex_rxns=False, **solver_parameters):
    """
    Gapfill a model using probabilistic weights
    :param default_penalties:
    :param model: cobra Model object, the model to be gapfilled
    :param universal_model: cobra Model object representing the database of reactions to choose from
    :param reaction_probabilities: reaction_probabilities dictionary
    :return:
    """

    if default_penalties is None:
        default_penalties = {'Universal': 1, 'Exchange': 100, 'Demand': 1, 'Reverse': 75}
    penalties = default_penalties
    reactions_to_remove = []
    for r in universal_model.reactions:
        if model.reactions.has_id(r.id):
            reactions_to_remove.append(r)
            penalties[r.id] = 0  # In the model
        elif r.id in reaction_probabilities:
            penalties[r.id] = max(0, 1 - reaction_probabilities.get_probability(r.id)) * (penalties[r.id] if r.id in penalties else 1)
    universal_model.remove_reactions(reactions_to_remove)
    return cobra.flux_analysis.growMatch(model, universal_model, penalties=penalties, dm_rxns=dm_rxns, ex_rxns=ex_rxns, **solver_parameters)


def build_universal_model(template_model_file, clean_exchange_reactions=False, compartments=[0]):
    """

    :param template_model_file: path to file where a template model is saved. e.g. templates/GramNegative.json
    :param compartments: list (usually integers) to be appended to all compounds/reactions indicating model compartments
                         e.g. if compartments = [0], rxn00001_c, rxn00001_c0 will be in the universal model, and will
                         have compounds of the form cpd00001_c, cpd00001_c0 respectively. This is used if your model has
                         multiple compartments or indicates the compartment as 0 e.g. rxn00001_c0 with cpd00001_c0
    :param clean_exchange_reactions: If true, re-formats exchange reactions to cobra form (see below)

    Cleaning Exchange Reactions:
    if a reaction in the template model uses external compounds to indicate exchange reactions, and the clean_exchange_reactions
    flag is set to True, then external compounds are removed to indicate exchange implicitly (CobraPy style)

    Example: "cpd00001_e0 <=> cpd00001_c0" is used to indicate exchange reaction of cpd00001. With the flag, this is
             changed to " <=> cpd00001_c0"
    :return:
    """
    universal = cobra.Model("Universal Reactions")
    template_model = _load_template_file(template_model_file)
    # Creates a dictionary of the reactions from the tab delimited database, storing their ID and the reaction string
    for template_rxn in template_model['reactions']:
        compartments = [str(i) for i in compartments]
        compartments.extend('')
        for suffix in compartments:
            _add_rxn_from_template_rxn(universal, template_rxn, clean_exchange_reactions=clean_exchange_reactions, compartment=suffix)
    return universal


def clean_exchange_reactions(model, regex='.*_e([0-9]*)$'):
    model = model.copy()
    compound_regex = re.compile(regex)
    mets_to_clean = [m for m in model.metabolites if compound_regex.match(m.id)]
    for m in mets_to_clean:
        m.remove_from_model(method='subtractive')
    return model


def _add_rxn_from_template_rxn(model, template_reaction, clean_exchange_reactions=False, compartment=''):
    rxn = cobra.Reaction(template_reaction['id'] + str(compartment))
    try:
        model.add_reaction(rxn)
    except Exception:  # Can't be more specific because this is exactly how it is specified in CobraPy
        return

    # Get equation as a string
    reactants = []
    products = []
    for compound in template_reaction['templateReactionReagents']:
        # coefficient is a string float representation. Convert to integer
        coef = int(float(compound['coefficient']))
        comp_str = compound['templatecompcompound_ref'].split('/')[-1] + str(compartment)
        if abs(coef) > 1:
            comp_str = '(' + str(abs(coef)) + ') ' + comp_str
        if coef < 0:
            reactants.append(comp_str)
        else:
            products.append(comp_str)
    if clean_exchange_reactions:
        reactants = [r for r in reactants if r.find('_e' + str(compartment)) < 0]
        products = [r for r in products if r.find('_e' + str(compartment)) < 0]

    rxn.reaction = ' + '.join(reactants) + '<=>' + ' + '.join(products)  # TODO adjust for if/not reversible
    rxn.name = template_reaction['name']
    # TODO ADD BOUNDS
    return rxn


def export_json(rxn_probs, filename):
    """
    Exports the given reaction probabilities into a JSON formatted file, saved at filename
    :param rxn_probs: reaction probabilities, as outputted by generate_reaction_probabilities
    :param filename: file name (str)
    :return: filename
    """
    with open(filename, 'w') as f:
        f.write(json.dumps(rxn_probs))
    return filename


def _load_template_file(template_file):
    # Create a dictionary from the json template file
    with open(template_file) as f:
        return json.loads(f.read())


def _complex_to_roles_dict(template):
    # Build a dictionary to look up roles in the template by ID.
    roles = dict()
    for index in range(len(template['roles'])):
        roles[template['roles'][index]['id']] = index

    # Create a dictionary to map a complex to a list of roles as defined in the template.
    complexes_to_roles = dict()
    for index in range(len(template['complexes'])):
        complex_id = template['complexes'][index]['id']
        if len(template['complexes'][index]['complexroles']) > 0:
            complexes_to_roles[complex_id] = list()
            for cr_index in range(len(template['complexes'][index]['complexroles'])):
                # A complex has a list of complexroles and each complexrole has a reference
                # to a role and each role has a name. Role ID is last element in reference.
                role_id = template['complexes'][index]['complexroles'][cr_index]['templaterole_ref'].split('/')[-1]
                # Mike Mundy used this statement, but it results in all probs being zero.
                # complexes_to_roles[complex_id].append(roleId)
                complexes_to_roles[complex_id].append(template['roles'][roles[role_id]]['name'])
    return complexes_to_roles


def _reactions_to_complexes_dict(template):
    # Create a dictionary to map a reaction to a list of complexes as defined in the template.
    reactions_to_complexes = dict()
    for index in range(len(template['reactions'])):
        reaction_id = template['reactions'][index]['id']
        if len(template['reactions'][index]['templatecomplex_refs']) > 0:
            reactions_to_complexes[reaction_id] = list()
            for complexRef in template['reactions'][index]['templatecomplex_refs']:
                # Complex ID is last element in reference.
                reactions_to_complexes[reaction_id].append(complexRef.split('/')[-1])
    return reactions_to_complexes


def _deserialize_cplx(complex_str, separator):
    complexes = complex_str.split(separator)
    result = []
    for complex in complexes:
        if len(complex) > 0:
            values = re.split('[\(\)\;]', complex)
            complex = str(values[0]).strip()
            probability = float(values[1])
            complex_type = str(values[2]).strip()
            result.append({'complex': complex, 'probability': probability, 'type': complex_type})
    return result


class ReactionProbabilities(object):
    def __init__(self, rxn_probs):
        self.data = dict([(r['reaction'], r) for r in rxn_probs]) if rxn_probs is not None else dict()

    def __str__(self):
        return self.to_json()

    def get_probability(self, reaction):
        """
        return the probability of a given reaction
        :param reaction:
        :return:
        """
        return self.__getitem__(reaction)

    def add_reaction(self, reaction_struct):
        """
        Add
        :param reaction_struct:
        :return:
        """
        self.data[reaction_struct['reaction']] = reaction_struct

    def to_json_file(self, path):
        """
        Serializes this object as a JSON stringrxn
        :return:
        """
        with open(path, 'w') as f:
            f.write(self.to_json())

    @staticmethod
    def from_json_file(path):
        """
        Deserialize a ReactionProbabilities from a JSON file
        :param path: File in JSON form of the data relevant to this object
        :return:
        """
        with open(path, 'r') as f:
            return ReactionProbabilities.from_json(f.read())

    def __getitem__(self, item):
        if item in self.data:
            return self.data[item]['probability']
        else:
            raise KeyError('No probability stored for reaction: ' + str(item))

    def __setitem__(self, key, value):
        value = float(value)  # throws an error if the value is not permitted
        self.data[key] = {'reaction': key, 'probability': value}

    def __contains__(self, item):
        return item in self.data

    def update(self, rxn_probs):
        """
        Updates the Reaction Probabilities
        :param rxn_probs:
        :return:
        """
        pass

    def to_json(self):
        return json.dumps([rxn[1] for rxn in self.data.items()] if self.data is not None else None)

    @staticmethod
    def from_json(json_str):
        data = json.loads(json_str)
        return ReactionProbabilities(data)



class FastaNotFoundError(Exception):
    pass
