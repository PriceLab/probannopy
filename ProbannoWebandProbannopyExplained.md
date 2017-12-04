# ProbannoWeb and Probannopy Explained

The tools and services provided in ProbannoWeb and Probannopy implement the probabalistic annotation and likelihood-based gap-filling algorithms. This guide is written wuth the assumption that the reader has a basic understanding of these algorithms and what they intend to accomplish. These algorithms, their uses, and motivations are described in full in the original paper on these methods ([Benedict et. al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4199484/)). A brief explanation of these algorithms is also included in our [publication on these tools](insertlink.com).

## Algorithm Interfaces Supported by Both ProbannoWeb and Probannopy

### Probabalistic Annotation

- **Input**: Proteome sequence in FASTA format
- **Output**: A list of reaction likelihoods

Probabalistic Annotation assigns likelihoods for each candidate reaction in a universal model which serves as a reaction database based on sequence homology between the argument FASTA sequence and a trusted set of annotations. A reaction likelihood consists of: 

- a reaction identifier 
- the calculated likelihood of that reaction occuring as a result of the metabolism encoded in the FASTA sequence
- the Gene-Protein-Reaction relationship associated with this reaction
- A list of associated complexes and their likelihoods

Below is an example:

```
{
	reaction: "rxn05250_c",
	type: "HASCOMPLEXES",
	gpr: "tr|Q6LZQ8|Q6LZQ8_METMP",
	probability: 0.5634711450747343,
	complexes: [
		{
			complex: "cpx00498",
			type: "CPLX_FULL",
			probability: 0.5635
		}
	]
}
```

### Likelihood-based Gap-filling
- **Input**: 
	- A metabolic model in a CobraPy compatible format
	- A list of reaction likelihoods
	- A choice of universal model to serve as a reaction database
- **Output**:	
   - A metabolic model gap-filled to using reactions from the universal model according to likelihood

Likelihood-based Gap-filling uses these assigned likelihoods to correspondingly adjust the penalty weights for reactions from the universal model in the Cobrapy gap-filling problem. The clearest use case for this is as an alternative step in workflows which already use Cobrapy for the gapfilling of metabolic models in ModelSEED indentifiers.

## Interfacing with ProbannoWeb and Probannopy

### Probannopy

Probannopy is a python package, and can be used in any python environment of version 2.6 or greater. Probannopy can be interacted with in an interactive python session, or imported and used in other python scripts and packages. Installation instructions can be found in the [GitHub repository for probannopy](https://github.com/PriceLab/probannopy)

##### Probabalistic Annotation Example

Each function in the probannopy package includes documentation in python documentation format. One way to read the documentation is by calling `help(probanno)` in an interactive python session. Here is the documentation for the function fresponsible for running probabilistic annotation:

```
generate_reaction_probabilities(fasta_file, template_model_file, genome_id=None)
        A function for generating reaction likelihoods for a given genome according to the Probabilistic Annotation
        algorithm as
        :param fasta_file: file name of a proteome sequence for an organism in FASTA format
        :param template_model_file: filename of a template model. Some are included, e.g. templates/GramNegative.json
        :param genome_id: (optional) genome id for the organism. Used in naming intermediate files
        :return: ReactionProbabilities object
```

As one can see above, the required parameters here are a path to a FASTA proteome sequence file and a template model file to serve as a reaction database. The package includes 3: `GramNegative.json`, `GramPositive.json`, `Microbial.json`.

Probannopy also supports retrieving this proteome sequence from UniProt by taxonomic idenitifier:

```
    get_fasta_by_id(proteome_id, output_file)
        Downloads a FASTA file for the proteome by organism ID
        :param proteome_id: ID of the organism, e.g. 267377
        :param output_file: path to a destination for the file
        :return: the name of the file where the FASTA is saved (should be output_file)
```

The output of this function is a path to the downloaded FASTA file, which can be used as an argument in `generate_reaction_probabilities`. With these two functions, we can run probabilistic annotation on the proteome sequence for an organism we have the taxonomic indentifier for:

```python
import probanno
proteome_id = '267377'  # methanococcus maripaludis
output_file = 'maripaludis.fasta'
fasta_file = probanno.get_fasta_by_id(proteome_id, output_file)
likelihoods = probanno.generate_reaction_probabilities('my_fasta_file.fasta', 'templates/GramNegative.json', genome_id='my_genome')
# Output is a ReactionProbabilities object, which has dictionary behavior for reaction_id -> likelihood
# Example:
print(likelihoods['rxn00383_c'])
```

Final output is a reaction probabilities object, which is essentially a dictionary of reaction identifiers to their probabilities. It also includes functions for serialization and deserialization to JSON, and has metadata accessible at `likelihoods.data`.


#### CobraPy Gap-filling interface

Probannopy includes a gapfilling interface similar to the cobrapy gap-filling interface in `cobra.flux_analysis.gapfill` interface. For those unfamiliar, the gap-filling function in its simplest form takes in an argument model to gapfill and a "universal" model which serves as a database of reactions which can be used in gap-filling the target model. Additionally, our likelihood based gap-filling function takes in a ReactionLikelihoods object, which we use to individually penalize reactions in the gap-filling formulation according to their calculated likelihood. Our fuction is fairly simple given these likelihoods, and simply calculates a per reaction penalty before calling the `cobra.flux_analysis.gapfill` funciton in `cobrapy` with these penalties.

Here is the interface:

```
probabilistic_gapfill(model, universal_model, reaction_probabilities, clean_exchange_rxns=True, default_penalties=None, dm_rxns=False, ex_rxns=False, **solver_parameters)
    Gapfill a model using probabilistic weights
    :param model: cobra Model object, the model to be gapfilled
    :param universal_model: cobra Model object representing the database of reactions to choose from
    :param reaction_probabilities: reaction_probabilities dictionary
    :param clean_exchange_rxns: If specified, model's with exterior metabolites (e.g. a_c0 -> a_e0) will be translated to ones without exterior metabolite
        which is the normal cobra form (e.g. a_c0 -> )
    :return: A list of solutions, where a solution is a single list of reactions that could be added to the model to gapfill it 
```

Here is an example usage:

```
# build universal model is a helper function for building a reaction database from a template file, a few of which are included in probannopy
universal_model = probanno.build_universal_model(temlate_file, clean_exchange_reactions=True)
reactions = probanno.probabilistic_gapfill(model, universal_model, likelihoods)[0]
```

**Note:** Although algorithmically we state our gap-filling tool as producing a gapfilled model as opposed to a list of solutions, we choose to output a list of solutions in the probannopy implementation in order to maintain consistency with the current cobrapy gapfilling interfaces. 

For comparison, here is the interface for `cobra.flux_analysis.gapfilling`:

```
gapfill(model, universal=None, lower_bound=0.05, penalties=None, demand_reactions=True, exchange_reactions=False, iterations=1)
    Perform gapfilling on a model.

    See documentation for the class GapFiller.

    Parameters
    ----------
    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model, None
        A universal model with reactions that can be used to complete the
        model. Only gapfill considering demand and exchange reactions if
        left missing.
    lower_bound : float
        The minimally accepted flux for the objective in the filled model.
    penalties : dict, None
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively.
    iterations : int
        The number of rounds of gapfilling to perform. For every iteration,
        the penalty for every used reaction increases linearly. This way,
        the algorithm is encouraged to search for alternative solutions
        which may include previously used reactions. I.e., with enough
        iterations pathways including 10 steps will eventually be reported
        even if the shortest pathway is a single reaction.
    exchange_reactions : bool
        Consider adding exchange (uptake) reactions for all metabolites
        in the model.
    demand_reactions : bool
        Consider adding demand reactions for all metabolites.

    Returns
    -------
    iterable
        list of lists with on set of reactions that completes the model per
        requested iteration.

    Examples
    --------
    >>> import cobra.test as ct
    >>> from cobra import Model
    >>> from cobra.flux_analysis import gapfill
    >>> model = ct.create_test_model("salmonella")
    >>> universal = Model('universal')
    >>> universal.add_reactions(model.reactions.GF6PTA.copy())
    >>> model.remove_reactions([model.reactions.GF6PTA])
    >>> gapfill(model, universal)
```

In our usage, likelihoods are used to create a dictionary of reaction-specific penalties, which is then passed to this function. For those familiar with the cobrapy package, our tool should provide a useful and convenient addition to a workflow.

### ProbannoWeb API

The ProbannoWeb API is a simple REST API interface through which users can run probablistic annotation and likelihood based gap-filling on our servers without having to install any probanno-specific software of their own. For those familiar with REST APIs and the HTTP system, our API should feel fairly simple and straight forward, and our full and complete API documentation is on SwaggerHub at [https://app.swaggerhub.com/apis/kingb12/ProbannoWeb/1.0.0](https://app.swaggerhub.com/apis/kingb12/ProbannoWeb/1.0.0). For others, we will provide some code examples for calling out APIs from a python environment, but a general understanding of REST would be advantageous and beyond the scope of this walk-through.

#### Sessions

One unique aspect of the ProbannoWeb API is the use of sessions and session identifiers. You'll notice that neither our web-site nor API use usernames/passwords or some other authentication scheme. This is for the simplicity of both our system and the user experience, and to emphasize that our service is intended for these specific algorithmic use cases as opposed to a general modeling resource and repository like such as [ModelSEED](http://modelseed.org/). Instead, to preserve state between requests and establishing a basic sense for identity, we use session identifiers, which are basically Universally Unique Identifer (UUID) keys which identify a session, it's jobs, models, and reaction likelihoods. 

The workflow begins with a user requesting a session. Examples here are shown using the `curl` command, but any languages tool supporting web requests (e.g. python `requests`) can be substituted.

Request:

```
curl -X GET "http://probannoweb.systemsbiology.net/api/session" -H "accept: application/json"
```

Response Body:

```
"f8ae14e6-304d-4771-8ae6-fbe6874968d1"
```

Here the user receives a UUID which can be treated as a 'key' to their sessions jobs, models, and reaction likelihoods. For each other kind of request, the `session` returned by this call is passed as a header parameter.


#### Jobs

Our core algorithm requests run asynchronously, meaning that we return a response for the request before guaranteeing that the calculation is done. This is because many of these calculations can run from several seconds to several minutes, and our servers can only handle so many at a time, and process them in a first-in first-out manner. Instead of returning results to a calculation request, we return a job, which encodes the status of the calculation, as well as what is being calculated.

Here is the anatomy of a Job:

```
{
  "jid": "795d7cad-00ca-4bb3-ab47-bd7cbf3afbc5",
  "job": "calculate_probanno",
  "sid": "f8ae14e6-304d-4771-8ae6-fbe6874968d1",
  "status": "Complete",
  "target": "267377"
}
```

The `jid` is the ID of the Job, and can be used to uniquely retrieve a job to check its status. `sid` is the session associated with this job. `job` is the kind of work being done in this Job. the `target` is the primary argument for the Job, and the `status` encodes the status in the running of the job, which will be one of `Not Started`, `Running`, `Complete` or `Failure`. For a probabilistic annotation of gapfilling, our endpoints return a Job. A Job (and its current status) can be retreived by Job ID at GET `/api/job`, which takes `job_id` as a query argument and `sesison` as a header parameter.

Request:

```
curl -X GET "http://probannoweb.systemsbiology.net/api/job?job_id=795d7cad-00ca-4bb3-ab47-bd7cbf3afbc5" -H "accept: application/json" -H "session: f8ae14e6-304d-4771-8ae6-fbe6874968d1"
``` 

Response:
```
{
  "jid": "795d7cad-00ca-4bb3-ab47-bd7cbf3afbc5",
  "job": "calculate_probanno",
  "sid": "f8ae14e6-304d-4771-8ae6-fbe6874968d1",
  "status": "Complete",
  "target": "267377"
}
```

#### Calculating Reaction Likelihoods

Calculating reaction likelihoods is done using the `/api/probanno/calculate` endpoint, using either a GET or PUT request.
We'll walk through the simplest case (GET), in which you are generating likelihoods for a proteome seequence available in Uniprot and you have a FASTA identifer. This endpoint takes a query parameter `fasta_id` and a header parameter `session`. Here is an example call:

Request:

```
curl -X GET "http://probannoweb.systemsbiology.net/api/probanno/calculate?fasta_id=267377" -H "accept: application/json" -H "session: f8ae14e6-304d-4771-8ae6-fbe6874968d1"
```

Response:

```
{
  "jid": "795d7cad-00ca-4bb3-ab47-bd7cbf3afbc5",
  "job": "calculate_probanno",
  "sid": "f8ae14e6-304d-4771-8ae6-fbe6874968d1",
  "status": "Complete",
  "target": "267377"
}
```

In the case where you are generating likelihoods for a FASTA sequence you have locally, you can use a PUT request on the same endpoint with the file included in the form data as `fasta`. More details available on [Swaggerhub](https://app.swaggerhub.com/apis/kingb12/ProbannoWeb/1.0.0).

When a job completes, the result is saved and can be retrieved by ID. In the case of reaction likelihoods, the ID is the `fasta_id`. This is either the uniprot taxonomic identifier for GET calculated likelihoods or the filename minus .fasta for PUT calculated likelihoods. You can always list your likelihoods using GET `/api/probanno/list`:

Request:

```
curl -X GET "http://probannoweb.systemsbiology.net/api/probanno/list" -H "accept: application/json" -H "session: f8ae14e6-304d-4771-8ae6-fbe6874968d1"
```

Response: 

```
[
  {
    "fasta_id": "267377",
    "name": "Methanococcus maripaludis (strain S2 / LL)"
  }
]
```
Which returns a list of likelihoods with their ID (`fasta_id`) and name.

You can then retrieve a reaction likelihoods by `fasta_id` at `/api/probanno`:

Request:

```
curl -X GET "http://probannoweb.systemsbiology.net/api/probanno?fasta_id=267377" -H "accept: application/json" -H "session: f8ae14e6-304d-4771-8ae6-fbe6874968d1"
```

Which returns the likelihoods in JSON format. The `/api/probanno/download` endpoint takes the same arguments but returns the likelihoods as a file download.

#### Gapfilling

To perform likelihood-based gapfilling, you must first upload a metabolic model. It is important this model uses ModelSEED identifiers, and is formatted in cobrapy JSON format. One can easily take a model in the cobrapy environment and export it to the proper JSON format using the `cobra.io.json.to_json` function in cobra. 

Model storage and retrieval occurs through the `/api/model` endpoint. GET is for retrieving a Model. POST is for uploading a new model, and PUT is for overwritting an existing one (by ID). To upload a model, one supplies a `model_id`, which is a string used to uniquely identify the model, a `file` containing the metabolic model in JSON format, and the `session` header parameter. Here is an example:

Request:

```
curl -X POST "http://probannoweb.systemsbiology.net/api/model" -H "accept: application/json" -H "session: f8ae14e6-304d-4771-8ae6-fbe6874968d1" -H "Content-Type: multipart/form-data" -F "model_id=my_model" -F "file=@maripaludis_model.json;type=application/json"
```

Which will simply return a 200 status if the upload is successful or a 400 indicating that some arguments were missing or illegal.

Next, one submits a gap-filling job to GET `/api/model/gapfill` with the following query arguments:
- `model_id`: ID of the uploaded model to gapfill
- `output_id`: output ID for the gap-filled metabolic model. The new model will be saved with this as its ID
- `fasta_id`: ID of the reaction likelihoods for this gap-filling job
- `template`: A string argument for choosing between our pre-supplied reaction database models (GramNegative, GramPositive, or Microbial)

as well as the `session` header parameter. Example: 

Request:

```
curl -X GET "http://probannoweb.systemsbiology.net/api/model/gapfill?model_id=my_model&output_id=my_gapfilled_model&fasta_id=267377&template=GramPositive" -H "accept: application/json" -H "session: f8ae14e6-304d-4771-8ae6-fbe6874968d1"
```

Response:

```
{
  "jid": "00ae756a-2f7d-4adf-86aa-5dad4b858131",
  "job": "gapfill_model",
  "sid": "f8ae14e6-304d-4771-8ae6-fbe6874968d1",
  "status": "Not Started",
  "target": "my_model"
}
```

Once the gap-filling is completed, one can retrieve the model at GET `/api/model?model_id=my_gapfilled_model` or as a file using the otherwise equivalent `/api/model/download?model_id=my_gapfilled_model`.

And that't the basics! More information can be found on [SwaggerHub](https://app.swaggerhub.com/apis/kingb12/ProbannoWeb/1.0.0), including other endpoints not described here for listing models, likelihoods, and jobs.

### ProbannoWeb Web Interface

ProbannoWeb is also available from a light-weight web site interface for users not familiar with python or REST, who would just like to try out the service without installation. The web-site provides basic functionality for calculating likelohoods, gapfilling models, listing jobs, and more from its primary home page at [http://probannoweb.systemsbiology.net/](http://probannoweb.systemsbiology.net/). One convenience is that sessions are managed entirely by the web-site, so users don't need to track or submit a session ID as they navigate the site. As such, users must accept cookies from the site (the session ID). 

## Troubleshooting

For the most part, users shouldn't experience any trouble with calculating reaction likelihoods provided the argument FASTA is a properly formatted amino acid sequence. Gap-filling, however, can prove more troublesome. 

One of the most common issues occurs when a **gap-filling job runs indefinitely or times out**. If you are using probannopy with a cobrapy environment, it would be good to make sure you are using a capable solver such as CPLEX or Gurobi. Glpk often won't be able to handle a problem of this size. Free academic licenses are often available for these stronger solvers.

If a gap-filling job times out or fails with infeasability errors, the problem is likely more nuanced. This often occurs when there is an identifier mis-match between the model and the reaction database. To resolve this:
- Make sure reactions in the model and reactions in the universal database model use the same identifiers.
- Make sure the reactions in the model and universal model properly format exchange reactions (e.g. in Cobra, it is common to represent metabolite intake as `( -> a_c0)` instead of `(a_e0 -> a_c0)`.
- Make sure reaction likeihood reactions appear in the universal model with the same identifiers.

A good rule of thumb is that if a model can be gapfilled using the `cobra.flux_analysis.gapfill` function with a universal model derived from one of the included templates, likelihood based gap-filling should also succeed. If you find this to not be the case, please get in contact with us using one of the contact methods listed below.

## Contact

If you have questions or are having issues getting started, or something is not working, please feel free to reach out. The best place to report an issue is on github at [https://github.com/PriceLab/probannopy/issues](https://github.com/PriceLab/probannopy/issues) for Probannopy and [https://github.com/kingb12/flask_probanno/issues](https://github.com/kingb12/flask_probanno/issues) for ProbannoWeb. You can also reach out directly via email at kingb12 'at' cs.washington.edu
