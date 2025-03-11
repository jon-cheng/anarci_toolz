![Version](https://img.shields.io/badge/anarci--toolz%20-%20version%200.1.1-brightgreen)

# anarci_toolz

[ANARCI](https://academic.oup.com/bioinformatics/article/32/2/298/1743894?login=true): Antigen Receptor Numbering and Receptor Classification is an antibody sequence numbering tool. 
Additionally, [AbNumber](https://abnumber.readthedocs.io/en/latest/) is a convenient wrapper of ANARCI that provides sequence annotations such as CDR and FWR definitions.

While the base [ANARCI tool](https://github.com/oxpig/ANARCI) and AbNumber are highly popular tools in the antibody engineering community, they may be somewhat difficult to use, especially for batch processing of large numbers of sequences.

Here, we present `anarci-toolz`, a performant convenient wrapper for ANARCI and AbNumber antibody sequence aligners based on ANARCI. The tool combines annotations from both tools in a consolidated output table. We provide improved usability of the ANARCI tool by providing support for tabular (`.csv` or DataFrame) inputs and outputs. The user only has to provide a table of amino acid sequences. This package utilizes Python `multiprocessing` to leverage some or all available CPU cores on an instance, which significantly reduces run times than base ANARCI without multiprocessing.

## Installation

```Bash
pip install git+https://github.com/jon-cheng/anarci_toolz.git
```

### Requirements

**Python 3.10**+

Additional dependencies outside of Python not included in the `pip install`:

- ANARCI
```Bash
conda install --yes -c bioconda anarci
```

- AbNumber
```Bash
conda install --yes -c bioconda abnumber
```

You should install `anarci-toolz` within the same conda environment or Docker image as ANARCI and AbNumber.


## Usage
### Input
At a minimum, you'll need this input:

- Antibody numbering scheme
- Directory of csv file(s) (in CLI mode) OR a Pandas Dataframe of AA sequences (in Python import mode) - files should contain at least one column containing AA sequences of a user-defined column name as input for the flag `seq_aa_header`
- Species

#### DNA Mode
Corresponding DNA sequences that accompany the amino acid sequence may also be input. `anarci-toolz` will linearly slice CDR and FR regions based on the ANARCI-derived amino acid result and apply that onto the DNA sequence.


### Arguments

| Command Prefix | Type | Required | Default | Choices | Description |
| --- | --- | --- | --- | --- | --- |
| `--scheme` | str | Yes | N/A | "imgt", "chothia", "kabat", "aho" | Numbering scheme: One of 'imgt', 'chothia', 'kabat', 'aho' |
| `--seq_aa_header` | str | Yes | N/A | N/A | The name of the column containing the amino acid sequences |
| `--seq_dna_header` | str | No | None | N/A | The name of the column containing the DNA sequences |
| `--allowed_species` | str | Yes | N/A | "human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "None" | Allowed species for germline assignment. Use "None" to allow all species or one or multiple of 'human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca' |
| `--retain_indices` | N/A (flag) | No | False | N/A | Whether to retain linear index columns in sequence annotation; defaults to False. Refers to the `_start` and `_end` suffixes that demarcate the boundaries of the variable regions and CDRs and frameworks on purely linear indexing, independent of numbering scheme. |
| `--skip_run_base_anarci` | N/A (flag) | No | False | N/A | Whether to skip running base ANARCI for detailed alignment quality metrics like e-value; reduces run-time; defaults to False |
| `--display_residue_view` | N/A (flag) | No | False | N/A | Whether to display wideform view showing all residue positions as columns; defaults to False |
| `--display_unformatted_residue_view` | N/A (flag) | No | False | N/A | Whether to override default formatted residue view that includes scheme prefix, only displaying the scheme residue positions as columns; defaults to False |


### Output
The schema of the output DataFrame:


| Field | Description | Required |
| --- | --- | --- |
| `scheme` | Input scheme | True |
| `passed_abnumber` | Whether AbNumber call was successful | True |
| `sequence_alignment_aa` | Extracted variable region sequence determined by AbNumber | True |
| `species` | Species determined by AbNumber | True |
| `chain_type` | One of H, K, L | True |
| `v_gene` | V-gene call | True |
| `j_gene` | J-gene call | True |
| `cdr1_aa` | CDR1 amino acid sequence | True |
| `cdr2_aa` | CDR2 amino acid sequence | True |
| `cdr3_aa` | CDR3 amino acid sequence | True |
| `fr1_aa` | FR1 amino acid sequence | True |
| `fr2_aa` | FR2 amino acid sequence | True |
| `fr3_aa` | FR3 amino acid sequence | True |
| `fr4_aa` | FR4 amino acid sequence | True |
| `cdr1_dna` | CDR1 DNA sequence | False |
| `cdr2_dna` | CDR2 DNA sequence | False |
| `cdr3_dna` | CDR3 DNA sequence | False |
| `fr1_dna` | FR1 DNA sequence | False |
| `fr2_dna` | FR2 DNA sequence | False |
| `fr3_dna` | FR3 DNA sequence | False |
| `fr4_dna` | FR4 DNA sequence | False |
| `passed_anarci` | Whether ANARCI call was successful | False |
| `variable_region_start_index` | Start index of the variable region | False |
| `variable_region_end_index` | End index of the variable region | False |
| `e_value` | E-value from ANARCI | False |
| `bitscore` | Bit score from ANARCI | False |
| `bias` | Bias from ANARCI | False |

### Use Cases


#### Importing the `run_anarci_tools` function

See also the `demo.ipynb` notebook.

```Python
import pandas as pd
from anarci_tools.pipeline import run_anarci_tools


df = pd.read_csv(
    "/path/to/myfiles/myfile.csv"
)
df_result = run_anarci_tools(
    df=df,
    scheme="imgt",
    allowed_species=["human"],
    seq_aa_header="sequence_aa",
)
```

#### Run Command Line tool

The Command Line tool can take multiple csv inputs, and will generate one or multiple csv table outputs into a subdirectory `anarci_annot`.

##### Example 1: Using required flags only
```Bash
anarci-toolz \
    --input "/path/to/myfiles/" \ 
    --scheme "imgt" \
    --seq_aa_header "sequence_aa" \
    --allowed_species "human" \
```

##### Example 2: Disable base ANARCI, saving run time
```Bash
anarci-toolz \
    --input "/path/to/myfiles/" \ 
    --scheme "imgt" \
    --seq_aa_header "sequence_aa" \
    --seq_dna_header "sequence_dna" \
    --allowed_species "human" \
    --skip_run_base_anarci
```

##### Example 3: Enable search of multiple species' databases
```Bash
anarci-toolz \
    --input "/path/to/myfiles/" \ 
    --scheme "imgt" \
    --seq_aa_header "sequence_aa" \
    --seq_dna_header "sequence_dna" \
    --allowed_species "human" "mouse" "rat" \
    --skip_run_base_anarci
```

##### Example 4: Enable search of ALL available species' databases
```Bash
anarci-toolz \
    --input "/path/to/myfiles/" \ 
    --scheme "imgt" \
    --seq_aa_header "sequence_aa" \
    --seq_dna_header "sequence_dna" \
    --allowed_species "human" "mouse" "rat" \
    --skip_run_base_anarci
```

##### Example 5: Run in DNA-slicing mode, enabling the ANARCI-generated regions to be sliced on corresponding input DNA sequences 
```Bash
anarci-toolz \
    --input "/path/to/myfiles/" \ 
    --scheme "imgt" \
    --seq_aa_header "sequence_aa" \
    --seq_dna_header "sequence_dna" \
    --allowed_species "human" \
```

##### Example 6: Display residue view
Display whole set of all possible residue positions in the sequence set in a wideform table:
```Bash
anarci-toolz \
    --input "/path/to/myfiles/" \ 
    --scheme "imgt" \
    --seq_aa_header "sequence_aa" \
    --allowed_species "human" \
    --display_residue_view \
```


## Addendum

### Requirement Note: Troubleshooting ANARCI installation
ANARCI has a dependency on HMMER. From HMMER's documentation: "HMMER requires Intel/AMD compatible machines, Apple OS/X Intel or ARM machines." For example, a direct install on Apple Mac Silicon machines will not work. If you are attempting to install ANARCI locally on a Mac using conda, you may run into this problem:

```
conda install --yes -c bioconda anarci
Retrieving notices: ...working... done
Channels:
 - bioconda
 - defaults
Platform: osx-arm64
Collecting package metadata (repodata.json): done
Solving environment: failed
 
LibMambaUnsatisfiableError: Encountered problems while solving:
  - nothing provides hmmer >=3.1 needed by anarci-2020.04.23-py_0
 
Could not solve for environment specs
The following package could not be installed
└─ anarci is not installable because it requires
   └─ hmmer >=3.1 , which does not exist (perhaps a missing channel).
```

This is due to ANARCI being incompatible with the Apple Silicon osx-arm64 architecture. ANARCI depends on HMMER. From HMMER's documentation: "HMMER requires Intel/AMD compatible machines, Apple OS/X Intel or ARM machines."

Here is one work-around: use Apple's [Rosetta](https://support.apple.com/en-us/102527) to emulate the Intel architecture osx-64 needed to run ANARCI. 

1. Install Rosetta:
```Bash
/usr/sbin/softwareupdate --install-rosetta --agree-to-license
```

2. Apply an environment variable: conda setting to run osx-64
```Bash
export CONDA_SUBDIR=osx-64
```

3. Create a new conda environment for running on osx-64
```Bash
conda create -n <env_name> python=3.10 # or whatever Python version you need
  
conda activate <env_name>
```

4. Try the ANARCI conda installation: 
```Bash
conda install --yes -c bioconda anarci
```

5. When you want to revert to your Mac's native M2 architecture, e.g. for things other than ANARCI, you can run: 
```Bash
unset CONDA_SUBDIR
```



### Installing ANARCI and AbNumber from source code

As an alternative to conda, you may install ANARCI and Abnumber from source code:

#### ANARCI
```Dockerfile
# Add to Dockerfile
RUN apt update && apt install -y \
    hmmer \
    && pip install biopython \
    && git clone https://github.com/oxpig/ANARCI.git \
    && cd ANARCI \
    && python3 setup.py install
ENV PATH="${PATH}:/ANARCI/bin"
```

#### AbNumber
```Dockerfile
# Add to Dockerfile
RUN git clone https://github.com/prihoda/AbNumber.git \
    && cd AbNumber \
    && python3 setup.py install
```
