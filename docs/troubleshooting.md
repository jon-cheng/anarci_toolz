# Troubleshooting: installing ANARCI and AbNumber without Docker

The [Docker route in the main README](../README.md#docker) is the recommended way to get ANARCI and AbNumber running, since it sidesteps everything below. This page is for cases where Docker isn't an option and you need a bare-metal or conda environment instead.

## Requirement Note: Troubleshooting ANARCI installation

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

## Installing ANARCI and AbNumber from source code

As an alternative to conda, you may install ANARCI and Abnumber from source code:

### ANARCI
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

### AbNumber
```Dockerfile
# Add to Dockerfile
RUN git clone https://github.com/prihoda/AbNumber.git \
    && cd AbNumber \
    && python3 setup.py install
```

Note: the repo's own `Dockerfile` follows this same approach, but installs AbNumber with `pip install --no-deps` rather than `python3 setup.py install` — the latter shells out to `easy_install` to resolve AbNumber's declared dependencies, which fails against modern sdist-less pandas releases.
