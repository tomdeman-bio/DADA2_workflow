# DADA2_workflow
A simple workflow for processing 16s sequence data with DADA2


### Create the DADA2 workflow environment, with all its dependencies, using Conda package manager and the supplied yml file
```bash

conda env create -f dada2_workflow_environment.yml
```

### Operate the Bash workflow script like this:
```bash

./dada2_pipeline_vC1_v6.sh -i raw_gzipped_fastqs -r name_of_sequencing_run_or_project
```
