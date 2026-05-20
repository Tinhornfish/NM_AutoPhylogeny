# AutoPhylogeny
**This script is for Biology/Zoology students of the University of Veterinary Medicine Budapest to automatize one of the more outdated homework assignments.**

## Setup

The workflow uses containers so that there are as few deppendencies as possible, however one has to install the **workflow manager** and a **container engine**.

### To start using the workflow, clone it to your machine using
`git clone https://github.com/Tinhornfish/NM_AutoPhylogeny.git`

### Then you need the workflow manager called **Nextflow**
A link to Nextflows website: https://www.nextflow.io/
To install nextflow, use:
```
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
# Confirm if nextflow is installed correctly:
nextflow info
```
### As for container engines you can use Docker or Apptainer (aka. Singularity)
https://www.docker.com/
https://apptainer.org/

## How to use
You just need to modify the first couple of lines in the autophylo.nf script to include the **species**and **genes** of your choice.
If you've given your desired species and genes, run the workflow with the following command:
```
nextflow run autophylo.nf -profile standard -resume
```
