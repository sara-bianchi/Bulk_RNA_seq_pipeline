# Bulk_RNA_seq_pipeline
Automated pipeline for the analysis of bulk gene expression data

INSTALLATION AND CONFIGURATION
All the packages required to run the pipeline have already been installed within a Docker image available on our repository as hedgelab/bulk_image:image4.
The use of Docker is recommended to ensure reproducibility of the results. In case you don't have Docker installed on your computer, you can install it following the instructions: https://docs.docker.com/get-started/get-docker/
Otherwise, it is possible to recreate the working environment following the installation commands in the provided Dockerfile (optimized for a bash shell).

Ensure that the RAM and CPU allocated to Docker are sufficient to run the analysis. If you wish to run only the analysis of already aligned counts 10 GB of RAM will be enough. Else, if you need to perform either trimming or alignment, at least 30 GB is required. In case you do not have an already-indexed genome, you will need to run the indexing step, which requires at least 100 GB of RAM. To set the proper RAM limits for your analysis some versions of Docker Desktop allow to directly edit these parameters in the settings, others require the use of an external settings file. In the latter case create a .wslconfig file in User/user_name folder of your computer like the one in the example below:
# Settings apply across all Linux distros running on WSL 2
[wsl2]
# Limits VM memory to use no more than 2 GB, this can be set as whole numbers using GB or MB
memory=100GB 
# Sets the VM to use two virtual processors (max 18 cores 36 threads)
processors=12

Once Docker is installed and properly configured, activate the docker engine direcly opening the Docker Desktop application or by running:
docker start
on your terminal
