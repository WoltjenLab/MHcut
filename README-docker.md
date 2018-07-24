## Docker

Docker simplifies the installation of tools by preparing an image of an OS with all the requirements and dependencies to run the tool.
No need for the user to install many different libraries and tools on their system, the docker image can be used directly.
The user still has to [install Docker](https://docs.docker.com/install/) though.

### Docker crash course

To run a command `command arg1 arg2 ...` within a docker container, a typical docker command looks like this:

```shell
docker run -v PATHL:PATHC -w PATHC imageName command arg1 arg2 ...
```

- `-v PATHL:PATHC` links the local folder `PATHL` to the folder `PATHC` in the container. Typically `PATHL` contains the input files and will be where we want the output files written.
- `-w PATHC` means that the working directory (where the command will be run) is `PATHC`, i.e. we want to run the command in the folder that we linked with the previous parameter (and that contains the input files).
- `imageName` is the name of the docker container. If not built manually, Docker will try to download it from Docker Hub or other platforms.

Hence for us, a dockerized command might look like that:

```shell
docker run -v `pwd`:/home -w /home jmonlong/mhcut python /root/MHcut.py -var clinvar-grch38-all-deletion.tsv -ref hg38.fa -out docker-test
```

We link the current folder (`` `pwd` ``) with a *home* folder in the container that we will use as working directory and run the python command.
The MHcut scripts are located in the `/root` folder of the container, hence the `python /root/MHcut.py`.

## Docker workflow for MHcut

*Note: this assumes that the reference genome file was downloaded and the input TSV file is ready.*

### Using BLAST

```shell
## Prepare for BLAST
docker run -v `pwd`:/home -w /home jmonlong/mhcut makeblastdb -in hg38.fa -dbtype nucl -title hg38
## Run MHcut
docker run -v `pwd`:/home -w /home jmonlong/mhcut python /root/MHcut.py -var clinvar-grch38-all-deletion.tsv -ref hg38.fa -out docker-test
```

## Using JellyFish

```shell
## Count 23-mers with JellyFish
docker run -v `pwd`:/home -w /home jmonlong/mhcut jellyfish count -m 23 -s 100M hg38.fa
## Run MHcut
docker run -v `pwd`:/home -w /home jmonlong/mhcut python /root/MHcut.py -var clinvar-grch38-all-deletion.tsv -ref hg38.fa -out docker-test -jf mer_counts.jf
```

*Note: `jellyfish count` takes ~1h to run on the human genome.*

## Temporary: build the image manually

Once MHcut becomes public, I will link the container to a docker database like [Docker Hub](https://hub.docker.com/).
Then, the user will be able to download the container directly (should be automatic the first time that a docker command is called).

In the meantime, we can build the container manually by running the following command in the folder containing the `Dockerfile`:

```shell
docker build -t jmonlong/mhcut .
```
