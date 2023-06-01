#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=interproscan
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.error
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bs66@mailbox.sc.edu



#goal of the script is to add putative functions of gene models to .gff file
#run the script in the desired output directory

#load correct environments
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

conda activate interproscan

# Set versions
version_major=5.61
version_minor=93.0
CONDA_PREFIX=/work/bs66/software/miniconda/envs/interproscan

# get the md5 of the databases
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5

# get the databases (with core because much faster to download)
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz

# checksum
md5sum -c interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5

# untar gz
tar xvzf interproscan-${version_major}-${version_minor}-64-bit.tar.gz

# remove the sample DB bundled by default
rm -rf $CONDA_PREFIX/share/InterProScan/data/

# copy the new db
cp -r interproscan-${version_major}-${version_minor}/data $CONDA_PREFIX/share/InterProScan/


