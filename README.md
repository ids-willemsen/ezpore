<img src="ezpore.png" alt="image info" width="35%">

# ezpore

Authors: Robbert van Himbeeck, Ids Willemsen

## About

**`ezpore` is still in the test phase, please use with caution.**

`ezpore` is a single-command pipeline to process bacterial (full 16S), fungal (full ITS) or Nematodal (full 18S) reads obtained by Nanopore sequencing. This pipeline is developed and intended for internal use by the Laboratory of Nematology (WUR). 

`ezpore` can perform following steps:

1) demultiplexing (dorado, barcoding kit EXP-NBD196)
2) filtering on length and quality (NanoFilt)
3) primer trimming (cutadapt)
4) ITS region extraction (ITSxpress, for fungal ITS)
5) cluster reads (vsearch)
6) read classification (emu)

<img src="ezpore_Diagram.png" alt="image info" width="70%">


Chimera detection and removal will be incorporated in the future.

This pipeline is inspired by the `decona` pipeline we use for nematode analysis.
More information about `decona` can be found here: https://github.com/Saskia-Oosterbroek/decona 





## Installation & prerequisites

`ezpore` is developed for Linux operating systems and will likely also work on other Unix-like OS (e.g. MacOS) but this has not been tested. 

Usage on Windows is not supported, however Windows Subsystem for Linux (WSL) can be used (see section "WSL installation instructions").

Be aware that running on windows takes way longer than on a Linux machine!

### ezpore installation instructions

To use `ezpore`, `conda` needs to be installed on your system. To install `conda` perform following steps:

1) 
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh #
```
2)
```
bash Miniconda3-latest-Linux-x86_64.sh #
```
navigate trough the interactive installation shell: Choose **yes** by “Do you wish the installer to initialize Miniconda3 by running conda init? [yes|no]" 
   
3) open a new terminal, (base) will appearing at the beginning of every rule.

To install `ezpore`:

1) download the whole repository, or clone it in your directory using:
```
git clone https://git.wur.nl/robbert.vanhimbeeck/ezpore
```
2) enter the ezpore directory by:
```
cd ezpore 
```
Check the content.

3) execute the `install.sh` file to setup and create the correct `conda` environment. 
Run the file by running:
```
./install.sh
```

After finishing you can check if the install worked by running 

```
conda activate ezpore
```

If you now see (ezpore) in stead of (base) in the front of the rule, the environment was succesfully created.

_note:_ 
Make sure that you've also downloaded/cloned the required reference databases. For 16S bacteria, the "16S_silva_full" folder can be downloaded and for ITS bacteria the "UNITE_ITS_emu" folder, both contains 3 files.
 


### WSL installation instructions

Please find below a small explaination of setting up WSL on windows:
1) install Ubuntu or WSL via the Microsoft Store (Windows 10)

If you never did this before, follow this tutorial, because Windows automatically generates resolv.conf file with the wrong nameserver.
https://stackoverflow.com/questions/62314789/no-internet-connection-on-wsl-ubuntu-windows-subsystem-for-linux

2) locate the file by running the following command:
```
sudo nano /etc/resolv.conf
```
You will see the following in the file:

```
#This file was automatically generated by WSL. To stop automatic generation of this file, add the following entry to /etc/resolv.conf
#[network]
#generateResolvConf = false
nameserver xxx.xx.xx
```
3) change the nameserver value to 8.8.8.8 and save the file. You should now be able to connect to the internet.

4) if you are able to connect to the internet now then you may also need to stop WSL from resetting this file when opening future terminals. You can do that by running these commands:
```
sudo rm /etc/resolv.conf
sudo bash -c 'echo "nameserver 8.8.8.8" > /etc/resolv.conf'
sudo bash -c 'echo "[network]" > /etc/wsl.conf'
sudo bash -c 'echo "generateResolvConf = false" >> /etc/wsl.conf'
sudo chattr +i /etc/resolv.conf
```
4) follow the steps of section "ezpore installation instructions".



## Usage

### Your environment and folder

Always make sure the `conda` environment `ezpore` is activated before you run the analysis. 
1) run to activate the environment 
```
conda activate ezpore
```
2) go and make a folder (e.g. ezpore_analysis_1) where you do the analysis by 
```
mkdir ezpore_analysis_1 && cd ezpore_analysis_1
```
This folder should contain following items:

1) a single input .fastq file containing basecalled data, e.g. calls.fastq from the basecalling (if demultiplexing is required) **OR** if your data is already demultiplexed, the data should be in a folder called "demultiplexed". **Note: If you have more than 96 samples, you should demultiplex the samples and rename the barcode fastq files to your sample names before running ezpore. You can then put the fastq files with your sample names in the "demultiplexed" folder and run `ezpore` without demultiplexing.**

2) the `settingsfile.txt` file, which contains the settings ezpore will use in your run.

3) the `barcode_files.txt` file, which specifies the barcodes that are included in your run. 
Change the content of the file to your own needs, so indicate witch barcodes you used! Don't change the formatting!


### `ezpore` usage and settingsfile
 
The `ezpore` pipeline is run using the `nem` command followed by the location of your `settingsfile.txt`. For example, if you run ezpore from the folder that contains your settingsfile the usage would be: 

```
nem settingsfile.txt
```

# The settingsfile

The settingsfile.txt contains all possible arguments that can be used by ezpore. Examples of how to use the settingsfile.txt for both demultiplexed and non-demultiplexed data, different organisms and other settings are given in each testrun folder. 

If demultiplexing is required (`-demultiplex TRUE`), the input file of `ezpore` should be a single .fastq file containing all data.

If your data is already demultiplexed (`-demultiplex FALSE`), this data should be present in the directory in a folder named "demultiplexed". 
The demultiplexed files within this folder, should be named as "barcodeXX.fastq", where XX has to be substituted by the correct barcode number. For example: barcode01.fastq or barcode21.fastq. (This feature is not tested yet). 

The `settingsfile.txt` takes following arguments:

| argument | description | input type | default value |
| -------- | ----------- | ------------  | ------------- |
|-demultiplex | demultiplexes the data using dorado | TRUE/FALSE | TRUE |
|-min | the minimum read length (in bp). Shorter reads are removed | integer | 100 | 
|-max |the maximum read length (in bp). Larger reads are removed | INTEGER | 10000 |
|-quality | the minimum average read quality to be retained. Reads with lower Q score are removed | INTEGER | 15 |
|-trim_primers | removes primers using cutadapt | TRUE/FALSE | TRUE |
|-primer_error_rate | the maximum allowed error rate for primer trimming. | UNIT INTERVAL[0-1} | 0.2 |  
|-min_abundance | the minimum relative abundance of an organism to be retained by emu | UNIT INTERVAL[0-1] | 0.0001 |
|-cluster_perc | the percentage identity to cluster on using vsearch. After clustering, consensus sequences are rereplicated for emu classification. If FALSE, no clustering is performed | UNIT INTERVAL[0-1} | FALSE |
|-rank | the taxonomic rank which emu uses to combine output of all files | not functional ATM | species | 
|-threads | the number of threads that emu uses for classification | INTEGER | 2 |
|-count_table | If TRUE, the read counts will be outputted. If FALSE, the relative abundances will be outputted. | TRUE/FALSE | FALSE|
|-group | the group of organisms: bacteria (16S_bac), nematodes (18S_nem) or fungi (ITS_fungi) | STRING | none |
|-barcode_file | .txt file containing the barcodes you used | .txt file | none |
|-input_file | the input file (.fastq) of the analysis | .fastq file | none | 


## Supported primers
At this moment, only the 27F (5'-AGAGTTTGATCMTGGCTCAG-3') and 1492R (5'-CGGTTACCTTGTTACGACTT-3') universal primers are supported for 16S bacterial analysis.

For 18S nematode primers, only ............

For fungi, only the ITS4ngsUni and ITS9MUNngs primers are supported.

If required, more primers can be added to the pipeline. Please open an issue if you would like to requests additional primer support.

## Test runs

It is strongly adviced to run the script with with both the 16S and ITS test dataset file as input, to see if everything is functioning properly.

### Test run for 16S

From the ezpore installation folder, go into the `testrun_bacteria` directory with
```
cd testrun_bacteria
```
This directory contains the 16S dataset which contains multiplexed data of barcode 57, 58, 59, and 60. The data is not yet demultiplexed. With demultiplexing, you will always detect a negligible amount of reads of non-used barcodes. We specify argument -barcode_file to retain only the included barcodes. Open a terminal in that directory, and run following command.

The `settingsfile.txt` located in the same folder contains our suggested settings to analyze this dataset. 

To start the testrun on the 16S dataset you can use the following command:

```
nem settingsfile.txt
```
### Test run for ITS
From the ezpore installation folder, go into the `testrun_bacteria` directory with
```
cd testrun_fungi
```

The `barcode58.fastq` file is already demultiplexed and contains only the data of barcode58. Therefore, there is no need to demultiplex the data (so `-demultiplex FALSE`). 

Please note!
As the data is already demultiplexed, the `barcode58.fastq` file is inside a folder called "demultiplexed". 
When running the script, ezpore will ask you to confirm this. Type 'yes' when asked. 

The `settingsfile.txt` located in the same folder contains our suggested settings to analyze this dataset. 

To start the testrun on the 16S dataset you can use the following command:
```
nem settingsfile.txt
```
Ignore the QIIME2 warning, because in this pipeline we are not using it.
###

If this completes without an error, you are ready to start analysing your own data!
 
### Results

Within your directory, a folder named "all_results_emu" will be created and your classified OTU tables will be stored there. Both a merged OTU table (emu-conbined-silva(-counts).tsv) and the separate barcode files will be given as output.

### Bugs and requests
If you encounter any bugs or you wish to request additional features, please open an issue on this GitLab page.
