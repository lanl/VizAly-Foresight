# Z-checker installer

 (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.

See COPYRIGHT in top-level directory.


Authors: Sheng Di, Hanqi Guo

## 3rd party libraries/tools

- gcc
- perl
- git
- texlive
- ghostscript(gsview) (z-checker-install.sh can install it automatically if missing)
- latexmk (z-checker-install.sh will install latexmk automatically if missing)
- gnuplot (z-checker-install.sh will install gnuplot automatically if missing)

## Testing/Installation method

z-checker-install.sh will download latexmk, gnuplot, Z-checker, ZFP, and SZ and install them one by one automatically, and then add the patches to let ZFP and SZ fit for Z-checker.

After installation, please download the two testing data sets, CESM-ATM and MD-simulation (exaalt). The two data sets are available only for the purpose of research of compression. Please ask for the data by contacting [sdi1@anl.gov]() if interested.

### Quick Start

Then, you are ready to conduct the compression checking.
You can generate compression results with SZ and ZFP using the following simple steps: 
(Note: you have to run z-checker-install.sh or z-checker-install2.sh to install the software before doing the following tests)

1. Create a new test-case, by executing "createNewZCCase.sh [test-case-name]". You need to replace [test-case-name] by a meaningful name.
   For example:
   [user@localhost z-checker-installer] ./createNewZCCase.sh CESM-ATM-tylor-data

2. Perform the checking by running the command "runZCCase.sh": runZCCase.sh [test-case-name] [data dir] [dimensions....].
   Example:
   [user@localhost z-checker-installer] ./runZCCase.sh CESM-ATM-tylor-data /home/shdi/CESM-testdata/1800x3600 3600 1800

Then, you can find the report generated in z-checker-installer/Z-checker/[test-case-name]/report.

### Step-by-step Checking

Unlike the above one-command checking, the following steps present the generation of compression results step by step.

1. Go to zfp/utils/, and then execute "zfp-zc-ratedistortion.sh [data directory] [dimension sizes....]". The compression results are stored in the compressionResults/ directory.
   For example, suppose the directory of CESM-ATM data set is here: /home/shdi/CESM-testdata/1800x3600, then the command is "zfp-zc-ratedistortion.sh /home/shdi/CESM-testdata/1800x3600 3600 1800". Note: the data files stored in the directory are also ending with .dat and the dimension sizes are the same (1800x3600) in this test-case.

2. Similarly, go to SZ/example/, and then generate compression results by SZ compressor as follows: "sz-zc-ratedistortion.sh [data directory] [dimension sizes....]". The compression results are stored in the compressionResults/ directory.
   As for the example CESM-ATM, the test command is "sz-zc-ratedistortion.sh /home/shdi/CESM-testdata/1800x3600 3600 1800".

3. Then, go to Z-checker/examples/ directory, and run the command "./analyzeDataProperty.sh [data directory] [dimension sizes....]" to generate the data properties based on the data sets. This step has nothing to do with the compressors. The data analysis results are stored in the dataProperties/ directory. 

4. Generate the figure files: run the command "./generateReport.sh" simply. The results of comparing different compressors (such as sz and zfp in this test-case) are stored in the directory called compareCompressors/.

### Create a new case

"createNewZCCase.sh [test-case-name]" allows you to create a new test-case.  This command will create a new workspace directory in Z-checker, SZ, and zfp respectively. The compression results will be put in those workspace directories to avoid bing messed with other test-cases.

For example, if you run the generateReport.sh in the directory ./Z-checker/examples, it is actually one test case, where the compression results and data analysis results will be put in the dataProperty/ and compressionResults/ under it.
For another test case with another set of data or application, you can create a new workspace directory by the script createNewZCCase.sh (which calls ./Z-checker/createNewCase.sh).

#### z-checker-update.sh

z-checker-update.sh can be used to update the repository (pull the new update from the server), so that you don't have to perform the update manually.

#### web installation

Web installation allows to install a web server on the local machine, such that you can visualize the data through a local webpage and other people can view the data/results via that page if public ip is provided. 
z-checker-web-install.sh
