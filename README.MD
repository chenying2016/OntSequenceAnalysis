# </a> Contents
* [Introduction](#S-introduction)
* [Installation](#S-installation)
* [Quick Start](#S-quick-start)
* [Before Running](#S-before-running)
* [A Description of Config File](#S-desc-config-file)
* [A Description of Ouput Format](#S-desc-output-format)
* [Contact](#S-contact)

# <a name="S-introduction"></a> Introduction
OntCns is an iterative consensus tool for Nanopore long noisy reads. It is a customized version of [MECAT](https://github.com/xiaochuanle/MECAT). OntCns is able to correct noisy Nanopore long noisy reads that with error rate as high as 30%.

# <a name="S-installation"></a> Installation
The current directory is `/home/sysu/projects`:
``` shell
pwd
/home/sysu/projects
```

* Build source code
``` shell
cd ontcns/src
make
```

* Add paths
``` shell
export PATH=$PATH:/home/sysu/projects/ontcns/Linux-amd64/bin
export PATH=$PATH:/home/sysu/projects/ontcns/shell_scripts
```

# <a name="S-quick-start"></a> Quick Start
### Input Data Format
Both `FASTA` and `FASTQ` formated files are accepted by OntCns.

### Config Files
Before running OntCns, we shall prepare two config files.
* Read names file: A file that contains all the `FASTA` or `FASTQ` files that to be corrected. ***Each file name must be given with full path***. Here is an example:
``` shell
cat /home/sysu/test/cns_test/ont_read_list.txt
/home/sysu/smrt/data/read1.fq
/home/sysu/smrt/data/read2.fq
```
* Config file: A file that gives all the parameters for running OntCns. You can find an example of config file in `/home/sysu/projects/ontcns/shell_scripts/ontcns_cfg.txt`:
``` shell
cat /home/sysu/projects/ontcns/shell_scripts/ontcns_cfg.txt
PROJECT=EcoliConsensus
THREADS=4
ONT_READ_LIST=/home/sysu/test/cns_test/ont_read_list.txt
OVLP_FAST_OPTIONS="-n 200 -z 20 -e 0.2"
OVLP_SENSITIVE_OPTIONS="-n 200 -z 10 -e 0.5"
CNS_FAST_OPTIONS="-c 20 -e 0.2"
CNS_SENSITIVE_OPTIONS="-c 20 -e 0.5"
NUM_ITER=3
```

### Running OntCns
The current directory is `/home/sysu/test/cns_test`:
``` shell
pwd
/home/sysu/test/cns_test
ls
ontcns_cfg.txt  ont_read_list.txt
```
The command to run OntCns is very simple:
``` shell
OntCns2Main.sh ontcns_cfg.txt
```

### Where is the Corrected Reads?
In the example above, the Corrected Reads is
``` shell
/home/sysu/test/cns_test/EcoliConsensus/WorkDirIter3/WorkDirIter3Cns.fasta
```
The (sub)sequences that failed to be corrected is
``` shell
/home/sysu/test/cns_test/EcoliConsensus/WorkDirIter3/WorkDirIter3UnCns.fasta
```

# <a name="S-before-running"></a> Before Assembly 
Nanopore raw noisy reads can be very long and so are the corrected ones. Please note that `mecat2canu` from [MECAT](https://github.com/xiaochuanle/MECAT) has hardcoded the maximal length of input reads to be about 60,000. So before running `mecat2canu` for assembly, please run `OntCns2SplitLongReads` to partition long reads into short fragments:
``` shell
OntCns2SplitLongReads 500 50000 WorkDirIter3Cns.fasta WorkDirIter3Cns50000.fasta
```
Type
``` shell
OntCns2SplitLongReads
```
to learn its usage.

# <a name="S-desc-config-file"></a> A Description of Config File
* `PROJECT` The project name. All the working files will be put in this directory.
* `THREADS` Number of CPU threads used by OntCns.
* `ONT_READ_LIST` The reads file name. ***It must be given with full path.***
* `OVLP_SENSITIVE_OPTIONS` The options for running `OntCns2CanFinder` in the first iteration. `OntCns2CanFinder` is a tool for detecting consensus candidates. Just type
``` shell
OntCns2CanFinder
```
to learn its options and descriptions.
* `OVLP_FAST_OPTIONS` The options for running `OntCns2CanFinder` in the remainder iterations.
* `CNS_SENSITIVE_OPTIONS` The options for running `OntCns2Consensus` in the first iteration. `OntCns2Consensus` is tool for correcting noisy reads. Just type
``` shell
OntCns2Consensus
```
to learn its options and descriptions.
* `CNS_FAST_OPTIONS` The options for running `OntCns2Consensus` in the remainder iterations.
* `Num_ITER` Number of iterations.

# <a name="S-desc-output-format"></a> A Description of Output Format
The corrected reads will be output in a `FASTA` format file. The header of a corrected read has the following form:
``` shell
>ontcns/A/B/C_D_E
```
where:
* `A` is the id of the input reads to OntCns. The first raw read is `read 0`; the second raw read is `read 1`, and so on.
* `B` is the length of `read A`.
* `C` is the left-most effective position.
*  `D` is the right-most effective position.
*  `E` is the length of the corrected read.

# <a name="S-contact"></a> Contact
* Chuan-Le Xiao: xiaochuanle@126.com
* Ying Chen: chenying2016@gmail.com

# <a name="S-update-history"></a> Update History
* 2018.1.4 Version 0.0.1.
