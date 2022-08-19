# MetaCompare2.0

MetaCompare2.0 is an updated version of MetaCompare. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine or server.

### System Requirements (*tested on linux Ubuntu 14.04*)

* git installed
* Python3 with `pandas` and `biopython` package installed
  * `pip install pandas`
  * `pip install biopython`
  * `pip install pprodigal`

* DIAMOND (https://github.com/bbuchfink/diamond)
* MMseqs2 (https://github.com/soedinglab/MMseqs2)

### Installing

**Step 1:** Change the current working directory to the location where you want the cloned `MetaCompare2.0` directory to be made.
**Step 2:** Clone the repository using git command
```
~$ git clone https://github.com/mrumi/MetaCompare2.0.git
```

**Step 3:** Make directory `BlastDB` under `MetaCompare2.0` directory and change the woring directory to it.

```
~$ cd MetaCompare2.0
~/MetaCompare2.0$ mkdir BlastDB
~/MetaCompare2.0$ cd BlastDB
~/MetaCompare2.0/BlastDB$
```

**Step 4:** Download the compressed Blast Database file from the web server and uncompress it.

```
~/MetaCompare2.0/BlastDB$ wget http://bench.cs.vt.edu/ftp/data/metacomp/BlastDB.tar.gz
~/MetaCompare2.0/BlastDB$ tar -zxvf BlastDB.tar.gz
```
*Note:* If you are encountering `ERROR: cannot verify bench.cs.vt.edu's certificate`, execute wget command with `--no-check-certificate` option (`wget http://bench.cs.vt.edu/ftp/data/metacomp/BlastDB.tar.gz --no-check-certificate`).

**Step 5:** Get back to working directory `MetaCompare` and run `metacmp.py`

```
~/MetaCompare2.0/BlastDB$ cd ..
~/MetaCompare2.0$ python metacompare.py
```

## Running MetaCompare2.0

### Preparing input files

MetaCompare2.0 requires one FASTA file as input and that is assembled contigs generated by a standard assembler. 

### Running

Suppose you have an assembled contigs file, `S1.fa` (*This file is already in your working directory*).

The following command runs MetaCompare2.0 with 128 threads.

```
~$ python metacompare.py -c S1.fa 
```
The output should be look like as follows:
```
Running Prodigal
Running blastn on ACLAME
Running blastx on CARD
Running blastn on PATRIC
Reading files...
Computing resistome risk score..
Resistome risk score: 38.64014990951873
```

You can see detailed description for command line options by using `-h` option.
```
~$ ./metacmp.py -h
```

## License

MIT License

Copyright (c) 2018 Min Oh (minoh@vt.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
