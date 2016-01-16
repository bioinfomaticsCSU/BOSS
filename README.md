# BOSS
A new algorithm for scaffolding

License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


Scaffolder: BOSS(Building Optimized Scaffold graph for Scaffolding)
=================

BOSS is more advantageous for scaffolding read library whose ratio of read number to genome size is large.

1)Installing.

BOSS is written C++ and therefore will require a machine with GNU C++ pre-installed.

Create a main directory (eg:BOSS). Copy all source code to this directory.

Type "g++ main.cpp scaffoldgraph.cpp scaffolding.cpp -o boss ./lp/liblpsolve55.a -lm -ldl -I include/ -L lib/ -lbamtools" 

2)Running.

Run command line: 

"boss contigs.fa bamfile read_length insert_size std min_weight min_number min_cov_time is_paired_end scaffold_file_name"

contigs.fa is the file includes contigs produced by one assembler.

bamfile is the prefix of two mapping files. Before scaffolding, users should using one mapping tool(bowtie2, bwa or bowtie) to map left read file and right read file to contigs respectively, and this will produce two bam files which should be named "bamfile_1" and "bamfile_2".

read_length is the length of read.

insert_size is the insert size of read library.

std is standard deviation of insert size (such as, std = 0.1*insert_size).

min_weight is one cutoff for removing suprious edgs in scaffold graph (such as, min_weight = 0.3).

min_number is the minimum number of links between contigs (such as, min_number = 3).

min_cov_time is the coverage time cutoff, when sequence region's coverage is greater than min_cov_time times larger than average coverage, it will be regarded as repetitive region(such as, min_cov_time = 4). 

is_paired_end is 0 or 1, 1 represents that read library is paired-end library, 0 represents that read library is mate-paired.

scaffold_file_name is the output file name, this file includes scaffolds produced by BOSS. 

If there two read libraries, the command line:
"boss contigs.fa bamfile1 read_length1 insert_size1 std1 min_weight1 min_number1 min_cov1 is_paired_end1 bamfile2 read_length2 insert_size2 std2 min_weight2 min_number2 min_cov2 is_paired_end2 scaffold_file_name"

3)Output.

scaffold_file_name is the file (fasta format) which includes scaffolds produced by BOSS.
