# BOSS(Building Optimized Scaffold graph for Scaffolding)

A novel algorithm is used for scaffolding contigs produced by any assembler. 

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


Scaffolder: BOSS
=================

1) Introduction

	BOSS is an scaffolder which aims to determine the orientations and orders of contigs. 
	The contigs can be produced by any assembler.

2) Before installing and running
	
	Before scaffolding, users should use one mapping tool to map left read library and right read libray to contigs.
	BOSS needs the bam files as input file, which can be produced by following commands:
	
	(a) Mapping tool: Bowtie2
	./bowtie2-build contigs.fa contigs
	./bowtie2 -x contigs left.fastq -S left.sam
	./bowtie2 -x contigs right.fastq -S right.sam
	./samtools view -Sb left.sam > left.bam
	./samtools view -Sb right.sam > right.bam
	
	(b) Mapping tool: Bowtie
	./bowtie-build contigs.fa contigs
	./bowtie -v 3 -q contigs left.fastq -S left.sam (or ./bowtie -v 0 -q contigs left.fastq -S left.sam)
	./bowtie -v 3 -q contigs right.fastq -S right.sam (or ./bowtie -v 0 -q contigs right.fastq -S right.sam)
	./samtools view -Sb left.sam > left.bam
	./samtools view -Sb right.sam > right.bam
	
	(c) Mapping tool: BWA
	./bwa index contigs.fa
	./bwa mem contigs.fa left.fastq > left.sam
	./bwa mem contigs.fa right.fastq > right.sam
	./samtools view -Sb left.sam > left.bam
	./samtools view -Sb right.sam > right.bam

3) Installing.

	BOSS is written by C++ and requires a machine with GNU C++ pre-installed.
	Create a main directory (eg:BOSS). Copy all source code to this directory.
	Type "make all" 

4) Running.

	Run command line: 
	"boss <contigs.fa> <bamfile_left.bam> <bamfile_right.bam> <read_length> <insert_size> <std> <min_weight> <min_number> <is_paired_end> <edge_weight_method> <scaffold_file_name>"
	<contigs.fa>: 
		The file includes contigs produced by one assembler.
	<bamfile_left.bam>:
		Before scaffolding, users should using one mapping tool(bowtie2, bwa or bowtie) to map left read library to contigs, and this will produce the bamfile_left.bam.
	<bamfile_right.bam>:
		Before scaffolding, users should using one mapping tool(bowtie2, bwa or bowtie) to map right read library to contigs, and this will produce the bamfile_right.bam.
	<read_length>: 
		The length of read.
	<insert_size>: 
		The insert size of read library.
	<std_percentage>: 
		The percentage of standard deviation to insert size, std = insert_size*std_percentage. In default, std_percentage = 0.07.
	<min_weight>: 
		One cutoff for removing suprious edgs in the scaffold graph. In default, min_weight = 0.2. If the coverage is larger than 100, the value of min_weight can set to be 0.3 or more large. Note that, min_weight is smaller than 1.
	<min_number>: 
		The minimum number of links between contigs. In default, min_number = 2. If the coverage is larger than 100, the value of min_number can set to be 3 or more large.
	<is_paired_end>: 
		It is equal to 0 or 1, 1 represents that read library is paired-end library, 0 represents that read library is mate-paired library.
	<edge_weight_method>: 
		It is equal to 0 or 1, 0 represents that the edge weight calculated by arithmetic mean, 1 represents that the edge weight calculated by geometric mean.
	<scaffold_file_name>: 
		The output file name, this file includes scaffolds produced by BOSS. 

5) Example:

	(a). If there is one read library:
	./boss contigs.fa left.bam right.bam 76 650 0.07 0.2 2 0 0 result
	This command will produce the scaffolding result: result_ScaffoldSet.fa
	
	(b). If there are two read libraries:
	./boss contigs.fa left1.bam right1.bam 76 650 0.07 0.2 2 0 0 left2.bam right2.bam 75 2700 0.07 0.2 2 0 0 result_com
	This command will produce the scaffolding result: result_com_ScaffoldSet.fa

6）Suggestion:

	When you want to enhance the accuracy of scaffolding results, please set large numbers for the parameters  <min_weight> and <min_number>;
	When the coverage number of paired read library is small (smaller than 30), you can try small numbers for the parameters <min_weight> and <min_number>;   
