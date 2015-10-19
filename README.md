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


Scaffolder
=================

1)Installing.

BOSS is written C++ and therefore will require a machine with GNU C++ pre-installed.

Create a main directory (eg:BOSS). Copy all source code to this directory.

Type "g++ main.cpp scaffoldgraph.cpp scaffolding.cpp -o boss ./lp/liblpsolve55.a -lm -ldl -I include/ -L lib/ -lbamtools" 

2)Running.

Run command line: 

"boss contigs.fa bamfile read_length insert_size std min_weight min_number min_cov is_paired_end out_put_path"

contigs.fa is the file includes contigs produced by one assembler.

read_length is the length of read.

insert_size is the insert size of read library.

std is standard deviation of insert size (default: 0.1*insert_size).

min_weight is one cutoff for removing suprious edgs in scaffold graph (default: 0.3).

min_number is the minimum number of links between contigs (default:5).

min_cov is the coverage cutoff, when sequence region's coverage is greater than min_cov, it will be regarded as repetitive region. 

is_paired_end is 0 or 1, 1 represents that read library is paired-end library, 0 represents that read library is mate-paired.

out_put_path is the path which is used for output. 


3)Output.

ScaffoldSet.fa is the file which includes scaffolds produced by LNSG.
