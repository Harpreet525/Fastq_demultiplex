## Description

The script is adapted from copilot and is used to extract the undetermined reads from the paired end illumina sequencing with unique dual indexes.

The process begins by extracting the index1 and index2 sequences from the undetermined paired-end FASTQ files. It then aligns index1 and the reverse complement of index2 with the user-provided reference sequences for index1 and index2.

Initially, it attempts direct alignment, allowing either a perfect match or up to two mismatches. If no match is found, it performs frameshift matching by shifting one base to the left and right with a tolerance of one mismatch, or by shifting two bases with no mismatches allowed.

When both index1 and index2 satisfy these conditions, the reads are saved to a new FASTQ file.

### Usage

Download or copy the script to your space.

Run the command below to compile the script. Needs only basic gcc libraries:

**g++ -std=c++17 -o output_name_of_file input_name_of_script.cpp -pthread**

Then you can run the compiled binary with the command below:

**./output_name_of_file tab_file.txt read1.fastq read2.fastq output_dir threads(int)**

Example:

./your_compiled_binary path_to_your_tab_file.txt path_to_your_read1.fastq path_to_your_read2.fastq path_to_your_output_directory 50

### Requirements

The tab file must consists of 3 columns without header. The 1st column will have your sample id, 2nd column will have sequence of index1 and 3rd column will have sequence of index2. The tab file can be obtained from your sample sheet.

The name of the output of fastq file will be same as the name of the sample id for the respective indexes.

col1        col2      col3

sample1   index1  index2

sample2   index1  index2

sample3   index1  index2


**No strict dependencies. Just basic c++ installation.**
