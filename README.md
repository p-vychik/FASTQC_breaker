# FastQC Breaker

### Hello!
This repository is for our Python project in Bioinformatics Institute. Here we created a program similar to [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool. 

### What is the FastQC?
FastQC is a tool for quality control of high throughput sequence raw data. This program takes `.fastq` files and creates html report with some characteristics of reads quality.  
Program makes:
- summary for input data
- graphics or tables for some data characteristics
- status (good, warning or fail) for each characteristic
- HTML report that contains all the above

### Files in this repository

- `fastqc_breaker.py` - it is the main script that is calling all of the functions and creates report;  
- `sequences.py` - class that parses `.fastq` file, creates 2D numpy arrays and make them suitable for subsequent work;  
- `reports.py` - file that contains functions for each report block. Each function makes picture or table and returns status;  
- `html_creator.py` - file that creates html report. There is template and pictures for it in the folder `HTML`.

### How to run the program on your computer:

The program is written and tested in Python 3.9

1. Clone this repository using `git clone` command
2. Install necessary requirements using next command:  
``` pip3 install -r requirements.txt ```
3. To run program use this:  
``` fastqc_breaker.py -i <fastq-file-name> -o <output-dir-name> ```  

    where:  
    `<fastq-file-name>` - input file (with path if it is nessesary)  
    `<output-dir-name>` - name of directory for output

The program creates a folder where you will find `report.html` and separate pictures for each graphic.

Our program was succesfully tested on the next OSs:
- Linux Mint 20.2
- Ubuntu 20.04
- Windows 10 64-bit (via PhyCharm 2021.2.2)

### Output of the program

Our program provides the next metrics:  

- Basic Statistics
- Per base sequence quality
- Per sequence quality scores
- Per base sequence content
- Per base GC content
- Per sequence GC content
- Per base N content
- Sequence Length Distribution
- Sequence Duplication Levels
- Overrepresented sequences

You can read more about all the characteristics (exept Per base GC content) in the [official documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) of the original FastQC project.  
Per base GC content plots out the GC content of each base position in a file. This module gives *warning* if the GC content of any base strays more than 5%  and *failure* if more than 10% from the mean GC content.

### Team members and contributions:
1. [Tatiana Kikalova](https://github.com/Tatiana-kik)
- Backbone of the project (templates for all the files and for html report)
- Basic Statistics
- Per Base Sequence Quality function
- Per Sequence Quality Scores function
2. [Elena Grigoreva](https://github.com/lengrigo)
- Per Base Sequence Content function
- Per Base GC Content function
- Per Sequence GC Content function
- README file
3. [Dmitrii Iliushchenko](https://github.com/DIliushchenko)
- Per Base N Content function
- Sequence Length Distribution function
- requirements.txt file
4. [Pavel Vychik](https://github.com/p-vychik)
- Sequence Duplication Levels
- Overrepresented sequences