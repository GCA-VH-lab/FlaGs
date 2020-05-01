
This readme.txt file stands for the python script descriptionCloud.py which
generates the wordcloud output. The script takes all the Flanking genes
descriptions, count them and generate a plot representing the 20 most common
descriptions. The size of the description text represents the abundance.
This is a easy way to get a better idea of what are the functions/names/origin
of the identified Flanking genes.

### Requirement:

In order to use this script, user needs to install following packages if not
already installed:

- numpy
- pandas
- pillow
- os
- wordcloud
- matplotlib

(use "pip install name_of_the_package")


### Input File:

In the FlaGs tool output directory,	a text (.txt) file with suffix “_outdesc.txt”
is the input for descriptionCloud.py

### Output:

descriptionCloud.py generates the output file with suffix "_wordcloud.png" in
the same directory

### Command:

python3 descriptionCloud.py -i input_outdesc.txt
