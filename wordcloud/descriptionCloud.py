__author__		= "Maxence Delannoy, Chayan Kumar Saha, Gemma C. Atkinson"
__copyright__	= "GNU General Public License v3.0"

import argparse
import numpy as np
import pandas as pd
from os import path
from PIL import Image
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator
import matplotlib.pyplot as plt

usage= ''' Description:  Use FlaGs generated description file with suffix "_outdesc.txt" as input and visualize them as wordcloud; Requirement= Python3, WordCloud (https://github.com/amueller/word_cloud) '''


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", help=" FlaGs generated description file with suffix '_outdesc.txt'. ")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.1.3')
args = parser.parse_args()
parser.parse_args()




all_word=""
count ={}
with open (args.input,"r") as Description:
	text = Description.read()
	line = text.split("\n")
	for l in line :
		if l != "":
			W=str.split(l,"\t")
			sep = str.split(W[2]," ")
			word=""
			for m in sep:
				if m!="MULTISPECIES:" and m!="":
					word=word + m + " "
			if word not in count :
				count[word]=1
			else:
				count[word]=count[word]+1
			all_word = all_word + word[:-1] + " "

wordcloud = WordCloud(max_font_size=100, max_words=20,background_color="white").generate_from_frequencies(count)
plt.imshow(wordcloud, interpolation="bilinear")
plt.axis("off")
pngOut=args.input.replace('.txt','_wordcloud.png')
plt.savefig(pngOut, format='png')
