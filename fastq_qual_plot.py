#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import sys
import logging
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


def read_fq(inputfq):
  q_list = []
  for record in SeqIO.parse(inputfq,'fastq'):
    list_qual = record.letter_annotations["phred_quality"]
    avg = sum(list_qual)/len(list_qual)
    q_list.append(avg)
  return q_list

def freq_stat(freq):
  freq_dict = {}
  for a in freq:
    if a not in freq_dict:
      freq_dict[a] = 1
    else:
      freq_dict[a] +=1
  return freq_dict

def freq_plot(qual_list,name):
  df = freq_stat(qual_list)
  qual_l = list(set(qual_list))
  qual_l.sort()
  x = [ i for i in qual_l]
  y = []
  for t in qual_l:
    ts = round((float(df[t])*100/len(qual_list)),2)
    y.append(ts)
  plt.bar(x, y, width=0.8, linewidth=0.5, color='blue', edgecolor="white")
  plt.xticks(fontsize=8)
  plt.yticks(fontsize=8)
  plt.xlim(20,40)
  plt.xlabel("Reads quality",fontsize=10,weight="bold")
  plt.ylabel("Percentage %",weight="bold", fontsize=10)
  plt.savefig("%s.quality_freq.png" % name)


def set_args(parser):
  parser.add_argument('-i','--input',metavar='FILE', type=str, required=True, help='Input file')
  parser.add_argument('-n','--name', metavar='STR',type=str, default='out', help='Name of png')
  return parser

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''

Description:
	fastq_len_q_plot.py	plot the quality frequency of reads in NGS data
Usage: 
	python fastq_len_q_plot.py fastq out_png_name

''')

    args = set_args(parser).parse_args()
    quality = (read_fq(args.input))
    freq_plot(quality,args.name)

if __name__ == "__main__":
  main()

