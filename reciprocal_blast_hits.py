#!/usr/bin/env python

#purpose: extract reciprocal best BLAST matches for a pair of datasets
#usage: ./reciprocal_blast_hits.py a_vs_b b_vs_a col_query col_match col_score sort_order out_file

#example, requires both blast hits attained highest bit score (12th column in blast's '-outfmt 6'):
# ./reciprocal_blast_hits.py a_vs_b.blastout b_vs_a.blastout 1 2 12 high a_b.hits.out

#example, requires both blast hits attained lowest evalue (11th column in -outfmt 6):
# ./reciprocal_blast_hits.py a_vs_b.blastout b_vs_a.blastout 1 2 11 low a_b.hits.out

import sys

def stop_err( msg ):
  sys.stderr.write("%s\n" % msg)
  sys.exit(1)
def get_col_index(col_str):
  if col_str[0]=="c":
    col_str = col_str[1:]
  return int(col_str)-1

def main():
#Parse Command Line
  try:
    a_vs_b, b_vs_a, c_query, c_match, c_score, sort_order, out_file = sys.argv[1:]
  except:
    stop_err("Expect 7 arguments: two input files, column settings, output file")


  want_highest = want_lowest = False
  if sort_order == "high":
    want_highest = True
  elif sort_order == "low":
    want_lowest = True
  else:
    stop_err("Sort order argument should be high or low")

  if out_file in [a_vs_b, b_vs_a]:
    stop_err("Output file would overwrite an input file")

  c_query = get_col_index(c_query)
  c_match = get_col_index(c_match)
  c_score = get_col_index(c_score)
  if len(set([c_query, c_match, c_score])) < 3:
    stop_err("Need three different column numbers!")

  best_a_vs_b = dict()
  for line in open(a_vs_b):
    if line.startswith("#"): continue
    parts = line.rstrip("\n").split("\t")
    a = parts[c_query]
    b = parts[c_match]
    score = float(parts[c_score])
    if (a not in best_a_vs_b) \
    or (want_highest and score > best_a_vs_b[a][1]) \
    or (want_lowest and score < best_a_vs_b[a][1]):
      best_a_vs_b[a] = (b, score, parts[c_score])
  b_short_list = set(b for (b,score, score_str) in best_a_vs_b.values())

  best_b_vs_a = dict()
  for line in open(b_vs_a):
    if line.startswith("#"): continue
    parts = line.rstrip("\n").split("\t")
    b = parts[c_query]
    a = parts[c_match]
    if a not in best_a_vs_b:
      continue
        #stop_err("The A-vs-B file does not have A-ID %r found in B-vs-A file" % a)
    if b not in b_short_list: continue
    score = float(parts[c_score])
    if (b not in best_b_vs_a) \
    or (want_highest and score > best_b_vs_a[b][1]) \
    or (want_lowest and score < best_b_vs_a[b][1]):
      best_b_vs_a[b] = (a, score, parts[c_score])
#TODO - Preserve order from A vs B?
  a_short_list = sorted(set(a for (a,score,score_str) in best_b_vs_a.values()))

  count = 0
  outfile = open(out_file, 'w')
  outfile.write("#A_id\tB_id\tA_vs_B\tB_vs_A\n")
  for a in a_short_list:
    b = best_a_vs_b[a][0]
    if b in best_b_vs_a and a == best_b_vs_a[b][0]:
      outfile.write("%s\t%s\t%s\t%s\n" % (a, b, best_a_vs_b[a][2], best_b_vs_a[b][2]))
      count += 1
  outfile.close()
#print "Done, %i RBH found" % count
if __name__ == '__main__':
  main()
