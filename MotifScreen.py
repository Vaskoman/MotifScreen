
## Read sequence file to be examined for motifs ####
file_path = input ("Enter absolute path to the file containing the sequence to be tested:\n")
# Define holder variable for sequence to be read from file:
seq = ""
# Open file containing sequence:
sequence = open (file_path, "r") # read the file
# Read sequence line-by-line and append if not FASTA title:
for line in sequence:
    li = line.strip() # remove trailing spaces and new lines
    if not li.startswith('>'):
        seq = seq + line.strip() # removes all trailing spaces and \n's
sequence.close()
seq = seq.upper()

## Create PSSM from a JASPAR file for motif ####
motif = input ("Enter absolute path to the file containing the jaspar TF motif:\n")
from Bio import motifs
# Open file and read:
with open(motif) as handle:
    m = motifs.read(handle, "jaspar")
# Get relative proportion of nucleotides from test sequence:
A = seq.count("A")/len(seq) # proportion of A
C = seq.count("C")/len(seq) # proportion of C
G = seq.count("G")/len(seq) # proportion of G
T = seq.count("T")/len(seq) # proportion of T
pcounts = {"A":A, "C":C, "G":G, "T":T}
# Generate Position Weight Matrix from jaspar:
pwm = m.counts.normalize(pseudocounts=
{"A":0.6,"C":0.4,"G":0.4,"T":0.6})
# Generate Position-Specific Scoring Matrix for motif Searching
pssm = pwm.log_odds(pcounts)

## Search for instances usine PSSMs ####
from Bio.Seq import Seq
import operator
test_seq = Seq(seq, m.alphabet)
# Forward:
print("MOTIF ALIGNMENT:\nNegative positions denote reverse strand:\n")
for position, score in pssm.search(test_seq, threshold=1.0):
    print("Position %d: score = %5.3f \n Sequence = %s" %(position,score,
    test_seq[position:position+len(pssm.consensus)]))
results = [(position, score, test_seq[position:position+len(pssm.consensus)])
for position,score in pssm.search(test_seq, threshold=1.0)]
results.sort(key = operator.itemgetter(1),reverse=True) # sort the list based on Score
# Reverse:
rpssm = pssm.reverse_complement() # generate reverse complement of motif
print("\n\nREVERSE COMPLEMENT:\n")
for position, score in rpssm.search(test_seq, threshold=1.0):
    print("Position %d: score = %5.3f \n Sequence = %s" %(position,score,
    test_seq[position:position+len(rpssm.consensus)]))
resultsRC = [(position, score, test_seq[position:position+len(rpssm.consensus)])
for position, score in rpssm.search(test_seq, threshold=1.0)]
resultsRC.sort(key= operator.itemgetter(1),reverse=True) # sort the list based on Score

## Create tables and save as csv files ####
csvtable = input ("Enter absolute path and filename of table followed by '.csv':\n")
import csv # import module
# Forward
with open(csvtable, "w") as out:
    csv_out = csv.writer(out)
    csv_out.writerow(['StartPosition','Score','Sequence'])
    #for x in results:
    #    csv_out.writerow(x)
    csv_out.writerows(results)
# Reverse
'''
with open('resultsRC.csv', "w") as out:
    csv_out = csv.writer(out)
    csv_out.writerow(['StartPosition','Score','Sequence'])
    #for x in results:
    #    csv_out.writerow(x)
    csv_out.writerows(resultsRC)
'''

'''
The negative positions refer to instances of the motif found on
the reverse strand of the test sequence, and follow the Python
convention on negative indices. Therefore, the instance of the motif
at pos is located at test_seq[pos:pos+len(m)] both for positive and
for negative values of pos.You may notice the threshold parameter,
here set arbitrarily to 3.0. This is in log2, so we are now looking only
for words, which are eight times more likely to occur under the motif model
than in the background. The default threshold is 0.0, which selects everything
that looks more like the motif than the background.
'''
