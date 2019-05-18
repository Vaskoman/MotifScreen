
## Read sequence file to be examined for motifs ####
# Define holder variable for sequence to be read from file:
seq = ""

# Open file containing sequence:
sequence = open ("testSeq.txt", "r") # read the file

# Read sequence line-by-line and append if not FASTA title:
for line in sequence:
    li = line.strip() # remove trailing spaces and new lines
    if not li.startswith('>'):
        seq = seq + line.strip() # removes all trailing spaces and \n's
sequence.close()
seq = seq.upper()


## Create PSSM from a JASPAR file for motif ####
from Bio import motifs
# Open file and read:
with open("VDR.jaspar") as handle:
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
print(pssm)

## Search for instances usine PSSMs ####

from Bio.Seq import Seq
test_seq = Seq(seq, m.alphabet)

for position, score in pssm.search(test_seq, threshold=1.0):
    print("Position %d: score = %5.3f" %(position,score))



'''The negative positions refer to instances of the motif found on
the reverse strand of the test sequence, and follow the Python
convention on negative indices. Therefore, the instance of the motif
at pos is located at test_seq[pos:pos+len(m)] both for positive and
for negative values of pos.'''
'''You may notice the threshold parameter, here set arbitrarily to 3.0.
This is in log2, so we are now looking only for words, which are eight
times more likely to occur under the motif model than in the background.
The default threshold is 0.0, which selects everything that looks more
like the motif than the background.'''

# Calculate the scores at all positions along the sequence:
'''print(pssm.calculate(test_seq)) # FORWARD STRAND only
# Same for reverse strand:
rpssm = pssm.reverse_complement()
print(rpssm.calculate(test_seq)) # for reverse strand now
## JASPAR motifs ####

print(pwm.anticonsensus)
print(pwm.degenerate_consensus)
print(m.degenerate_consensus)
rpwm = pwm.reverse_complement()
print(rpwm)



print("%4.2f" % pssm.max)
print("%4.2f" % pssm.min)

mean = pssm.mean(background)
std = pssm.std(background)
print("mean = %0.2f, standard deviation = %0.2f" %(mean, std))'''
''' The mean is equal to relative entropy and is a measure for the
information content of the motif compared to background '''

## Search for a motif in a sequence ####
'''from Bio.Seq import Seq
test_seq=Seq("TACACTGCATTACAACCCAAGCATTA", m.alphabet)
print(len(test_seq))
# Searching for exact matches:
instances = [Seq("TACAA"),
Seq("TACGC"),
Seq("TACAC"),
Seq("TACCC"),
Seq("AACCC"),
Seq("AATGC"),
Seq("AATGC")]
m = motifs.create(instances)
for pos, seq in m.instances.search(test_seq):
    print("%i %s" %(pos, seq))

r = m.reverse_complement()
print(r)
for pos, seq in r.instances.search(test_seq):
    print("%i %s" %(pos, seq))

# Searching for matches using PSSM score:

'''
# This program screens for VDREs
# It requests and input sequence from user
# Returns the location of the VDRE within sequence
# Returns the VDRE itself
