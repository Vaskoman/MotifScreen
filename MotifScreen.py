from Bio import motifs

## Reading from JASPAR "Binding sites information" ####
arnt = motifs.read(open("Arnt.sites"), "sites")
print(arnt.instances[:3])
for instance in arnt.instances:
    print (instance)
print(arnt.counts)

## JASPAR motifs ####
fh = open("VDR.jaspar")
for m in motifs.parse (fh, "jaspar"):
    print(m)
pwm = m.counts.normalize(pseudocounts={"A":0.6, "C":0.4,
"G":0.4, "T":0.6})
print(pwm)
print(pwm.consensus)
print(pwm.anticonsensus)
print(pwm.degenerate_consensus)
print(m.degenerate_consensus)
rpwm = pwm.reverse_complement()
print(rpwm)

## Search for a motif in a sequence ####
from Bio.Seq import Seq
test_seq=Seq("TACACTGCATTACAACCCAAGCATTA", m.alphabet)
print(len(test_seq))

instances = [Seq("TACAA"),
Seq("TACGC"),
Seq("TACAC"),
Seq("TACCC"),
Seq("AACCC"),
Seq("AATGC"),
Seq("AATGC")]
m = motifs.create(instances)
print(m)
for pos, seq in m.instances.search(test_seq):
    print("%i %s" %(pos, seq))

r = m.reverse_complement()
print(r)
for pos, seq in r.instances.search(test_seq):
    print("%i %s" %(pos, seq))
# This program screens for VDREs
# It requests and input sequence from user
# Returns the location of the VDRE within sequence
# Returns the VDRE itself
