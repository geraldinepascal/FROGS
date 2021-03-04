import sys
sys.path.append("/Users/moussa/FROGS_moussa/tools/lib")
from frogsUtils import *
from frogsSequenceIO import * #Pour parser les Sequences(fasta et fastq)
from frogsBiom import BiomIO


input_fasta = "test.fasta"
output_fasta = "mon_fichier_renamed.fasta"

FH_input = FastaIO( input_fasta )
# pour écrire dans un fichier au format fasta
FH_output = FastaIO( output_fasta, "w" )

for record in FH_input:
    # récupération du nnom de la séquence 
    print record.id
    # tu peux modifier l'identifiant, par exemple en ajoutant Moussa
    record.id = record.id + "Moussa"
    # récupération de la description
    print record.description
    # récupération de la séquence
    print record.string
    
    # pour écrire la sequence (identifiant, description, et sequence) dans ton fichier de sortie)
    FH_output.write(record)
    

