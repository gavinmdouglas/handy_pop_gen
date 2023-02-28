FASTA=$1
WORKING_DIR=$2
HYPHY_DIR="/home/gdouglas/local/prg/hyphy-analyses"

STARTING_DIR=$PWD

FASTA_FILE=$(basename $FASTA)

BASENAME=$(basename $FASTA_FILE .fna)
BASENAME=$(basename $BASENAME .fa)
BASENAME=$(basename $BASENAME .fasta)
BASENAME=$(basename $BASENAME .fas)

mkdir $WORKING_DIR
cp $FASTA $WORKING_DIR/$FASTA_FILE
cd $WORKING_DIR 

hyphy $HYPHY_DIR/codon-msa/pre-msa.bf --input ./$FASTA_FILE &> pre-msa_log.txt

PROTEIN_INTERMEDIATE=$FASTA_FILE"_protein.fas"
NUCL_INTERMEDIATE=$FASTA_FILE"_nuc.fas"
PROTEIN_MSA=$FASTA_FILE"_protein.msa"
FINAL_MSA="$BASENAME.msa.fna"

muscle -in $PROTEIN_INTERMEDIATE -out $PROTEIN_MSA  2> muscle_log.txt

hyphy $HYPHY_DIR/codon-msa/post-msa.bf --protein-msa $PROTEIN_MSA --nucleotide-sequences $NUCL_INTERMEDIATE --output $FINAL_MSA &> post-msa_log.txt

cd $STARTING_DIR

