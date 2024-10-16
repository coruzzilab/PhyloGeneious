#!/bin/sh

FORCE=1
ALIAS="OIDrun"
FILEDIR=""
PATHTOA=`pwd`

while getopts 'a:d:s:o:fh' opt; do
  case "$opt" in
    a)
      ALIAS="$OPTARG"
      ;;

    d)
      FILEDIR="$OPTARG"
      ;;

    s)
      SPECIESFILE="$OPTARG"
      ;;
   
    o)
      OUTGROUP="$OPTARG"
      ;;

    f)
      FORCE=0
      ;;

    ?|h)
      echo
      echo "$(basename $0) -d PATH/TO/SEQUENCE_FILES [-a RUN_DIR_ALIAS -s SPECIESFILE -o OUTGROUP -f]"
      echo
      echo "A script for auto-building an OrthologID run directory."
      echo -e "-d\tPath to directory containing sequence files \n-a\tName to give run directory [default: OIDrun]\n-s\tTab-delimited file where first column contains the species \n\tidentifiers used in the sequence files and second column contains the \n\tspecies names\n-o\tList of outgroup species identifiers [ie. \"Species1 Species2 Species3\"]\n-f\tForce overwrite of run directory\n-h\tDisplay help\n"
      echo "Note: Sequence IDs should contain the species identifier (eg. >SpeciesID#GeneID) and sequence files names should only contain the species identifier (eg. SpeciesID.faa)"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

if [[ ${#FILEDIR} -eq 0 ]]; then
    echo "Error: Sequence file directory must be specified."
    exit 1
fi

if [[ -d $ALIAS ]]; then
    if [[ $FORCE -eq 0 ]]; then
        echo "Overwriting run directory $ALIAS setup"
        rm -r $ALIAS/blastdb
    else
        echo -e "Error: directory $ALIAS already exists.\nTo overwrite use the -f argument."
        exit 1
    fi
else
    mkdir $ALIAS
fi


if [[ ! -d $FILEDIR ]]; then
    echo "$FILEDIR not found. Exiting."
    exit 1
else
    NUMFASTA=`ls $FILEDIR | grep ".seq$\|.fa$\|.faa$\|.fasta$\|.pep" | wc -l`
#    NUMFASTA=`ls $FILEDIR/*.seq $FILEDIR/*.fa $FILEDIR/*.faa $FILEDIR/*.fasta $FILEDIR/*.pep | wc -l`
    if [[ $NUMFASTA -eq 0 ]]; then
        echo "No sequence files found in $FILEDIR. Setup incomplete."
        exit 1
    fi
    echo "$NUMFASTA sequence files found."
    mkdir $ALIAS/blastdb
    for FILE in `ls $FILEDIR | grep ".seq$\|.fa$\|.faa$\|.fasta$\|.pep"`; do
        ln -s $FILEDIR/$FILE $ALIAS/blastdb/$FILE
    done
    #check sequence files for invalid characters
    COUNT=`grep -c '[^A-Za-z0-9#._>\*]\|-' $ALIAS/blastdb/* | awk -F":" '{sum+=$2;} END{print sum;}'`
    if [[ $COUNT -gt 0 ]]; then
        echo "Invalid characters found in the following sequence files: please fix before running."
        grep -c '[^A-Za-z0-9#._>\*]\|-' $ALIAS/blastdb/* | awk -F":" '$2>0'
    fi
fi

#ln -s $PATHTOA/$FILEDIR $ALIAS/blastdb
cp /scratch/cgsb/coruzzi/bigplant_v4_share/testdata/config $ALIAS #need to generalize home directory
cp /scratch/cgsb/coruzzi/bigplant_v4_share/testdata/setoid.sh $ALIAS
cp /scratch/cgsb/coruzzi/bigplant_v4_share/testdata/procfiles.txt $ALIAS
sed -i "s|\(OID_USER_DIR=\).*|\1${PATHTOA}/$ALIAS|" $ALIAS/setoid.sh
sed -i "s/\(INGROUP=\).*/\1/" $ALIAS/config
sed -i "s/\(OUTGROUP=\).*/\1/" $ALIAS/config

unspec_species_error () { echo "Species identifiers not specified. Please adjust $ALIAS/config before running OrthologID."; }

if [[ ${#SPECIESFILE} -ne 0 && ${#OUTGROUP} -ne 0 ]]; then
    if [[ ! -f $SPECIESFILE ]]; then
        echo "Error: Species file not found."
        unspec_species_error
        break
    fi
    cp $SPECIESFILE $ALIAS/species.txt
    awk '{print $1}' $SPECIESFILE | grep -q '[^A-Za-z0-9#._]\|-'
    if [[ $? -eq 0 ]]; then
        echo "Error: Illegal characters in species identifiers in $SPECIESFILE."
        echo "Only letter, digits, ., and _ are allowed."
        unspec_species_error
    else
        TEMP=`echo $OUTGROUP | sed "s/ /\\\\\\|/g"`
        INGROUP=`awk '{print $1}' $SPECIESFILE | grep -v "${TEMP}" | tr "\n" " "`
        sed -i "s/\(INGROUP=\).*/\1$INGROUP/" $ALIAS/config
        sed -i "s/\(OUTGROUP=\).*/\1$OUTGROUP/" $ALIAS/config
    fi
elif [[ ${#SPECIESFILE} -eq 0 && ${#OUTGROUP} -eq 0 ]]; then
    unspec_species_error
else
    echo "Error: Species file and outgroup identifiers must be specified together. Setup incomplete."
    unspec_species_error
fi

echo -e "Directory $PATHTOA/$ALIAS setup complete. Please double-check that your sequence files and $ALIAS/config match specifications before running.";

