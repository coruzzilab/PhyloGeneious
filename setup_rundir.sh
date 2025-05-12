#!/bin/sh

FORCE=1
ALIAS="PGrun"
FILEDIR=""
PATHTOA=`pwd`
IMAGE=0
image_path=""

while getopts 'a:d:s:o:fh-:' opt; do
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

    -)
      OPTION=$(echo $OPTARG |cut -d"=" -f1)
      case "$OPTION" in
        docker)
          IMAGE=1
          image_path="$(echo $OPTARG |cut -d"=" -f2)"
          if [[ $image_path == "" ]]; then echo "Error: Image name/path not specified"; exit 1; fi
          ;;

        singularity)
          IMAGE=2
          image_path="$(echo $OPTARG |cut -d"=" -f2)"
          if [[ $image_path == "" ]]; then echo "Error: Image name/path not specified"; exit 1; fi
          if [[ "${image_path::1}" != "/" ]]; then
            image_path="$PATHTOA/$image_path"
          fi
          if [[ "${image_path: -1}" == "/" ]]; then
            image_path=${image_path::-1}
          fi
          ;;

        *)
          echo "Error: Unrecognized argument"
          echo
          echo "$(basename $0) -d PATH/TO/SEQUENCE_FILES [-a RUN_DIR_ALIAS -s SPECIESFILE -o OUTGROUP -f --docker=IMAGE_NAME --singularity=PATH/TO/IMAGE]"
          echo
          echo "A script for auto-building a PhyloGeneious run directory."
          echo -e "-d\tPath to directory containing sequence files \n-a\tName to give run directory [default: PGrun]\n-s\tTab-delimited file where first column contains the species \n\tidentifiers used in the sequence files and second column contains the \n\tspecies names\n-o\tList of outgroup species identifiers [ie. \"Species1 Species2 Species3\"]\n-f\tForce overwrite of run directory\n--docker\tUse provided docker image as run enviroment\n--singularity\tUse provided singularity image as run enviroment\n-h\tDisplay help\n"
          echo "Note: Sequence IDs should contain the species identifier (eg. >SpeciesID#GeneID) and sequence files names should only contain the species identifier (eg. SpeciesID.faa)"
          exit 1
          ;;
      esac
      ;;

    ?|h)
      echo
      echo "$(basename $0) -d PATH/TO/SEQUENCE_FILES [-a RUN_DIR_ALIAS -s SPECIESFILE -o OUTGROUP -f --docker=IMAGE_NAME --singularity=PATH/TO/IMAGE]"
      echo
      echo "A script for auto-building a PhyloGeneious run directory."
      echo -e "-d\tPath to directory containing sequence files \n-a\tName to give run directory [default: PGrun]\n-s\tTab-delimited file where first column contains the species \n\tidentifiers used in the sequence files and second column contains the \n\tspecies names\n-o\tList of outgroup species identifiers [ie. \"Species1 Species2 Species3\"]\n-f\tForce overwrite of run directory\n--docker\tUse provided docker image as run enviroment\n--singularity\tUse provided singularity image as run enviroment\n-h\tDisplay help\n"
      echo "Note: Sequence IDs should contain the species identifier (eg. >SpeciesID#GeneID) and sequence files names should only contain the species identifier (eg. SpeciesID.faa)"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

if [[ ${#FILEDIR} -eq 0 ]]; then
    echo "Error: Sequence file directory must be specified."
    exit 1
else
    if [[ "${FILEDIR::1}" != "/" ]]; then
        FILEDIR="$PATHTOA/$FILEDIR"
    fi
    if [[ "${FILEDIR: -1}" == "/" ]]; then
        FILEDIR=${FILEDIR::-1}
    fi
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
    if [[ -f $SPECIESFILE ]]; then
        echo "Only using species in species files."
        for ID in `awk -F"\t" '{print $1}' $SPECIESFILE`; do
            for FILE in `ls $FILEDIR | grep "^${ID}\." | grep ".seq$\|.fa$\|.faa$\|.fasta$\|.pep"`; do
                ln -s $FILEDIR/$FILE $ALIAS/blastdb/$FILE
        done
    else
        for FILE in `ls $FILEDIR | grep ".seq$\|.fa$\|.faa$\|.fasta$\|.pep"`; do
            ln -s $FILEDIR/$FILE $ALIAS/blastdb/$FILE
        done
    fi
    #check sequence files for invalid characters
    COUNT=`grep -c '[^A-Za-z0-9#._>\*]\|-' $ALIAS/blastdb/* | awk -F":" '{sum+=$2;} END{print sum;}'`
    if [[ $COUNT -gt 0 ]]; then
        echo "Invalid characters found in the following sequence files: please fix before running."
        grep -c '[^A-Za-z0-9#._>\*]\|-' $ALIAS/blastdb/* | awk -F":" '$2>0'
    fi
fi

#ln -s $PATHTOA/$FILEDIR $ALIAS/blastdb
cp OID_HOME/testdata/config $ALIAS #need to generalize home directory
cp OID_HOME/testdata/run.sh $ALIAS
cp OID_HOME/testdata/procfiles.txt $ALIAS
sed -i "s|\(OID_USER_DIR=\).*|\1${PATHTOA}/$ALIAS|" $ALIAS/run.sh
sed -i "s/\(INGROUP=\).*/\1/" $ALIAS/config
sed -i "s/\(OUTGROUP=\).*/\1/" $ALIAS/config

unspec_species_error () { echo "Species identifiers not specified. Please adjust $ALIAS/config before running PhyloGeneious."; }

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

image_setup_error () { echo  "Failed to set up environment container for pipeline run."; }

if [[ $IMAGE -gt 0 ]]; then
    echo Using image $image_path
    if [[ $IMAGE -eq 1 ]]; then #docker
        cp OID_HOME/docker.bash $ALIAS
        sed -i "s|RUNDIR|${PATHTOA}/$ALIAS|" $ALIAS/docker.bash
        docker_path=$(which docker)
        if [[ $? -gt 0 ]]; then echo "Error: Docker not found on system."; image_setup_error; exit 1; fi
        sed -i "s|DOCKER|${docker_path}|" $ALIAS/docker.bash
        sed -i "s|phylogeneious:latest|${image_path}|" $ALIAS/docker.bash
        sed -i "s|\(ENV_WRAPPER=\).*|\1${PATHTOA}/$ALIAS/docker.bash|" $ALIAS/run.sh
    fi
    if [[ $IMAGE -eq 2 ]]; then #singularity
        cp OID_HOME/singularity.bash $ALIAS
        sed -i "s|RUNDIR|${PATHTOA}/$ALIAS|" $ALIAS/singularity.bash
        singularity_path=$(which singularity)
        if [[ $? -gt 0 ]]; then echo "Error: Singularity not found on system."; image_setup_error; exit 1; fi
        sed -i "s|SINGULARITY|${singularity_path}|" $ALIAS/singularity.bash
        sed -i "s|IMAGEPATH|${image_path}|" $ALIAS/singularity.bash
        sed -i "s|\(ENV_WRAPPER=\).*|\1${PATHTOA}/$ALIAS/singularity.bash|" $ALIAS/run.sh
    fi
fi

echo -e "Directory $PATHTOA/$ALIAS setup complete. Please double-check that your sequence files and $ALIAS/config match specifications before running.";

