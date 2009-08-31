#!/bin/bash

EXE_FILE=/home/magania/.yp/finder/pi_finder
META_ELIST=/work/tlaloc-clued0/magania/elists/aadst/meta_elist_last500
OUT_DIR=/prj_root/2677/ckm_write/magania/pi
LOCK_DIR=/prj_root/2677/ckm_write/magania/pi_lock
RECOVER=pi.root
RECOVER1=pi_elist

ROOT_DIR=/prj_root/2677/ckm_write/magania/local/root


# if metaelist no existe salir

function get_mirror {
 MIRROR0=elsanto-clued0
 MIRROR1=tlaloc-clued0
 MIRROR2=cinvescuatro-clued0

 CMIRROR=MIRROR$[$RANDOM%3]

 MIRROR=${!CMIRROR}

 if [ -n "${PBS_JOBID:+x}" ]; then
  kbatch
 else
  klist -s || kinit 
 fi

 echo "Mirror: " $MIRROR
}

function get_myid_workdir {
if [ -n "${PBS_JOBID:+x}" ]; then
 BASE_DIR=/scratch
 if echo $PBS_JOBID | grep -q timber ; then BASE_DIR=/batch; fi
 WORKDIR=$BASE_DIR/$PBS_JOBID
 MYID=$PBS_JOBID
else
 ALEATORIO=$RANDOM
 WORKDIR=/tmp/$ALEATORIO
 MYID=$ALEATORIO.`hostname`
 mkdir $WORKDIR
 echo Doing local $WORKDIR
fi
}

function set_workspace {
 . $ROOT_DIR/bin/thisroot.sh
 get_mirror
# scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P14.tar.gz p14.tar.gz
# scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P17.tar.gz p17.tar.gz
# scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P20.tar.gz p20.tar.gz
 scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P21.tar.gz p21.tar.gz
# tar xzf p14.tar.gz 
# tar xzf p17.tar.gz
# tar xzf p20.tar.gz
 tar xzf p21.tar.gz
# scp $MIRROR.fnal.gov:${EXE_FILE}_P14 P14/aa_do 
# scp $MIRROR.fnal.gov:${EXE_FILE}_P17 P17/aa_do
# scp $MIRROR.fnal.gov:${EXE_FILE}_P20 P20/aa_do
 scp $MIRROR.fnal.gov:${EXE_FILE} P21/aa_do
 scp $MIRROR.fnal.gov:${META_ELIST} meta_elist

 echo "Workspace set."
}

function get_dataset_elist_nlock {
n_meta_elist=`cat meta_elist | wc -l`
n_choose=0

while [ true ];do
 #echo "ls -1 $OUT_DIR | grep $RECOVER | cut -d. -f2 | sort -n -u | awk '{while (n!=$1){print n++};n++ } END{while(n<'$n_meta_elist'){print n++}}'"
 MISSING=`ls -1 $OUT_DIR | grep $RECOVER | cut -d. -f3 | sort -n -u | awk '{while (n!=$1){print n++};n++ } END{while(n<'$n_meta_elist'){print n++}}'`
 # echo $MISSING
 n_missing=`echo $MISSING | wc -w`
 echo Total $n_meta_elist Missing $n_missing

 n_locks=`ls -1 $LOCK_DIR | wc -l`
 if [ "$n_missing" -lt "$n_locks" ] ; then  
    break
 fi

 x_choose=$[$RANDOM%$n_missing+1]
 echo $x_choose
 n_choose=`echo $MISSING | cut -d' ' -f$x_choose`
# echo $n_choose

 e_choose=`awk '$1=='$n_choose meta_elist`
 DATASET=`echo $e_choose | cut -f2 -d' '`
  ELIST=`echo $e_choose | cut -f3 -d' '`

  if ls $OUT_DIR/$RECOVER.$n_choose.* 2> /dev/null ; then
    continue
  fi

  if [ -f $LOCK_DIR/lock.$n_choose ]; then
    continue
  fi

  sleep $[$RANDOM%10]

  if [ -f $LOCK_DIR/lock.$n_choose ]; then
    continue
  fi

  touch $LOCK_DIR/lock.$n_choose

  NLOCK=$n_choose
  echo "Choose $NLOCK $DATASET $ELIST"
  return 0
done
NLOCK=X
return 1
}

function copy_elist {
 get_mirror
 echo Copy ...
 cp $ELIST aadst
 echo aadst > file
}

function do_task {
 echo "Working ..."
 ./aa_do -i file > aa.out
}

function send_result {
 mv $RECOVER $OUT_DIR/$RECOVER.$n_choose.$MYID
 mv $RECOVER1 $OUT_DIR/$RECOVER1.$n_choose.$MYID
 mv aa.out $OUT_DIR/out.$n_choose.$MYID
 rm $LOCK_DIR/lock.$NLOCK
}


#    ------------- MAIN -----------

get_myid_workdir
cd $WORKDIR
set_workspace
while get_dataset_elist_nlock ; do
  echo "dataset: $DATASET"
  cd $DATASET
  copy_elist
  do_task
  send_result
  cd $WORKDIR
done
