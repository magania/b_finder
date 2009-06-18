#!/bin/bash

EXE_FILE=/home/magania/Bs/finder/bs_finder
META_ELIST=/work/elsanto-clued0/magania/elists/bs/meta_elist
OUT_DIR=/prj_root/2677/ckm_write/magania/bs
LOCK_DIR=/prj_root/5002/magania/lock_bs
RECOVER=bs.root
RECOVER1=bs_elist

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
 scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P14.tar.gz p14.tar.gz
 scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P17.tar.gz p17.tar.gz
 scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P20.tar.gz p20.tar.gz
 scp $MIRROR.fnal.gov:/work/$MIRROR/magania/tar/AA_P21.tar.gz p21.tar.gz
 tar xzf p14.tar.gz 
 tar xzf p17.tar.gz
 tar xzf p20.tar.gz
 tar xzf p21.tar.gz
 scp $MIRROR.fnal.gov:${EXE_FILE}_P14 P14/aa_do 
 scp $MIRROR.fnal.gov:${EXE_FILE}_P17 P17/aa_do
 scp $MIRROR.fnal.gov:${EXE_FILE}_P20 P20/aa_do
 scp $MIRROR.fnal.gov:${EXE_FILE}_P21 P21/aa_do
 scp $MIRROR.fnal.gov:${META_ELIST} meta_elist

 echo "Workspace set."
}

function get_dataset_elist_nlock {
n_meta_elist=`cat meta_elist | wc -l`
n_choose=0

while [ "$n_choose" -lt "$n_meta_elist" ];do
 let N_INTENTOS=N_INTENTOS+1

 if [ $N_INTENTOS -lt 50 ]; then
  n_choose=$[$RANDOM%$n_meta_elist]
 else
  n_choose=$[$N_INTENTOS-50]
 fi

 if [ "$n_choose" -ge "$n_meta_elist" ] ; then 
    break
 fi

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
 scp $MIRROR.fnal.gov:/work/$MIRROR/magania/elists/$ELIST elist
}

function do_task {
 echo "Working ..."
 ./aa_do ed elist > aa.out
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
