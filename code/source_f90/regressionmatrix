#!/bin/bash

NAME=$TELEMANNAME
if [ -z $NAME ]  ; then
	NAME=essai
	echo "no name given - prefix changed to $NAME"
fi


cd $TELEMANROOT/samples
for SAMPLE in *
do
echo "regression test for" $SAMPLE

cd $TELEMANROOT/samples_ref/
for RESULT1 in $SAMPLE* 
do
for RESULT2 in $SAMPLE* 
do

	echo -n  $RESULT1 $RESULT2
 cat $RESULT1 $RESULT1 > reference
 cat $RESULT1 $RESULT2 > try
 rm -f result1$SAMPLE
 rm -f reference.bz2
 rm -f try.bz2
 bzip2 reference try
# stat  -f"%z"  reference.bz2  > refs
# stat  -f"%z"  try.bz2 > trys
stat  --format=%s  reference.bz2  > refs
 stat  --format=%s  try.bz2 > trys

 echo '-' > plus
 echo '(' > leftp
 echo ')' > rightp
 echo '/' > over
 cat leftp refs plus trys rightp over refs | tr -d '\n' > bc
 echo -e '\n' >> bc
echo -n  "mutual change "
 bc -l  < bc 
rm -f bc *.bz2 leftp refs plus trys rightp over refs
done
done
done
