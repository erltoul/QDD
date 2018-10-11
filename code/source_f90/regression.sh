#!/bin/bash

NAME=$TELEMANNAME
if [ -z $NAME ]  ; then
	NAME=essai
	echo "no name given - prefix changed to $NAME"
fi


if [ -z $1 ] ; then
    SAMPLE=C2
else
    SAMPLE=$1
fi

echo "regression test for" $SAMPLE

cd $TELEMANROOT/samples/$SAMPLE
rm -f psp* 
rm -f  2st-* pescel* pMP* pminpos* pposion* pgemion* pvelion* pkinenion* pr*
rm -f ptempion* penerclu* pforce* plforce* projforce* pfrcel* pforce* 
rm -f  pstatmom* infosp* penerstat* 

 $TELEMANROOT/code/telemanexec > result
	cat Time >> ../../samples_ref/Time_$SAMPLE
	if [[ $(find result -type f -size +200c 2>/dev/null) ]] 
	then
# 		sed '/time/d' result >  result_t
		rm -f  result
		rm -f for006*
		rm -f  pstat*
		rm -f *dat
		rm -f gmon.out
		tail -2 Time
		rm -f  Time
		rm -f save*  rsave*

		cd $TELEMANROOT/samples
		tar -cf result2$SAMPLE $SAMPLE
	else
		echo "result file to small - probable error"
		exit 114
	 fi
