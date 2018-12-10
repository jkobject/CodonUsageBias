 for i in `ls -d nonlinfit*.txt`
do
sed "/.[A-Z] ,[0-9]* , [A-Z]./d" $i > tmp
mv tmp $i
#echo $i
done
