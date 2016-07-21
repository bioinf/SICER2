for f in /home/dima/PycharmProjects/BAMplay/*
do
 echo "Processing $f"
 /usr/bin/python2.7 /home/dima/PycharmProjects/BAMplay/broadPeaks.py $f -w 400
 # do something on $f
done
