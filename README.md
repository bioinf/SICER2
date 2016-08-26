# SICER2

Our re-implementation of SICER broad peak calling algorithm


### Install latest version:

~~~
curl https://raw.githubusercontent.com/bioinf/SICER2/master/SICER2 > SICER2
chmod +x SICER2
~~~

### Usage:

~~~
./SICER2 -t [input bam-file] [Optional arguments]
~~~

#### Optional arguments:
~~~
  -t,  --input      Path to `input.bam` file
  -o,  --output     Output file name
  -l,  --log        Output log file name
  -w,  --window     Window size (bp). Default: 200
  -f,  --fragment   Fragment size (bp). Default: 250
  -e,  --gms        Proportion of effective genome length; has to be in [0.0, 1.0]. Default: auto
  -g,  --gap        Gap size shows how many bases could be skipped. Default: 200
  -p,  --pvalue     P-value; has to be in [0.0, 1.0]. Default: 0.2
  -x,  --threshold  Island score threshold. Default: 0
~~~

#### Examples:

~~~
./SICER2 input.bam
./SICER2 -t input.bam -w 200 -g 600 -e 0.77
~~~
