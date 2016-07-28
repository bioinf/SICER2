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
  -w,  --window     Window size (bp).  DEFAULT: 200
  -e,  --gms        Proportion of effective genome length; has to be in [0.0, 1.0]  DEFAULT: auto
  -g,  --gap        Gap size shows how many bases could be skipped  DEFAULT: 200
  -p,  --pvalue     P-value; has to be in [0.0, 1.0]  DEFAULT: 0.1
  -x,  --threshold  Island score threshold  DEFAULT: 0
~~~

#### Examples:

~~~
./SICER2 input.bam
./SICER2 -t input.bam -c control.bam -w 100 -g 100 -e 0.845
~~~
