#!/bin/bash
# cd /Molly/igor/SICER2/ && ./src/@compile.sh && chmod +x ./SICER2

appname="SICER2"

# Тупо цинично сливаем все питоновы файлы в один
rm -f $appname
cat src/init.py >> $appname
cat src/pybam.py >> $appname
cat src/island_threshold.py >> $appname
cat src/species.py >> $appname
cat src/functions.py >> $appname
cat src/application.py >> $appname
