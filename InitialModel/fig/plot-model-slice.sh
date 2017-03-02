#!/bin/bash
model=$1
range=-130/-120/200/350
#scale=16/-9
scale=16/-6
awk '{print $1,$3}' testmoho.dat     | psxy -R$range -JX$scale -W8/255/0/0  -P -K -Y5 > $model.ps 
awk '{print $1,$3}' testslabtopo.dat | psxy -R$range -JX$scale -W8/25/255/0 -P -K -O >> $model.ps 
awk '{print $1,$3}' testbasement.dat | psxy -R$range -JX$sacle -W8/25/5/250 -P -K -O >> $model.ps 
psxy $model -R$range -JX$scale -Sp0.08 -Ba2:"Lon (degree)":/a10f5:"Depth (km)":WSne -P -K -O >> $model.ps
echo "$model.ps"
