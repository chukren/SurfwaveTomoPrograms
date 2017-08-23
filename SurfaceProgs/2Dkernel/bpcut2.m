* change cut1 and cut2 for different periods for example, from 20 to 50 we set
* cut1=700 and cut2=900

echo on

 setbb cut1 500
 setbb cut2 900
 setbb ctb1 0
 setbb ctb2 1600
 evaluate to cut1 %cut1 - 50
 evaluate to cut2 %cut2 + 50
    cut %cut1 %cut2
    evaluate to totlen %cut2 - %cut1
    evaluate to tprfrac 50 / %totlen

     r $fn
     p1
     taper w %tprfrac
     w temp2
    cuterr fillz
    cut %ctb1 %ctb2
    r temp2
    cuterr u
     w $fn1
     cut off
     r $fn1
     p1
echo off

sc rm  temp?



