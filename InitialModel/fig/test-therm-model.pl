#!/usr/bin/perl
use POSIX;
#
#
#
@age = (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);

foreach $age(@age){
    open(out, "> ageinput.xx");
    print out "$age\n";
    close(out);
    qx{../src/testtherm < ageinput.xx > tmodel.$age.dat };
} 

#=== plot figures ===
$T_RANGE = "0/1429/0/230";
$V_RANGE = "3.8/4.8/0/320";
$Q_RANGE = "0.001/0.1/0/230";

$SCALE   = "X12/-18";
$LOGSCALE = "X12l/-18";

$PST = "temp-age.ps";
$PSV = "vs-age.ps";
$PSQ = "Q-age.ps";

qx{psbasemap -R$T_RANGE -J$SCALE -K -P -Y3 -Ba200:"Temperature (C)":/a40:"Depth (km)":WNes > $PST};
qx{psbasemap -R$V_RANGE -J$SCALE -K -P -Y3 -Ba0.2g0.05:"Shearwave speed (km/s)":/a40g10:"Depth (km)":WNes > $PSV};
#qx{psbasemap -R$Q_RANGE -J$SCALE -K -P -Y3 -Ba100:"Attenuation":/a40:"Depth (km)":WNes > $PSQ};
qx{psbasemap -R$Q_RANGE -J$LOGSCALE -K -P -Y3 -Ba1f3pg3p:"Attenuation (1/Q)":/a40g10:"Depth (km)":WNes > $PSQ};

foreach $age(@age){
    open(out1, "> xx.tmp");
    open(out2, "> yy.tmp");
    open(out3, "> zz.tmp");

    open(in, "< tmodel.$age.dat");
    @line=<in>; close(in);
    foreach $line(@line){
        ($dummy,$depth,$dummy,$temp,$vs,$Q) = split(/\s+/,$line);
        $Qinv=1.0/$Q;
        print out1 "$temp $depth\n";
        print out2 "$vs $depth\n";
        print out3 "$Qinv $depth\n";
    }
    close(out1);   
    close(out2);   
    close(out3);   

    $color = $age*25.5;
    qx{psxy xx.tmp -J$SCALE -R$T_RANGE -W8/$color/0/0 -K -P -O >> $PST};
    qx{psxy yy.tmp -J$SCALE -R$V_RANGE -W8/$color/0/0 -K -P -O >> $PSV};
    qx{psxy zz.tmp -J$LOGSCALE -R$Q_RANGE -W8/$color/0/0 -K -P -O >> $PSQ};
}

open(out4, "> rr.tmp");
open(in, "< ../mantlemodel.dat");
@line=<in>; close(in);

$count = 0;
$depthbot = 20.;
foreach $line(@line){
    ($thick,$dummy,$dummy,$vs) = split(/\s+/,$line);
    $count = $count + 1;
    if ($count > 1){
        $depthtop = $depthbot;
        $depthbot = $depthtop + $thick;
        print out4 "$vs $depthtop\n";
        print out4 "$vs $depthbot\n";
    }
}
close(out4);

qx{psxy rr.tmp -J$SCALE -R$V_RANGE -W8/0 -K -P -O >> $PSV};

# clean up
#qx{rm xx.tmp yy.tmp zz.tmp rr.tmp};

print "$PST\n";
print "$PSV\n";
print "$PSQ\n";
