echo on
*  WARNING Before using this program, create directories named
*  "work1,work2,work3" for each event

* 1. add event information to SAC headers
* 2. change KZTIME to origin time of the event 
* 3.  correct stations to all look like STS2 response
* 4.  decimate records to ten samples/s and correct to common time window
* 5.  divide stations into geographical groupings with overlap
* 6.  equalize station name length and simplify file names by origin time

setbb EVLA $EVLA
setbb EVLO $EVLO
setbb EVDP $EVDP
* enter year,Julian day, hr, min, sec (integer) and msec (3 digit integer,
* i.e. seconds of 13.03 should be entered as 13 and 030 msecs
* HR, MIN, SEC should be entered as two digits, i.e.  06 05 03 etc. if
* less than 10
setbb YEAR $YEAR JDAY $JDAY HOUR $HOUR MIN $MIN SEC $SEC MSEC $MSEC

echo off

readerr badfile ignore nofiles ignore memory delete
setbb LSTNM 'ERROR' 

cut off

do STA list ISA 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    transfer from polezero subtype ../sts1.pz to polezero subtype ../sts2.pz
    taper
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    echo on
      fft amph
      writesp amph temp
      cp ../corISAtoSTS1 .
      corISAtoSTS1
     readsp amph tempf
     ifft
     w over %fn%
    echo off
  endif
enddo

do STA list CMB  

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    transfer from polezero subtype ../sts1.pz to polezero subtype ../sts2.pz
    taper
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    echo on
      fft amph
      writesp amph temp
      cp ../corCMBtoSTS1 .
      corCMBtoSTS1
     readsp amph tempf
     ifft
     w over %fn%
    echo off
  endif
enddo


do STA list OSI 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    rmean
    taper width 0.02
    transfer from polezero subtype ../sts1.pz to polezero subtype ../sts2.pz
    taper
    if %YEAR EQ 2005
      if %JDAY LE 030
        setbb suffix '0.BHZ.sac.dec'
      endif
      if %JDAY GT 030
        setbb suffix '1.BHZ.sac.dec'
      endif
    endif
    if %YEAR LT 2005
        setbb suffix '0.BHZ.sac.dec'
    endif
    if %YEAR EQ 2006
        setbb suffix '1.BHZ.sac.dec'
    endif
    if %YEAR EQ 2007
      if %JDAY LT 031
        setbb suffix '1.BHZ.sac.dec'
      endif
      if %JDAY GE 031
        setbb suffix '2.BHZ.sac.dec'
      endif
    endif
    if %YEAR GE 2008
        setbb suffix '2.BHZ.sac.dec'
    endif    
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    echo on
    if %suffix eq '1.BHZ.sac.dec'
      fft amph
      writesp amph temp
      cp ../corOSI1toOSI2 .
      corOSI1toOSI2
     readsp amph tempf
     ifft
     w over %fn%
    endif
    echo off
  endif
enddo


do STA list PFO GSC GLA 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    transfer from polezero subtype ../sts1.pz to polezero subtype ../sts2.pz
    taper
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo
  
do STA list SNCC HOPS

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    transfer from polezero subtype ../sts1.pz to polezero subtype ../sts2.pz
    taper
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list Q03C Q04C LAVA R05C R06C Q07A R07C R04C S05C S06C S04C BNLO T05C

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list T06C HAST U04C U05C V03C V04C V05C HELL

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list SCZ2

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for station gain varying with time
    if %YEAR LE 2005
       div 1.28
    endif
    if %YEAR EQ 2006
      if %JDAY LE 045
        div 1.28
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list Q08A R08A S08C EDW2

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    setbb fn2 work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    w %fn2%
  endif
enddo

do STA list SDP MPP SBC VES ARV MPI LGU CIA FMP DEC CVS BDM

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list SMM 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for station gain varying with time
    if %YEAR LE 2006
       div 1.26
    endif
    if %YEAR EQ 2007
      if %JDAY LE 215
        div 1.26
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list RCT

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for station gain varying with time
    if %YEAR LE 2005
       div 1.34
    endif
    if %YEAR EQ 2006
      if %JDAY LE 160
        div 1.34
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list MNRC MCCM POTR WENL JRSC PACP

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo


do STA list TIN CWC LRL BFS BCC GOR

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    setbb fn2 work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    w %fn2%
  endif
enddo

do STA list EDW.

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix 'EDW0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%suffix%
    setbb fn2 work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    w %fn2%
  endif
enddo

do STA list Q09A R09A S09A Q10A R10A S10A U10A Q11A R11A S11A T11A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list U11A V11A R12A S12A T12A U12A V12A W12A Y12C NEE2

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list GRA MPM TUQ RRX HEC BBR DNR BEL

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list PDM

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR GE 2006
       div 1.27
    endif
*  correct for varying station gain with time
    if %YEAR EQ 2005
      if %JDAY GE 116
        div 1.27
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list SHO

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for varying station gain with time
    if %YEAR LE 2004
       div 1.29
    endif
    if %YEAR EQ 2005
      if %JDAY LE 193
        div 1.29
      endif
      if %JDAY GE 194
         div 0.96
      endif
    endif
    if %YEAR GE 2006
       div 0.96
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list GMR

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for varying station gain with time
    if %YEAR EQ 2006
      if %JDAY GE 207
         if %JDAY LE 317
            div 1.27
         endif
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list BC3

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for varying station gain with time
    if %YEAR GE 2008
       div 1.26
    endif
    if %YEAR EQ 2007
      if %JDAY GE 127
        div 1.26
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list ADO

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*  correct for varying station gain with time
    if %YEAR LE 2006
       div 1.27
    endif
    if %YEAR EQ 2007
      if %JDAY LE 207
        div 1.27
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list NEE.

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix 'NEE0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list Q13A R13A S13A T13A U13A V13A W13A X13A Y13A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    setbb fn2 work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    w %fn2%
  endif
enddo

do STA list Q14A R14A S14A T14A U14A V14A W14A X14A Y14A Z14A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list T15A V15A W15A X15A Y15A Z15A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list Q15A R15A S15A U15A U16A S17A T17A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    transfer from polezero subtype ../trillium240.pz to polezero subtype ../sts2.pz
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list Q16A R16A S16A T16A X16A Y16A Z16A 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list R17A U17A X17A Y17A Z17A Y18A Z18A 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list W16A 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR GE 2009
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
    endif    
    if %YEAR EQ 2008
      if %JDAY LE 086
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
      endif
      if %JDAY GE 317
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2007
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list W17A 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR GE 2009
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
    endif    
    if %YEAR EQ 2008
      if %JDAY LE 085
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
      endif
      if %JDAY GE 321
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2007
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list 217A 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR GE 2008
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
    endif    
    if %YEAR EQ 2007
      if %JDAY GE 197
        transfer from polezero subtype ../Trillium240.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list MLAC 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR EQ 2006
      if %JDAY LE 213
        transfer from polezero subtype ../CMG1T.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2005
        transfer from polezero subtype ../CMG1T.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list DAN 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR EQ 2006
      if %JDAY GE 116
        transfer from polezero subtype ../CMG40T.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR GE 2007
        transfer from polezero subtype ../CMG40T.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list TPNV 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR EQ 2005
      if %JDAY LE 292
        transfer from polezero subtype ../USNN.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2004
        transfer from polezero subtype ../USNN.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list FUR 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR EQ 2005
      if %JDAY LE 193
        transfer from polezero subtype ../CMG40T.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2005
        transfer from polezero subtype ../CMG40T.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list MOD 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR EQ 2006
      if %JDAY LE 200
        transfer from polezero subtype ../STS1.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2005
        transfer from polezero subtype ../STS1.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list WUAZ 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    if %YEAR EQ 2004
      if %JDAY LE 263
        transfer from polezero subtype ../USNN.pz to polezero subtype ../sts2.pz
        taper
      endif
    endif
    if %YEAR LE 2003
        transfer from polezero subtype ../USNN.pz to polezero subtype ../sts2.pz
        taper
    endif
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo
echo on
do STA list TUC 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
*    transfer from polezero subtype ../../sts1.pz to polezero subtype ../../sts2.pz
    transfer from polezero subtype ../sts1.pz to polezero subtype ../sts2.pz
    taper
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo
  

do STA list SCI2 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo


do STA list 109C

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    setbb fn2 work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn2 ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    w %fn2%
  endif
enddo

do STA list DVT SWS

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list MONP 112A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo


do STA list 113A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    setbb fn2 work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn2 ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    w %fn2%
  endif
enddo

do STA list MONP2

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix 'MON2.BHZ.sac.dec'
    setbb fn2 work2/%YEAR%.%JDAY%.%HOUR%.%MIN%.%suffix%
*    setbb fn2 ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%suffix%
    rmean
    taper width 0.02
    w %fn2%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn2%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list 114A 115A 116A 117A 118A 216A 218A 318A

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work3/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*s    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list NE70 NE71 NE72 NE73 NE80

 r *$STA$*BHZ*
 
 if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '.BHZ.sac.dec'
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
*    setbb fn ../work4/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    rmean
    taper width 0.02
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
  endif
enddo

do STA list KCC 

  r *$STA$*BHZ*

  if &1,KSTNM ne %LSTNM
    setbb STN $STA$
    rmean
    ch evla %EVLA evlo %EVLO evdp %EVDP
    ch O gmt %YEAR %JDAY %HOUR %MIN %SEC %MSEC
    setbb otime &1,O
    evaluate to cotime %otime * (-1.0)
    ch ALLT %cotime IZTYPE IO
    setbb suffix '0.BHZ.sac.dec'
    rmean
    taper width 0.02
    if %YEAR EQ 2006
      if %JDAY LE 178
        transfer from polezero subtype ../STS1.pz to polezero subtype ../sts2.pz
        setbb suffix '1.BHZ.sac.dec'
        taper
      endif
    endif
    if %YEAR LE 2005
        transfer from polezero subtype ../STS1.pz to polezero subtype ../sts2.pz
        setbb suffix '1.BHZ.sac.dec'
        taper
    endif
    if %YEAR EQ 2007
      if %JDAY GE 307
        transfer from polezero subtype ../STS1.pz to polezero subtype ../sts2.pz
        setbb suffix '2.BHZ.sac.dec'
        taper
      endif
    endif
    if %YEAR GE 2008
        transfer from polezero subtype ../STS1.pz to polezero subtype ../sts2.pz
        setbb suffix '2.BHZ.sac.dec'
        taper
    endif    
    setbb fn work1/%YEAR%.%JDAY%.%HOUR%.%MIN%.%STN%%suffix%
    w %fn%
    cuterr fillz
    cut O 0 O 5400
*  this will make all records start within one sample time of origin
    r %fn%
    ch LOVROK true
    ch LCALDA true
    cut off
    cuterr u
*  presuming all records at 40 samples/s, decimate by 10
*  if records at 20 samples/s, decimate by 5
    if &1,DELTA eq 5.000000e-02
     dec 2
     w over
    endif
    if &1,DELTA eq 2.500000e-02
     dec 4
     w over
    endif
    if %suffix eq '1.BHZ.sac.dec'
      fft amph
      writesp amph temp
      cp ../corKCC1toKCC0 .
      corKCC1toKCC0
      readsp amph tempf
      ifft
      w over %fn%
    endif
  endif
  endif
enddo

readerr badfile warning nofiles fatal memory delete
cuterr usebe
cut off
