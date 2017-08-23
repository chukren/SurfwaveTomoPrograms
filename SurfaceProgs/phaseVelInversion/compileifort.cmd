#!/bin/bash
ifort -132 -check all -fpe0 -f77rtl -o srchwave591.JdF.ifort srchwave591.JdF.f
ifort -132 -check all -fpe0 -f77rtl -o rdsetupsimul.1.JdF.ifort rdsetupsimul.1.JdF.f -lm /Users/yruan/src/sac-64/lib/libsacio.a
