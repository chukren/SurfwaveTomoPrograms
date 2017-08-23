#!/bin/bash
g77 -o srchwave591.JdF srchwave591.JdF.f
g77 -o rdsetupsimul.1.JdF rdsetupsimul.1.JdF.f -lm /usr/local/sac/lib/libsacio.a
