#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../liosolo/liosolo
fi

$LIOBIN -i carotenox.profile.in -b DZVP -c caroteno.xyz -v 
