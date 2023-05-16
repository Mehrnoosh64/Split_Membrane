#!/bin/bash

MAINDIR=/YOUR_PATH

SUBFOLDER=/Split_Slipids_2020.ff

cd $MAINDIR/$SUBFOLDER

echo -n "Write the value for Radius: "
read RADIUS

echo -n "Write the value for Potential: "
read POTENTIAL

for ITPs in PA.itp PC.itp PE.itp PS.itp SM.itp; 
do

    awk -v OFS=" "  -v RADIUS="$RADIUS" -v POTENTIAL="$POTENTIAL" '{ if ( $2 == "RADIUS" )
        {$3=RADIUS}
        if ( $2 == "POTENTIAL" )
        {$3=POTENTIAL}
        print $0
        }' $ITPs > file.tmp && mv file.tmp $ITPs

done

echo -n "Write the value for rRadius: "
read rRADIUS

echo -n "Write the value for rPotential: "
read rPOTENTIAL

for ITPs in PA.itp PC.itp PE.itp PS.itp SM.itp; 
do

    awk -v OFS=" "  -v rRADIUS="$rRADIUS" -v rPOTENTIAL="$rPOTENTIAL" '{ if ( $2 == "rRADIUS" )
        {$3=rRADIUS}
        if ( $2 == "rPOTENTIAL" )
        {$3=rPOTENTIAL}
        print $0
        }' $ITPs > file.tmp && mv file.tmp $ITPs

done

