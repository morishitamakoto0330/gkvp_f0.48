#!/bin/sh

tail -n 58 ../log/*.000000.0.log.001 > wrk.txt

awk '
  {
    rows[NR]=$0
  }
  END{
    for (i= 3;i<=12;i++){print rows[ i] >> "elt_coarse.dat"}

    for (i= 6;i<= 7;i++){print rows[ i] >> "elt_medium.dat"}
    for (i=16;i<=33;i++){print rows[ i] >> "elt_medium.dat"}
                         print rows[12] >> "elt_medium.dat"

    for (i= 6;i<= 7;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=16;i<=22;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=37;i<=48;i++){print rows[ i] >> "elt_fine.dat"}
                         print rows[27] >> "elt_fine.dat"
    for (i=49;i<=51;i++){print rows[ i] >> "elt_fine.dat"}
                         print rows[29] >> "elt_fine.dat"
    for (i=52;i<=54;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=31;i<=32;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=55;i<=57;i++){print rows[ i] >> "elt_fine.dat"}
                         print rows[12] >> "elt_fine.dat"
  }' wrk.txt

mv *.dat ./data

rm wrk.txt

