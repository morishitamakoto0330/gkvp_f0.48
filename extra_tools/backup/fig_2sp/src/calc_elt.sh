#!/bin/sh

tail -n 67 ../log/*.000000.0.log.001 > wrk.txt

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
    for (i=16;i<=18;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=37;i<=45;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=20;i<=22;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=46;i<=57;i++){print rows[ i] >> "elt_fine.dat"}
                         print rows[27] >> "elt_fine.dat"
    for (i=58;i<=60;i++){print rows[ i] >> "elt_fine.dat"}
                         print rows[29] >> "elt_fine.dat"
    for (i=61;i<=63;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=31;i<=32;i++){print rows[ i] >> "elt_fine.dat"}
    for (i=64;i<=66;i++){print rows[ i] >> "elt_fine.dat"}
                         print rows[12] >> "elt_fine.dat"
  }' wrk.txt

mv *.dat ./data

rm wrk.txt

