set style data dots
set nokey
set xrange [0: 4.33512]
set yrange [ 10.64903 : 19.17198]
set arrow from  0.54940,  10.64903 to  0.54940,  19.17198 nohead
set arrow from  1.40476,  10.64903 to  1.40476,  19.17198 nohead
set arrow from  1.69101,  10.64903 to  1.69101,  19.17198 nohead
set arrow from  1.97730,  10.64903 to  1.97730,  19.17198 nohead
set arrow from  3.00657,  10.64903 to  3.00657,  19.17198 nohead
set xtics ("G"  0.00000,"X"  0.54940,"M"  1.40476,"Y"  1.69101,"K"  1.97730,"F"  3.00657,"G"  4.33512)
 plot "NbP_band.dat"
