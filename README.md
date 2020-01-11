# plane-averaged-charge-density-difference
plane-averaged charge density difference along z-axis
copyright: Dr. Yaqiong Su, Xiamen
contact me: Y-Q.Su@tue.nl; yqsu1989@gmail.com; 398168203@qq.com
How to perform:
species A's CHGCAR is named as CHGCAR0, species B's CHGCAR is named as CHGCAR1, and species AB's CHGCAR is named as CHGCAR2. And then put CHGCAR0, CHGCAR1 and CHGCAR2 under the same directory, and then only direct run this script, then you will get a file named 'zchgdiff.dat', and then you can use this file to make the figure
CHGDIFF = CHGCAR2 - CHGCAR1 - CHGCAR0
