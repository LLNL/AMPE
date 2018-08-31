for y in `ls *.py`;
do
sed "s/THE U.S. DEPARTMENT/UT BATTELLE, LLC,\n\/\/ THE U.S. DEPARTMENT/g" $y > temp;
mv temp $y;
sed "s/Security, LLC./Security, LLC and\n\/\/ UT-Battelle, LLC./g" $y > temp;
mv temp $y;
sed "s/Produced at the Lawrence Livermore National Laboratory/Produced at the Lawrence Livermore National Laboratory and\n\/\/ the Oak Ridge National Laboratory/g" $y > temp;
mv temp $y;
done

