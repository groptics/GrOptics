# run this script in the GrOptics directory to change the version listing
# in all relevant files

#  Charlie Duke
#  10MAY2012
#  Grinnell College

#  command line parameters are:
#     1. current version line
#     2. new version line
#     3. old date line
#     4. new data line
# e.g. changeVersion.sh VERSION2.1 VERSION2.2 1MARCH2012 10May2012

echo "old version:  $1         new version:   $2"
echo "old date:     $3         new date:      $4"

echo -n "   IF THIS IS CORRECT, HIT ENTER TO CONTINUE OR CTRL/C TO STOP  "
read reply
echo "  starting sed to make the file edits "

sed  -i "s/$1/$2/g" include/*.h
sed  -i "s/$3/$4/g" include/*.h

sed  -i "s/$1/$2/g" src/*.cpp
sed  -i "s/$3/$4/g" src/*.cpp

sed  -i "s/$1/$2/g" Makefile
sed  -i "s/$3/$4/g" Makefile

sed  -i "s/$1/$2/g" Makefile.common
sed  -i "s/$3/$4/g" Makefile.common

sed  -i "s/$1/$2/g" Config/*.pilot
sed  -i "s/$3/$4/g" Config/*.pilot

sed  -i "s/$1/$2/g" Config/*.cfg
sed  -i "s/$3/$4/g" Config/*.cfg

echo "BE SURE TO MANUALLY CHANGE VERSION IN GDEFINITION.H"
