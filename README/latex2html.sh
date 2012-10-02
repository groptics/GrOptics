# script to convert the users guide to a single page html document 
# requires plastex package (my installation is from Debian)
#   Charlie Duke
#   Grinnell College
#   September 28, 2012

# this script produces both a single html file:
#         GrOptics/README/UsersGuideGrOptics.html
# and linked chapter html files:
#         GrOptics/README/html/index.html

plastex  --dir=./ --renderer=XHTML --dir=./html UsersGuideGrOptics.tex
rm html/*.jhm html/*.xml html/*.hs html/*.paux html/*.hh* 
rm -rf html/icons html/images html/styles

plastex --split-level=-4 --dir=./ --filename=UsersGuideGrOptics.html UsersGuideGrOptics.tex

rm *.jhm *.xml *.hs *.paux *.hh* 
rm -rf icons images styles

# lastly, produce a pdf file from UsersGuideGrOptics.tex

latex UsersGuideGrOptics
dvipdf UsersGuideGrOptics.dvi
rm UsersGuideGrOptics.aux UsersGuideGrOptics.dvi UsersGuideGrOptics.log