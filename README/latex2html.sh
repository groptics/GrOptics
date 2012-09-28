# script to convert the users guide to a single page html document 
# requires plastex package (my installation is from Debian)
#   Charlie Duke
#   Grinnell College
#   September 28, 2012
plastex --split-level=-2 --dir=./ --filename=UsersGuideGrOptics.html UsersGuideGrOptics.tex

rm -rf UsersGuideGrOptics
rm *.jhm *.xml *.hs *.paux *.hh* 
rm -rf icons images styles