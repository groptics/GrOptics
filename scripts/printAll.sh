# print command for your system
# run from include directory

printCommand="a2ps -2 -Pscisec2"

#add linenumbers command for your printer
linenumbers="-C"

filenames=(GArrayTel GDCGeometry GDCRayTracer \
GDCTelescopeFactory GDCTelescope GDefinition GGeometryBase \
GOrderedGrid GPilot GRayTracerBase GReadDCStdBase GReadDCStdGrISU \
ReadPhotonBase GReadPhotonGrISU GReadSCStd GRootDCNavigator \
GRootDCNavigatorLinkDef GRootWriter GRootWriterLinkDef \
GSCTelescopeFactory GSCTelescope SimulateOptics ATelescopeFactory \
GTelescope GUtilityFuncts )

linenumbers="-C"

for i in "${filenames[@]}"; do
    echo "$printCommand $linenumbers ${i}.h"
    echo "$printCommand $linenumbers ../src/${i}.cpp"
done

echo "a2ps -2 -Pscisec2 $linenumbers ../src/grOptics.cpp"

for i in "${filenames[@]}"; do
    $printCommand $linenumbers ${i}.h
    $printCommand $linenumbers ../src/${i}.cpp
done

"a2ps -2 -Pscisec2 $linenumbers ../src/grOptics.cpp"
