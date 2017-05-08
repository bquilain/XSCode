#!/bin/sh -x
uname=`uname`
for base in $@
do
   if [ -f ${CLASSDIR}/${base}.h -a -f ${CLASSDIR}/${base}.cc ] ; then
    awk 'BEGIN {
                print "#ifdef __CINT__";
                        print "#pragma link off all globals;";
                        print "#pragma link off all classes;";
                        print "#pragma link off all functions;";
                        print "";
        }
    {
                if ($1 == "class" && match($5, "^T")) {
                                print "#pragma link C++ class ", $2, ";";
                        }
        }
        END {
                print "#endif";
                }' ${CLASSDIR}/${base}.h > ${base}LinkDef.h
        rootcint -f ${base}Dict.cc -c ${CLASSDIR}/${base}.h ${base}LinkDef.h
        case $uname in
                *BSD )
                c++ `root-config --cflags` -fpic -fPIC -c ${base}Dict.cc -o ${base}Dict.So
                c++ `root-config --cflags` -fpic -fPIC -c ${CLASSDIR}/${base}.cc -o ${base}.So
                c++ -shared -Wl,-x -o $base.so `lorder ${base}Dict.So ${base}.So | tsort -q` `root-config --libs`
                ;;
                Linux )
                c++ ${INCDIRS} ${COPTFLAGS} `root-config --cflags` -fpic -fPIC -c ${base}Dict.cc -o ${base}Dict.o
                c++ ${INCDIRS} ${COPTFLAGS} `root-config --cflags` -fpic -fPIC -c ${CLASSDIR}/${base}.cc -o ${base}.o
                c++ -shared -Wl,-x -o ${base}.so ${base}Dict.o ${base}.o `root-config --libs`
                ;;
        esac
   fi
done
