#!/bin/bash

NUM_UP=$1
ECMS=$2
JobText_SaveDir=/besfs/groups/cal/dedx/$USER/bes/pid_eff_check/run/jobs_text/data/$ECMS
ROOT_SaveDir=/besfs/groups/cal/dedx/$USER/bes/pid_eff_check/run/rootfile/data/$ECMS
mkdir -p $JobText_SaveDir
mkdir -p $ROOT_SaveDir
rm -r $JobText_SaveDir/*txt
rm -r $ROOT_SaveDir/*.root

for ((num = 0; num <= $NUM_UP; num++))
do
    file_list=data_$ECMS\_track_sel_5G-${num}.txt
    rootfile=track_sel_data_$ECMS-${num}.root
    jobOptions=jobOptions_track_sel_data_$ECMS-${num}.txt
    echo "#include \"\$ROOTIOROOT/share/jobOptions_ReadRec.txt\"                                   "  > ${JobText_SaveDir}/${jobOptions}
    echo "#include \"\$VERTEXFITROOT/share/jobOptions_VertexDbSvc.txt\"                            " >> ${JobText_SaveDir}/${jobOptions}
    echo "#include \"\$MAGNETICFIELDROOT/share/MagneticField.txt\"                                 " >> ${JobText_SaveDir}/${jobOptions}
    echo "#include \"\$ABSCORROOT/share/jobOptions_AbsCor.txt\"                                    " >> ${JobText_SaveDir}/${jobOptions}
    echo "#include \"\$MEASUREDECMSSVCROOT/share/anaOptions.txt\"                                  " >> ${JobText_SaveDir}/${jobOptions}
    echo "#include \"/besfs/groups/cal/dedx/$USER/bes/pid_eff_check/run/samples/$ECMS/$file_list\" " >> ${JobText_SaveDir}/${jobOptions}
    echo "                                                                                         " >> ${JobText_SaveDir}/${jobOptions}
    echo "ApplicationMgr.DLLs += {\"TrackSelAlg\"};                                                " >> ${JobText_SaveDir}/${jobOptions}
    echo "ApplicationMgr.TopAlg += { \"TrackSelAlg\" };                                            " >> ${JobText_SaveDir}/${jobOptions}
    echo "                                                                                         " >> ${JobText_SaveDir}/${jobOptions}
    echo "// Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )            " >> ${JobText_SaveDir}/${jobOptions}
    echo "MessageSvc.OutputLevel = 5;                                                              " >> ${JobText_SaveDir}/${jobOptions}
    echo "                                                                                         " >> ${JobText_SaveDir}/${jobOptions}
    echo "// Number of events to be processed (default is 10)                                      " >> ${JobText_SaveDir}/${jobOptions}
    echo "ApplicationMgr.EvtMax = -1;                                                              " >> ${JobText_SaveDir}/${jobOptions}
    echo "                                                                                         " >> ${JobText_SaveDir}/${jobOptions}
    echo "ApplicationMgr.HistogramPersistency = \"ROOT\";                                          " >> ${JobText_SaveDir}/${jobOptions}
    echo "NTupleSvc.Output = {\"FILE1 DATAFILE='$ROOT_SaveDir/$rootfile' OPT='NEW' TYP='ROOT'\"};      " >> ${JobText_SaveDir}/${jobOptions}
done
