#!/usr/bin/env bash

# Main driver to submit jobs 
# Author SHI Xin <shixin@ihep.ac.cn>
# Modified by JING Maoqiang <jingmq@ihep.ac.cn>
# Created [2019-12-11 Dec 14:56]

usage() {
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-9s  %-40s"  "0.1"       "[run on data sample for pid efficieny]"
    printf "\n\n" 
}

usage_0_1() {
    printf "\n\t%-9s  %-40s"  ""          ""   
    printf "\n\t%-9s  %-40s"  "0.1.1"     "Split data sample with each group 5G"
    printf "\n\t%-9s  %-40s"  "0.1.2"     "Generate Condor jobs on data ---- 1"
    printf "\n\t%-9s  %-40s"  "0.1.3"     "Test for data"
    printf "\n\t%-9s  %-40s"  "0.1.4"     "Submit Condor jobs on data ---- 2"
    printf "\n\t%-9s  %-40s"  "0.1.5"     "Synthesize data root files"
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

usage_0_2() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.2.1"     ""
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

usage_0_3() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.3.1"     ""
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

if [[ $# -eq 0 ]]; then
    usage
    echo "Please enter your option: "
    read option
else
    option=$1
fi

sub_0_1() {

case $option in
    
    # --------------------------------------------------------------------------
    #  Data  
    # --------------------------------------------------------------------------

    0.1.1) echo "Split data sample with each group 5G ..."
	       ./python/get_samples.py /besfs4/offline/data/705-1/jpsi/round12/dst run/samples/data_pid_eff.txt 5G
	       # made 13914 groups
	       ;;

    0.1.2) echo "Generate Condor jobs on data ---- 1..." 
	       cd run/pid_eff/gen_script
	       ./make_jobOption_file_Data_pid_eff.sh 50
	       cd /besfs/groups/cal/dedx/$USER/bes/pid_eff_check/run/jobs_text/data
	       cp -r jobOptions_pid_eff_data-1.txt jobOptions_pid_eff_data-0.txt
           sed -i "s/ApplicationMgr\.EvtMax = -1/ApplicationMgr\.EvtMax = 5000/g" jobOptions_pid_eff_data-0.txt
	       ;;

    0.1.3) echo "Test for data" 
           echo "have you changed test number?(yes / no)
           ./run/jobs_text/data/jobOptions_pid_eff_data-0.txt"
           read opt
           if [ $opt == "yes" ]
               then
               echo "now in yes"  
               cd run/jobs_text/data
               boss.exe jobOptions_pid_eff_data-0.txt
           else
               echo "Default value is 'no', please change test number."
           fi
           ;;

    0.1.4) echo "Submit Condor jobs on data ---- 2..." 
	       cd run/jobs_text/data
           find . -name "*.out.*" | xargs rm
           find . -name "*.err.*" | xargs rm
		   boss.condor -g physics -n 50 jobOptions_pid_eff_data-%{ProcId}.txt
	       ;;

    0.1.5) echo "Synthesize data root files..."
           cd run/rootfile/data
           rm -rf pid_eff_data.root
           hadd pid_eff_data.root *.root
	       ;;

esac

}

sub_0_2() {

case $option in
    
    # --------------------------------------------------------------------------
    #   
    # --------------------------------------------------------------------------

    0.2.1) echo "..."
	       ;;

esac

}

sub_0_3() {

case $option in
    
    # --------------------------------------------------------------------------
    #   
    # --------------------------------------------------------------------------

    0.3.1) echo "..."
	       ;;

esac

}
case $option in
    
    # --------------------------------------------------------------------------
    #  Data  
    # --------------------------------------------------------------------------

    0.1) echo "Running on data sample..."
         usage_0_1 
         echo "Please enter your option: " 
         read option  
         sub_0_1 option 
	     ;;

    0.1.*) echo "Running on data sample..."
           sub_0_1 option  
           ;;  
        
    # --------------------------------------------------------------------------
    #  
    # --------------------------------------------------------------------------

    0.2) echo "..."
         usage_0_2 
         echo "Please enter your option: " 
         read option  
         sub_0_2 option 
	     ;;

    0.2.*) echo "..."
           sub_0_2 option  
           ;;  

    # --------------------------------------------------------------------------
    #  
    # --------------------------------------------------------------------------

    0.3) echo "..."
         usage_0_3 
         echo "Please enter your option: " 
         read option  
         sub_0_3 option 
	     ;;

    0.3.*) echo "..."
           sub_0_3 option  
           ;;  

esac
