#!/bin/bash
#
# LoadLeveler submission script, MPI/GPU job with Intel MPI
# CRIHAN v 1.00 - March 2013  
# support@crihan.fr

# Job name
# @ job_name = impi.run
# Batch output file
# @ output = $(job_name).o$(jobid)
# Batch error file
# @ error  = $(job_name).e$(jobid)

# Job type
# @ job_type = MPICH

# GPU resources
# @ requirements = (Feature == "gpu")

# Job time (hh:mm:ss)
# @ wall_clock_limit = 2:00:00

# -------------------
# Nodes number 
# @ node = 4 
#
# MPI tasks per node
# @ tasks_per_node = 4
#
# MPI task maximum memory (mb, gb)
# @ data_limit = 8000mb
# -------------------

export BASE=/home/2011007/fcalva02/pwteleman_benchcrihan/;
# Input files directory
# @ cri_initialdir  = /home/2011007/fcalva02/pwteleman_benchcrihan/samples/pormgnospin
# Output files directory
# @ cri_finaldir =  /home/2011007/fcalva02/pwteleman_benchcrihan/samples/pormgnospin

# @ notification = complete
# User e-mail address
# @ notify_user = prenom.nom@adresse.fr 

# @ queue

#cd $LOCAL_WORK_DIR  

# User commands 

# ----------------------------------------
# Select I_MPI_FABRICS (default : "shm:tmi") 

export I_MPI_FABRICS=shm:tmi
# export I_MPI_FABRICS=tmi
# export I_MPI_FABRICS=shm:dapl
# export I_MPI_FABRICS=dapl
# ----------------------------------------

# Source most recent installed CUDA/OpenCL (and CAPS Compilers if needed) environment(s)

# MPI code execution (binary linked with Intel MPI)
mpirun.Ompi $BASE/code/essai.par > logg
#$BASE/code/essai.seq > logg

# Move output files to $LOCAL_SPOOL_DIR
# (before automatic copy to "cri_finaldir")
mv *.res *.dat *.log $LOCAL_SPOOL_DIR
