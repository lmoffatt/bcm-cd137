#!/bin/bash
#
# Script para correr trabajo serial
#

# Opciones SGE
#$ -V
# Exporta las variables de entorno
#$ -N bcm-r=1.0
# El nombre del job
#$ -d /home/lmoffatt/bcm-cd137/runs/radius=1.0
# Cambia al directorio runs

# Comando para correr el programa, tal cual lo llamaríamos desde la línea de comandos
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lmoffatt.inquimae/opt/gmp-5.1.1/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lmoffatt.inquimae/opt/mpfr-3.1.2/lib:/home/lmoffatt.inquimae/opt/mpc-1.0.1/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lmoffatt.inquimae/opt/gcc-4.8.0/lib64/
export C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH
export OBJC_INCLUDE_PATH=$C_INCLUDE_PATH
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
export PATH=/opt/gcc-4.8.0/bin:$PATH
/home/lmoffatt.inquimae/bcm-cd137/hpc-cluster/BCMdraft optimize 2000 1.0
