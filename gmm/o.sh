#!/bin/bash
#
# Script para correr trabajo serial
#

# Opciones SGE
#PBS -V
# Exporta las variables de entorno
#PBS -N bcm-t
# El nombre del job
#PBS -d /home/lmoffatt/bcm-cd137/runs/dario
# Cambia al directorio runs


# Comando para correr el programa, tal cual lo llamaríamos desde la línea de comandos
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lmoffatt/opt/gmp-5.1.1/lib:/home/lmoffatt/opt/mpfr-3.1.2/lib:/home/lmoffatt/opt/mpc-1.0.1/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lmoffatt/opt/gcc-4.8.0/lib64/
export C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH
export OBJC_INCLUDE_PATH=$C_INCLUDE_PATH
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
export PATH=/opt/gcc-4.8.0/bin:$PATH
/home/lmoffatt/bcm-cd137/gmm/BCMdraft optimize 5 1
