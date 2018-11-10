def call(solverConfig) {
    """
        export INITIAL_SURFACE_SNIPPET = ${solverConfig.INITIAL_SURFACE_SNIPPET}
        
        source /etc/profile
        source /usr/share/Modules/init/bash

        module load gcc/7.2.0 cmake/3.11.1
        module unload galois
        module list

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/local/lib

        echo "Local libraries"
        LOCAL_LIB=$HOME/usr/local/lib
        ls -la $LOCAL_LIB

        echo "Local includes"
        LOCAL_INC=$HOME/usr/local/include
        ls -la $LOCAL_INC

        echo "INITIAL_SURFACE_SNIPPET"

        echo "Inserting INITIAL_SURFACE_SNIPPET snippet"
        perl -pe 's/##INITIAL_SURFACE_SNIPPET##/`echo "$ENV{'INITIAL_SURFACE_SNIPPET'}"`/e' -i ./src/problems/cahn_hilliard/ch_2d.hpp
        cat ./src/problems/cahn_hilliard/ch_2d.hpp


        cmake . \
        -DBLAS_LIBRARIES=$LOCAL_LIB/libblas.so \
        -DLAPACK_LIBRARIES=$LOCAL_LIB/liblapack.so \
        -DBOOST_INCLUDEDIR=/opt/boost/include \
        -DGalois_INCLUDE_DIRS=$LOCAL_INC/Galois \
        -DGalois_LIBRARIES=$LOCAL_LIB \
        && make
    """
}