def call() {
    '''
        source /etc/profile
        source /usr/share/Modules/init/bash

        module load gcc/7.2.0 cmake/3.11.1
        module unload galois
        module list

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/local/lib

        CDIR=$(pwd)

        # Install galois
        cd $CDIR
        wget http://iss.ices.utexas.edu/projects/galois/downloads/Galois-2.2.1.tar.gz
        rm -rf Galois-2.2.1 || true
        tar xzvf Galois-2.2.1.tar.gz
        cd Galois-2.2.1/build
        mkdir release
        cd release
        cmake -DBoost_INCLUDE_DIR=/opt/boost/include -DSKIP_COMPILE_APPS=ON ../..
        make && \
        make install DESTDIR=~/


        # Add lapack & blas
        wget https://dl.dropboxusercontent.com/s/ay1v02bijizv2dc/lb.zip?dl=0 -O lb.zip
        unzip lb.zip \
        && mv *.so $HOME/usr/local/lib/ \
        && rm lb.zip
    '''
}