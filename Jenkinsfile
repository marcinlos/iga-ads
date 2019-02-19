pipeline {
    agent { 
        node { 
            label 'node2'
        }
    }
    environment {
        GIT_URL = "${scm.userRemoteConfigs[0].url}"
        BUILD_NUMBER = "${env.BUILD_NUMBER}"
    }

    parameters {
        string(
            name: 'DIMENSION',
            defaultValue: "2",
            description: 'The dimension count of the problem. Choose from {2,3}.'
        )
        string(
            name: 'ELEMENTS',
            defaultValue: "25",
            description: 'The numer of elements along each dimension (25 means 25x25 elements in mesh)'
        )
        string(
            name: 'SPLINEDEG',
            defaultValue: "4",
            description: 'The order of B-Splines used.'
        )
        string(
            name: 'ORDER',
            defaultValue: "2",
            description: 'The order.'
        )
        string(
            name: 'STEPS',
            defaultValue: "10000",
            description: 'The number of steps to simulate'
        )
        string(
            name: 'DELTA',
            defaultValue: "0.001",
            description: 'Time step length.'
        )
        string(
            name: 'SCHEME',
            description: 'The scheme to be used. It has to comply with the format \"s | a(s-1) ... a(0) | b(s) b(s-1) ... b(0)\" or the pre-built values \"AM-0", ..., "AM-3", "BDF-1", ..., "BDF-4\""',
            defaultValue:  "1 | -1 | 0.5 0.5"
        )
        text(
            name: 'INITIAL_FORMULA_SNIPPET',
            description: 'CPP code snippet which should return a value for (x,y,i) for 2D or (x,y,z,i) for 3D problem. Parameter "i" is the step number (not that this is a not time). The size of the time step can be accessed like in the example beneath.',
            defaultValue: """
                constexpr double k = 2 * M_PI * M_PI;
                double t = i * steps.dt;
                double e = std::exp(-k * t); 
                return e * std::sin(x * M_PI) * std::sin(y * M_PI);
            """
        )
        string(
            name: 'EMAIL_RECIPIENTS',
            description: 'Comma-separated recipients of the email notifications (will be sent after success or failure)',
            defaultValue: ""
        )
    }

    stages {

        stage('Checkout') {
            steps {
                echo "Checking out ${GIT_URL} on branch ${GIT_BRANCH}"
                git branch: "${GIT_BRANCH}", url: "${GIT_URL}"
            }
        }
        
        stage('Install dependencies') {
            steps {
                sh '''#!/bin/bash
                    source /etc/profile
                    source /usr/share/Modules/init/bash
                    module load gcc/7.2.0 cmake/3.11.1
                    module unload galois
                    module list
                    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/local/lib
                    export CDIR=$(pwd)
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
                    wget https://dl.dropboxusercontent.com/s/ay1v02bijizv2dc/lb.zip?dl=0 -O lb.zip
                    unzip lb.zip \
                        && mv *.so $HOME/usr/local/lib/ \
                        && rm lb.zip
                '''
            }
        }

        stage('Initialise') {
            steps {
                sh '''#!/bin/bash
cat > initformula<<EOL
${INITIAL_FORMULA_SNIPPET}
EOL
sed -i '/##INITSTART##/,/##INITEND##/!b;//!d;/##INITSTART##/r initformula' src/problems/multistep/multistep*d.hpp
                '''
            }
        }

        stage('Compile') {
            steps {
                sh '''#!/bin/bash
                    source /etc/profile
                    source /usr/share/Modules/init/bash
                    module load gcc/7.2.0 cmake/3.11.1
                    module unload galois
                    module list
                    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/local/lib
                    echo "Local libraries"
                    LOCAL_LIB=$HOME/usr/local/lib
                    LOCAL_INC=$HOME/usr/local/include
                    
                    cmake . \
                    -DUSE_GALOIS=ON \
                    -DBLAS_LIBRARIES=$LOCAL_LIB/libblas.so \
                    -DLAPACK_LIBRARIES=$LOCAL_LIB/liblapack.so \
                    -DBOOST_INCLUDEDIR=/opt/boost/include \
                    -DGalois_INCLUDE_DIRS=$LOCAL_INC/Galois \
                    -DGalois_LIBRARIES=$LOCAL_LIB \
                    && make
                '''
            }
        }

        stage('Run') {
            steps {
                sh '''#!/bin/bash
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
                    
                    ./multistep \
                    ${DIMENSION} \
                    ${SPLINEDEG} \
                    ${ELEMENTS} \
                    "${SCHEME}" \
                    ${ORDER} \
                    ${STEPS} \
                    ${DELTA} \
                    > errors.out
                '''
                stash name: 'results', includes: '*.data,*.out'
            }
        }

        stage('Process results') {
            steps {
                sh '''#!/bin/bash
                    echo 'Producing charts'
                    gnuplot plot_norm

                    echo 'Creating movie'
                    gnuplot plot
                    ./movie

                    echo "Compressing results\n"
                    zip data.zip *.data
                    zip images.zip *.png
                    zip movies.zip *.mp4

                    RESULTS_DIR=/home/proj/jenkins_pub/pub/multistep-$BUILD_NUMBER/

                    mkdir $RESULTS_DIR
                    cp errors.png $RESULTS_DIR/
                    cp images.zip $RESULTS_DIR/
                    cp movies.zip $RESULTS_DIR/
                    cp data.zip $RESULTS_DIR/
                    chmod -R 777 $RESULTS_DIR
                '''
            }
        }

    }

    post {
        always {
            cleanWs()
        }
        success {  
             mail(
                 subject: "Simulations no. ${BUILD_NUMBER} are complete",
                 to: "${EMAIL_RECIPIENTS}",
                 from: 'jenkins@a2s.agh.edu.pl',
                 body: """
                    Download the results from:
                    <ul>
                        <li>
                        <a href="https://jenkins.a2s.agh.edu.pl/pub/multistep-${BUILD_NUMBER}/movies.zip">Movies</a>
                        </li>
                        <li>
                        <a href="https://jenkins.a2s.agh.edu.pl/pub/multistep-${BUILD_NUMBER}/images.zip">Images</a>
                        </li>
                    </ul>
                 """,
                 charset: 'UTF-8',
                 mimeType: 'text/html',
                 replyTo: 'gurgul.grzegorz@gmail.com'
             )
         }  
         failure {  
             mail(
                 subject: "Simulations no. ${BUILD_NUMBER} failed",
                 to: "${EMAIL_RECIPIENTS}",
                 from: 'jenkins@a2s.agh.edu.pl',
                 body: """
                    The simulations failed.
                    Please see the <a href="https://jenkins.a2s.agh.edu.pl/job/IGA-ADS-SCAN/job/multistep-atari/${BUILD_NUMBER}/console">logs</a>.
                    Reply to this e-mail if you want to report an issue.
                 """,
                 charset: 'UTF-8',
                 mimeType: 'text/html',
                 replyTo: 'gurgul.grzegorz@gmail.com'
             )
         }  
    }
    
}