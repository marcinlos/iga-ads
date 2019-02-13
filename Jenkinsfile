pipeline {

    environment {
        GIT_URL = "${scm.userRemoteConfigs[0].url}"
        GIT_BRANCH = "multistep"
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

        stage('Initialise') {
            steps {
                sh '''#!/bin/bash
                    cat > initformula <<'EOL'
                    ${INITIAL_FORMULA_SNIPPET}
                    EOL
                    sed -i '/##INITSTART##/,/##INITEND##/!b;//!d;/##INITSTART##/r initformula' src/problems/multistep/multistep*d.hpp
                '''
            }
        }

        stage('Compile') {
            steps {
                sh '''#!/bin/bash
                    cmake -DUSE_GALOIS=ON . && make
                '''
            }
        }

        stage('Run') {
            steps {
                sh '''#!/bin/bash
                    ./multistep \
                    ${DIMENSION} \
                    ${SPLINEDEG} \
                    ${ELEMENTS} \
                    ${SCHEME} \
                    ${ORDER} \
                    ${STEPS} \
                    ${DELTA} \
                    > errors.sum
                '''
                stash name: 'results', includes: '*.data,*.sum'
            }
        }

        stage('Process results') {
            steps {
                sh '''#!/bin/bash
                    gnuplot plot
                    ./movie
                    echo "Compressing results\n"
                    zip data.zip *.data
                    zip images.zip *.png
                    zip movies.zip *.mp4

                    RESULTS_DIR=/usr/share/nginx/html/sim-$BUILD_NUMBER/

                    mkdir $RESULTS_DIR
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
                        <a href="https://jenkins.a2s.agh.edu.pl/pub/ch-${BUILD_NUMBER}/movies.zip">Movies</a>
                        </li>
                        <li>
                        <a href="https://jenkins.a2s.agh.edu.pl/pub/ch-${BUILD_NUMBER}/images.zip">Images</a>
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
                    Please see the <a href="https://jenkins.a2s.agh.edu.pl/job/IGA-ADS-SCAN/job/cahn-hilliard/${BUILD_NUMBER}/console">logs</a>.
                    Reply to this e-mail if you want to report an issue.
                 """,
                 charset: 'UTF-8',
                 mimeType: 'text/html',
                 replyTo: 'gurgul.grzegorz@gmail.com'
             )
         }  
    }
    
}