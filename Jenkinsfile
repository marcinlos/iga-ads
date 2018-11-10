static final def PHASE_SEPARATION_SIMULATION = "PHASE_SEPARATION"
static final def CUSTOM_SIMULATION = "PHASE_SEPARATION"
static final def PREDEFINED_SIMULATION_TYPES = [
    PHASE_SEPARATION_SIMULATION,
    CUSTOM_SIMULATION
]

def problem
def solverConfig

pipeline {

    agent { 
        node { 
            label 'node2'
        }
    }

    options {
        retry(2)
        timeout time: 24, unit:'HOURS'
    }

    environment {
        GIT_URL = "${scm.userRemoteConfigs[0].url}"
        GIT_BRANCH = "cahn-hilliard"
        BUILD_NUMBER = "${env.BUILD_NUMBER}"
    }

    parameters {
        choice(
            name: 'SIMULATION_TYPE',
            description: 'The type of the simulation to run. Predefined simulations contain sensible defaults.',
            choices: PREDEFINED_SIMULATION_TYPES
        )
    }

    stages {

        stage("Configure") {
            steps {
                echo "Loading defaults for ${params.SIMULATION_TYPE}"
                script {
                    problem = load(problemScriptFor(params.SIMULATION_TYPE))
                }
            }
        }

        stage("Define problem") {
            options {
                timeout(time:15, unit:'MINUTES')
            }

            steps {
                script {
                    SOLVER_ENV = input(
                        message: 'Customize defaults',
                        parameters: [
                            string(
                                name: 'ORDER',
                                defaultValue: "${problem.defaultOrder()}"
                            ),
                            string(
                                name: 'ELEMENTS',
                                defaultValue: "${problem.defaultElements()}"
                            ),
                            string(
                                name: 'STEPS',
                                defaultValue: "${problem.defaultSteps()}"
                            ),
                            string(
                                name: 'DELTA',
                                defaultValue: "${problem.defaultDelta()}"
                            ),
                            string(
                                name: 'MOBILITY_FORMULA',
                                defaultValue: "${problem.defaultMobilityFormula()}"
                            ),
                            string(
                                name: 'CHEMICAL_POTENTIAL_FORMULA',
                                defaultValue: "${problem.defaultChemicalPotentialFormula()}"
                            ),
                            text(
                                name: 'INITIAL_SURFACE_SNIPPET',
                                description: 'CPP code snippet which should return a value in x,y (x,y are double inputs)',
                                defaultValue: "${problem.defaultInitialSurfaceSnippet()}"
                            ),
                        ]
                    )
                }
                echo "Your parameters ${SOLVER_ENV}"
            }
        }

        stage('Checkout') {
            steps {
                echo "Checking out ${GIT_URL} on branch ${GIT_BRANCH}"
                git branch: "${GIT_BRANCH}", url: "${GIT_URL}"
            }
        }

        stage('Install dependencies') {
            steps {
                script {
                    def installDependenciesScript = load('jenkins/installDependencies.groovy')()
                }
                bash installDependenciesScript
            }
        }

        stage('Build') {
            steps {
                bash load('jenkins/build.groovy')(solverConfig)
            }
        }

        stage('Run') {
            steps {
                bash "./cahn_hilliard ${solverConfig.ORDER} ${solverConfig.ELEMENTS} ${solverConfig.STEPS} ${solverConfig.DELTA} ${solverConfig.MOBILITY_FORMULA} ${solverConfig.CHEMICAL_POTENTIAL_FORMULA}"
            }
        }

    }
    
}

@NonCPS
def problemScriptFor(String type) {
    def simulationDir = 'jenkins/simulations'
    switch(type) {
        case { it == env.PHASE_SEPARATION_SIMULATION}:
            return "${simulationDir}/phaseSeparation.groovy"
        default:
            return "${simulationDir}/custom.groovy"
    }
}
