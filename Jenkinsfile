enum SimulationType {
    PHASE_SEPARATION,
    CUSTOM
]

def problemScriptFor(SimulationType type) {
    def simulationDir = 'jenkins/simulations'
    switch(type) {
        case PHASE_SEPARATION:
            return "${simulationDir}/phaseSeparation.groovy"
        default:
            return "${simulationDir}/custom.groovy"
    }
}

def problem
def solverConfig

pipeline {

    agent { 
        node { 
            abel 'node2'
        }
    }

    options {
        retry(2)
        timeout time: 20, unit:'MINUTES'
    }

    environment {
        GIT_URL = scm.userRemoteConfigs[0].url
        GIT_BRANCH = "cahn-hilliard"
        BUILD_NUMBER = "${env.BUILD_NUMBER}"
    }

    parameters {
        choice(
            name: 'SIMULATION_TYPE',
            description: 'The type of the simulation to run. Predefined simulations contain sensible defaults.',
            defaultValue: SimulationType.PHASE_SEPARATION
        )
    }

    stages {
        stage("Define problem") {
            steps {
                script {
                    problem = load(problemScriptFor(params.SIMULATION_TYPE))
                }

                script {
                    timeout(time:1, unit:'HOURS') {
                        solverConfig = input(
                            message: 'Customize defaults',
                            parameters: [
                                string(
                                    name: 'ORDER',
                                    defaultValue: problem.defaultOrder()
                                ),
                                string(
                                    name: 'ELEMENTS',
                                    defaultValue: problem.defaultElements()
                                ),
                                string(
                                    name: 'STEPS',
                                    defaultValue: problem.defaultSteps()
                                ),
                                string(
                                    name: 'DELTA',
                                    defaultValue: problem.defaultDelta()
                                ),
                                string(
                                    name: 'MOBILITY_FORMULA',
                                    defaultValue: problem.defaultMobilityFormula()
                                ),
                                string(
                                    name: 'CHEMICAL_POTENTIAL_FORMULA',
                                    defaultValue: problem.defaultChemicalPotentialFormula()
                                ),
                                text(
                                    name: 'INITIAL_SURFACE_SNIPPET',
                                    description: 'CPP code snippet which should return a value in x,y (x,y are double inputs)',
                                    defaultValue: problem.defaultInitialSurfaceSnippet()
                                )
                            ]
                        )
                    }
                }
            }
        }

        stage('Checkout') {
            steps {
                git branch: "${params.GIT_BRANCH}", url: "${GIT_URL}"
            }
        }

        stage('Install dependencies') {
            steps {
                bash load('jenkins/installDependencies.groovy')()
            }
        }

        stage('Build') {
            steps {
                bash load('jenkins/build.groovy')(solverConfig)
            }
        }

        stage('Run') {
            steps {
                bash ./cahn_hilliard solverConfig.ORDER solverConfig.ELEMENTS solverConfig.STEPS $DELTA solverConfig.MOBILITY_FORMULA solverConfig.CHEMICAL_POTENTIAL_FORMULA
            }
        }

    }
    
}