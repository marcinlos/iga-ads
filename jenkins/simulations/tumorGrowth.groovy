def defaultOrder() {
    "3"
}

def defaultElements() {
    "50"
}

def defaultSteps() {
    "200"
}

def defaultDelta() {
    "0.0000000001"
}

def defaultChemicalPotentialFormula() {
    "4*(x^3-6*x^2+2*x)"
}

def defaultMobilityFormula() {
    "200"
}

def defaultInitialSurfaceSnippet() {
    """
    if ( (x-0.2)*(x-0.2)+(y-0.2)*(y-0.2)<=0.1*0.1 || (x-0.8)*(x-0.8)+(y-0.8)*(y-0.8)<=0.1*0.1 ) {
        return 0.8;
    } else {
        return 0.1;
    }
    """
}
return this
