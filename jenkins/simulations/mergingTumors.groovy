def defaultOrder() {
    "3"
}

def defaultElements() {
    "50"
}

def defaultSteps() {
    "250000"
}

def defaultDelta() {
    "0.0000000001"
}

def defaultChemicalPotentialFormula() {
    "4*(x^3-6*x^2+2*x)"
}

def defaultMobilityFormula() {
    "800*x*(1-x)"
}

def defaultInitialSurfaceSnippet() {
    """
    if ( (x-0.65)*(x-0.65)+(y-0.65)*(y-0.65)<=0.15*0.15 || (x-0.38)*(x-0.38)+(y-0.38)*(y-0.38)<=0.2*0.2 ) {
        return 0.8;
    } else {
        return 0.1;
    }
    """
}
return this