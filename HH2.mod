TITLE Hippocampal HH channels

COMMENT

    Fast Na+ and K+ currents responsible for action potentials
    Iterative equations

    Equations modified by Traub, for Hippocampal Pyramidal cells, in:
    Traub & Miles, Neuronal Networks of the Hippocampus, Cambridge, 1991

    range variable vtraub adjust threshold
    Written by Alain Destexhe, Salk Institute, Aug 1992

    Modified Oct 96 for compatibility with Windows: trap low values of arguments

    2017-02-07 Parameters modified to be consistent with Sohal & Huguenard 2003
    2017-02-10 Removed ena, eka, celsius from the PARAMETER block
    2017-03-05 Moved some variables from PARAMETER to ASSIGNED
    2017-03-07 Initialize to m_inf, h_inf, n_inf instead of 0
    2018-05-09 Changed tabs to spaces

ENDCOMMENT

NEURON {
    SUFFIX hh2
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    RANGE vtraub, gnabar, gkbar
    GLOBAL tadj
    RANGE m_inf, h_inf, n_inf
    RANGE tau_m, tau_h, tau_n
    RANGE m_exp, h_exp, n_exp
    RANGE m, h, n
    RANGE ina, ik
}


UNITS {
    (mA)    = (milliamp)
    (mV)    = (millivolt)
}

PARAMETER {
    : RANGE variables whose values are specified in hoc
    ena             (mV)        : Na+ reversal potential (mV)
    ek              (mV)        : K+ reversal potential (mV)
    vtraub  = -63   (mV)        : action potential threshold (mV)
    gnabar  = .05   (mho/cm2)   : Na+ conductance (S/cm2), 
                                :   Sohal & Huguenard 2003
    gkbar   = .005  (mho/cm2)   : K+ conductance (S/cm2), 
                                :   Sohal & Huguenard 2003
}


ASSIGNED {
    : Variables that are assigned outside the mod file
    v               (mV)
    celsius         (degC)
    dt              (ms)

    : GLOBAL variables that are assigned in the INITIAL block
    tadj            (1)         : temperature adjustion to time constants

    : RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
    m_inf           (1)
    h_inf           (1)
    n_inf           (1)
    tau_m           (ms)
    tau_h           (ms)
    tau_n           (ms)
    m_exp           (1)
    h_exp           (1)
    n_exp           (1)

    : RANGE variables that are assigned in the BREAKPOINT block
    ina             (mA/cm2)
    ik              (mA/cm2)
}

STATE {
    m h n
}

INITIAL {
:
:  Q10 was assumed to be 3 for both currents
:
    : Calculate constants from temperature
    tadj = 3.0 ^ ((celsius - 36 (degC))/ 10  (degC))

    : Initialize variables
    evaluate_fct(v)
    m = m_inf
    h = h_inf
    n = n_inf
}

BREAKPOINT {
    : Update gating variables m, h & n from voltage v
    SOLVE states METHOD cnexp

    : Calculate currents from m, h, n, ena, ek, v
    ina = gnabar * m*m*m*h * (v - ena)
    ik  = gkbar * n*n*n*n * (v - ek)
}


DERIVATIVE states {   : exact Hodgkin-Huxley equations
    evaluate_fct(v)
    m' = (m_inf - m) / tau_m
    h' = (h_inf - h) / tau_h
    n' = (n_inf - n) / tau_n
}

UNITSOFF

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

    v2 = v - vtraub : convert to traub convention

    a = 0.32 * vtrap(13-v2, 4)
    b = 0.28 * vtrap(v2-40, 5)
    tau_m = 1 / (a + b) / tadj
    m_inf = a / (a + b)

    a = 0.128 * Exp((17-v2)/18)
    b = 4 / ( 1 + Exp((40-v2)/5) )
    tau_h = 1 / (a + b) / tadj
    h_inf = a / (a + b)

    a = 0.032 * vtrap(15-v2, 5)
    b = 0.5 * Exp((10-v2)/40)
    tau_n = 1 / (a + b) / tadj
    n_inf = a / (a + b)

    m_exp = 1 - Exp(-dt/tau_m)
    h_exp = 1 - Exp(-dt/tau_h)
    n_exp = 1 - Exp(-dt/tau_n)
}

FUNCTION vtrap(x,y) {
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    } else {
        vtrap = x/(Exp(x/y)-1)
    }
}

FUNCTION Exp(x) {
    if (x < -100) {
        Exp = 0
    } else {
        Exp = exp(x)
    }
} 

UNITSON


COMMENT
OLD CODE:

ENDCOMMENT

