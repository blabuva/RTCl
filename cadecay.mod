TITLE Fast mechanism for submembranal Ca++ concentration (cai)

COMMENT

    Takes into account:
       - increase of cai due to calcium currents
       - extrusion of calcium with a simple first order equation

    This mechanism is compatible with the calcium pump "cad" and has the 
    same name and parameters; however the parameters specific to the pump
    are dummy here.

    Parameters:
       - depth: depth of the shell just beneath the membran (in um)
       - cainf: equilibrium concentration of calcium (2e-4 mM)
       - taur: time constant of calcium extrusion (must be fast)
       - kt,kd: dummy parameters

    Written by Alain Destexhe, Salk Institute, 1995

    2017-02-07 Modified the Faraday constant
    2017-02-07 Parameters modified to be consistent with Sohal & Huguenard 2003
    2017-03-07 Initialize cai from hoc
    2018-05-09 Changed tabs to spaces

ENDCOMMENT

NEURON {
    SUFFIX cad
    USEION ca READ ica, cai WRITE cai
    RANGE depth, taur, cainf
    RANGE drive_channel
}

UNITS {
    (molar) = (1/liter)        : moles do not appear in units
    (mM)    = (millimolar)
    (um)    = (micron)
    (mA)    = (milliamp)
    (msM)    = (ms mM)
    FARADAY    = (faraday) (coulombs)
}

PARAMETER {
    : RANGE variables whose values are specified in hoc
    depth    = .1        (um)        : depth of shell for Ca++ (um)
    taur    = 24        (ms)        : Ca++ removal time constant (ms), Sohal & Huguenard 2003
    cainf    = 2.4e-4    (mM)        : steady state intracellular [Ca++] (mM)
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    ica            (mA/cm2)

    : RANGE variables that are assigned in the DERIVATIVE block
    drive_channel        (mM/ms)
}
    
STATE {
    cai            (mM)        : intracellular [Ca++] (mM)
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
    
    : Calculate flux due to ica
    drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
    if (drive_channel <= 0.) { drive_channel = 0. }    : since the pump is outward, do not let cai decrease from ica

    : Calculate change in calcium concentration
    cai' = drive_channel + (cainf-cai)/taur
}

COMMENT
OLD CODE:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}    : not necessary in NEURON

    FARADAY = 96489        (coul)        : moles do not appear in units
:    FARADAY = 96.489    (k-coul)    : moles do not appear in units
    taur    = 5    (ms)        : rate of calcium removal, Dextexhe et al 1996
    cainf    = 2e-4    (mM)        : steady state concentration of calcium inside the cell
CONSTANT {
    FARADAY = 96485.3329    (coul)    : moles do not appear in units
}

    RANGE depth, taur, cainf, kt, kd
    kt    = 0        (mM/ms)        : dummy
    kd    = 0        (mM)        : dummy

    if (drive_channel <= 0.) { drive_channel = 0. }    : cannot pump inward

INITIAL {
    : Initialize variables
    cai = cainf
}

ENDCOMMENT
