TITLE Fast mechanism for submembranal Cl- concentration, version 0 (cld)

COMMENT

    Takes into account:
      - increase of cli due to chloride currents
      - extrusion of chloride modelled with a simple first order equation

    Parameters:
      - depth: depth of the shell just beneath the membrane (in um)
      - clinf: equilibrium concentration of chloride (10 mM)
      - tauKCC2: time constant of chloride extrusion (must be fast)

    2017-02-07 Adapted from cadecay.mod by Alain Destexhe, Salk Institute, 1995
    2017-02-22 This is not used; cldecay2.mod is used instead
    2017-03-07 Initialize cli from hoc
    2017-03-12 This is now an option
    2017-03-15 Made drive_channel & drive_extrusion RANGE variables
    2017-03-15 Forced drive_extrusion to be negative instead
    2017-03-31 Changed units of tauKCC2 to seconds
    2018-05-09 Changed tabs to spaces and set column width at 80
    2018-05-09 Improved comments

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX cld
    USEION cl READ icl, cli WRITE cli VALENCE -1
    RANGE depth, tauKCC2, clinf
    RANGE drive_channel, drive_extrusion
}

UNITS {
    (molar) = (1/liter)     : moles do not appear in units
    (mM)    = (millimolar)
    (um)    = (micron)
    (mA)    = (milliamp)
    (msM)   = (ms mM)
    FARADAY = (faraday) (coulombs)
}

PARAMETER {
    : RANGE variables whose values are specified in hoc
    depth   = .1    (um)    : depth of shell for Cl- (um)
    tauKCC2 = 30    (s)     : Cl- removal time constant (s), Peter's value
                            :   (Jedlicka et al 2011 used 3 s)
    clinf   = 8     (mM)    : steady state intracellular [Cl-] (mM)
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    icl             (mA/cm2)

    : RANGE variables that are assigned in the DERIVATIVE block
    drive_channel   (mM/ms) : driving Cl- flux (mM/ms) due to channel opening
    drive_extrusion (mM/ms) : driving Cl- flux (mM/ms) due to leak/extrusion
}
    
STATE {
    cli             (mM)    : intracellular [Cl-] (mM)
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
    : Calculate driving flux due to channel opening
    drive_channel = - (10000) * icl / ((-1) * FARADAY * depth)

    : Calculate driving flux due to leak/extrusion
    drive_extrusion = (clinf - cli) / (tauKCC2 * (1000))

    : Calculate change in chloride concentration
    cli' = drive_channel + (clinf - cli) / (tauKCC2 * (1000))
}

COMMENT
OLD CODE:

ENDCOMMENT
