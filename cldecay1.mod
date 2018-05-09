TITLE Fast mechanism for submembranal Cl- concentration (cli)

COMMENT

    Takes into account:
       - chloride ion accumulation by chloride pump (Lineweaver-Burke equation) 
          and chloride leak
        - radial diffusion
       - longitudinal diffusion

    Diffusion model is modified from Ca diffusion model in Hines & Carnevale: 
    Expanding NEURON with NMODL, Neural Computation 12: 839-851, 2000 (Example 8)

    2017-03-12 Modified from cldecay2.mod, but now uses a simple exponential decay for KCC2 as in cldif.mod
    2017-03-15 Made drive_channel & drive_extrusion RANGE variables
    2017-03-15 Added cli1
    2017-03-31 Changed units of tauKCC2 to seconds
    2018-05-09 Changed tabs to spaces

ENDCOMMENT

NEURON {
    SUFFIX cld1
    USEION cl READ icl, cli WRITE cli VALENCE -1
    RANGE depth, DCl, tauKCC2, clinf
    GLOBAL vrat                : vrat must be GLOBAL, i.e., same across sections
    RANGE leak
    RANGE drive_channel, drive_extrusion, cli, cli1
}

DEFINE Nannuli 2

UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    FARADAY = (faraday) (10000 coulomb)
    PI = (pi) (1)
}

PARAMETER {
    : RANGE variables whose values are specified in hoc
    depth    = 0.02        (1)        : relative depth of shell (to diameter) for Cl-
    DCl    = 2         (um2/ms)     : Cl- diffusion coefficient (um2/ms), Brumback & Staley 2008
                        :     also Kuner & Augustine 2000, Neuron 27: 447
    tauKCC2    = 30        (s)        : Cl- removal time constant (s), Peter's value (Jedlicka et al 2011 used 3 s)
    clinf    = 8        (mM)        : steady state intracellular [Cl-] (mM)
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    diam             (um)        : diameter is defined in hoc, always in um
    icl             (mA/cm2)    : chloride current is written by gabaaCl

    : GLOBAL variables that are assigned in the INITIAL block
    vrat[Nannuli]        (1)        : numeric value of vrat[i] equals the volume
                        : of annulus i of a 1 um diameter cylinder
                        : multiply by diam^2 to get volume per unit length

    : RANGE variables that are assigned in the KINETIC block
    drive_channel        (um2 mM/ms)    : driving Cl- flux (um2 mM/ms) due to channel opening
    cli             (mM)        : [Cl-] at outermost annulus just inside the membrane (mM)

    : RANGE variables that are assigned in the BREAKPOINT block
    drive_extrusion        (um2 mM/ms)    : driving Cl- flux (um2 mM/ms) due to KCC2
    cli1             (mM)        : [Cl-] at 2nd outermost annulus just inside the membrane (mM)
}

STATE {
    cl[Nannuli]        (mM)     <1e-10>    : cl[0] is equivalent to cli
                        : cl[] are very small, so specify absolute tolerance
}


LOCAL factors_done                : LOCAL variables are shared across sections but not visible in hoc
                        : LOCAL variables are initialized at 0

INITIAL {
    : Calculate vrat & frat
    if (factors_done == 0) {        : flag becomes 1 in the first segment
                        :     to avoid unnecessary recalculation of vrat & frat
        factors()
        factors_done = 1        : Note: vrat must GLOBAL, otherwise all subsequent segments will have vrat = 0 
    }

    : Initialize variables
    FROM i=0 TO Nannuli-1 {
        cl[i] = cli            : initialize chloride concentration in all annuli to be the same
    }
}

LOCAL frat[Nannuli]                : scales the rate constants for model geometry

PROCEDURE factors() {
    LOCAL r, hth

    r = 1/2                    : start at edge of a cylinder with diameter 1 um
    hth = (r-depth)/(2*(Nannuli-1)-1)    : half thickness of all other annuli

    vrat[0] = PI*(r-depth/2)*2*depth    : volume of outermost annulus
    frat[0] = 2*PI*r/depth            : circumference/depth (this is not used)
    r = r - depth                : outer radius of second outermost annulus
    frat[1] = 2*PI*r/(depth+hth)        : surface area per unit length in between annuli/distance between centers 
    r = r - hth                : center radius for second outermost annulus
    vrat[1] = PI*(r+hth/2)*2*hth        : volume of outer half of second outermost annulus
    if (Nannuli > 2) {
        FROM i=1 TO Nannuli-2 {
            vrat[i] = vrat[i] + PI*(r-hth/2)*2*hth    : add volume of inner half of this annulus
            r = r - hth                : outer radius of next annulus
            frat[i+1] = 2*PI*r/(2*hth)        : circumference in between annuli/distance between centers 
            r = r - hth                : center radius for next annulus
            vrat[i+1] = PI*(r+hth/2)*2*hth        : volume of outer half of next annulus
        }
    }
}

LOCAL vol_fac                    : can't define LOCAL in KINETIC block or use in COMPARTMENT statement

BREAKPOINT {                     : this is the main execution block, makes everything consistent with time t
    SOLVE state METHOD sparse 

    : Calculate driving flux due to KCC2
    vol_fac = diam*diam*vrat[0]                : volume factor for unit correction
    drive_extrusion = ((clinf - cli)/(tauKCC2*(1000)))*vol_fac    : extrusion modeled by simple exponential    

    : For other mechanisms, cli1 is the concentration at the second outermost annulus
    cli1 = cl[1]
                
}

KINETIC state {
    COMPARTMENT i, diam*diam*vrat[i] {cl}            : index, volume[index] {state}
    LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}    : diffusion between segments
    vol_fac = diam*diam*vrat[0]                : volume factor for unit correction

    : Equilibrium at the outermost annulus
    drive_channel = icl*PI*diam/FARADAY            : calculate driving flux due to icl
    ~ cl[0] << (drive_channel + ((clinf - cl[0])/(tauKCC2*(1000)))*vol_fac) : Cl- accumulation & extrusion
                                        : extrusion modeled by simple exponential

    : Equilibrium at other annuli
    FROM i=0 TO Nannuli-2 {
        ~ cl[i] <-> cl[i+1]    (DCl*frat[i+1], DCl*frat[i+1])    : Cl- diffusion
    }

    : For other mechanisms, cli is the concentration at the outermost annulus
    cli = cl[0]
}

COMMENT
OLD CODE:
    : ~ cl[0] << (icl*PI*diam/FARADAY + ((clinf - cl[0])/tauKCC2)*vol_fac) : positive icl is Cl- influx 
    : if (drive_channel <= 0.) { drive_channel = 0. }    : since the pump is outward, do not let cli decrease from icl

: 2017-03-15 Forced drive_extrusion to be negative instead
    if (drive_extrusion > 0.) { drive_extrusion = 0. }    : KCC2 pumps Cl- outward only

    tauKCC2    = 30000        (ms)        : Cl- removal time constant (ms), Peter's value (Jedlicka et al 2011 used 3000 ms)

ENDCOMMENT
