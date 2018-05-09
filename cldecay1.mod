TITLE Fast mechanism for submembranal Cl- concentration (cli)

COMMENT

    Takes into account:
       - chloride ion accumulation by chloride pump (Lineweaver-Burke equation) 
          and chloride leak
        - radial diffusion
       - longitudinal diffusion

    Diffusion model is modified from Ca diffusion model in Hines & Carnevale: 
    Expanding NEURON with NMODL, Neural Computation 12: 839-851, 2000 
        (Example 8)

    2017-03-12 Modified from cldecay2.mod, but now uses a simple 
                exponential decay for KCC2 as in cldif.mod
    2017-03-15 Made drive_channel & drive_extrusion RANGE variables
    2017-03-15 Added cli1
    2017-03-31 Changed units of tauKCC2 to seconds
    2018-05-09 Changed tabs to spaces and set column width at 80
    2018-05-09 Improved comments

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
    (mM)    = (millimolar)
    (um)    = (micron)
    (mA)    = (milliamp)
    FARADAY = (faraday) (10000 coulomb)
    PI      = (pi) (1)
}

PARAMETER {
    : RANGE variables whose values are specified in hoc
    depth   = 0.02  (1)     : relative depth of shell (to diameter) for Cl-
    DCl     = 2     (um2/ms): Cl- diffusion coefficient (um2/ms), 
                            :   Brumback & Staley 2008
                            :   also Kuner & Augustine 2000, Neuron 27: 447
    tauKCC2 = 30    (s)     : Cl- removal time constant (s), Peter's value 
                            :   (Jedlicka et al 2011 used 3 s)
    clinf   = 8     (mM)    : steady state intracellular [Cl-] (mM), 
                            :   Peter's value
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    diam            (um)    : diameter is defined in hoc, always in um
    icl             (mA/cm2): chloride current is written by gabaaCl

    : GLOBAL variables that are assigned in the INITIAL block
    vrat[Nannuli]   (1)     : numeric value of vrat[i] equals the volume
                            : of annulus i of a 1 um diameter cylinder
                            : multiply by diam^2 to get volume per unit length

    : RANGE variables that are assigned in the KINETIC block
    drive_channel   (um2 mM/ms) : driving Cl- flux (um2 mM/ms) 
                                :   due to channel opening
    cli             (mM)        : [Cl-] at outermost annulus just 
                                :   inside the membrane (mM)

    : RANGE variables that are assigned in the BREAKPOINT block
    drive_extrusion (um2 mM/ms) : driving Cl- flux (um2 mM/ms) due to KCC2
    cli1            (mM)        : [Cl-] at 2nd outermost annulus just 
                                :   inside the membrane (mM)
}

STATE {
    cl[Nannuli]     (mM)    <1e-10> : cl[0] is equivalent to cli
                                    : cl[] are very small, 
                                    :   so specify absolute tolerance
}

: Note: LOCAL variables are shared across sections but not visible in hoc
:       LOCAL variables are initialized at 0
LOCAL factorsDone      

INITIAL {
    : Calculate vrat & frat
    if (factorsDone == 0) {    
        : flag becomes 1 in the first segment
        :   to avoid unnecessary recalculation of vrat & frat
        : Note: vrat must GLOBAL, otherwise all subsequent segments 
        :   will have vrat = 0 
        factors()
        factorsDone = 1
    }

    : Initialize variables
    FROM i=0 TO Nannuli-1 {
        : Initialize chloride concentration in all annuli to be the same
        cl[i] = cli
    }
}

LOCAL frat[Nannuli]             : scales the rate constants for model geometry

PROCEDURE factors() {
    LOCAL r, hth

    : Start at edge of a cylinder with diameter 1 um
    r = 1/2                     

    : Compute half thickness of all other annuli
    hth = (r-depth)/(2*(Nannuli-1)-1)

    : Compute volume of outermost annulus
    vrat[0] = PI*(r-depth/2)*2*depth

    : Compute circumference/depth (this is not used)
    frat[0] = 2*PI*r/depth

    : Compute outer radius of second outermost annulus
    r = r - depth

    : Compute surface area per unit length in between annuli,
    :   i.e., distance between centers 
    frat[1] = 2*PI*r/(depth+hth)

    : Compute center radius for second outermost annulus
    r = r - hth

    : Compute volume of outer half of second outermost annulus
    vrat[1] = PI*(r+hth/2)*2*hth
    if (Nannuli > 2) {
        FROM i=1 TO Nannuli-2 {
            : Add volume of inner half of this annulus
            vrat[i] = vrat[i] + PI*(r-hth/2)*2*hth

            : Compute outer radius of next annulus
            r = r - hth

            : Compute circumference in between annuli/distance between centers 
            frat[i+1] = 2*PI*r/(2*hth)

            : Compute center radius for next annulus
            r = r - hth

            : Compute volume of outer half of next annulus
            vrat[i+1] = PI*(r+hth/2)*2*hth
        }
    }
}

: Note: can't define LOCAL in KINETIC block or use in COMPARTMENT statement
LOCAL volFactor                 : volume factor for unit correction

BREAKPOINT {
    : Note: This is the main execution block, 
    :   making everything consistent with time t
    SOLVE state METHOD sparse 

    : Compute driving flux due to KCC2
    : Compute volume factor for unit correction
    volFactor = diam * diam * vrat[0]                

    : Compute extrusion modeled with simple exponential
    drive_extrusion = ((clinf - cli)/(tauKCC2 * (1000))) * volFactor

    : For other mechanisms, cli1 is the concentration 
    :   at the second outermost annulus
    cli1 = cl[1]
}

KINETIC state {
    : The compartments are index, volume[index] {state}
    COMPARTMENT i, diam*diam*vrat[i] {cl}

    : There is diffusion between segments
    LONGITUDINAL_DIFFUSION i, DCl*diam*diam*vrat[i] {cl}

    : Compute volume factor for unit correction
    volFactor = diam * diam * vrat[0]

    : Equilibrium at the outermost annulus
    : Compute driving flux due to channel opening
    drive_channel = icl * PI * diam / FARADAY

    : Compute Cl- accumulation & extrusion with
    :   extrusion modeled with simple exponential
    ~ cl[0] << (drive_channel + ((clinf - cl[0])/(tauKCC2*(1000))) * volFactor) 
                                        

    : Equilibrium at other annuli
    FROM i=0 TO Nannuli-2 {
        : Compute Cl- radial diffusion
        ~ cl[i] <-> cl[i+1]    (DCl*frat[i+1], DCl*frat[i+1])
    }

    : For other mechanisms, cli is the concentration at the outermost annulus
    cli = cl[0]
}

COMMENT
OLD CODE:

ENDCOMMENT
