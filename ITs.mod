TITLE Low threshold calcium current

COMMENT

    Ca++ current responsible for low threshold spikes (LTS)
    RETICULAR THALAMUS
    Differential equations

    Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
    The kinetics is described by standard equations (NOT GHK)
    using a m2h format, according to the voltage-clamp data
    (whole cell patch clamp) of Huguenard & Prince, J Neurosci.
    12: 3804-3817, 1992.

     - Kinetics adapted to fit the T-channel of reticular neuron
     - Time constant tau_h refitted from experimental data
     - shift parameter for screening charge

    Model described in detail in:
      Destexhe, A., Contreras, D., Steriade, M., Sejnowski, T.J. and
      Huguenard, J.R.  In vivo, in vitro and computational analysis of
      dendritic calcium currents in thalamic reticular neurons.
      Journal of Neuroscience 16: 169-185, 1996.

    Written by Alain Destexhe, Salk Institute, Sept 18, 1992
    Modified by Vikaas Sohal to use constant-field equation, August, 1997

    2017-02-07 Modified the Faraday constant
    2017-03-04 Changed the integration method from euler to cnexp
    2017-03-05 Moved some variables from PARAMETER to ASSIGNED
    2017-03-05 Made coeff a RANGE variable (so that it is thread safe)
    2017-03-05 Fixed units for the GHK equation when v is small
    2017-03-07 Use permeability pcabar directly, make it a RANGE variable
    2017-03-07 Renamed P to gpc
    2018-05-09 Changed tabs to spaces

ENDCOMMENT

NEURON {
    SUFFIX Ts
    USEION ca READ cai, cao WRITE ica
    GLOBAL qm, qh
    RANGE pcabar
    GLOBAL phi_m, phi_h, coeff
    RANGE gpc
    RANGE m_inf, h_inf, tau_m, tau_h
    RANGE m, h
    RANGE ica
}

UNITS {
    (molar) = (1/liter)
    (mV) =    (millivolt)
    (mA) =    (milliamp)
    (mM) =    (millimolar)
    FARADAY    = (faraday) (coulombs)
    R    = (k-mole) (joule/degC)
}

PARAMETER {
    : GLOBAL variables whose values are specified in hoc
    qm    = 2.5    (1)        : q10 for T channel activation
    qh     = 2.5    (1)        : q10 for T channel inactivation

    : RANGE variables whose values are specified in hoc
    pcabar    = 1e-4  (cm/s)        : Ca++ permeability (cm/s) for T channels
}

ASSIGNED {
    : Variables that are assigned outside the mod file
    v        (mV)
    celsius        (degC)
    cai        (mM)        : intracellular [Ca++] (mM)
    cao        (mM)        : exftracellular [Ca++] (mM)

    : GLOBAL variables that are assigned in the INITIAL block
    phi_m        (1)        : temperature adjustion to tau_m
    phi_h        (1)        : temperature adjustion to tau_h
    coeff        (1/mV)        : zF/RT, which is about 1/(13 mV) at 34 degC

    : RANGE variables that are assigned in the INITIAL block
    gpc        (mho/cm2 mM)    : conductance per unit concentration (mho/cm2 mM) (permeability times (zF)^2/RT)
                    : for pcabar = 1e-4 cm/s, we have gpc = 0.0015 mho/cm2 mM at 34 degC
                    : which for cao = 2 mM means a conductance of about g = 0.003 mho/cm2

    : RANGE variables that are assigned in the INITIAL & DERIVATIVE blocks
     m_inf        (1)
    h_inf        (1)
    tau_m        (ms)
    tau_h        (ms)

    : RANGE variables that are assigned in the BREAKPOINT block
    ica        (mA/cm2)    : calcium current (mA/cm2) through the T channel
}

STATE {
    m        (1)        : activation gate
    h        (1)        : inactivation gate
}

INITIAL {
:
: Activation functions and kinetics were obtained from
: Huguenard & Prince, and were at 23-25 deg.
: Transformation to celsius using Q10 (qm & qh)
:

    : Calculate constants from temperature
    phi_m = qm ^ ((celsius-24 (degC))/10 (degC))
    phi_h = qh ^ ((celsius-24 (degC))/10 (degC))
    coeff = (1e-3)*2*FARADAY/(R*(273.15 (degC) + celsius))            : conversion factor (1e-3 V/mV)

    : Calculate gpc from temperature & pcabar
    gpc = pcabar*(1e-6)*((2*FARADAY)^2)/(R*(273.15 (degC) + celsius))    : conversion factor (1e-6 ML/cm3)

    : Initialize variables
    ica = 0                : no current in the beginning %%% TO EXAMINE
    evaluate_fct(v)            : calculate m_inf & h_inf
    m = m_inf
    h = h_inf
}

BREAKPOINT { LOCAL tmp
    : Update gating variables m & h from voltage v
    SOLVE castate METHOD cnexp

    : Calculate calcium current ica from m, h, cai, cao, v
    if ((v > .1) || (v < -1)) {                    : %%% why not v < -.1?
        tmp = exp(-coeff*v)
        ica = m*m*h * gpc * (cai - cao*tmp)/(1 - tmp) * v
    } else {
        tmp = exp(-coeff*v)
        ica = m*m*h * gpc * (cai - cao*tmp) * (1/coeff)
    }
}

DERIVATIVE castate {
    : Update m_inf, h_inf, tau_m, tau_h from v
    evaluate_fct(v)

    : Update m & h from m_inf, h_inf, tau_m, tau_h
    m' = (m_inf - m) / tau_m
    h' = (h_inf - h) / tau_h
}

PROCEDURE evaluate_fct(v (mV)) { 
:
: Time constants were obtained from J. Huguenard
:

UNITSOFF
    m_inf = 1.0 / ( 1 + exp(-(v+50)/7.4) )
    h_inf = 1.0 / ( 1 + exp((v+78)/5.0) )

    tau_m = ( 3 + 1.0 / ( exp((v+25)/10) + exp(-(v+100)/15) ) ) / phi_m
    tau_h = ( 85 + 1.0 / ( exp((v+46)/4) + exp(-(v+405)/50) ) ) / phi_h
UNITSON
}

COMMENT
OLD CODE:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}    : not necessary in NEURON

    gcabar    = .003    (mho/cm2)

    FARADAY = 96520 (coul)
    celsius    = 37    (degC)
    FARADAY = 96485.3329    (coul)        : moles do not appear in units
    R = 8.314 (V coul / degC)

    cai    = 2.4e-4 (mM)        : adjusted for eca=120 mV    : this is ignored and set by hoc
    cao    = 2    (mM)                        : this is ignored and set by hoc

     GLOBAL P, qm, qh, coeff, p
    GLOBAL P, qm, qh

    SOLVE castate METHOD euler                : "NOT THREAD SAFE"

     P = P/(.0029 (cm2))
    P    = .0067    (mho/cm2 mM)        : permeability times (zF)^2/RT

        ica = m*m*h*P*R*(celsius+273 (degC))*(cai-cao*tmp)    : %%% Units are wrong, to understand

    coeff = (.001)*2*FARADAY/(R*(273 (degC) + celsius))
     p    = 1.8e-8    (cm3/s)    : permeability of T channels %%% Wrong units ???
    P = p*4*(FARADAY^2)/(R*(273 (degC) + celsius))            : permeability times (zF)^2/RT


ENDCOMMENT


