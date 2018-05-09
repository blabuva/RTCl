TITLE Pulse train generator

COMMENT

    ipulse2.mod
    Generates a train of current pulses
    User specifies dur (pulse duration), per (period, i.e. interval between pulse onsets),
    and number of pulses.
    Ensures that period is longer than pulse duration.
    2/6/2002 NTC
    http://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=137

    2018-05-09 Changed tabs to spaces

ENDCOMMENT

NEURON {
    POINT_PROCESS Ipulse2
    ELECTRODE_CURRENT i
    RANGE del, dur, amp, per, num, i, on
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    del    (ms)        : stimulation delay (ms)
    dur    (ms) <0, 1e9>    : current pulse duration (ms)
    amp    (nA)        : current pulse amplitude (nA)
    per    (ms) <0, 1e9>    : current pulse period (ms), i.e. interval between pulse onsets
    num    (1)        : number of current pulses
}

ASSIGNED {
    ival    (nA)
    i    (nA)
    on    (1)        : whether current pulse is on
    tally    (1)        : how many more current pulses to deliver
}

INITIAL {
    if (dur >= per) {
        per = dur + 1 (ms)
        printf("per must be longer than dur\n")
        UNITSOFF
        printf("per has been increased to %g ms\n", per)
        UNITSON
    }
    i = 0
    ival = 0
    tally = num
    if (tally > 0) {
        net_send(del, 1)
        on = 0
        tally = tally - 1
    }
}

BREAKPOINT {
    i = ival
}

NET_RECEIVE (w) {
    if (flag == 1) {                : ignore any but self-events with flag == 1
        if (on == 0) {
            : turn it on
            ival = amp
            on = 1
            : prepare to turn it off
            net_send(dur, 1)
        } else {
            : turn it off
            ival = 0
            on = 0
            if (tally > 0) {
                : prepare to turn it on again
                net_send(per - dur, 1)
                tally = tally - 1
            }
        }
    }
}

COMMENT
OLD CODE:
    printf("event sent!\n")
    printf("flag == %d\n", flag)

ENDCOMMENT
