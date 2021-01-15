CAUTION (12/29/2016 updated)
----------------------------

Currently, this is ALPHA version i.e. just developed and only tested
internally. DO NOT use it if you don't want to run a risk. I plan to
replace my old plugin std_ErpCalc() once this is tested well.

Features
--------

-   Compatible with STUDY.design i.e. supports up to 2-way factorial
    design.
-   Uses STUDY's parameters to limit time ranges to plot.
-   Shows average ERP plus standard error (envelope) of a cluster for
    user-specified combination of conditions.
-   Supports EEGLAB pop_eegfiltnew()'s default Hamming low-pass filter
    (transient bandwidth == (upper pass-band edge \[Hz\])/3).
-   Re-defines baseline period.
-   Uses manually selected time windows to performs either 1) mean
    amplitude test, or 2) positive/negative peak amplitude and latency
    tests. The window selected is graphically shown in the plots.
-   Uses EEGLAB's statcond() function to perform either 1) permutation
    test, or 2) bootstrap test with 100,000 iterations. The results are
    shown in left bottom.
-   Reports how many unique subjects and number of ICs per each group in
    right bottom.
-   Reports and stores the individual subject's mean values. This would
    be useful for exporting data to SPSS to perform correlation
    analysis, for example.
    -   Note that there are missing subjects and redundant subjects in
        the output list. This group-level inconsistency derives from
        ICA. To address the issue of redundant ICs for some subjects,
        you should either 1) compute mean across the multiple entries,
        or 2) select the one that is closer to the cluster mean.
-   Uses the same amplitude and latency scale across the ERP plots.
-   For 2x2 design, reports simple effects i.e. (A-B) and (C-D)
    separately for your convenience. These are NOT main effects.
    -   By the way you do not need to perform ANOVA first to perform
        simple effect test; see [this
        paper](http://beheco.oxfordjournals.org/content/19/3/690.short).
-   Can make the background white for screen capturing and easy
    touch-up.

![Screenshot-std_erpstudio().png](images/Screenshot-std_erpstudio().png)

This page was written by Makoto Miyakoshi.