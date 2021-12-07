# DahuHunting

This repository is part of the Eurocrypt submission 63.
It provides python modules to manipulate Boolean Functions and Rotational Symmetric Functions.

The BF module allows to work on Boolean Function objects, defined by their Truth Table (TT), Algebraic Normal Form (ANF) and Walsh Spectrum (WS). One of its representation can be arbitrarily modified (e.g. the set_TT method) and the other ones updated accordingly (e.g. update_ANF or update_WS). Once all representations are up-to-date, the resiliency and algebraic imminity of the function can then be verified (is_resilient or is_algebraic_immune).

The RSF module allows to work on Rotational Symmetric Functions, defined by their Simplified Truth Table (STT), Simplified Algebraic Normal Form (SANF) and Simplified Walsh Spectrum (SWS). The methods are similar to those of the BF class.

The find_BF and find_RSF modules provide some functions to search, more or less exhaustively, Boolean Functions or Rotational Symmetric Functions with a specified resiliency and/or algebraic immunity.

The file example.py replays somes results of the submission using the above mentionned modules.


# Current limitations

The repository is gradually updated when our original code is considered "proper enough" to be released.
All results of the submission can be retrieved using the current released version, except for Algorithm 2 of Section 7.2 and the exhaustive search of dahus in six variables of Section 7.1, the code of which is not robust enough yet. The examples given can still be verified with the BF or RSF classes.

A few transitions, such as updating WS to TT has not been implemented yet, since they were nor necessary for our purpose.
