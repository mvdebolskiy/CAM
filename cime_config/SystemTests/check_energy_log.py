#!/usr/bin/env python3

"""
Given two log files, compare the diagnostics from check_energy.
Report any differences.
Treat the first file as the 'standard' for reporting of relative differences.
"""

import argparse
import gzip
import os
import re
import sys

_DEFAULT_EPS = 1.0e-14

class LogEntry:
    """Class to hold the data from a single energy diagnostic entry from a
    CAM log file."""

    __entry_re = re.compile(r"nstep, te")
    __nstep_diff_msg = "nstep differs"
    __table_spc = 24
    __log_entries = {"nstep" : "nstep",
                     "te_input" : "TE input", "te_output" : "TE output",
                     "heat" : "Heat", "p_surf" : "Surf. Pres.", "p_top" : "Pres. top"}

    def __init__(self, input_line, linenum=None, eps=_DEFAULT_EPS):
        """Initialize the log entry directly from a line of a CAM log file.
        If not None, <linenum> is the input line number."""
        if self.is_log_entry(input_line):
            tokens = [x for x in input_line.strip().split(" ") if x]
            if len(tokens) != len(self.__log_entries) + 2:
                raise ValueError(f"Bad atm log input line, '{input_line.strip()}'")
            # end if
            self.__linenum = linenum
            self.__nstep = int(tokens[2])
            self.__te_input = float(tokens[3])
            self.__te_output = float(tokens[4])
            self.__heat = float(tokens[5])
            self.__p_surf = float(tokens[6])
            self.__p_top = float(tokens[7])
            self.__eps = eps
        else:
            raise ValueError(f"Input line ({input_line}) is not a valid entry")
        # end if

    def compare(self, log_entry):
        """Compare <self> against <log_entry>.
        If the entries are equivalent, return an empty list
        If the entries differ, return a list of difference messages."""
        errs = []
        if self.__nstep != log_entry.nstep:
            errs.append(self.__nstep_diff_msg)
            errs.append(f"File 1: nstep = {self.__nstep}")
            errs.append(f"File 2: nstep = {log_entry.__nstep}")
        else:
            diff = self._find_val_diffs(self.__te_input,
                                        log_entry.__te_input, "te_input")
            if diff:
                errs.append(diff)
            # end if
            diff = self._find_val_diffs(self.__te_output,
                                        log_entry.__te_output, "te_output")
            if diff:
                errs.append(diff)
            # end if
            diff = self._find_val_diffs(self.__heat,
                                        log_entry.__heat, "heat")
            if diff:
                errs.append(diff)
            # end if
            diff = self._find_val_diffs(self.__p_surf,
                                        log_entry.__p_surf, "p_surf")
            if diff:
                errs.append(diff)
            # end if
            diff = self._find_val_diffs(self.__p_top,
                                        log_entry.__p_top, "p_top")
            if diff:
                errs.append(diff)
            # end if
        # end if

        return errs

    def _find_val_diffs(self, val1, val2, key):
        """Compare two real values. If they are equal, return None.
        If the two values differ, return a list:
           Value Header
           File 1 value
           File 2 value
           relative difference
        <key> is the type of value being compared
        """
        err = None
        if val1 == val2:
            return err
        # end if
        if val1 > self.__eps:
            diff = abs(val2 - val1) / val1
        else:
            diff = abs(val2 - val1)
        # end if
        f1val = f"{val1:18.16e}"
        f2val = f"{val1:18.16e}"
        diffval = f"{diff:5.3e}"
        err = [' '*6+f"{self.__log_entries[key]}"+' '*(self.__table_spc-len(self.__log_entries[key])-6),
                f1val+' '*(self.__table_spc-len(f1val)), f2val+' '*(self.__table_spc-len(f2val)),
                diffval+' '*(self.__table_spc-len(diffval))]
        return err

    @property
    def nstep(self):
        """Return the step number of this log entry"""
        return self.__nstep

    @property
    def linenum(self):
        """Return the source line number of this log entry"""
        return self.__linenum

    @classmethod
    def is_log_entry(cls, input_line):
        """Return True iff <input_line> is a CAM log file energy diagnostic"""
        return cls.__entry_re.search(input_line.strip()) is not None

    @classmethod
    def nstep_diff_msg(cls):
        """Return the special nstep difference error message string"""
        return cls.__nstep_diff_msg

###############################################################################

def gather_energy_diags_from_log(pathname):
    """Gather the energy logs from <pathname>."""
    if not os.path.exists(pathname):
        raise ValueError(f"input file, '{pathname}', does not exist")
    # end if
    log_entries = {}
    linenum = 0
    try:
        infile = gzip.open(pathname, "rt")
    except gzip.BadGzipFile as exc:
        infile = open(pathname, "rt")
    # end try

    for line in infile:
        linenum += 1
        if LogEntry.is_log_entry(line):
            new_entry = LogEntry(line, linenum=linenum)
            if new_entry.nstep in log_entries:
                emsg = f"Duplicate nstep {new_entry.nstep} on line {linenum}."
                emsg += f"\nOriginal entry on line {log_entries[new_entry.nstep].linenum}."
                raise ValueError(emsg)
            # end if
            log_entries[new_entry.nstep] = new_entry
        # end if (no else, just ignore other lines)
    # end for

    # Don't leave without this
    infile.close()

    return log_entries

###############################################################################

def compare_log_files(logfile1, logfile2):
    """Compare the diagnostics from check_energy of <logfile1> vs. <logfile2>
    Report any differences.
    Treat <logfile1> as the 'standard' for reporting of relative differences.
    Return True iff the files are considered equivalent
    """

    errs = []
    if not os.path.exists(logfile1):
        raise ValueError(f"logfile1, '{logfile1}', does not exist")
    # end if
    if not os.path.exists(logfile2):
        raise ValueError(f"logfile2, '{logfile2}', does not exist")
    # end if

    log_entries1 = gather_energy_diags_from_log(logfile1)
    log_entries2 = gather_energy_diags_from_log(logfile2)

    # We go through the entries from logile2 since it should have a subset
    # of the entries from logfile1.
    # It is an error for an entry in logfile2 to not be represented on logfile1
    for nstep in log_entries2:
        if nstep in log_entries1:
            diffs = log_entries1[nstep].compare(log_entries2[nstep])
            if diffs:
                if diffs[0] == LogEntry.nstep_diff_msg:
                    errs.extend(diffs)
                else:
                    errs.append(f"\nnstep {nstep}:")
                    errs.append(' '*7+f"{''.join([x[0] for x in diffs])}")
                    errs.append(f"File  1: {''.join([x[1] for x in diffs])}")
                    errs.append(f"File  2: {''.join([x[2] for x in diffs])}")
                    errs.append(f"RelDiff: {''.join([x[3] for x in diffs])}")
                # end if
            # end if
        else:
            errs.append(f"Log entry for nstep = {nstep} missing on File 2")
        # end if
    # end for
    if errs:
        print(f"Energy diagnostics differ:")
        print(f"File 1: {logfile1}\nFile 2: {logfile2}:")
        print("\n".join(errs))
    else:
        print(f"{logfile1} and {logfile2} have equivalent energy diagnostics:")
    # end if

###############################################################################

def _main_func():
    """Parse the two files to check from the command line
    and call the main compare function.
    Raise an exception if the input line does not contain two items.
    The return value is True if the files are equivalent"""
    if len(sys.argv) != 3:
        raise ValueError(f"Usage: {sys.argv[0]} <pathname1> <pathname2>")
    # end if
    logfile1 = sys.argv[1]
    logfile2 = sys.argv[2]
    compOK = compare_log_files(logfile1, logfile2)
    return compOK

###############################################################################

if __name__ == "__main__":
    _main_func()
