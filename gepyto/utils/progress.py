"""
Utilities to display task progress.
"""

# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


from __future__ import division, print_function


__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import sys

try:
    import progressbar
    HAS_PROGRESSBAR = True
except ImportError:
    HAS_PROGRESSBAR = False


def flush_print(s, end="\n"):
    sys.stdout.write("{}{}".format(s, end))
    sys.stdout.flush()


class Progress(object):
    def __init__(self, maxval, step=0.1):
        self.started = False
        self.maxval = maxval
        self.step = step
        self.next_milestone = self.step
        if HAS_PROGRESSBAR:
            self._bar = progressbar.ProgressBar(maxval=self.maxval)
        else:
            self._bar = None

    def update(self, val):
        if not self.started:
            self.started = True
            if self._bar is not None:
                self._bar.start()
            else:
                flush_print("0% ", end="")

        if self._bar is not None:
            return self._bar.update(val)

        # Check if we're at a new milestone.
        status = val / self.maxval
        if status > self.next_milestone:
            flush_print("{}% ".format(int(round(status * 100))), end="")
            self.next_milestone += self.step

    def finish(self):
        if self._bar is not None:
            return self._bar.finish()
        else:
            flush_print("done!")
