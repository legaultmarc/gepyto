
# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"

import unittest

from ..structures import region


class TestRegion(unittest.TestCase):

    def setUp(self):
        pass

    def test_segment_case1(self):
        """Vanilla overlap."""
        li = [
            region._Segment("3", 1, 5),
            region._Segment("3", 7, 20),
            region._Segment("3", 12, 25),
            region._Segment("3", 24, 31),
            region._Segment("3", 35, 42)
        ]

        merged = region._Segment.merge_segments(li)

        merged_ans = [
            region._Segment("3", 1, 5),
            region._Segment("3", 7, 31),
            region._Segment("3", 35, 42),
        ]
        self.assertEqual(merged, merged_ans)

    def test_segment_case2(self):
        """One large chunk overlapping two smaller in the middle."""
        li = [
            region._Segment("3", 1, 2),
            region._Segment("3", 4, 15),
            region._Segment(3, 5, 6),
            region._Segment("3", 8, 12),
            region._Segment("3", 18, 21)
        ]

        merged = region._Segment.merge_segments(li)

        merged_ans = [
            region._Segment("3", 1, 2),
            region._Segment("3", 4, 15),
            region._Segment("3", 18, 21),
        ]
        self.assertEqual(merged, merged_ans)

    def test_segment_case3(self):
        """Overlap on the first one."""
        li = [
            region._Segment("3", 1, 5),
            region._Segment("3", 2, 11),
            region._Segment("3", 13, 21),
        ]

        merged = region._Segment.merge_segments(li)

        merged_ans = [
            region._Segment("3", 1, 11),
            region._Segment("3", 13, 21),
        ]
        self.assertEqual(merged, merged_ans)

    def test_segment_case4(self):
        """Overlapping border."""
        li = [
            region._Segment("3", 1, 5),
            region._Segment("3", 7, 11),
            region._Segment("3", 7, 15),
            region._Segment("3", 21, 26),
        ]

        merged = region._Segment.merge_segments(li)

        merged_ans = [
            region._Segment("3", 1, 5),
            region._Segment("3", 7, 15),
            region._Segment("3", 21, 26),
        ]
        self.assertEqual(merged, merged_ans)

    def test_segment_case5(self):
        """Different chromosomes raises Exception."""
        li = [
            region._Segment("3", 1, 5),
            region._Segment("5", 7, 11),
            region._Segment("3", 21, 26),
        ]

        self.assertRaises(Exception, region._Segment.merge_segments, li)

    def test_segment_case6(self):
        """Unsorted list raises Exception."""
        li = [
            region._Segment("3", 1, 5),
            region._Segment("3", 7, 11),
            region._Segment("3", 21, 26),
            region._Segment("3", 7, 15),
        ]

        self.assertRaises(Exception, region._Segment.merge_segments, li)

    def test_two_overlapping_segments(self):
        li = [
            region._Segment("3", 1, 5),
            region._Segment("3", 2, 6),
        ]

        merged = region._Segment.merge_segments(li)

        merged_ans = [
            region._Segment("3", 1, 6),
        ]
        self.assertEqual(merged, merged_ans)

    def test_two_not_overlapping_segments(self):
        li = [
            region._Segment("3", 1, 5),
            region._Segment("3", 10, 20),
        ]

        merged = region._Segment.merge_segments(li)

        merged_ans = [
            region._Segment("3", 1, 5),
            region._Segment("3", 10, 20),
        ]
        self.assertEqual(merged, merged_ans)
