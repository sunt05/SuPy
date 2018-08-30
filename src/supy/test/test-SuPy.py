from unittest import TestCase

import supy as sp

class TestSuPy(TestCase):
    def test_is_string(self):
        s = sp.list_var_output
        self.assertTrue(isinstance(s[0], basestring))
