# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram, TDPourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO, TDPourbaixEntry


class TestPourbaixDiagram(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self._entries = entries
        self._pd = PourbaixDiagram(entries)
        self.list_of_stable_entries = ["ZnO(s)", "Zn[2+]", "ZnO2(s)", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"]

        
    def test_pourbaix_diagram(self):
        self.assertEqual(len(self._pd.facets), 6, "Incorrect number of facets")
        for entry in self._pd.stable_entries:
            self.assertIn(entry.name, self.list_of_stable_entries, "List of stable entries does not match")

class TestTDPourbaixDiagram(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, nontd_entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                             "test_entries.csv"))
        td_entries = []
        for nontd_entry in nontd_entries:
            entropy = 0.0
            temperature = 298.15
            correction = 0.0
            if nontd_entry.name == "ZnHO[+]":
                entropy = 0.1
            td_entry = TDPourbaixEntry(entry=nontd_entry.raw_entry,
                                       gibbs_energy=nontd_entry.energy,
                                       entropy=entropy,
                                       temperature=temperature,
                                       correction=correction)
            td_entries.append(td_entry)
        self._entries = td_entries
        self._pd = TDPourbaixDiagram(entries=td_entries,
                                     plot_type="T_pH",
                                     const_pot=3.0)
        self.list_of_stable_entries = ["ZnHO[+]", "ZnO(s)", "Zn[2+]", "ZnO2(s)", "ZnHO2[-]", "ZnO2[2-]", "Zn(s)"]

    def test_td_pourbaix_diagram(self):
        print("Number of factes:", len(self._pd.facets))
        print("stable entries:")
        for entry in self._pd.stable_entries:
            print(entry.name)


if __name__ == '__main__':
    unittest.main()
