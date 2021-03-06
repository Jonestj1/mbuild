__author__ = 'jonestj1'

import mbuild as mb

from mbuild.lib.atoms import H
from mbuild.lib.moieties import CH3
from mbuild.lib.moieties import CH2
from mbuild.lib.moieties import Silane

class PegMonomer(mb.Compound):
    def __init__(self):
        super(PegMonomer, self).__init__()

        mb.load('peg_monomer.pdb', compound=self, relative_to_module=self.__module__)
        mb.translate(self, -self[0].pos)

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.translate(self['down'], [0, -0.07, 0])

        self.add(mb.Port(anchor=self[6]), 'up')
        mb.translate(self['up'], [0, 0.364, 0])
