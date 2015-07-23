__author__ = 'jonestj1'

import mbuild as mb

from mbuild.components.atoms.H import H
from mbuild.components.small_groups.ch3 import CH3
from mbuild.components.small_groups.silane import Silane
from mbuild.components.surfaces.betacristobalite import Betacristobalite


class PegMonomer(mb.Compound):
    def __init__(self):
        super(PegMonomer, self).__init__()

        mb.load('peg_monomer.pdb', compound=self, relative_to_module=self.__module__)
        mb.translate(self, -self.C[0])

        self.add(mb.Port(anchor=self.C[0]), 'down')
        mb.translate(self.down, [0, -0.07, 0])


        self.add(mb.Port(anchor=self.O[0]), 'up')
        mb.translate(self.up, [0, 0.364, 0])


class PegPolymer(mb.Compound):
    def __init__(self, n=3, cap_front=True, cap_end=True):

        if n < 1:
            raise Exception('n must be 1 or more')
        super(PegPolymer, self).__init__()


        if not cap_front:
            n += 1
        if not cap_end:
            n += 1

        chain = mb.Polymer(PegMonomer(), n=n-1, port_labels=('down', 'up'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(CH3(), 'methyl_front')
            mb.equivalence_transform(
                self.chain, self.chain.up, self.methyl_front.up)
        else:
            self.add(chain.up, 'up', containment=False)

        if cap_end:
            self.add(CH3(), 'methyl_end')
            mb.equivalence_transform(
                self.methyl_end, self.methyl_end.up, self.chain.down)
        else:
            self.add(chain.down, 'down', containment=False)


class PegPolymerSilane(mb.Compound):

    def __init__(self, chain_length):
        super(PegPolymerSilane, self).__init__()

        peg_polymer = PegPolymer(chain_length, cap_end=False)
        self.add(peg_polymer, 'peg_polymer')
        silane = Silane()
        self.add(silane, 'silane')
        mb.equivalence_transform(
            self.peg_polymer, self.peg_polymer.down, self.silane.up)

        self.add(silane.down, 'down', containment=False)


class PegMonolayer(mb.Compound):

    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=5):
        super(PegMonolayer, self).__init__()

        surface = Betacristobalite(area_per_port=.25)

        tc = mb.TiledCompound(surface, n_tiles=(tile_x, tile_y, 1))
        self.add(tc, 'tiled_surface')
        peg_polymer_silane = PegPolymerSilane(chain_length)
        hydrogen = H()

        peg_polymer_silanes, hydrogens = mb.apply_mask(host=self.tiled_surface,
                                                       guest=peg_polymer_silane,
                                                       mask=mask, backfill=hydrogen)

        self.add(peg_polymer_silanes)
        self.add(hydrogens)