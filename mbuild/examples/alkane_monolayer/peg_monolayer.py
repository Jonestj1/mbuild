__author__ = 'jonestj1'

import mbuild as mb

from mbuild.lib.atoms import H
from mbuild.lib.moieties import CH3
from mbuild.lib.moieties import CH2
from mbuild.lib.moieties import Silane
from mbuild.lib.surfaces import Betacristobalite
from mbuild.lib.bulk_materials import AmorphousSilica

class PegMonomer(mb.Compound):
    def __init__(self):
        super(PegMonomer, self).__init__()

        mb.load('peg_monomer.pdb', compound=self, relative_to_module=self.__module__)
        mb.translate(self, -self[0].pos)

        self.add(mb.Port(anchor=self[0]), 'down')
        mb.translate(self['down'], [0, -0.07, 0])


        self.add(mb.Port(anchor=self[6]), 'up')
        mb.translate(self['up'], [0, 0.364, 0])


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
                self['chain'], self['chain']['up'], self['methyl_front']['up'])
        else:
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(CH3(), 'methyl_end')
            mb.equivalence_transform(
                self['methyl_end'], self['methyl_end']['up'], self['chain']['down'])
        else:
            self.add(chain['down'], 'down', containment=False)


class PegPolymerSilane(mb.Compound):

    def __init__(self, chain_length):
        super(PegPolymerSilane, self).__init__()

        peg_polymer = PegPolymer(int((chain_length-1)/3), cap_end=False)
        self.add(peg_polymer, 'peg_polymer')
        silane = Silane()
        self.add(silane, 'silane')
        mb.equivalence_transform(
            self['peg_polymer'], self['peg_polymer']['down'], self['silane']['up'])

        self.add(silane['down'], 'down', containment=False)


class PegMonolayer(mb.Monolayer):

    def __init__(self, pattern, surface=Betacristobalite(), tile_x=1, tile_y=1, chain_length=10):
        # bulk denotes surface structure
        # if bulk is 'crystalline':
        #     surface = Betacristobalite()
        # elif bulk is 'amorphous':
        #     surface = AmorphousSilica()
        peg_polymer_silane = PegPolymerSilane(chain_length)
        hydrogen = H()

        super(PegMonolayer, self).__init__(surface, peg_polymer_silane, backfill=hydrogen, 
                                           pattern=pattern, tile_x=tile_x, 
                                           tile_y=tile_y)


class PegPolymerOxyCap(mb.Compound):
    def __init__(self, n=3, cap_front=True, cap_end=True):

        if n < 1:
            raise Exception('n must be 1 or more')
        super(PegPolymerOxyCap, self).__init__()


        if not cap_front:
            n += 1
        if not cap_end:
            n += 1

        chain = mb.Polymer(PegMonomer(), n=n-1, port_labels=('down', 'up'))
        self.add(chain, 'chain')

        if cap_front:
            self.add(H(), 'oxy_front')
            mb.equivalence_transform(
                self['chain'], self['chain']['up'], self['oxy_front']['up'])
        else:
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(CH3(), 'methyl_end')
            mb.equivalence_transform(
                self['methyl_end'], self['methyl_end']['up'], self['chain']['down'])
        else:
            self.add(chain['down'], 'down', containment=False)

class PegPolymerSilaneOxyCap(mb.Compound):

    def __init__(self, chain_length):
        super(PegPolymerSilaneOxyCap, self).__init__()

        peg_polymer = PegPolymerOxyCap(int((chain_length-1)/3), cap_end=False)
        self.add(peg_polymer, 'peg_polymer')
        
        carbon = CH2()
        self.add(carbon, 'ch2')
        mb.equivalence_transform(
            self['peg_polymer'], self['peg_polymer']['down'], self['ch2']['down'])

        silane = Silane() 
        self.add(silane, 'silane')
        mb.equivalence_transform(
            self['silane'], self['silane']['up'], self['ch2']['up'])

        self.add(silane['down'], 'down', containment=False)

class PegOxyCapMonolayer(mb.Monolayer):

    def __init__(self, pattern, surface=Betacristobalite(), tile_x=1, tile_y=1, chain_length=10):
        # bulk denotes surface structure
        # if bulk is 'crystalline':
        #     surface = Betacristobalite()
        # elif bulk is 'amorphous':
        #     surface = AmorphousSilica()
        peg_polymer_silane_oxy_cap = PegPolymerSilaneOxyCap(chain_length)
        hydrogen = H()

        super(PegOxyCapMonolayer, self).__init__(surface, peg_polymer_silane_oxy_cap, 
                                                 backfill=hydrogen, pattern=pattern, 
                                                 tile_x=tile_x, tile_y=tile_y)
