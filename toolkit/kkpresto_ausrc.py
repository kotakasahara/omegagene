#!/usr/bin/python

import re
import kkpdb
import kkpresto
from optparse import OptionParser

class PrestoAtomSelector(object):
    def __init__(self, molid, res_id,
                 res_name, atom_name):
        super(PrestoAtomSelector, self).__init__()
        self.molid = molid
        self.res_id = res_id
        self.res_name = res_name
        self.atom_name = atom_name
        self.atom_id = -1
        return 
    def get_text(self):
        line = "%4d%5d %-4s%4s"%(self.molid,
                                 self.res_id,
                                 self.res_name,
                                 self.atom_name)
        return line
    def set_atom_id(self, tpl):
        ## tpl is an instance of kkpresto.TPL
        mol_type_id = 0
        mol_id_in_type = 0
        i = 0
        for mol in tpl.mols:
            for i_mol in range(mol.mol_num):
                if i == self.molid: break
                mol_id_in_type += 1
                i += 1
            if i == self.molid: break
            mol_type_id += 1
        ## tpl.mols[mol_type_id]
        ## mol_id_in_type -th instance of the mol

        mol = tpl.mols[mol_type_id]
        for i, atom in enumerate(mol.atoms):
            if atom.res_id == self.res_id and \
                    atom.res_name == self.res_name and \
                    atom.atom_name == self.atom_name:
                atom_id = mol.head_atom_id + mol_id_in_type * len(mol.atoms) + i
                self.atom_id = atom_id
                return atom_id
        return -1

class PrestoAUSReactCoord(object):
    def __init__(self):
        return

class PrestoAUSReactCoordDist(PrestoAUSReactCoord):
    def __init__(self):
        super(PrestoAUSReactCoordDist, self).__init__()

        ## atom_groups[0] ... list of PrestoAtomSelector
        ## atom_groups[1] ... list of PrestoAtomSelector
        self.atom_groups = ([], [])
        self.atom_ids = ([], [])
        self.grp_names = ("aus_group_a", "aus_group_b")
        return 
    def add_atom(self, group_id,
                 molid, res_id, res_name, atom_name):
        assert(group_id < 2)
        self.atom_groups[group_id].append(
            PrestoAtomSelector(molid, res_id, res_name, atom_name)
            )
    #def get_text(self):
    #    line = "%4d%5d %-4s%4s%6d%5d %-4s%4s"%(self.molid[0],
    #                                           self.res_id[0],
    #                                           self.res_name[0],
    #                                           self.atom_name[0],
    #                                           self.molid[1],
    #                                           self.res_id[1],
    #                                           self.res_name[1],
    #                                           self.atom_name[1])
    #    return line
    def set_atom_ids(self, tpl):
        for grp_id, grp in enumerate(self.atom_groups):
            for at in grp:
                at_id = at.set_atom_id(tpl)
                if at_id < 0:
                    print "Atom not found:"
                    print at.get_text()
                else:
                    self.atom_ids[grp_id].append(at_id)
        return
    def get_group_text(self):
        text = ""
        for grp_id, atoms in enumerate(self.atom_ids):
            text += self.grp_names[grp_id] + " " 
            text += " ".join([str(x+1) for x in atoms])
            text += "\n"
        return text
    
class PrestoAUSReactCoordReader(kkpresto.PrestoAsciiReader):
    def __init__(self, fn):
        super(PrestoAUSReactCoordReader, self).__init__(fn)
    def read(self):
        rcd = PrestoAUSReactCoordDist()
        self.open()
        read_mode = 0
        aus_type = 0
        add_group = 0
        line = "DUM"
        while 1:
            line = self.readline_comment()
            if not line : break
            terms = line.strip().split()
            if not line: break
            #print line
            if terms[0] == "AUSAMP>":
                if terms[1] == "TYPE":
                    read_mode = 1                
                elif terms[1] == "LIST":
                    read_mode = 2
                elif terms[1] == "STOP":
                    read_mode = 3
                else:
                    sys.stderr.write("Invalid AUSAMP> (it should be TYPE, LIST, or STOP): \n")
                    sys.stderr.write(line + "\n")
            elif read_mode == 1:
                ##"AUSTYPE 2"
                if terms[0] == "AUSTYPE":
                    austype = int(terms[1])
                    if not austype in (1,2):
                        sys.stderr.write("Invalid AUSTYPE (1 or 2): \n")
                        sys.stderr.write(line + "\n")
                        sys.exit(1)
                else:
                    sys.stderr.write("Invalid file format. AUSTYPE should be specified: \n")
                    sys.stderr.write(line + "\n")
                    sys.exit(1)
            elif read_mode == 2:
                if terms[0] == "GROUP":
                    if terms[1] == "A":
                        add_group = 0
                    elif terms[1] == "B":
                        add_group = 1
                    else:
                        sys.stderr.write("Invalid GROUP (it should be A or B): \n")
                        sys.stderr.write(line + "\n")
                        sys.exit(1)
                elif terms[0] == "END":
                    read_mode = 3
                else:
                    molid = int(terms[0])
                    res_id = int(terms[1])
                    res_name = terms[2]
                    atom_name = terms[3]
                    rcd.add_atom(add_group,
                                 molid, res_id, res_name, atom_name)
            elif read_mode == 3:
                pass
        self.close()
        return rcd

def get_options():
    p = OptionParser()
    p.add_option('-i', '--aus', dest='fn_i',
                 help="File name for input AUS group definition file")
    p.add_option('-p', '--pdb', dest='fn_pdb',
                 help="File name for pdb file")
    p.add_option('-o', '--out', dest='fn_out',
                 help="File name for output, atom group definition file")
    p.add_option('-t', '--tpl', dest='fn_tpl',
                 help="File name for topology file")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

def _main():
    opts, args = get_options()
    reader = PrestoAUSReactCoordReader(opts.fn_i)
    tpl = kkpresto.TPLReader(opts.fn_tpl).read_tpl()
    reactcoord = reader.read()
    reactcoord.set_atom_ids(tpl)
    
    fo = open(opts.fn_out, "w")
    fo.write(reactcoord.get_group_text())
    fo.close()

if __name__ == "__main__":
    _main()
