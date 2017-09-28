#!/usr/bin/python2.7

import kkkit
import re

class AtomGroupsReader(kkkit.FileI):
    def __init__(self, fn):
        super(AtomGroupsReader, self).__init__(fn)
        self.groups = [[]]
        self.names = ["all"]
    def parse_unit(self, term):
        at = re.compile("^(\d+)$")
        at_range = re.compile("^(\d+)\-(\d+)$")
        m1 = at.match(term)
        m2 = at_range.match(term)
        val = []
        if m1:
            val = [int(m1.group())]
        elif m2:
            val = [x for x in range(int(m2.group(1)), int(m2.group(2))+1)]
        else:
            msg = "[ " + term + " ]"
            kkkit.error_msg("011_002",msg)
        return val
    def read_groups(self):
        self.open()
        for orig_line in self.f:
            line = kkkit.eliminate_comment(orig_line).strip()
            line2 = line.replace(","," ")
            
            terms = line2.split()

            if len(line) == 0: continue
            g_name = terms[0]
            grp = set()
            for x in terms[1:]:
                for x_val in self.parse_unit(x):
                    grp.add(x_val)
            
            #self.groups[len(self.groups)+1] = (g_name, sorted(list(grp)))
            #self.groups[g_name] = sorted(list(grp))
            self.groups.append(sorted(list(grp)))
            self.names.append(g_name)
            #self.groups[terms[0]] = [int(x) for x in terms[1:]]
            print "Group-ID: " + str(len(self.groups)-1) + " n_atoms: " + str(len(self.groups[-1])) + " " + self.names[-1]
        self.close()
        return self.groups, self.names
