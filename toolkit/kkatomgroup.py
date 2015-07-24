#!/usr/bin/python2.7

import kkkit

class AtomGroupsReader(kkkit.FileI):
    def __init__(self, fn):
        super(AtomGroupsReader, self).__init__(fn)
        self.groups = {}
    def read_groups(self):
        self.open()
        for orig_line in self.f:
            line = kkkit.eliminate_comment(orig_line).strip()
            terms = line.split()
            if len(line) == 0: continue
            self.groups[terms[0]] = [int(x) for x in terms[1:]]
        self.close()
        return self.groups
