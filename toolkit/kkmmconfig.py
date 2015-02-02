import kkkit
import re
import sys

class Config(object):
    def __init__(self):
        self.config = {}
        self.type_def = {
            "fn-i-tpl":ConfigReader.STR,
            "fn-i-initial-pdb":ConfigReader.STR,
            "fn-i-restart":ConfigReader.STR,
            "fn-i-shake":ConfigReader.STR,
            "cell-x":ConfigReader.FLOAT,
            "cell-y":ConfigReader.FLOAT,
            "cell-z":ConfigReader.FLOAT,
            "cell-center-x":ConfigReader.FLOAT,
            "cell-center-y":ConfigReader.FLOAT,
            "cell-center-z":ConfigReader.FLOAT,
            "cell-origin-x":ConfigReader.FLOAT,
            "cell-origin-y":ConfigReader.FLOAT,
            "cell-origin-z":ConfigReader.FLOAT,

            "integrator":ConfigReader.STR,

            "cutoff":ConfigReader.FLOAT,
            "n-steps":ConfigReader.INT,
            "dt":ConfigReader.FLOAT,

            "electrostatic":ConfigReader.STR,
            "ele-ewaldalpha":ConfigReader.FLOAT,

            "print-interval-coord":ConfigReader.INT,
            "print-interval-velo":ConfigReader.INT,
            "print-interval-log":ConfigReader.INT,
            "print-interval-energy":ConfigReader.INT,
            "print-interval-energyflow":ConfigReader.INT,
            "fn-o-coord":ConfigReader.STR,
            "fn-o-log":ConfigReader.STR,
            "fn-o-energy":ConfigReader.STR,
            "fn-o-energyflow":ConfigReader.STR,

            "thermostat":ConfigReader.STR,
            "temperature":ConfigReader.FLOAT,

            "barostat":ConfigReader.STR,
            "pressure":ConfigReader.FLOAT,

            "center-of-motion":ConfigReader.STR,
            "particle-cluster-shake":ConfigReader.INT,
            }
        return
    def get_val(self, key):
        if key in self.config:
            return self.config[key]
        else:
            return None
    def set_val(self, key, val):
        self.config[key] = val

class ConfigReader(kkkit.FileI):
    FLOAT = "float"
    INT = "int"
    STR = "str"
    DEBUG = True
    def __init__(self, fn):
        super(ConfigReader, self).__init__(fn)

    def read_config(self):
        config = Config()
        self.open()
        line = self.read_line_com()
        while line:
            terms = re.compile("\s+").split(line.strip())
            if len(terms) < 2:
                line = self.read_line_com()
                continue
            if terms[0][0:2] != "--": 
                sys.stderr.write("Invalid description : " + line)
                line = self.read_line_com()
                continue
            key = terms[0][2:]
            val = terms[1]
            if ConfigReader.DEBUG: print key + " : " + val
            if key in config.type_def:
                val_type = config.type_def[key]
                if val_type == ConfigReader.FLOAT:
                    val = float(val)
                elif val_type == ConfigReader.INT:
                    val = int(val)
                config.set_val(key,val)
            else:
                sys.stderr.write("Invalid key : " + key + "\n")
            line = self.read_line_com()
        self.close()
        return config


