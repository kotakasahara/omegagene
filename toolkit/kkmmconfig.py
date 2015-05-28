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
            #"particle-cluster-shake":ConfigReader.INT,
            
            "fn-i-ttp-v-mcmd-inp":ConfigReader.STR,
            "fn-i-ttp-v-mcmd-initial":ConfigReader.STR,
            "ttp-v-mcmd-initial-vs":ConfigReader.INT,
            "ttp-v-mcmd-seed":ConfigReader.INT,

            "fn-i-atom-groups":ConfigReader.STR,
            "fn-i-dist-restraint":ConfigReader.STR,
            "mol-settle":ConfigReader.ARRAY_STR
            }
        for key, k_type in self.type_def.items():
            if is_array_type(k_type):
                self.config[key] = []
        return

    def get_val(self, key):
        if key in self.config:
            return self.config[key]
        else:
            return None
    def set_val(self, key, val):
        if is_array_type(self.type_def[key]):
            self.config[key].extend(val)
        else:
            self.config[key] = val[0]
            

def is_array_type(in_t):
    if in_t in (ConfigReader.ARRAY_STR,
                ConfigReader.ARRAY_FLOAT,
                ConfigReader.ARRAY_INT):
        return True
    else: return False

class ConfigReader(kkkit.FileI):
    FLOAT = "float"
    INT = "int"
    STR = "str"
    MULTI = "multi"
    ARRAY_STR = "array_str"
    ARRAY_INT = "array_int"
    ARRAY_FLOAT = "array_float"

    DEBUG = True
    def __init__(self, fn):
        super(ConfigReader, self).__init__(fn)
    def is_float_type(self, in_t):
        if in_t in (ConfigReader.FLOAT,
                    ConfigReader.ARRAY_FLOAT):
            return True
        else: return False
    def is_int_type(self, in_t):
        if in_t in (ConfigReader.INT,
                    ConfigReader.ARRAY_INT):
            return True
        else: return False
    def is_str_type(self, in_t):
        if in_t in (ConfigReader.STR,
                    ConfigReader.ARRAY_STR):
            return True
        else: return False
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
            vals = terms[1:]
            if ConfigReader.DEBUG: print key + " : " + ", ".join(vals)
            vals_conv = []
            if not key in config.type_def:
                sys.stderr.write("Invalid key : " + key + "\n")
            for val in vals:
                val_c = val
                val_type = config.type_def[key]
                if self.is_float_type(val_type):
                    val_c = float(val)
                elif self.is_int_type(val_type):
                    val_c = int(val)
                vals_conv.append(val_c)

            config.set_val(key,vals_conv)

            line = self.read_line_com()
        self.close()
        return config


