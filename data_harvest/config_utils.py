from os import path
import json

def get_config(CONFIG_FILE):

    if path.isfile(CONFIG_FILE) == False:
        raise IOError("No config file found at given location: "+CONFIG_FILE)

    f = open(CONFIG_FILE,'r')

    config_dict = json.load(f)

    f.close()

    for key, value in config_dict.items():
        if '[' in value and ']' in value and ',' in value:
            entries = value.replace('[','').replace(']','').split(',')
            l = []
            for e in entries:
                try:
                    l.append(float(e))
                except:
                    l.append(e)
            config_dict[key] = l

    return config_dict
