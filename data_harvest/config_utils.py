from os import path
import json

def get_config(CONFIG_FILE):

    if path.isfile(CONFIG_FILE) == False:
        raise IOError("No config file found at given location: "+CONFIG_FILE)

    f = open(CONFIG_FILE,'r')

    config_dict = json.load(f)

    f.close()

    for key, value in config_dict.items():
        if '[' in value and ']' in value:
            if ',' in value:
                entries = value.replace('[','').replace(']','').split(',')
            else:
                entries = [value.replace('[','').replace(']','')]
            l = []
            for e in entries:
                try:
                    l.append(float(e))
                except:
                    l.append(e.strip(' '))
            config_dict[key] = l

    config_dict = parse_booleans(config_dict)

    return config_dict

def parse_booleans(config_dict):

    for key, value in config_dict.items():
        if str(value).lower() in ['true', 'false']:
            if "true" in str(value).lower():
                value = True
            elif "false" in str(value).lower():
                value = False
            config_dict[key] = value

    return config_dict
