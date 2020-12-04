import json
import csv
from os import path
from sys import argv
from astropy.coordinates import SkyCoord
from astropy import units as u

def parse_event_data():

    params = get_args()

    events_per_field = json.loads(open(params['events_per_field_file'],'r').read())
    events = json.loads(open(params['events_file'],'r').read())

    events = assign_field_per_event(events_per_field, events)

    output_csv(events, params)

def output_csv(events, params):

    with open(params['output_file'], 'w') as f:
        event_order = list(events.keys())
        event_order.sort()

        for event in event_order:
            entry = event+', '
            for i,item in enumerate(events[event]):
                if i < len(events[event])-1:
                    entry += str(item)+', '
                else:
                    entry += str(item)+'\n'
            f.write(entry)


def assign_field_per_event(events_per_field, events):

    for field in range(1,21,1):
        for event in events_per_field[str(field)]:
            data = events[event]
            c = SkyCoord(data[0], data[1], unit=(u.deg, u.deg))
            (ra, dec) = str(c.to_string('hmsdms', sep=':')).split()
            data.append(field)
            data.append(ra)
            data.append(dec)
            events[event] = data

    return events

def get_args():
    params = {}
    if len(argv) == 1:
        params['events_per_field_file'] = input('Please enter the path to the JSON file of events per field: ')
        params['events_file'] = input('Please enter the path to the JSON file of events: ')
        params['output_file'] = input('Please enter the path to the CSV output file: ')
    else:
        params['events_per_field_file'] = argv[1]
        params['events_file'] = argv[2]
        params['output_file'] = argv[3]

    return params

if __name__ == '__main__':
    parse_event_data()
