import crossmatch_to_ipactable
import json
from os import path

class CodeArguments():
    def __init__(self):
        self.field_name = 'ROME-FIELD-01'
        self.qid = 1
        self.output_dir = './test_output'

def test_load_starcounts():
    """
    Star count files contain JSON-format information on the number of valid datapoints in the
    lightcurves of each star in a quadrant.  The files come in sets of 4 because it is easiest
    to generate them on a per quadrant basis, but the read function needs to combine this
    information into a single dictionary.
    """
    args = CodeArguments()

    # Generate test data.  The star_counts dictionary contains the expected results.
    nstars = 100
    nimages = 10
    star_counts = {}
    for qid in range(1,5,1):
        quad_counts = {}
        for j in range(0, nstars, 1):
            quad_counts[j+1] = {
                'gp': nimages,
                'rp': nimages,
                'ip': nimages
            }
        star_counts.update(quad_counts)

        file_path = path.join(args.output_dir, args.field_name + '_starcounts_Q' + str(qid) + '.json')

        json_data = json.dumps(quad_counts)
        with open(file_path, 'w') as write_file:
            write_file.write(json_data)
            write_file.close()

    # Call test function
    test_result = crossmatch_to_ipactable.load_starcounts(args)

    assert(star_counts == test_result)