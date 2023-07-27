from pyDANDIA import aws_utils
from pyDANDIA import automatic_pipeline
from pyDANDIA import logs
import glob
from os import path

def upload_all_lightcurves():
    config = automatic_pipeline.get_config()

    log = logs.start_pipeline_log(config['log_dir'], 'upload_lcs_aws')

    red_dir_list = glob.glob(path.join(config['data_red_dir'],'*'))
    for red_dir in red_dir_list:
        lc_dir = path.join(red_dir, 'lc')
        if path.isdir(red_dir) and path.isdir(lc_dir):
            lc_list = glob.glob(path.join(lc_dir, '*.dat'))
            if len(lc_list) > 0:
                for lc_file in lc_list:
                    aws_utils.upload_lightcurve_aws(config, lc_file, log=log)

    logs.close_log(log)

if __name__ == '__main__':
    upload_all_lightcurves()
