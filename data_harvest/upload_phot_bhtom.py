import argparse
import requests

URLROOT = 'https://bh-tom2.astrolabs.pl/api/token-auth/'

def upload_single_image_photometry():
    """
    Function to upload the photometry file for a single image to BHTOM
    """

    ur = {'target': target_pk, 'data_product_type': 'photometry', 'groups': target_groups}
    file_data = {'file': (payload['file_path'], open(payload['file_path'],'rb'))}
    dataupload_url = concat_urls(config['tom_url'],config['dataproducts_endpoint'])
    response = requests.post(dataupload_url, data=ur, files=file_data, auth=login)

    if log!= None:
        log.info('Uploaded lightcurve file to TOM at URL: '+repr(response.url))
        log.info('with response: '+repr(response.text))
