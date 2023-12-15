import argparse
from os import path
import json
from pyDANDIA import crossmatch
from pyDANDIA import metadata
import recombine_image_stamps
import h5py
import numpy as np

PRI_REF = 'lsc-doma-1m0-05-fa15_ip'

def run_extraction(args):
    """Driver function to extract image sections around a set of microlensing targets"""

    # Load information on known microlensing events within the field
    known_events = load_event_list(args)

    # Select high-magnification events with calibrated magnitude detections in the SDSS-i band
    # primary reference dataset
    selected_events = select_events(known_events)

    # Load the crossmatch table for this field.  This provides information on the star identifiers in the
    # primary reference dataset, as well as the available imagesets
    xmatch = crossmatch.CrossMatchTable()
    xmatch.load(args.crossmatch_file, log=None)

    # Load the pyDANDIA metadata file for the primary reference dataset for this field, and use it to
    # annotate the selected microlensing events with the (x,y) coordinates source in the field of view
    meta = metadata.MetaData()
    pri_ref_dir = path.join(args.red_dir, args.field_name + '_' + PRI_REF)
    meta.load_all_metadata(metadata_directory=pri_ref_dir,
                           metadata_name='pyDANDIA_metadata.fits')
    selected_events = fetch_xy_for_events(args, selected_events, xmatch, meta, pri_ref_dir)

    # Calculate the boundaries of thumbnail stamps around all microlensing events, weeding out from the
    # event set those where the thumbnail boundaries would extend over the image boundaries.  This ensures
    # that the thumbnail images are all the same size, which makes life easier later on
    selected_events = calc_thumbnail_boundaries(selected_events, meta)

    # Extract thumbnail images around all selected events for those images within t0 +/- 0.5tE
    # Record which images in the lightcurve these are for later reference
    selected_events = extract_thumbnail_images(args, selected_events, xmatch)

    # Output a summary of selected events
    output_catalog(args, selected_events)

def output_catalog(args, selected_events):
    """
    Function to output the JSON summary of the selected events
    """
    output_file = path.join(args.output_dir, 'selected_events.json')
    print(selected_events)

    json_data = json.dumps(selected_events, indent=4)

    with open(output_file, 'w') as write_file:
        write_file.write(json_data)
        write_file.close()

def output_thumbnails(args, event_name, thumbnails):
    """
    Function to output the sets of thumbnail images for each event as an HDF5 file
    """

    output_path = path.join(args.output_dir, event_name + '_thumbs.hdf5')

    with h5py.File(output_path, "w") as f:
        for image_name, data in thumbnails.items():
            dset = f.create_dataset(image_name,
                                    data.shape,
                                    dtype='float64',
                                    data=data)
        f.close()

def extract_thumbnail_images(args, selected_events, xmatch):
    """Function to work out which images from a field dataset took place during each microlensing event,
    and extract the pixel data for the thumbnails from the corresponding difference images"""

    for event_name, event_data in selected_events.items():
        tmin = event_data['target_data']['t0'] - 0.5*event_data['target_data']['tE']
        tmax = event_data['target_data']['t0'] + 0.5*event_data['target_data']['tE']

        idx1 = np.where(xmatch.images['hjd'] >= tmin)[0]
        idx2 = np.where(xmatch.images['hjd'] <= tmax)[0]
        image_idx = list(set(idx1).intersection(set(idx2)))

        if len(image_idx) > 0:
            print('Extracting ' + str(len(image_idx)) + ' thumbnail difference images for ' + event_name)
            output_path = path.join(args.output_dir, event_name + '_thumbs.hdf5')

            with h5py.File(output_path, "w") as f:

                for i in image_idx:
                    red_dir = path.join(args.red_dir, xmatch.images['dataset_code'][i])
                    stamps_dir = path.join(red_dir, 'diffim', xmatch.images['filename'][i])

                    # Stamps may be missing in the case of images that were flagged by the pipeline's
                    # quality control.  If so, we skip the attempt to extract a thumbnail
                    if path.isdir(stamps_dir):
                        # Load the metadata for the reduction corresponding to this image and extract the stamps
                        # table
                        red_meta = metadata.MetaData()
                        red_meta.load_all_metadata(metadata_directory=red_dir,
                                               metadata_name='pyDANDIA_metadata.fits')
                        stamps = recombine_image_stamps.parse_stamps_table(red_meta)

                        # Recombine the stamps for the appropriate differenced image
                        full_image = recombine_image_stamps.stamps_to_fullframe_image(stamps, stamps_dir, 'diff_stamp')

                        # Extract the thumbnail around the target event
                        box_boundaries = event_data['rome_star']['box_boundaries']
                        thumb_image = full_image[
                                    box_boundaries['ymin']:box_boundaries['ymax'],
                                    box_boundaries['xmin']:box_boundaries['xmax']
                        ]

                        # Add this dataset to the HDF5 file
                        dset = f.create_dataset(xmatch.images['filename'][i],
                                                thumb_image.shape,
                                                dtype='float64',
                                                data=thumb_image)

                        print(' -> ' + xmatch.images['filename'][i])

            f.close()

        event_data['thumbnails'] = xmatch.images['filename'][image_idx].data

    return selected_events

def calc_thumbnail_boundaries(selected_events, meta):
    """Function to calculate the boundaries of the thumbnail image around each.  This function filters
    out any events where the thumbnail boundaries would exceed the image boundaries."""

    # The boundaries of the images used in the reduction are stored in the metadata
    # from the primary reference dataset, as well as the pixel scale of the data
    image_boundaries = {
        'xmin': meta.reduction_parameters[1]['IMAGEX1'][0],
        'xmax': meta.reduction_parameters[1]['IMAGEX2'][0],
        'ymin': meta.reduction_parameters[1]['IMAGEY1'][0],
        'ymax': meta.reduction_parameters[1]['IMAGEY2'][0]}
    pixscale = meta.reduction_parameters[1]['PIX_SCALE'][0]

    # The thumbnails are set to a size of 30 arcsec square; npix represents the half-width of the box:
    npix = int(round(15.0 / pixscale,0))

    # Review the dictionary of events and calculate the thumbnail pixel boundaries,
    # eliminating any events that are too close to the edge of the frame to have a full
    # thumbnail in order to keep the thumbnail sizes consistent
    revised_events = {}
    for event_name, event_data in selected_events.items():
        x_center = int(round(event_data['rome_star']['x'],0))
        y_center = int(round(event_data['rome_star']['y'],0))
        box_boundaries = {
            'x_center': x_center,
            'y_center': y_center,
            'xmin': x_center - npix,
            'xmax': x_center + npix,
            'ymin': y_center - npix,
            'ymax': y_center + npix
        }
        if check_image_boundaries(box_boundaries, image_boundaries ):
            event_data['rome_star']['box_boundaries'] = box_boundaries
            revised_events[event_name] = event_data

    print('Selected ' + str(len(revised_events)) + ' events with thumbnail boxes within the image FoV')

    return revised_events

def check_image_boundaries(box_boundaries, image_boundaries ):
    """Function to verify whether or not a box lies fully within the image field of view"""

    result = (box_boundaries['xmin'] >= image_boundaries['xmin']) \
            and (box_boundaries['xmin'] < image_boundaries['xmax']) \
            and (box_boundaries['xmax'] > image_boundaries['xmin']) \
            and (box_boundaries['xmax'] <= image_boundaries['xmax']) \
            and (box_boundaries['ymin'] >= image_boundaries['ymin']) \
            and (box_boundaries['ymin'] < image_boundaries['ymax']) \
            and (box_boundaries['ymax'] > image_boundaries['ymin']) \
            and (box_boundaries['ymax'] <= image_boundaries['ymax'])

    return result

def fetch_xy_for_events(args, selected_events, xmatch, meta, pri_ref_dir):
    """Function to identify the pixel coordinates of all selected events within the primary reference
    for the given field"""

    # Primary reference dataset identifier
    pri_ref = path.basename(pri_ref_dir)

    for event_name, event_data in selected_events.items():
        field_id = event_data['rome_star']['field_id']
        field_idx = field_id - 1

        # Fetch the index for this object in the field primary reference from the crossmatch field index
        # Note this should almost always be the same as the field_id because of the selection cuts
        # we're doing here, but we do this step to be sure
        pri_ref_id = xmatch.field_index[pri_ref + '_index'][field_idx]
        pri_ref_idx = pri_ref_id - 1
        event_data['rome_star'][pri_ref + '_index'] = pri_ref_id

        # Fetch the pixel coordinates of the star from the metadata and append the information
        # to the selected_events dictionary
        event_data['rome_star']['x'] = meta.star_catalog[1]['x'][pri_ref_idx]
        event_data['rome_star']['y'] = meta.star_catalog[1]['y'][pri_ref_idx]

        selected_events[event_name] = event_data

    print('Extracted the pixel coordinates of the selected events')

    return selected_events

def select_events(known_events):
    """Function to select microlensing events that meet the criteria:
    u0 <=0.2
    cal_i_mag_lsc_doma != None
    """

    selected_events = {}
    for event_name, event_data in known_events.items():

        # MOA events provide amax rather than u0 like all other surveys, so
        # perform an approximate conversion
        if 'u0' in event_data['target_data'].keys():
            u0 = event_data['target_data']['u0']
        elif 'amax' in event_data['target_data'].keys():
            u0 = 1.0 / event_data['target_data']['amax']
        else:
            u0 = None

        if u0:
            if u0 <= 0.2:

                # Events may have multiple matches in the ROME catalog, and the input file
                # includes all hits within 2 arcsec of the target.  We assume the nearest star
                # is the true target
                min_sep = 30.0
                rome_star = None
                for star in event_data['rome_stars']:
                    if star['separation_deg'] <= min_sep:
                        min_sep = star['separation_deg']
                        rome_star = star

                if rome_star and rome_star['cal_i_mag_lsc_doma'] > 0.0 \
                        and rome_star['cal_i_mag_lsc_doma'] <= 21.0:
                    selected_events[event_name] = {
                                                    'target_data': event_data['target_data'],
                                                    'rome_star': rome_star
                    }

    print('Of the ' + str(len(known_events)) + ' events, '
        + str(len(selected_events)) + ' met the initial selection critiera')

    return selected_events

def load_event_list(args):
    """Function to load a list of known events within a given field in JSON format"""

    if not path.isfile(args.events_file):
        raise IOError('Cannot find input list of targets at '
                        + args.events_file)

    with open(args.events_file, "r") as read_file:
        data = json.load(read_file)

    return data

def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('events_file', help='Path to JSON file of known events from the current field')
    parser.add_argument('crossmatch_file', help='Path to the crossmatch file for the current field')
    parser.add_argument('red_dir', help='Path to the reduction directory of the primary reference dataset')
    parser.add_argument('field_name', help='Name of the current field')
    parser.add_argument('output_dir', help='Path to the output directory')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    run_extraction(args)