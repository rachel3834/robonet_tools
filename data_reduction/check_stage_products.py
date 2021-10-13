from os import path
from sys import argv
import glob

def review_field_reduced_datasets(top_dir, field_id):

    datasets = glob.glob(path.join(top_dir, field_id+'*'))

    log = open(path.join(top_dir,'logs','data_products_report.txt'),'w')

    for red_dir in datasets:
        count_reduction_products(red_dir, log)
        count_kernel_stamps(red_dir, log)

    log.write('\n\n')

    for red_dir in datasets:
        find_missing_frames(red_dir, log)

    log.close()

def count_reduction_products(red_dir, log):

    n_input_images = count_dir_fits_products(path.join(red_dir,'data'))
    n_resampled = count_dir_fits_products(path.join(red_dir,'resampled'))
    p_resampled = round( (n_resampled/n_input_images)*100.0, 1 )
    n_kernel = count_dir_fits_products(path.join(red_dir,'kernel'))
    p_kernel = round( (n_kernel/n_input_images)*100.0, 1 )
    n_diffim = count_dir_fits_products(path.join(red_dir,'diffim'))
    p_diffim = round( (n_diffim/n_input_images)*100.0, 1 )

    report = path.basename(red_dir)+': Ninput_images='+str(n_input_images)+\
                                 ' Nresampled='+str(n_resampled)+'('+str(p_resampled)+'%) '+\
                                 ' Nkernel='+str(n_kernel)+'('+str(p_kernel)+'%) '+\
                                 ' Ndiffim='+str(n_diffim)+'('+str(p_diffim)+'%)'
    log.write(report+'\n')
    print(report)

def find_missing_frames(red_dir, log):

    log.write('============================================================\n')
    log.write('Missing frame report for '+path.basename(red_dir)+'\n')

    input_images = list_dir_fits_products(path.join(red_dir,'data'))
    resampled_images = list_dir_fits_products(path.join(red_dir,'resampled'))
    kernel_images = list_dir_fits_products(path.join(red_dir,'kernel'))
    diff_images = list_dir_fits_products(path.join(red_dir,'diffim'))

    log.write('Images missing from resampled data products: ')
    missing = diff_image_lists(input_images, resampled_images)
    record_missing_images(missing, log)

    log.write('Images missing from kernel data products: ')
    missing = diff_image_lists(input_images, kernel_images)
    record_missing_images(missing, log)

    log.write('Images missing from diffim data products: ')
    missing = diff_image_lists(input_images, diff_images)
    record_missing_images(missing, log)

def count_dir_fits_products(dir_path):
    file_list = glob.glob(path.join(dir_path,'*.fits'))
    return float(len(file_list))

def count_kernel_stamps(dir_path, log):
    image_list = glob.glob(path.join(dir_path,'kernel', '*.fits'))

    for image in image_list:
        kernels_npy = glob.glob(path.join(image,'kernel_stamp_?.npy'))
        kernels_fits = glob.glob(path.join(image,'kernel_stamp_?.fits'))
        ukernels_fits = glob.glob(path.join(image,'kernel_err_stamp_?.fits'))

        log.write(path.basename(image)+': Nkernels_npy='+str(len(kernels_npy))+\
                                        ' Nkernel_fits='+str(len(kernels_fits))+\
                                        ' Nkernel_err_fits='+str(len(ukernels_fits))+'\n')

def list_dir_fits_products(dir_path):
    file_list = glob.glob(path.join(dir_path,'*.fits'))
    image_list = [path.basename(image_path) for image_path in file_list]
    image_list.sort()
    return image_list

def diff_image_lists(data_list, red_data_list):
    missing_entries = set(data_list) - set(red_data_list)
    return missing_entries

def record_missing_images(missing_entries, log):
    for entry in missing_entries:
        log.write(entry+'\n')
    log.write('\n')

if __name__ == '__main__':
    if len(argv) == 1:
        top_dir = input('Please enter the path to the top level directory: ')
        field_id = input('Please enter the name of the field: ')
    else:
        top_dir = argv[1]
        field_id = argv[2]

    review_field_reduced_datasets(top_dir, field_id)
