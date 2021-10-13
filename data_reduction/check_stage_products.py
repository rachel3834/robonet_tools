from os import path
from sys import argv
import glob

def review_field_reduced_datasets(top_dir, field_id):

    datasets = glob.glob(path.join(top_dir, field_id+'*'))

    log = open(path.join(top_dir,'logs','data_products_report.txt'),'w')

    kernel_stamp_data = {}
    for red_dir in datasets:
        count_reduction_products(red_dir, log)
        kernel_stamp_data = count_kernel_stamps(red_dir, kernel_stamp_data)
        stamp_data = kernel_stamp_data[path.basename(red_dir)]
        log.write('Percentage of images with expected number of kernel stamps: '+str(stamp_data['%OK'])+'\n')

    log.write('\n\n')

    for red_dir in datasets:
        find_missing_frames(red_dir, log)

    log.write('\n\n')

    record_kernel_stamp_data(kernel_stamp_data, log)

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

def count_kernel_stamps(dir_path, kernel_stamp_data):
    image_list = glob.glob(path.join(dir_path,'kernel', '*.fits'))

    stamps_data = {}
    NOK = 0
    for image in image_list:
        kernels_npy = glob.glob(path.join(image,'kernel_stamp_*.npy'))
        kernels_fits = glob.glob(path.join(image,'kernel_stamp_*.fits'))
        ukernels_fits = glob.glob(path.join(image,'kernel_err_stamp_*.fits'))

        if len(kernels_npy) == 16 and len(kernels_fits) == 16 and \
                len(ukernels_fits) == 16:
            status = 'OK'
            NOK += 1
        else:
            status = 'Stamps missing'

        stamps_data[path.basename(image)] = {'Nkernels_npy': len(kernels_npy),
                                             'Nkernel_fits': len(kernels_fits),
                                             'Nkernel_err_fits': len(ukernels_fits),
                                             'status': status}
    if len(image_list) > 0:
        stamps_data['%OK'] = round( (float(NOK) / float(len(image_list)))*100.0, 1)
    else:
        stamps_data['%OK'] = 0.0

    kernel_stamp_data[path.basename(dir_path)] = stamps_data

    return kernel_stamp_data

def record_kernel_stamp_data(kernel_stamp_data, log):

    for dataset,stamps_data in kernel_stamp_data.items():
        log.write('Kernel stamps listing for '+dataset+'\n')

        for key, data in stamps_data.items():
            if key == '%OK':
                log.write('Percentage of images with expected number of kernel stamps: '+str(data)+'\n')
            else:
                log.write(key+': Nkernels_npy='+str(data['Nkernels_npy'])+\
                                        ' Nkernel_fits='+str(data['Nkernel_fits'])+\
                                        ' Nkernel_err_fits='+str(data['Nkernel_err_fits'])+\
                                        ' '+data['status']+'\n')

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
