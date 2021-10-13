from os import path
from sys import argv
import glob

def review_field_reduced_datasets(top_dir, field_id):

    datasets = glob.glob(path.join(top_dir, field_id+'*'))

    log = open(path.join(top_dir,'logs','data_products_report.txt'),'w')

    for red_dir in datasets:
        report = check_reduction_products(red_dir)
        log.write(report+'\n')
        print(report)

    log.close()

def check_reduction_products(red_dir):

    n_input_images = count_dir_fits_products(path.join(red_dir,'data'))
    n_resampled = count_dir_fits_products(path.join(red_dir,'resampled'))
    p_resampled = round( (n_resampled/n_input_images)*100.0, 0 )
    n_kernel = count_dir_fits_products(path.join(red_dir,'kernel'))
    p_kernel = round( (n_kernel/n_input_images)*100.0, 0 )
    n_diffim = count_dir_fits_products(path.join(red_dir,'diffim'))
    p_diffim = round( (n_diffim/n_input_images)*100.0, 0 )

    report = path.basename(red_dir)+': Ninput_images='+str(n_input_images)+\
                                 ' Nresampled='+str(n_resampled)+'('+str(p_resampled)+'%) '+\
                                 ' Nkernel='+str(n_kernel)+'('+str(p_kernel)+'%) '+\
                                 ' Ndiffim='+str(n_diffim)+'('+str(p_diffim)+'%)'
    return report

def count_dir_fits_products(dir_path):
    file_list = glob.glob(path.join(dir_path,'*.fits'))
    return float(len(file_list))


if __name__ == '__main__':
    if len(argv) == 1:
        top_dir = input('Please enter the path to the top level directory: ')
        field_id = input('Please enter the name of the field: ')
    else:
        top_dir = argv[1]
        field_id = argv[2]

    review_field_reduced_datasets(top_dir, field_id)
