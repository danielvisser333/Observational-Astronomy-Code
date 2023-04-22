import os
import numpy as np
import astropy.io.fits as fits
import ccdproc


# Print a 2D aray with header
def table(data_in, header):
    for header_element in header:
        print(header_element + " ", end='')
    print()
    for column in data_in:
        for element in column:
            print(str(element) + " ", end="")
        print()


def filter_header_data(name, data, mean):
    return [name, data["FILTER"], data["EXPTIME"], mean]


# Add the PEDESTAL from the header to the data
def add_pedestal(f_in, f_out):
    data, header = fits.getdata(f_in, header=True)
    pedestal = 1.0 * header['PEDESTAL']
    data_new = data + pedestal
    header['PEDESTAL'] = 0
    fits.writeto(f_out, data_new, overwrite=True)


def inverse_average(data):
    return 1 / data.mean()


def create_masterflat(files, f_out):
    print("Creating masterflat " + f_out)
    temp_files = [os.path.basename(filename) for filename in files]
    for i in range(len(files)):
        add_pedestal(files[i], temp_files[i])
    masterflat = ccdproc.combine(temp_files, method='average', unit='adu',
                                 scale=inverse_average)
    masterflat.write(f_out, overwrite=True)
    return masterflat


# Start
print("Observational Astronomy Data")

flat_g_file_names = ["./data/flat_g/pic%d.fit" % i for i in range(1, 11)]
flat_r_file_names = ["./data/flat_r/pic%d.fit" % i for i in range(1, 11)]
sloan_g_file_names = ["./data/sloan_g/pic%d.fit" % i for i in range(1, 11)]
sloan_r_file_names = ["./data/sloan_r/pic%d.fit" % i for i in range(1, 11)]

if not os.path.isdir('./processing'):
    os.makedirs('./processing')

# Print the data from the FITS files
flat_g_image_headers = []
flat_r_image_headers = []
sloan_g_image_headers = []
sloan_r_image_headers = []

for i in range(10):
    flat_g_data, flat_g_header = fits.getdata(flat_g_file_names[i], header=True)
    flat_r_data, flat_r_header = fits.getdata(flat_r_file_names[i], header=True)
    sloan_g_data, sloan_g_header = fits.getdata(sloan_g_file_names[i], header=True)
    sloan_r_data, sloan_r_header = fits.getdata(sloan_r_file_names[i], header=True)
    flat_g_image_headers.append(filter_header_data(flat_g_file_names[i], flat_g_header, np.mean(flat_g_data)))
    flat_r_image_headers.append(filter_header_data(flat_r_file_names[i], flat_r_header, np.mean(flat_r_data)))
    sloan_g_image_headers.append(filter_header_data(sloan_g_file_names[i], sloan_g_header, np.mean(sloan_g_data)))
    sloan_r_image_headers.append(filter_header_data(sloan_r_file_names[i], sloan_r_header, np.mean(sloan_g_data)))

table(flat_g_image_headers, ["FLAT NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(flat_r_image_headers, ["FLAT NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(sloan_g_image_headers, ["NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(sloan_r_image_headers, ["NAME", "FILTER", "EXPOSURE", "AVERAGE"])

# Create the master-flats
flat_g_master = create_masterflat(flat_g_file_names, './processing/flat_g_master.fit')
flat_r_master = create_masterflat(flat_r_file_names, './processing/flat_r_master.fit')
sloan_g_master = create_masterflat(sloan_g_file_names, './processing/sloan_g_master.fit')
sloan_r_master = create_masterflat(sloan_r_file_names, './processing/sloan_r_master.fit')
