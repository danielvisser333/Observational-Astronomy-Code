import os
import numpy as np
import astropy.io.fits as fits
import ccdproc
import astroalign as aa
import math


# Print a 2D array with header
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
    temp_files = ["./processing/" + os.path.basename(filename) for filename in files]
    for i in range(len(files)):
        add_pedestal(files[i], temp_files[i])
    masterflat = ccdproc.combine(temp_files, method='average', unit='adu',
                                 scale=inverse_average)
    masterflat.write(f_out, overwrite=True)
    return masterflat


# Normalize the data from the files
def flatfield_normalisation(filenames_in, filenames_out, masterflat):
    masterflat = fits.getdata(masterflat)
    for j in range(len(filenames_in)):
        file_in = filenames_in[j]
        file_out = filenames_out[j]
        print("Processing file: %s -> %s" % (file_in, file_out))
        data, header = fits.getdata(file_in, header=True)

        # - Add back the PEDESTAL to the image data and divide by the masterflat
        pedestal = header["PEDESTAL"]
        datap = data + pedestal
        dataf = datap / masterflat
        header['PEDESTAL'] = 0

        # Save the flat-fielded science image together with the updated header
        fits.writeto(file_out, dataf, header, overwrite=True)


# Compute and save image overlaid onto another one
def align_image(source, target, output):
    targetdata = fits.getdata(target)
    sourcedata = fits.getdata(source)
    t_data, (source_pos_array, target_pos_array) = aa.find_transform(sourcedata, targetdata)
    print(f"Image {source}, translation {t_data.translation}, rotation {t_data.rotation}, scale {t_data.scale}")
    assert (abs(t_data.rotation * 180 / math.pi) < 0.01)  # Ensure the images properly align
    assert (abs(t_data.scale - 1) < 0.0002)
    data_translated, footprint = aa.apply_transform(
        t_data, np.int32(sourcedata), np.int32(targetdata))
    fits.writeto(output, data_translated, overwrite=True)


def combine_images(files_in, file_out):
    aligned_image_paths = ["./processing/aligned%d.fit" % i for i in range(len(files_in))]
    for i in range(len(files_in)):  # Align all images
        align_image(files_in[i], files_in[0], aligned_image_paths[i])
    combined_image = ccdproc.combine(aligned_image_paths, method="average", unit="adu")
    print(f"Average of image {file_out}, {np.array(combined_image).mean()}")
    combined_image.write('./processing/combined_temporary.fit', overwrite=True)
    data, header = fits.getdata('./processing/combined_temporary.fit', 0, header=True)  # Only take the PRIMARY image
    fits.writeto(file_out, data, header, overwrite=True)
    print(f"Created new combined image {file_out}")


# Start
print("Observational Astronomy Data")

# Flat files are the calibration images
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

# Display the data
table(flat_g_image_headers, ["FLAT NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(flat_r_image_headers, ["FLAT NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(sloan_g_image_headers, ["NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(sloan_r_image_headers, ["NAME", "FILTER", "EXPOSURE", "AVERAGE"])

# Create the master-flats
flat_g_master = create_masterflat(flat_g_file_names, './processing/flat_g_master.fit')
flat_r_master = create_masterflat(flat_r_file_names, './processing/flat_r_master.fit')

print(
    f"Averages of masterflats: {np.mean(flat_g_master)}, {np.mean(flat_r_master)}")
print(
    f"Standarddeviations of masterflats: {np.std(flat_g_master)}, {np.std(flat_r_master)}")

# Correct the images (normalise)
flat_g_file_names_out = ["./processing/picC%d.fit" % i for i in range(1, 11)]
flat_r_file_names_out = ["./processing/picC%d.fit" % i for i in range(1, 11)]

flatfield_normalisation(flat_g_file_names, flat_g_file_names_out, './processing/flat_g_master.fit')
flatfield_normalisation(flat_r_file_names, flat_r_file_names_out, './processing/flat_r_master.fit')

# Align and combine the images

combine_images(sloan_g_file_names, "./processing/sloan_g_combined.fit")
combine_images(sloan_r_file_names, "./processing/sloan_r_combined.fit")
