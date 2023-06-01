import os
import numpy as np
import astropy.io.fits as fits
import ccdproc
import astroalign as aa
import math
import astropy.io.ascii as asc
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import cmd_read


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


def match_star_magnitudes(urat_name):
    pcoo_element = None
    for element in pcoo:
        pcoo_element = element
        if element[2] == urat_name:
            star_urat = element
            break
    star_g_data = None
    for star in g_data:
        if np.linalg.norm(np.array(star[1], star[2]) - np.array(star_urat[0], star_urat[1])) < 3:
            star_g_data = star
            star_r_data = r_data[star_g_data[0] - 1]
            break
    if star_g_data is None:
        return []
    return urat_name, star_g_data[5], star_g_data[6], star_r_data[5], star_r_data[6], pcoo_element


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

# # Display the data
table(flat_g_image_headers, ["FLAT NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(flat_r_image_headers, ["FLAT NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(sloan_g_image_headers, ["NAME", "FILTER", "EXPOSURE", "AVERAGE"])
table(sloan_r_image_headers, ["NAME", "FILTER", "EXPOSURE", "AVERAGE"])
#
# # Create the master-flats
flat_g_master = create_masterflat(flat_g_file_names, './processing/flat_g_master.fit')
flat_r_master = create_masterflat(flat_r_file_names, './processing/flat_r_master.fit')

print(
    f"Averages of masterflats: {np.mean(flat_g_master)}, {np.mean(flat_r_master)}")
print(
    f"Standarddeviations of masterflats: {np.std(flat_g_master)}, {np.std(flat_r_master)}")
#
# # Correct the images (normalise)
flat_g_file_names_out = ["./processing/picC%d.fit" % i for i in range(1, 11)]
flat_r_file_names_out = ["./processing/picC%d.fit" % i for i in range(1, 11)]

flatfield_normalisation(flat_g_file_names, flat_g_file_names_out, './processing/flat_g_master.fit')
flatfield_normalisation(flat_r_file_names, flat_r_file_names_out, './processing/flat_r_master.fit')
#
# # Align and combine the images
#
combine_images(sloan_g_file_names, "./processing/sloan_g_combined.fit")
combine_images(sloan_r_file_names, "./processing/sloan_r_combined.fit")
#
os.system("sex ./processing/sloan_g_combined.fit -CATALOG_NAME ./processing/g.cat")
os.system("sex ./processing/sloan_g_combined.fit,./processing/sloan_r_combined.fit -CATALOG_NAME ./processing/r.cat")
#
g_data = asc.read('./processing/g.cat')
r_data = asc.read('./processing/r.cat')

# # Load the data from the URAT catalogue around M67
name, right_ascension, declination = [], [], []
with open("./data/urat.txt") as f:
    for line in f:
        line_data = line.split()
        name.append(line_data[0])
        right_ascension.append(float(line_data[1]))
        declination.append(float(line_data[2]))

w = WCS("./data/new-image.fits")
xy = w.all_world2pix(right_ascension, declination, 0)

pcoo = []
for _x, _y, _name in zip(xy[0], xy[1], name):
    pcoo.append((_x, _y, _name))

urat = []
urat_with_instrumental_magnitudes = []
pcoo_new = []
with open("./data/urat.txt") as f:
    for line in f:
        lsplit = line.split()
        urat.append(lsplit)
for star in urat:
    matched = match_star_magnitudes(star[0])
    if matched == [] or matched[2] > 0.015 or matched[4] > 0.015:
        continue
    if len(star) != 7:
        continue
    new_star = star
    for data in matched[1:]:
        new_star.append(str(data))
    print(f"Matched {new_star[0]}")
    urat_with_instrumental_magnitudes.append(new_star)
    pcoo_new.append(matched[-1])

g_offsets = []
r_offsets = []
for star in urat_with_instrumental_magnitudes:
    g_offsets.append(float(star[3]) - float(star[7]))
    r_offsets.append(float(star[5]) - float(star[9]))
g_offsets = np.array(g_offsets)
r_offsets = np.array(r_offsets)
g_offset = np.mean(g_offsets)
r_offset = np.mean(r_offsets)
print("g offset:", g_offset)
print("r offset:", r_offset)

for i in range(len(g_data)):  # Add the zero point value to the source extractor data
    g_data[i]["MAG_APER"] += g_offset
    r_data[i]["MAG_APER"] += r_offset

g_r = []
g = []
for i in range(len(g_data)):
    if g_data[i]["MAGERR_APER"] > 1 or r_data[i]["MAGERR_APER"] > 1:
        print(g_data[i], r_data[i])
    #    continue  # I had quite a lot of really inaccurate data, this gets rid of it
    if g_data[i]["X_IMAGE"] > 2000 and g_data[i]["Y_IMAGE"] < 220:
        continue  # This is the area of the field that had ice on it
    if r_data[i]["MAG_APER"] > 20:
        continue
    g_r.append(g_data[i]["MAG_APER"] - r_data[i]["MAG_APER"])
    g.append(g_data[i]["MAG_APER"])
g_r = np.array(g_r)
g = np.array(g)

B_V = (g_r+0.23)/1.09
v = g - 0.6 * B_V + 0.12

MVH, BVH = np.loadtxt('./data/Hyades.txt', usecols=(1, 2), unpack=True)
MVP, BVP = np.loadtxt('./data/Pleiades.txt', usecols=(1, 2), unpack=True)

# Plot the apparent magnitude diagram
plt.plot(g_r, g, '.', label="M67")
plt.legend()
plt.xlim(-0.3, 1.8)
plt.ylim(np.max(g) + 1, np.min(g) - 1)
plt.xlabel(r'$g-r$')
plt.ylabel(r'$M_g$')
plt.title("CMD of M67 with apparent magnitudes")
plt.minorticks_on()
plt.show()

V_absolute_offset = -9.7  # This variable has to be redetermined for every cluster
# Plot the hyades and pleiades for calibration of V_absolute_offset
plt.plot(BVH, MVH, '.', label='Hyades')
plt.plot(BVP, MVP, '.', label='Pleiades')
plt.plot(B_V, g + V_absolute_offset, '.', label="M67")
plt.legend()
plt.xlim(-0.3, 1.8)
plt.ylim(12, -4)
plt.xlabel(r'$B-V$')
plt.ylabel(r'$M_V$')
plt.title("CMD of M67 with absolute magnitudes")
plt.minorticks_on()
plt.show()

iso = cmd_read.ISOCMD('./data/isochrones.cmd')
isochrones = []
ages = []
for i in reversed(range(len(iso.isocmds))):
    if i % 4 == 0 and i > 20:
        isochrones.append(iso.isocmds[i])
        ages.append(iso.ages[i])

for i in range(len(isochrones)):
    B = isochrones[i]['Bessell_B']
    V = isochrones[i]['Bessell_V']
    plt.plot(B-V, V)

plt.scatter(B_V, g + V_absolute_offset)
legend = list(map(lambda x: f"10^{x:.1f}",ages))

legend.append("Measured data")
plt.legend(legend)
#plt.xlim(0.2,1.3)
#plt.ylim(2,9)
plt.show()

D = 10 * 10 ** (-V_absolute_offset / 5)
offset_max = -10
offset_min = -9
distance_min = 10 * 10 ** (-offset_min / 5)
distance_max = 10 * 10 ** (-offset_max / 5)

# print(f"FWHM in g and r filters: {np.mean(fwhm_g)} and {np.mean(fwhm_r)}")
print(
    f"The distance to M36 is determined at {D}pc (+{distance_max - D}/-{D - distance_min}).")

r_filter_average_exposures = [231.9452680717721,230.39620414842662,233.99146127854004,230.8995351530245,230.8995351530245,233.7834886193657,233.94982367779753,234.49782868216516]
g_filter_average_exposures = [147.4994852287601,143.47954525362337,145.71146496952542,147.51761104762073,149.45035904716238,147.74415561933466,148.68092629874025,147.73124161054378]

sky_brightness_g = -2.5 * \
    np.log10(np.mean(g_filter_average_exposures)) + g_offset
sky_brightness_r = -2.5 * \
    np.log10(np.mean(r_filter_average_exposures)) + r_offset
print(f"Brightnesses. {sky_brightness_g}, {sky_brightness_r}")
