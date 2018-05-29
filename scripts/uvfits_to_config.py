import argparse
import pyuvsim.simsetup

# Take a uvfits file, and save sim parameters as a yaml file.

parser = argparse.ArgumentParser(description=("Generate basic simulation parameters from uvfits."))

parser.add_argument('file_in', metavar='<FILE>', type=str, nargs='+')
parser.add_argument('-p', '--config_filename', default=None)
parser.add_argument('-t', '--telescope_config_path', default='')
parser.add_argument('-l', '--layout_csv_path', default='')

args = parser.parse_args()

pyuvsim.simsetup.uvfits_to_config_file(args.file_in[0], config_filename=args.config_filename,
                                       keep_path=False, telescope_config_path=args.telescope_config_path,
                                       layout_csv_path=args.layout_csv_path, catalog='mock', path='.')
