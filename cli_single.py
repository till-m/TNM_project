"""Command-line wrapper around the rDCM_from_fMRI function.

Expected function call:
python main.py path/to/file01.nii.gz [path/to/file02.nii.gz ...] [--scheme (harvox | yeo | aal)] [--manual]
"""
import argparse
from process import rDCM_from_fMRI

parser = argparse.ArgumentParser()
parser.add_argument('filepaths',
                    help='Path(s) to .nii or .nii.gz file(s).',
                    nargs='+')
parser.add_argument('--scheme',
                    default='yeo',
                    help='Which parcellation scheme to use.')
parser.add_argument('--manual',
                    dest='manual',
                    action='store_true',
                    help='Run matlab code manually.')
parser.set_defaults(manual=False)

args = parser.parse_args()

filepaths = args.filepaths
scheme_name = args.scheme
manual = args.manual
rDCM_from_fMRI(filepaths, meta=None, scheme=scheme_name, manual=manual)
