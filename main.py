import yaml  # pip install pyyaml
from process import batch_rDCM_from_fMRI

with open('config.yml') as config_file:
    config = yaml.safe_load(config_file)
scheme_name = config['scheme']
manual = config['manual']
data = config['data']
batch_rDCM_from_fMRI(scheme_name, manual, data)
