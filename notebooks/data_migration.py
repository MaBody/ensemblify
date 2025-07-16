import os 
import shutil
from pathlib import Path
from tqdm import tqdm
import re

# --------

names = ["ClustalO", "Mafft_FFT-NS-2", "Mafft_L-INS-i", 
"Mafft_G-INS-i", "Muscle3", "Muscle5"]
names_new = names_new = ["ClustalO", "MAFFT_FFT-NS-2", 
"MAFFT_L-INS-i", "MAFFT_G-INS-i", "Muscle3", "Muscle5"]
name_map = dict(zip(names, names_new))

regex_default = re.compile(r"^msa\.[0-9]+\..*\.fasta$")
regex_muscle5 = re.compile(r"^[a-z]{3,4}\.[0-9]+\.msa\.fasta$")

# --------

base_dir = Path("/hits/fast/cme/bodynems/data/output/treebase_v2/")
print(len(os.listdir(base_dir)))
counts = dict(zip(names, [0]*len(names)))
missing = []
for dataset in tqdm(os.listdir(base_dir)):
    data_dir = base_dir / dataset
    if os.path.isdir(data_dir):
        count = 0
        for name in names:
            ens_dir = data_dir / name
            regex = regex_muscle5 if name == "Muscle5" else regex_default
            count += len(list(filter(lambda msa_file: re.match(regex, msa_file), os.listdir(ens_dir))))
        if count != 48:
            missing.append(data_dir)
            print(data_dir)
            # for msa_name in filter(lambda msa_file: re.match(regex, msa_file), os.listdir(ens_dir)):    
            #     counts[name] += 1


# --------

dest_dir = Path("/hits/fast/cme/bodynems/data/paper/treebase_v1")
for dataset in tqdm(os.listdir(base_dir)):
    data_dir = base_dir / dataset
    if os.path.isdir(data_dir) & (data_dir not in missing):
        os.makedirs(dest_dir / dataset / "ensemble", exist_ok=True)
        
        if len(os.listdir(dest_dir / dataset / "ensemble")) == 48:
            continue
        
        for name in names:
            ens_dir = data_dir / name
            regex = regex_muscle5 if name == "Muscle5" else regex_default
            msa_files = sorted(os.listdir(ens_dir))
            
            for i, msa_name in enumerate(filter(lambda msa_file: re.match(regex, msa_file), msa_files)):
                dest_path = dest_dir / dataset / "ensemble" / ".".join([name_map[name], str(i), "fasta"])
                shutil.copy(ens_dir / msa_name, dest_path)
        
        