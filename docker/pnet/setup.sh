# install dependencies
apt-get update
apt-get install -y --no-install-recommends unzip

# clone the PNET git repository and change `pipeline/one_split.py` and
# `train/run_me.py` to allow selection of random seeds
git clone https://github.com/marakeby/pnet_prostate_paper.git
git -C pnet_prostate_paper checkout -b pnet-robustness 2b16264
git -C pnet_prostate_paper apply ../patch_seeds.diff

# setup conda environment
conda env create --name pnet_env --file=environment_pnet.yml
echo "conda activate pnet_env" >> ~/.bashrc
echo "export PYTHONPATH=/app/pnet_prostate_paper:$PYTHONPATH" >> ~/.bashrc

# download PNET data
source /opt/conda/etc/profile.d/conda.sh
conda activate pnet_env
gdown 17nssbdUylkyQY1ebtxsIw5UzTAd0zxWb
unzip _database.zip "_database/*" -d pnet_prostate_paper
rm _database.zip
