# clone the DTox git repository and change `code/dtox.py` and
# `code/dtox_learning.py` to allow selection of random seeds
git clone https://github.com/EpistasisLab/DTox.git
git -C DTox checkout -b dtox-robustness 10c909b
git -C DTox apply ../patch_seeds.diff

# setup conda environment
conda env create --name dtox_env --file=environment_dtox.yml
echo "conda activate dtox_env" >> ~/.bashrc
echo "export PYTHONPATH=/app/DTox/code:$PYTHONPATH" >> ~/.bashrc
