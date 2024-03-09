cd workflow/envs
rm -rf scenicplus
mamba create -n scenicPlus2 python=3.8 -y
bash -c ". $HOME/.bashrc 
    mamba activate scenicPlus2
    git clone https://github.com/aertslab/scenicplus
    cd scenicplus
    pip install -e .
    ln -s $(conda env list | grep scenicPlus2 | awk '{{print $NF}}')"