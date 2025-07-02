'''

conda create -n conda_gcc_env -c conda-forge gcc_linux-64
conda activate conda_gcc_env

cd build

cmake -DCMAKE_BUILD_TYPE=Release ../src


'''