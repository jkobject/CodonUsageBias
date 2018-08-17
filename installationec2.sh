#installation on EC2
## do iterm install things
## installing rsub for sublime
sudo wget -O /usr/local/bin/rsub https://raw.github.com/aurora/rmate/master/rmate
sudo chmod +x /usr/local/bin/rsub
sudo yum install R
unset pip pip3
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
## some prompts will appear
conda install pip
conda install jupyter notebook
git clone https://github.com/jkobject/PyCUB.git
cd PyCUB
pip install -r requirements.txt

git clone https://github.com/DmitryUlyanov/Multicore-TSNE.git
cd Multicore-TSNE/
pip install .
#if you find any problem with rpy2, either you don't have R installed, either you have but it is not the latest version
#if the problem is about gcc compiler :
ln -s /usr/lib/gcc/$os/$version/libgomp.spec /usr/lib64/libgomp.spec
ln -s /usr/lib/gcc/$os/$version/libgomp.a /usr/lib64/libgomp.a
ln -s /usr/lib64/libgomp.so.1.0.0 /usr/lib64/libgomp.so
#
cd ../..
mkdir ssl
cd ssl
sudo openssl req -x509 -nodes -days 365 -newkey rsa:1024 -keyout "cert.key" -out "cert.pem" -batch
ipython
from IPython.lib import passwd
-> passwd()
exit
vi ~/.jupyter/jupyter_notebook_config.py
c = get_config()  # Get the config object.
c.NotebookApp.certfile = u'/home/ubuntu/ssl/cert.pem' # path to the certificate we generated
c.NotebookApp.keyfile = u'/home/ubuntu/ssl/cert.key' # path to the certificate key we generated
c.IPKernelApp.pylab = 'inline'  # in-line figure when using Matplotlib
c.NotebookApp.ip = '*'  # Serve notebooks locally.
c.NotebookApp.open_browser = False  # Do not open a browser window by default when using notebooks.
c.NotebookApp.password = 'sha1:fc216:3a35a98ed980b9...'

## you can add more memory with 
swapon -s
lsblk
sudo mkswap /dev/$disk
sudo swapon /dev/$disk
# and add to /etc/fstab
/dev/$disk       none    swap    sw  0       0

## to load stuff from the server
scp -i "~/Documents/Projects/amzkey/AWS-2.pem" ec2-user@ec2-52-213-244-82.eu-west-1.compute.amazonaws.com:PyCUB/CodonUsageBias/utils/$FILE utils/$FILE

## to continue the work while disconnected 
disown -h %PAUSEDJOBNUMBER









