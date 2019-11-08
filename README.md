# SIF tools

The public github repository was created to support a fluorescence application workshop at AGU 2018: https://agu.confex.com/agu/fm18/meetingapp.cgi/Session/54189

The main purpose of the examples and tools provided here is to help users not yet familiar with SIF data to learn about its use, the pitfalls and what programming tools are necessary. Some script programming experience is necessary to work with satellite data anyhow and we opted to used open source software, using either python and R. Jupyter notebooks for python (as well as Julia, etc) provide an ideal way to walk you through some of the code with explanations and example datasets. Note that most of the coding was done shortly before the workshop, so things mights be in a state of flux. Also, we used the GPLv3 public license for the code provided here, it would be good to make this a community resource and we would highly appreciate it if any imporvements can be shared with us so that updates will be available to the growing SIF community.

Requirements: You will need the following python libraries for sure:

h5py, netCDF4, numpy, matplotlib

----

### Download data:
#### ftp://fluo.gps.caltech.edu/data/AGU_workshop/ 
This provides a quick link to tar'ed directories (careful, large files!). You can download these and extract them on your computer.

Traditional download: Paths for all OCO-2 and TROPOMI data: ftp://fluo.gps.caltech.edu/data/<satellite>. The scripts will only focus on April through October 2018 (overlapping time period for TROPOMI and OCO-2)

----

For the python examples, please ensure that you have python3 and Jupyter notebooks installed: A basic intro is here:
https://jupyter.readthedocs.io/en/latest/install.html

----

### Building the Docker Container for Running Notebooks

navigate to the repo's dir e.g.

`cd SIF_tools`

Build a Docker image called <image name> from the Dockerfile in this repo 

`sudo docker build -t <image name> .`

Now create a container named <container name>  from that image, expose port 8888 for jupyter, and then mount your local user's home dir to /host so you can see all your code. If that port is already in use you’ll see an error that looks like: ’Bind for 0.0.0.0:8888 failed: port is already allocated.’

`sudo docker run --name <container name> -p 8888:8888 -v ~/:/host -it <image name>`

To stop the container run

`sudo docker stop <container name>` 

Or if you are in the container you can just type `exit`

To restart it run 

`sudo docker start <container name>`

And to attach the command line of this container run to your terminal run

`sudo docker attach <container name>`

Now you’re ready to develop! You can start a jupyter notebook server from the Docker container command line with:

	jupyter notebook --allow-root /host --ip 0.0.0.0

You’ll now be able to access your jupyter notebooks by navigating your browser to the address 127.0.0.1:<port> where the port is the one you specified above (e.g. 8888). The Docker container had your regular directories mounted at /host and because we’ve told the jupyter notebook with the above command to use /host as the root directory you should see all your files just like you would on your own computer but now you can access them from Docker.

### Using Conda to Run Notebooks

Conda is slightly easier but sometimes less reproducible and scalable. If you have conda installed you can create a conda environment to run these notebooks with:

`conda env create -f /tmp/environment.yml`

Then start that env by running:

`source activate earthml`


![NIR pic of Oak tree in front of Linde center](http://web.gps.caltech.edu/~cfranken/linde_small.jpg)
